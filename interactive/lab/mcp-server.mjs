import {McpServer} from '@modelcontextprotocol/sdk/server/mcp.js';
import {StdioServerTransport} from '@modelcontextprotocol/sdk/server/stdio.js';
import {z} from 'zod';
import {chromium} from 'playwright';
import {createHash} from 'node:crypto';
import {access,mkdtemp,readFile,realpath,stat,writeFile} from 'node:fs/promises';
import {openAsBlob} from 'node:fs';
import {dirname,isAbsolute,join,relative,resolve,sep} from 'node:path';
import {fileURLToPath,pathToFileURL} from 'node:url';
import {boundary,compare} from './science.js';
import {readCase} from './import-case.js';
import {alignmentTsv,inspectQuery,parseSam,reportAlignment,validateAlignmentCase,MAX_REPLAY_BYTES} from './alignment.js';
import {openIndexedBam} from './bam-input.js';

const own=dirname(fileURLToPath(import.meta.url));
const [workspaceArgument]=process.argv.slice(2);
if(!workspaceArgument)throw new Error('Pass the authorized workspace directory as the first argument.');
const workspace=await realpath(resolve(workspaceArgument));
if(!(await stat(workspace)).isDirectory())throw new Error('The workspace must be a directory.');
const tryApp=join(own,'CFTR_OBSERVABILITY_LAB.html');
const app=await access(tryApp).then(()=>tryApp,()=>join(own,'..','index.html'));
await access(app);
const {audit}=JSON.parse(await readFile(join(own,'data.json'),'utf8'));
const hash=bytes=>createHash('sha256').update(bytes).digest('hex');
const server=new McpServer({name:'read-observability-lab',version:'0.3.0'});
const output=result=>({content:[{type:'text',text:JSON.stringify(result)}],structuredContent:result});
const int=z.number().int().min(0).max(2_000_000);
const coord=z.number().int().min(0).max(2**31-1);
const interval=z.tuple([int,int]);
const fileInput={
  case_path:z.string().max(1000).optional(),
  case_name:z.string().max(240).optional(),
  case_text:z.string().max(6_000_000).optional()
};

async function localFile(path){
  if(typeof path!=='string'||isAbsolute(path))throw new Error('Use a workspace-relative file path.');
  const lexical=resolve(workspace,path);
  const lexicalRelative=relative(workspace,lexical);
  if(!lexicalRelative||isAbsolute(lexicalRelative)||lexicalRelative==='..'||lexicalRelative.startsWith('..'+sep))
    throw new Error('The case file must resolve inside the authorized workspace.');
  const candidate=await realpath(lexical);
  const rel=relative(workspace,candidate);
  if(!rel||isAbsolute(rel)||rel==='..'||rel.startsWith('..'+sep))
    throw new Error('The case file must resolve inside the authorized workspace.');
  if(!(await stat(candidate)).isFile())throw new Error('Choose a regular workspace file.');
  return candidate;
}

async function openCase({case_path,case_name,case_text}){
  if(case_path!==undefined){
    if(case_text!==undefined||case_name!==undefined||isAbsolute(case_path))
      throw new Error('Use one workspace-relative case_path or a case_name with case_text.');
    const candidate=await localFile(case_path);
    const info=await stat(candidate);
    if(!info.isFile()||info.size>6_000_000)throw new Error('Use a regular case file up to 6 MB.');
    const bytes=await readFile(candidate);
    if(bytes.length>6_000_000)throw new Error('Case file exceeded the 6 MB limit.');
    return {name:case_path.split(/[/\\]/).at(-1),body:bytes.toString('utf8'),sha256:hash(bytes)};
  }
  if(case_name===undefined||case_text===undefined||Buffer.byteLength(case_text,'utf8')>6_000_000)
    throw new Error('Supply case_name and case_text, or a workspace-relative case_path.');
  return {name:case_name,body:case_text,sha256:hash(case_text)};
}

function tool(name,description,inputSchema,handler){
  server.registerTool(name,{description,inputSchema},async args=>{
    try{return output(await handler(args));}
    catch(error){return {isError:true,content:[{type:'text',text:error.message}]};}
  });
}

tool('lab_boundary','Calculate finite-read anchor containment on zero-based half-open intervals.',
  {start:int,end:int,shift:z.number().int().min(-1_000_000).max(1_000_000),read_end:int},a=>
    ({inputs:a,result:boundary(a.start,a.end,a.shift,a.read_end)}));

tool('lab_compare_spans','Extract and compare two intervals on an archived record or an uploaded DNA case.',
  {...fileInput,archive_read:z.enum(audit.cases.map(r=>r.id)).optional(),a:interval,b:interval,
    reverse_complement:z.boolean().default(false)},async args=>{
    const hasFile=args.case_path!==undefined||args.case_text!==undefined||args.case_name!==undefined;
    if(Number(hasFile)+Number(Boolean(args.archive_read))!==1)throw new Error('Choose exactly one archive read or input case.');
    let record,sha256;
    if(args.archive_read)record=audit.cases.find(r=>r.id===args.archive_read);
    else{
      const input=await openCase(args);
      sha256=input.sha256;
      const normalized=readCase(input.name,input.body);
      if(normalized.kind!=='sequence')throw new Error('Span comparison requires a sequence case.');
      record={id:normalized.name,length:normalized.sequence.length,
        domain:[0,normalized.sequence.length],segment:normalized.sequence};
    }
    return {record:record.id,source:sha256?{sha256}:{archivedCommit:audit.commit},
      result:compare(record,args.a,args.b,args.reverse_complement)};
  });

tool('lab_normalize_case','Read a supported case and report its format, size and SHA-256.',
  fileInput,async args=>{
    const input=await openCase(args),result=readCase(input.name,input.body);
    return {sha256:input.sha256,format:result.format,kind:result.kind,
      ...(result.kind==='structure'?{atoms:result.atoms.length,bonds:result.bonds.length,
        chains:[...new Set(result.atoms.map(a=>a.chain))]}:{bases:result.sequence.length})};
  });

tool('lab_alignment_locus','Audit a SAM slice, an alignment case JSON, or a locus in local BAM + BAI/CSI. Save exact CIGAR and depth evidence.',
  {...fileInput,bam_path:z.string().max(1000).optional(),index_path:z.string().max(1000).optional(),
    reference:z.string().max(1000).optional(),locus_start:coord.optional(),locus_end:coord.optional(),
    alignment_index:int.default(0),query_index:coord.default(0),original_orientation:z.boolean().default(false),
    coverage_start:coord.optional(),coverage_end:coord.optional(),min_mapq:z.number().int().min(0).max(254).default(0),
    include_non_primary:z.boolean().default(false)},async args=>{
    let input;
    if(args.bam_path!==undefined||args.index_path!==undefined){
      if(!args.bam_path||!args.index_path||!args.reference||args.locus_start===undefined||args.locus_end===undefined||
        args.case_path!==undefined||args.case_text!==undefined||args.case_name!==undefined)
        throw new Error('Supply either a BAM + index + reference locus or one SAM/JSON case.');
      const blob=async path=>Object.defineProperty(await openAsBlob(await localFile(path)),'name',
        {value:path.split(/[/\\]/).at(-1)});
      const [bam,index]=await Promise.all([blob(args.bam_path),blob(args.index_path)]);
      const reader=await openIndexedBam(bam,index);
      input=await reader.query(args.reference,args.locus_start,args.locus_end);
    }else{
      const raw=await openCase(args);
      input=raw.name.toLowerCase().endsWith('.sam')?parseSam(raw.name,raw.body,{inputSHA256:raw.sha256}):
        validateAlignmentCase(JSON.parse(raw.body));
    }
    const settings={index:args.alignment_index,at:args.query_index,original:args.original_orientation,
      lo:args.coverage_start,hi:args.coverage_end,minMapq:args.min_mapq,includeNonPrimary:args.include_non_primary};
    const {read,coverage}=reportAlignment(input,settings);
    const position=inspectQuery(read,settings.at,settings.original);
    const {runs,...summary}=coverage;
    settings.lo=coverage.window[0];settings.hi=coverage.window[1];
    const report={schema:'read-observability-result/1',appVersion:'0.3.0',kind:'alignment_locus_audit',
      portableInput:input,inputs:settings,result:{read,position,coverage}};
    const files={'alignment-replay.json':JSON.stringify(report,null,2)+'\n','alignment-evidence.tsv':alignmentTsv(read,coverage)};
    const manifest={appVersion:report.appVersion,source:input.source,
      files:Object.fromEntries(Object.entries(files).map(([name,contents])=>[name,hash(contents)]))};
    const out=await mkdtemp(join(workspace,'alignment-evidence-'));
    await Promise.all(Object.entries({...files,'manifest.json':JSON.stringify(manifest,null,2)+'\n'})
      .map(([name,body])=>writeFile(join(out,name),body,{flag:'wx'})));
    return {folder:relative(workspace,out),manifest,read:read.name,position,coverage:summary};
  });

tool('lab_render','Replay a lab session or open a case, then save a zoom-preserving PNG, replay JSON and hash manifest in the authorized workspace.',
  {...fileInput,preset:z.enum(['edge','span','same-label','atom']).optional(),
    replay:z.object({schema:z.literal('read-observability-result/1')}).passthrough().optional(),
    scale:z.union([z.literal(1),z.literal(2),z.literal(4)]).default(2)},async args=>{
    const hasFile=args.case_path!==undefined||args.case_text!==undefined||args.case_name!==undefined;
    if(Number(hasFile)+Number(Boolean(args.preset))+Number(Boolean(args.replay))!==1)
      throw new Error('Choose exactly one preset, replay or case input.');
    const caseFile=hasFile?await openCase(args):null;
    if(args.replay&&Buffer.byteLength(JSON.stringify(args.replay),'utf8')>MAX_REPLAY_BYTES)
      throw new Error('The replay exceeds the 24 MB limit.');
    const browser=await chromium.launch({channel:'chrome',headless:true,args:['--enable-unsafe-swiftshader']});
    try{
      const page=await browser.newPage({viewport:{width:1440,height:960}});
      await page.route(/^https?:/,route=>route.abort());
      await page.goto(pathToFileURL(app).href,{waitUntil:'load'});
      await page.waitForFunction(()=>window.readObservabilityLab&&document.querySelector('#metrics').children.length>0);
      const {result,png,tsv}=await page.evaluate(async ({preset,replay,input,scale})=>{
        const lab=window.readObservabilityLab;
        const result=preset?await lab.preset(preset):replay?await lab.replay(replay):await lab.open(input.name,input.body);
        return {result,png:lab.renderPng(scale),tsv:lab.alignmentTsv()};
      },{preset:args.preset,replay:args.replay,input:caseFile,scale:args.scale});
      const bytes=Buffer.from(png.dataUrl.split(',')[1],'base64');
      if(bytes.length<1024||bytes.readUInt32BE(0)!==0x89504e47)throw new Error('Browser did not return a valid PNG.');
      const out=await mkdtemp(join(workspace,'read-observability-'));
      const resultJson=JSON.stringify(result,null,2)+'\n';
      const files={'scene.png':bytes,'replay.json':resultJson,...(tsv?{'alignment-evidence.tsv':tsv}:{})};
      const manifest={appVersion:result.appVersion,width:png.width,height:png.height,scale:args.scale,
        files:Object.fromEntries(Object.entries(files).map(([name,body])=>[name,hash(body)])),
        inputSHA256:caseFile?.sha256??null};
      await Promise.all(Object.entries({...files,'manifest.json':JSON.stringify(manifest,null,2)+'\n'})
        .map(([name,body])=>writeFile(join(out,name),body,{flag:'wx'})));
      const summary=result.kind==='alignment_locus_audit'?
        {read:result.result.read,position:result.result.position,covered:result.result.coverage.covered}:
        result.result;
      return {folder:relative(workspace,out),...manifest,kind:result.kind,result:summary};
    }finally{await browser.close();}
  });

await server.connect(new StdioServerTransport());
