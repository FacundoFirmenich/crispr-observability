import {createHash} from 'node:crypto';
import {createReadStream, openAsBlob} from 'node:fs';
import {open, realpath, stat} from 'node:fs/promises';
import {basename, isAbsolute, relative, resolve, sep} from 'node:path';
import {z} from 'zod';
import {boundary, compare, revcomp} from './science.js';
import {readCase} from './import-case.js';
import {openIndexedBam} from './bam-input.js';
import {alignmentTsv, inspectQuery, MAX_CASE_BYTES, MAX_LOCUS_WIDTH, parseSam,
  reportAlignment, validateAlignmentCase} from './alignment.js';

const id = z.string().regex(/^[A-Za-z][A-Za-z0-9_.-]{0,63}$/);
const path = z.string().min(1).max(500).regex(/^[^:\\\x00-\x1f]+$/)
  .refine(s => s.split('/').every(p => p && p !== '.' && p !== '..'), 'Use a relative path without dot segments.');
const sha = z.string().regex(/^[a-f0-9]{64}$/);
const coord = z.number().int().min(0).max(2147483647);
const span = z.tuple([coord, coord]).refine(([a,b]) => a < b, 'End must exceed start.');
const file = z.strictObject({path, sha256:sha.optional()});
const single = type => z.strictObject({id,type:z.literal(type),file});
const dataset = z.discriminatedUnion('type', [
  single('sequence'), single('structure'), single('alignment'),
  z.strictObject({id,type:z.literal('indexed-bam'),bam:file,index:file})
]);
const window = z.strictObject({ref:z.string().min(1).max(1000),span})
  .refine(x=>x.span[1]-x.span[0]<=MAX_LOCUS_WIDTH, 'Locus exceeds the 50,000-base limit.');
const operand = z.strictObject({input:id,span,reverse_complement:z.boolean().default(false)});
const job = z.discriminatedUnion('op', [
  z.strictObject({id,op:z.literal('inspect'),input:id}),
  z.strictObject({id,op:z.literal('anchor-containment'),anchor:span,shift:z.number().int().min(-2147483647).max(2147483647),read_end:coord}),
  z.strictObject({id,op:z.literal('compare-sequence-spans'),a:operand,b:operand}),
  z.strictObject({id,op:z.literal('atom-distance'),input:id,atoms:z.tuple([coord,coord])}),
  z.strictObject({id,op:z.literal('alignment-locus'),input:id,locus:window,
    alignment_index:coord.default(0),query_index:coord.optional(),original_orientation:z.boolean().default(false),
    min_mapq:z.number().int().min(0).max(254).default(0),include_non_primary:z.boolean().default(false)})
]);
export const projectSchema=z.strictObject({
  schema:z.literal('observability-project/1'),project_id:id,
  datasets:z.array(dataset).max(32),jobs:z.array(job).min(1).max(64)
});

const digest=bytes=>createHash('sha256').update(bytes).digest('hex');
function canonical(v){
  if(Array.isArray(v))return v.map(canonical);
  if(v!==null&&typeof v==='object')
    return Object.fromEntries(Object.keys(v).sort().map(k=>[k,canonical(v[k])]));
  return v;
}
export function validateProject(input){
  const plan=projectSchema.parse(input);
  const types=new Map(),jobIds=new Set();
  for(const d of plan.datasets){
    if(types.has(d.id))throw new Error('Duplicate dataset '+d.id);
    types.set(d.id,d.type);
  }
  function requireDataset(id,allowed){
    const type=types.get(id);
    if(!type)throw new Error('Unknown dataset '+id);
    if(!allowed.includes(type))throw new Error('Dataset '+id+' has incompatible type '+type);
  }
  for(const j of plan.jobs){
    if(jobIds.has(j.id))throw new Error('Duplicate operation ID '+j.id);
    jobIds.add(j.id);
    if(j.op==='inspect')requireDataset(j.input,['sequence','structure','alignment','indexed-bam']);
    if(j.op==='atom-distance')requireDataset(j.input,['structure']);
    if(j.op==='alignment-locus')requireDataset(j.input,['alignment','indexed-bam']);
    if(j.op==='compare-sequence-spans'){
      requireDataset(j.a.input,['sequence']);requireDataset(j.b.input,['sequence']);
    }
    if(j.op==='anchor-containment')boundary(...j.anchor,j.shift,j.read_end);
  }
  return plan;
}
function within(root,candidate){
  const sub=relative(root,candidate);
  return Boolean(sub)&&!isAbsolute(sub)&&sub!=='..'&&!sub.startsWith('..'+sep);
}
export async function workspaceFile(root,requested){
  const safe=path.parse(requested);
  const base=await realpath(root);
  const lexical=resolve(base,safe);
  if(!within(base,lexical))throw new Error('Input outside the authorized workspace');
  const physical=await realpath(lexical);
  if(!within(base,physical)||(await stat(physical)).isFile()===false)
    throw new Error('Input outside the authorized workspace or not a regular file');
  return physical;
}
function checkDigest(spec,actual){
  if(spec.sha256&&actual!==spec.sha256)throw new Error('SHA-256 mismatch for '+spec.path);
}
async function boundedBytes(path,limit){
  const file=await open(path,'r');
  try{
    const info=await file.stat();
    if(info.size>limit)throw new Error('Input exceeds '+limit+' bytes');
    const buffer=Buffer.allocUnsafe(info.size+1);
    let offset=0;
    while(offset<buffer.length){
      const {bytesRead}=await file.read(buffer,offset,buffer.length-offset,null);
      if(!bytesRead)break;
      offset+=bytesRead;
    }
    if(offset!==info.size)throw new Error('Input size changed while reading');
    return buffer.subarray(0,offset);
  }finally{await file.close();}
}
export async function readManifest(root,path){
  const f=await workspaceFile(root,path);
  return JSON.parse(await boundedBytes(f,1024*1024));
}
async function indexedFile(root,spec){
  const full=await workspaceFile(root,spec.path);
  const blob=Object.defineProperty(await openAsBlob(full),'name',{value:basename(full)});
  let sha256=null;
  if(spec.sha256){
    const h=createHash('sha256');
    for await(const block of createReadStream(full))h.update(block);
    sha256=h.digest('hex');checkDigest(spec,sha256);
  }
  return {blob,source:{path:spec.path,bytes:blob.size,sha256}};
}
async function loadDataset(root,d){
  if(d.type==='indexed-bam'){
    const bam=await indexedFile(root,d.bam);
    const index=await indexedFile(root,d.index);
    const reader=await openIndexedBam(bam.blob,index.blob);
    return {type:d.type,reader,source:{bam:bam.source,index:index.source},
      summary:{references:Object.keys(reader.refLengths).length}};
  }
  const physical=await workspaceFile(root,d.file.path);
  const bytes=await boundedBytes(physical,MAX_CASE_BYTES);
  const sha256=digest(bytes);
  checkDigest(d.file,sha256);
  const text=new TextDecoder('utf-8',{fatal:true}).decode(bytes);
  const source={path:d.file.path,bytes:bytes.length,sha256};
  let data;
  if(d.type==='alignment'){
    data=basename(physical).toLowerCase().endsWith('.sam')?
      parseSam(basename(physical),text,{inputSHA256:sha256}):validateAlignmentCase(JSON.parse(text));
  }else{
    data=readCase(basename(physical),text);
    if(data.kind!==d.type)throw new Error('Parsed '+data.kind+' but declared '+d.type);
  }
  const summary=d.type==='sequence'?{format:data.format,bases:data.sequence.length}:
    d.type==='structure'?{format:data.format,units:data.units,atoms:data.atoms.length,bonds:data.bonds.length}:
      {format:data.format,records:data.records.length};
  return {type:d.type,data,source,summary};
}
function sequenceOperand(dataset,operand){
  const seq=dataset.data.sequence;
  const [a,b]=operand.span;
  if(b>seq.length)throw new Error('Span outside sequence '+operand.input);
  const bases=seq.slice(a,b);
  return {input:operand.input,span:operand.span,orientation:operand.reverse_complement?'reverse-complement':'forward',
    length:b-a,sequence:operand.reverse_complement?revcomp(bases):bases};
}
async function compute(j,get){
  if(j.op==='anchor-containment')return boundary(...j.anchor,j.shift,j.read_end);
  if(j.op==='compare-sequence-spans'){
    const a=sequenceOperand(await get(j.a.input),j.a);
    const b=sequenceOperand(await get(j.b.input),j.b);
    let overlap=null;
    if(a.input===b.input){
      const seq=(await get(a.input)).data.sequence;
      overlap=compare({domain:[0,seq.length],length:seq.length,segment:seq},a.span,b.span).overlap;
    }
    return {a,b,identicalString:a.sequence===b.sequence,sharedCoordinateSpace:a.input===b.input,overlap};
  }
  const d=await get(j.input);
  if(j.op==='inspect')return d.summary;
  if(j.op==='atom-distance'){
    const [a,b]=j.atoms.map(i=>d.data.atoms[i]);
    if(!a||!b)throw new Error('Atom index outside normalized atom list');
    return {units:'angstrom',indices:j.atoms,atoms:[a,b],
      distance:Math.hypot(...a.xyz.map((x,i)=>x-b.xyz[i]))};
  }
  const [lo,hi]=j.locus.span;
  const data=d.reader?await d.reader.query(j.locus.ref,lo,hi):d.data;
  const {read,coverage}=reportAlignment(data,{index:j.alignment_index,lo,hi,minMapq:j.min_mapq,
    includeNonPrimary:j.include_non_primary});
  if(read.ref!==j.locus.ref)throw new Error('Selected alignment is on '+read.ref+', not '+j.locus.ref);
  const query=j.query_index===undefined?null:inspectQuery(read,j.query_index,j.original_orientation);
  const {name,flag,mapq,ref,start,referenceEnd,cigar,ordinal}=read;
  return {source:canonical(data.source),read:{name,flag,mapq,ref,start,referenceEnd,cigar,ordinal},
    query,coverage,evidence_tsv:alignmentTsv(read,coverage)};
}
export async function runProject(input,root){
  const plan=validateProject(input);
  const definitions=new Map(plan.datasets.map(x=>[x.id,x]));
  const loaded=new Map();
  const get=async id=>{
    if(!loaded.has(id)){
      try{loaded.set(id,await loadDataset(root,definitions.get(id)));}
      catch(error){throw new Error('Dataset '+id+': '+error.message,{cause:error});}
    }
    return loaded.get(id);
  };
  const results=[];
  for(const j of plan.jobs){
    try{results.push({id:j.id,op:j.op,parameters:j,result:await compute(j,get)});}
    catch(error){throw new Error('Operation '+j.id+': '+error.message,{cause:error});}
  }
  return {schema:'observability-project-report/1',engine_version:'0.4.0',project_id:plan.project_id,
    manifestSHA256:digest(JSON.stringify(canonical(plan))),
    datasets:plan.datasets.map(x=>{
      const d=loaded.get(x.id);
      return {id:x.id,type:x.type,status:d?'used':'not-read',...(d?{source:d.source,summary:d.summary}:{})};
    }),results};
}
