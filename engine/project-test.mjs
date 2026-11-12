import assert from 'node:assert/strict';
import {readFile,writeFile,mkdir,copyFile} from 'node:fs/promises';
import {join} from 'node:path';
import {execFile} from 'node:child_process';
import {promisify} from 'node:util';
import {fileURLToPath} from 'node:url';
import {Client} from '@modelcontextprotocol/sdk/client/index.js';
import {StdioClientTransport} from '@modelcontextprotocol/sdk/client/stdio.js';
import {runProject,validateProject} from './project-engine.mjs';

const [root]=process.argv.slice(2);
const file=(path,sha256)=>({path,...(sha256?{sha256}:{})});
const node=fileURLToPath(new URL('.',import.meta.url));
const run=promisify(execFile);
const oracle=JSON.parse(await readFile(join(root,'ALIGNMENT_FIXTURE_ORACLE.json'),'utf8'));
const [hg,chm]=oracle.checks;
const a=hg.records[0].sequence.slice(0,42);
const b=chm.records[0].sequence.slice(0,42);
assert.match(a,/^[ACGTN]{42}$/);
assert.match(b,/^[ACGTN]{42}$/);
await writeFile(join(root,'hg.fa'),'>HG00096 first 42 stored bases\n'+a+'\n');
await writeFile(join(root,'chm.fa'),'>CHM1 first 42 stored bases\n'+b+'\n');
const guideInput={
  schema:'observability-project/1',project_id:'public-cross-input-check',
  datasets:[
    {id:'hg-sequence',type:'sequence',file:file('hg.fa')},
    {id:'chm-sequence',type:'sequence',file:file('chm.fa')},
    {id:'structure',type:'structure',file:file('qa-1BNA.cif')},
    {id:'sam',type:'alignment',file:file(hg.sam)},
    {id:'bam',type:'indexed-bam',bam:file(hg.bam,hg.bamSHA256),index:file(hg.bam+'.bai')},
    {id:'unused',type:'sequence',file:file('never-open-this-file.fa')}
  ],
  jobs:[
    {id:'compare',op:'compare-sequence-spans',a:{input:'hg-sequence',span:[0,12]},
      b:{input:'chm-sequence',span:[0,12],reverse_complement:true}},
    {id:'same',op:'compare-sequence-spans',a:{input:'hg-sequence',span:[1,7]},
      b:{input:'hg-sequence',span:[5,10]}},
    {id:'sam-depth',op:'alignment-locus',input:'sam',locus:{ref:hg.ref,span:[hg.lo,hg.hi]},
      min_mapq:0,query_index:0,original_orientation:true},
    {id:'bam-depth',op:'alignment-locus',input:'bam',locus:{ref:hg.ref,span:[hg.lo,hg.hi]},
      min_mapq:0,query_index:0,original_orientation:true},
    {id:'structure-distance',op:'atom-distance',input:'structure',atoms:[0,1]},
    {id:'threshold',op:'anchor-containment',anchor:[210,222],shift:3,read_end:224},
    {id:'inspect',op:'inspect',input:'hg-sequence'}
  ]
};
const out=await runProject(guideInput,root);
const [cross,same,sam,bam,structure,anchor]=out.results.map(x=>x.result);
assert.equal(cross.overlap,null);
assert.equal(cross.a.sequence,a.slice(0,12));
const reverse=[...b.slice(0,12)].reverse().map(c=>({A:'T',C:'G',G:'C',T:'A',N:'N'})[c]).join('');
assert.equal(cross.b.sequence,reverse);
assert.equal(same.overlap,2);
assert.deepEqual(sam.coverage,bam.coverage);
assert.equal(sam.coverage.covered,67);
assert.equal(sam.query.reference,10065);
assert.equal(bam.query.reference,10065);
assert.equal(sam.evidence_tsv,bam.evidence_tsv);
assert.equal(structure.units,'angstrom');
assert.ok(structure.distance>1 && structure.distance<4);
assert.equal(anchor.observed,11);
assert.equal(out.results.at(-1).result.bases,42);
assert.equal(out.datasets.at(-1).status,'not-read');
assert.ok(out.datasets.filter(x=>x.status==='used').every(x=>x.source));

const rerun=await runProject(JSON.parse(JSON.stringify(guideInput)),root);
assert.deepEqual(out,rerun);
const second={
  schema:'observability-project/1',project_id:'second-study',
  datasets:[{id:'hg-sequence',type:'sequence',file:file('chm.fa')}],
  jobs:[{id:'compare',op:'compare-sequence-spans',
    a:{input:'hg-sequence',span:[0,6]},b:{input:'hg-sequence',span:[2,8]}}]
};
const secondRoot=join(root,'second-study');
await mkdir(secondRoot,{recursive:true});
await copyFile(join(root,'chm.fa'),join(secondRoot,'chm.fa'));
await writeFile(join(secondRoot,'project.json'),JSON.stringify(second));
const secondOut=await runProject(second,secondRoot);
assert.equal(secondOut.project_id,'second-study');
assert.equal(secondOut.results[0].result.a.sequence,b.slice(0,6));
assert.notEqual(secondOut.manifestSHA256,out.manifestSHA256);
assert.equal((await runProject(guideInput,root)).results[0].result.a.sequence,a.slice(0,12));
const clone=()=>structuredClone(guideInput);
let bad=clone();bad.jobs[0].a.input='unknown';
assert.throws(()=>validateProject(bad),/Unknown dataset/);
bad=clone();bad.jobs[0].a.input='structure';
assert.throws(()=>validateProject(bad),/incompatible type/);
bad=clone();bad.jobs[1].id=bad.jobs[0].id;
assert.throws(()=>validateProject(bad),/Duplicate operation/);
bad=clone();bad.datasets[0].file.sha256='0'.repeat(64);
await assert.rejects(runProject(bad,root),/SHA-256 mismatch/);
bad=clone();bad.jobs[2].locus.ref='missing-reference';
await assert.rejects(runProject(bad,root),/not missing-reference/);
bad=clone();bad.datasets[0].file.path='../outside.fa';
await assert.rejects(runProject(bad,root),/without dot segments/);

const manifestPath=join(root,'project.example.json');
await writeFile(manifestPath,JSON.stringify(guideInput,null,2)+'\n');
const cli=await run(process.execPath,[join(node,'project-cli.mjs'),'run',manifestPath]);
assert.deepEqual(JSON.parse(cli.stdout),out);

const client=new Client({name:'multi-input-smoke',version:'0.4.0'});
const transport=new StdioClientTransport({command:process.execPath,
  args:[join(node,'project-mcp.mjs'),root],stderr:'inherit'});
await client.connect(transport);
try{
  const {tools}=await client.listTools();
  assert.deepEqual(tools.map(x=>x.name),['lab_project_run']);
  const result=await client.callTool({name:'lab_project_run',arguments:{manifest_path:'project.example.json'}});
  assert.ok(!result.isError, result.content?.[0]?.text);
  assert.deepEqual(result.structuredContent,out);
  const nested=await client.callTool({name:'lab_project_run',arguments:{manifest_path:'second-study/project.json'}});
  assert.ok(!nested.isError, nested.content?.[0]?.text);
  assert.deepEqual(nested.structuredContent,secondOut);
}finally{await client.close();}
const proof={engineVersion:out.engine_version,actualInputs:out.datasets.filter(x=>x.status==='used').length,
  declaredInputs:out.datasets.length,operations:out.results.length,manifestSHA256:out.manifestSHA256,
  independentProjects:true,sameTsv:true,crossDatasetOverlap:cross.overlap,sameDatasetOverlap:same.overlap,
  coveredReferenceBases:sam.coverage.covered,originalIndexZeroReference:sam.query.reference,structureAngstrom:structure.distance,
  exactApiCliMcp:true,stableRerun:true,rejections:['unknown dataset','incompatible type','duplicate operation','checksum mismatch','reference mismatch','path escape']};
await writeFile(join(root,'PROJECT_V04_CHECKS.json'),JSON.stringify(proof,null,2)+'\n');
await writeFile(join(root,'PROJECT_V04_RESULT.json'),JSON.stringify(out,null,2)+'\n');
console.log(JSON.stringify(proof));
