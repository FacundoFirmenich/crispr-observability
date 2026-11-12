import assert from 'node:assert/strict';
import {writeFile,readFile} from 'node:fs/promises';
import {join} from 'node:path';
import {runProject,validateProject} from './project-engine.mjs';
import {parseSam} from './alignment.js';

const [root]=process.argv.slice(2);
const locus={ref:'2',span:[90,160]};
const junctions=[[104,114],[118,138]];
const expected={
  'complete':['observed','observed'],
  'split-mates':['not-observed-in-input','not-observed-in-input'],
  'split-alternatives':['not-observed-in-input','not-observed-in-input'],
  'missing-nh':['observed','indeterminate-nh'],
  'multi':['observed','not-observed-in-input'],
  'secondary-only':['observed','not-observed-in-input'],
  'deletion':['not-observed-in-input','not-observed-in-input'],
  'competing':['observed','indeterminate-primary'],
  'absent':['not-observed-in-input','not-observed-in-input'],
};
const inputs=[
  {id:'sam',type:'alignment',file:{path:'junction-fixture.sam'}},
  {id:'bam',type:'indexed-bam',bam:{path:'junction-fixture.bam'},index:{path:'junction-fixture.bam.bai'}}
];
const job=(input,read_name)=>({
  id:input+'-'+read_name,op:'junction-evidence',input,locus,read_name,mate:'first',junctions
});
const manifest={
  schema:'observability-project/1',project_id:'junction-regression',
  datasets:inputs,
  jobs:inputs.flatMap(({id})=>Object.keys(expected).map(name=>job(id,name)))
};
const report=await runProject(manifest,root);
const byId=Object.fromEntries(report.results.map(x=>[x.id,x.result]));
for(const [name,states] of Object.entries(expected)){
  for(const id of ['sam','bam']){
    const {endpoints}=byId[id+'-'+name];
    assert.deepEqual([endpoints.exactAnyReported,endpoints.exactUniquePrimary],states);
  }
  const normalized=({records,endpoints,localRecordCount})=>({
    records:records.map(({ordinal,...r})=>r),endpoints,localRecordCount
  });
  assert.deepEqual(normalized(byId['sam-'+name]),normalized(byId['bam-'+name]));
}
const [{mapqRaw,mappingQuality,alignmentTagsAndMate,skips}]=byId['bam-complete'].records;
assert.equal(mapqRaw,255);
assert.equal(mappingQuality,null);
assert.deepEqual(alignmentTagsAndMate,{nh:1,hi:1,rg:'g1',mateRef:'1',mateStart:200,templateLength:-52});
assert.deepEqual(skips.map(r=>r.flankingPhred),[[40,40],[40,40]]);
assert.deepEqual(byId['bam-missing-nh'].records[0].skips[0].flankingPhred,[null,null]);
assert.equal(byId['sam-secondary-only'].endpoints.exactUniqueReported,'observed');
assert.equal(byId['bam-secondary-only'].endpoints.exactUniqueReported,'observed');
assert.equal(byId['sam-missing-nh'].endpoints.exactUniqueReported,'indeterminate-nh');
assert.equal(byId['sam-deletion'].records[0].skips.length,0);
assert.equal(byId['sam-split-mates'].records.length,2);
assert.equal(byId['bam-split-alternatives'].records.length,2);
assert.equal(byId['bam-absent'].localRecordCount,0);

const ambiguous=job('sam','ambiguous');
await assert.rejects(runProject({...manifest,jobs:[ambiguous]},root),/multiple read groups/);
const explicit=await runProject({...manifest,jobs:[{...ambiguous,read_group:'g2'}]},root);
assert.equal(explicit.results[0].result.localRecordCount,1);
assert.throws(()=>validateProject({...manifest,jobs:[{...ambiguous,junctions:[[104,114],[104,114]]}]}),/nonoverlapping/);
assert.throws(()=>validateProject({...manifest,jobs:[{...ambiguous,junctions:[[104,114],[155,165]]}]}),/inside the locus/);
const emptyJob={...job('bam','complete'),locus:{ref:'2',span:[390,450]},junctions:[[404,414],[418,438]]};
const empty=await runProject({...manifest,jobs:[emptyJob]},root);
assert.equal(empty.results[0].result.localRecordCount,0);
assert.equal(empty.results[0].result.endpoints.exactAnyReported,'not-observed-in-input');
const sam=await readFile(join(root,'junction-fixture.sam'),'utf8');
assert.throws(()=>parseSam('bad.sam',sam.replace('NH:i:1','NH:Z:1')),/Invalid type/);
assert.throws(()=>parseSam('bad.sam',sam.replace('NH:i:1','NH:i:1\tNH:i:1')),/duplicate NH/);
const checks={version:report.engine_version,comparisons:Object.keys(expected).length,
  samBamParity:true,mateAndAlternativeFalseRescuesRejected:true,unknownMapqAndNhPreserved:true,
  readGroupAmbiguityRejected:true,coordinateErrorsRejected:true,emptyBamLocusHandled:true,
  scope:'Controlled alignment software regression, not biological fragment replay'};
await writeFile(join(root,'PROJECT_V05_JUNCTION_CHECKS.json'),JSON.stringify(checks,null,2)+'\n');
await writeFile(join(root,'PROJECT_V05_JUNCTION_RESULT.json'),JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify(checks));
