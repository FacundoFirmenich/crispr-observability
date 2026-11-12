import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {parseSam,reportAlignment,inspectQuery,parseCigar,coverage} from './alignment.js';
import {openIndexedBam} from './bam-input.js';

const {checks,source}=JSON.parse(await readFile('ALIGNMENT_FIXTURE_ORACLE.json','utf8'));
const proof=[];
for(const fixture of checks){
  const local=async name=>new File([await readFile(name)],name);
  const reader=await openIndexedBam(await local(fixture.bam),await local(fixture.bam+'.bai'));
  const indexed=await reader.query(fixture.ref,fixture.lo,fixture.hi);
  assert.equal(indexed.source.bamSHA256,fixture.bamSHA256);
  const sam=parseSam(fixture.sam,await readFile(fixture.sam,'utf8'));
  const normalize=rows=>rows.map(({ordinal,...r})=>r).sort((a,b)=>JSON.stringify(a).localeCompare(JSON.stringify(b)));
  assert.deepEqual(normalize(indexed.records),normalize(fixture.records));
  assert.deepEqual(normalize(sam.records),normalize(fixture.records));
  const inputs={index:0,lo:fixture.lo,hi:fixture.hi,minMapq:fixture.threshold};
  const i=reportAlignment(indexed,inputs),s=reportAlignment(sam,inputs);
  assert.deepEqual(i.coverage,s.coverage);
  const depth=i.coverage.runs.flatMap(x=>Array(x.end-x.start).fill(x.depth));
  assert.deepEqual(depth,fixture.depth);
  assert.throws(()=>reportAlignment(indexed,{...inputs,hi:fixture.hi+1}),/beyond the indexed/);
  proof.push({bam:fixture.bam,ref:fixture.ref,window:[fixture.lo,fixture.hi],records:indexed.records.length,
    minimumMapq:fixture.threshold,covered:i.coverage.covered,basesCompared:depth.length,peak:i.coverage.maxDepth,
    sourceSHA256:indexed.source.bamSHA256,sliceSHA256:indexed.source.sliceSHA256});
}

const fixture=checks[0],sam=parseSam(fixture.sam,await readFile(fixture.sam,'utf8'));
const {read,coverage:unfiltered}=reportAlignment(sam,{index:0,lo:fixture.lo,hi:fixture.hi,minMapq:0});
assert.equal(unfiltered.covered,67);
assert.equal(read.reverse,true);
assert.equal(inspectQuery(read,0).op,'S');
assert.equal(inspectQuery(read,0).reference,null);
assert.equal(inspectQuery(read,183).reference,fixture.lo);
assert.equal(inspectQuery(read,0,true).reference,fixture.hi-1);
assert.equal(inspectQuery(read,249,true).reference,null);
assert.throws(()=>parseCigar('3M2S4M',100,9),/terminal/);
assert.throws(()=>parseCigar('2M',100,1),/SEQ disagree/);
const mixed=parseCigar('4H3S2=1I2X1D2N3M2S2H',10,13);
assert.deepEqual(mixed.ops.map(r=>[r.op,r.query,r.reference]),[
  ['H',[0,0],[10,10]],['S',[0,3],[10,10]],['=',[3,5],[10,12]],['I',[5,6],[12,12]],
  ['X',[6,8],[12,14]],['D',[8,8],[14,15]],['N',[8,8],[15,17]],['M',[8,11],[17,20]],
  ['S',[11,13],[20,20]],['H',[13,13],[20,20]]
]);
const base={name:'coordinate-contract',flag:0,mapq:40,ref:'r',start:10,
  cigar:'4H3S2=1I2X1D2N3M2S2H',sequence:'ACGTACGTACGTA',quality:'*'};
const cov=coverage([base,{...base,flag:256},{...base,mapq:255}],'r',10,20);
assert.equal(cov.used,1);assert.equal(cov.sumDepth,7);assert.equal(cov.covered,7);
assert.deepEqual(cov.excluded,{unknownMapq:1,nonPrimary:1,lowMapq:0});
console.log(JSON.stringify({source,independentOracle:'Python stdlib BGZF sequential binary decode and per-base accumulation',
  exactCrossFormat:true,reverseClip:true,mapqSweep:{minimum0:67,minimum20:0},
  handCheckedCigarContract:true,checks:proof},null,2));
