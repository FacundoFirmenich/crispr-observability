import {validateAlignment} from './alignment.js';

const alignedOps=new Set(['M','=','X']);
const key=span=>span.join(':');

function readEnd(flag){
  if(!(flag&1))return 'unpaired';
  if(Boolean(flag&64)===Boolean(flag&128))return 'unspecified';
  return flag&64?'first':'last';
}

function adjacentQuality(read,block,query){
  if(!block||!alignedOps.has(block.op)||read.quality==='*')return null;
  return read.quality.charCodeAt(query)-33;
}

function observedJunctions(read){
  return read.ops.flatMap((block,i)=>{
    if(block.op!=='N')return [];
    const [q]=block.query;
    return [{span:block.reference,storedQueryBoundary:q,
      flankingPhred:[
        adjacentQuality(read,read.ops[i-1],q-1),
        adjacentQuality(read,read.ops[i+1],q)
      ]}];
  });
}

export function junctionEvidence(input,job){
  const {locus,read_name,read_group,mate,junctions}=job;
  const [start,end]=locus.span;
  const candidates=input.records.map(validateAlignment).filter(r=>
    r.ref===locus.ref&&r.start<end&&r.referenceEnd>start&&r.name===read_name&&
    (read_group===undefined||r.evidence?.rg===read_group));
  const groups=new Set(candidates.map(r=>r.evidence?.rg??null));
  if(groups.size>1)throw new Error('Read name spans multiple read groups. Specify read_group.');
  const required=junctions.map(key);
  const records=candidates.map(r=>{
    const skips=observedJunctions(r),present=new Set(skips.map(x=>key(x.span)));
    const exact=required.every(k=>present.has(k));
    const segment=readEnd(r.flag);
    return {
      name:r.name,mate:segment,flag:r.flag,primary:!(r.flag&0x900),
      ref:r.ref,start:r.start,end:r.referenceEnd,cigar:r.cigar,ordinal:r.ordinal,
      mapqRaw:r.mapq,mappingQuality:r.mapq===255?null:r.mapq,
      alignmentTagsAndMate:r.evidence??null,
      skips,matchedJunctions:junctions.filter(span=>present.has(key(span))),
      allRequiredInRecord:exact,eligibleCarrier:mate===segment
    };
  });
  const witnesses=records.filter(r=>r.allRequiredInRecord&&r.eligibleCarrier);
  const primary=witnesses.filter(r=>r.primary);
  const competingPrimary=records.filter(r=>r.eligibleCarrier&&r.primary).length>1;
  const unique=primary.some(r=>r.alignmentTagsAndMate?.nh===1);
  const unresolved=primary.some(r=>!(r.alignmentTagsAndMate?.nh>=1));
  return {
    scope:input.source?.scope??'supplied alignment records',
    reference:locus.ref,window:locus.span,requestedJunctions:junctions,
    readName:read_name,readGroup:read_group??[...groups][0]??null,mate,
    localRecordCount:records.length,
    endpoints:{
      exactAnyReported:witnesses.length>0?'observed':'not-observed-in-input',
      exactUniquePrimary:primary.length&&competingPrimary?'indeterminate-primary':unique?'observed':unresolved?'indeterminate-nh':'not-observed-in-input'
    },
    records
  };
}
