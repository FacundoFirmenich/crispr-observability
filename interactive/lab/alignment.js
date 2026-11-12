// Coordinate adjudication for an existing SAM/BAM alignment. SAM input SEQ is
// in reference traversal order for FLAG 0x10. Every interval below is half-open.
const consumesRead=new Set(['M','I','S','=','X']);
const consumesReference=new Set(['M','D','N','=','X']);
const complement={A:'T',C:'G',G:'C',T:'A',U:'A',M:'K',R:'Y',W:'W',S:'S',Y:'R',K:'M',
  V:'B',H:'D',D:'H',B:'V',N:'N','=':'='};
export const MAX_READS=2500;
export const MAX_CASE_BYTES=6_000_000;
export const MAX_REPLAY_BYTES=24_000_000;
export const MAX_LOCUS_WIDTH=50_000;
const int=(x,min,max)=>Number.isInteger(x)&&x>=min&&x<=max;

export function parseCigar(cigar,referenceStart,sequenceLength){
  if(!int(referenceStart,0,2**31-1)||typeof cigar!=='string'||cigar.length>50_000)
    throw new Error('Invalid alignment start or CIGAR.');
  const ops=[];
  let reference=referenceStart,query=0,used=0;
  for(const match of cigar.matchAll(/([1-9]\d*)([MIDNSHP=X])/g)){
    const n=Number(match[1]),op=match[2];
    if(match.index!==used||!int(n,1,2**28-1)||ops.length>=20000)throw new Error('Malformed or oversized CIGAR.');
    const readEnd=query+(consumesRead.has(op)?n:0);
    const refEnd=reference+(consumesReference.has(op)?n:0);
    if(readEnd>2**31-1||refEnd>2**31-1)throw new Error('CIGAR overflows the coordinate range.');
    ops.push({op,length:n,query:[query,readEnd],reference:[reference,refEnd]});
    query=readEnd;reference=refEnd;used+=match[0].length;
  }
  if(!ops.length||used!==cigar.length)throw new Error('Malformed or empty CIGAR.');
  if(sequenceLength!==null&&query!==sequenceLength)throw new Error('CIGAR query length and SEQ disagree.');
  const first=ops.findIndex(x=>x.op!=='H'&&x.op!=='S');
  const last=ops.findLastIndex(x=>x.op!=='H'&&x.op!=='S');
  for(const [i,block] of ops.entries()){
    if(block.op==='H'&&i!==0&&i!==ops.length-1)throw new Error('Hard clipping must be terminal.');
    if(block.op==='S'&&i>first&&i<last)throw new Error('Soft clipping must be terminal, inside any hard clipping.');
  }
  const hardStart=ops[0].op==='H'?ops[0].length:0,hardEnd=ops.at(-1).op==='H'?ops.at(-1).length:0;
  return {ops,queryLength:query,referenceEnd:reference,hardStart,hardEnd};
}

export function validateAlignment(r){
  if(typeof r?.name!=='string'||!r.name||r.name.length>254||
    !int(r.flag,0,65535)||!int(r.mapq,0,255)||typeof r.ref!=='string'||r.ref==='*'||!r.ref||
    !int(r.start,0,2**31-1)||typeof r.cigar!=='string'||r.cigar==='*'||(r.flag&4))
    throw new Error('Expected a mapped alignment with valid identity, FLAG, MAPQ, reference, POS and CIGAR.');
  const sequence=r.sequence==='*'?null:r.sequence;
  if(sequence!==null&&(typeof sequence!=='string'||!(/^[=ACMGRSVTWYHKDBNU]+$/).test(sequence)))
    throw new Error('Unsupported SEQ alphabet.');
  const quality=r.quality==='*'?null:r.quality;
  if(quality!==null&&(sequence===null||typeof quality!=='string'||quality.length!==sequence.length||!/^[!-~]*$/.test(quality)))
    throw new Error('Malformed QUAL or QUAL length does not match SEQ.');
  const derived=parseCigar(r.cigar,r.start,sequence?.length??null);
  return {name:r.name,flag:r.flag,mapq:r.mapq,ref:r.ref,start:r.start,cigar:r.cigar,
    sequence:sequence??'*',quality:quality??'*',ordinal:int(r.ordinal,0,2**31-1)?r.ordinal:null,
    reverse:Boolean(r.flag&16),...derived};
}

function rawRecord(r) {
  const {name,flag,mapq,ref,start,cigar,sequence,quality,ordinal}=r;
  return {name,flag,mapq,ref,start,cigar,sequence,quality,ordinal};
}

export function parseSam(name,text,source={}){
  const records=[],refLengths={};
  let inputRecords=0,unmapped=0,lineNo=0;
  for(const line of text.split(/\r?\n/)){
    lineNo++;
    if(!line)continue;
    if(/^@(HD|SQ|RG|PG|CO)(\t|$)/.test(line)){
      if(line.startsWith('@SQ\t')){
        const parts=line.split('\t');
        const id=parts.find(x=>x.startsWith('SN:'))?.slice(3),length=Number(parts.find(x=>x.startsWith('LN:'))?.slice(3));
        if(id&&int(length,1,2**31-1))Object.defineProperty(refLengths,id,{value:length,enumerable:true,configurable:true});
      }
      continue;
    }
    inputRecords++;
    const f=line.split('\t');
    if(f.length<11||!/^\d+$/.test(f[1])||!/^\d+$/.test(f[3])||!/^\d+$/.test(f[4])||
      !int(Number(f[1]),0,65535)||!int(Number(f[4]),0,255))
      throw new Error(`SAM line ${lineNo} has invalid mandatory fields.`);
    if(Number(f[1])&4){unmapped++;continue;}
    try{
      if(records.length>=MAX_READS)throw new Error('Too many mapped records. Load a smaller locus slice.');
      const r=validateAlignment({name:f[0],flag:Number(f[1]),ref:f[2],start:Number(f[3])-1,mapq:Number(f[4]),
        cigar:f[5],sequence:f[9].toUpperCase(),quality:f[10],ordinal:inputRecords-1});
      if(Object.hasOwn(refLengths,r.ref)&&r.referenceEnd>refLengths[r.ref])throw new Error('Alignment extends past @SQ LN.');
      records.push(rawRecord(r));
    }catch(error){throw new Error(`SAM line ${lineNo}: ${error.message}`);}
  }
  return validateAlignmentCase({schema:'observability-alignment/1',kind:'alignment',format:'SAM',name,records,
    source:{...source,scope:'complete input slice',inputRecords,unmappedExcluded:unmapped,refLengths}});
}

export function validateAlignmentCase(input){
  if(input?.schema!=='observability-alignment/1'||input.kind!=='alignment'||
    typeof input.name!=='string'||!input.name.trim()||input.name.length>240||
    !['SAM','BAM'].includes(input.format)||!Array.isArray(input.records)||!input.records.length||
    input.records.length>MAX_READS)
    throw new Error('Expected a supported alignment case with 1–2,500 mapped records.');
  const records=input.records.map(r=>rawRecord(validateAlignment(r)));
  const region=input.source?.region;
  if(region&&(!int(region.start,0,2**31-1)||!int(region.end,region.start+1,2**31-1)||
      region.end-region.start>MAX_LOCUS_WIDTH||typeof region.ref!=='string'||
      records.some(r=>r.ref!==region.ref)))
    throw new Error('Invalid or inconsistent indexed locus.');
  return {schema:input.schema,kind:input.kind,format:input.format,name:input.name,records,
    source:input.source&&typeof input.source==='object'&&!Array.isArray(input.source)?input.source:{}};
}

export function inspectQuery(r,index,original=false){
  if(!int(index,0,r.queryLength-1))throw new Error('Read index is outside the stored query.');
  const query=original&&r.reverse?r.queryLength-1-index:index;
  const op=r.ops.find(block=>query>=block.query[0]&&query<block.query[1]);
  const letter=r.sequence==='*'?null:r.sequence[query];
  const aligned=consumesReference.has(op.op);
  const beforeHard=r.reverse?r.hardEnd:r.hardStart;
  const originalQuery=r.reverse?r.queryLength-1-query:query;
  return {query,originalQuery,unclippedOriginalIndex:beforeHard+originalQuery,op:op.op,
    reference:aligned?op.reference[0]+query-op.query[0]:null,
    insertionAnchor:op.op==='I'?op.reference[0]:null,
    base:original&&r.reverse&&letter?complement[letter]:letter,
    quality:r.quality==='*'?null:r.quality.charCodeAt(query)-33};
}

export function coverage(records,ref,lo,hi,{minMapq=0,includeNonPrimary=false}={}){
  if(!int(lo,0,2**31-1)||!int(hi,lo+1,2**31-1)||hi-lo>MAX_LOCUS_WIDTH||!int(minMapq,0,254))
    throw new Error('Choose a non-empty reference window of at most 50,000 bases and MAPQ from 0 to 254.');
  const delta=new Int32Array(hi-lo+1);
  let considered=0,used=0,unknownMapq=0,nonPrimary=0,lowMapq=0;
  for(const row of records){
    if(row.ref!==ref)continue;
    const r=row.ops?row:validateAlignment(row);
    if(r.start>=hi||r.referenceEnd<=lo)continue;
    considered++;
    if(!includeNonPrimary&&(r.flag&(256|2048))){nonPrimary++;continue;}
    if(r.mapq===255){unknownMapq++;continue;}
    if(r.mapq<minMapq){lowMapq++;continue;}
    used++;
    for(const b of r.ops){
      if(!consumesRead.has(b.op)||!consumesReference.has(b.op))continue;
      const start=Math.max(lo,b.reference[0]),end=Math.min(hi,b.reference[1]);
      if(start<end){delta[start-lo]++;delta[end-lo]--;}
    }
  }
  let depth=0,covered=0,maxDepth=0,sumDepth=0;
  const runs=[];
  for(let i=0;i<hi-lo;i++){
    depth+=delta[i];
    if(depth)covered++;
    maxDepth=Math.max(depth,maxDepth);sumDepth+=depth;
    const prev=runs.at(-1);
    if(prev?.depth===depth)prev.end=lo+i+1;
    else runs.push({start:lo+i,end:lo+i+1,depth});
  }
  return {ref,window:[lo,hi],considered,used,excluded:{unknownMapq,nonPrimary,lowMapq},
    covered,maxDepth,sumDepth,runs,filters:{minMapq,includeNonPrimary}};
}

export function reportAlignment(input,settings){
  const data=validateAlignmentCase(input);
  if(!int(settings?.index,0,data.records.length-1))throw new Error('Select an alignment in the imported case.');
  const read=validateAlignment(data.records[settings.index]);
  if(settings.minMapq!==undefined&&!int(settings.minMapq,0,254))throw new Error('MAPQ threshold must be 0–254.');
  const region=data.source.region;
  const lo=settings.lo??Math.max(read.start,region?.start??read.start);
  const hi=settings.hi??Math.min(region?.end??Infinity,lo+1000,Math.max(lo+1,read.referenceEnd));
  if(region&&(lo<region.start||hi>region.end))
    throw new Error('Coverage cannot extend beyond the indexed source locus.');
  const knownLength=Object.hasOwn(data.source.refLengths??{},read.ref)?data.source.refLengths[read.ref]:null;
  if(knownLength!==null&&hi>knownLength)throw new Error('Coverage exceeds the declared reference length.');
  return {read,coverage:coverage(data.records,read.ref,lo,hi,settings)};
}

export function alignmentTsv(read,depth){
  const header=['kind','ref','start_0','end_exclusive','stored_query_start_0','stored_query_end_exclusive','operation','depth'];
  const ref=read.ref.replace(/[\t\r\n]/g,' ').replace(/^[=+\-@]/,"'$&");
  return [header,...read.ops.map(b=>['cigar',ref,...b.reference,...b.query,b.length+b.op,'']),
    ...depth.runs.map(b=>['coverage',ref,b.start,b.end,'','','',b.depth])].map(row=>row.join('\t')).join('\n')+'\n';
}
