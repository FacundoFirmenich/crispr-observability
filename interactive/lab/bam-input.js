import {BamFile} from '@gmod/bam';
import {BlobFile} from 'generic-filehandle2';
import {MAX_CASE_BYTES,MAX_LOCUS_WIDTH,MAX_READS,validateAlignmentCase} from './alignment.js';

const maxRangeBytes=24*1024*1024;
class BoundedBam extends BlobFile {
  budget=maxRangeBytes;
  async read(length,position,options){
    if(!Number.isSafeInteger(length)||length<0||length>this.budget)
      throw new Error('This locus exceeds the 24 MB indexed-read budget. Narrow the query.');
    this.budget-=length;
    return super.read(length,position,options);
  }
  async readFile(){
    throw new Error('Indexed BAM loading does not read the complete file.');
  }
}

export async function sha256(bytes){
  const digest=await globalThis.crypto.subtle.digest('SHA-256',bytes);
  return Array.from(new Uint8Array(digest),n=>n.toString(16).padStart(2,'0')).join('');
}

export async function openIndexedBam(bamBlob,indexBlob,signal){
  if(!bamBlob||!indexBlob||indexBlob.size>32*1024*1024)
    throw new Error('Choose a BAM and its BAI or CSI index (index up to 32 MB).');
  const isCsi=indexBlob.name.toLowerCase().endsWith('.csi');
  if(!isCsi&&!indexBlob.name.toLowerCase().endsWith('.bai'))throw new Error('Choose a BAI or CSI index.');
  const bamHandle=new BoundedBam(bamBlob);
  const indexHandle=new BlobFile(indexBlob);
  const bam=new BamFile({bamFilehandle:bamHandle,[isCsi?'csiFilehandle':'baiFilehandle']:indexHandle,
    maxCacheBytes:4*1024*1024});
  const header=await bam.getHeaderText({signal});
  const refLengths=Object.create(null);
  for(const line of header.split(/\r?\n/)){
    if(!line.startsWith('@SQ\t'))continue;
    const parts=line.split('\t');
    const id=parts.find(x=>x.startsWith('SN:'))?.slice(3),length=Number(parts.find(x=>x.startsWith('LN:'))?.slice(3));
    if(id&&Number.isSafeInteger(length)&&length>0)refLengths[id]=length;
  }
  if(!Object.keys(refLengths).length)throw new Error('BAM header has no reference dictionary.');
  let identity;
  return {
    refLengths,
    async query(ref,start,end,signal){
      if(!Object.hasOwn(refLengths,ref)||!Number.isInteger(start)||!Number.isInteger(end)||
        start<0||end<=start||end>refLengths[ref]||end-start>MAX_LOCUS_WIDTH)
        throw new Error('Choose a reference interval within the BAM header and at most 50,000 bases wide.');
      bamHandle.budget=maxRangeBytes;
      const records=await bam.getRecordsForRange(ref,start,end,{signal});
      const mapped=records.filter(r=>!r.isSegmentUnmapped());
      if(!mapped.length)throw new Error('No mapped records overlap that locus.');
      if(mapped.length>MAX_READS)throw new Error('The locus has more than 2,500 mapped records. Choose a narrower locus.');
      let estimatedBytes=0;
      const rows=mapped.map((r,ordinal)=>{
        const sequence=r.seq||'*',cigar=r.CIGAR||'*';
        estimatedBytes+=2*sequence.length+cigar.length+160;
        if(estimatedBytes>MAX_CASE_BYTES)throw new Error('This locus is too large to replay locally. Narrow the interval.');
        const q=r.qual;
        const quality=q?.length&&q[0]!==255?Array.from(q,n=>String.fromCharCode(n+33)).join(''):'*';
        return {name:r.name,flag:r.flags,mapq:r.mq??255,ref,start:r.start,cigar,sequence,quality,ordinal};
      });
      identity??=Promise.all([
        bamBlob.size<=MAX_CASE_BYTES?bamBlob.arrayBuffer().then(sha256):null,
        indexBlob.arrayBuffer().then(sha256)
      ]);
      const [bamSHA256,indexSHA256]=await identity;
      signal?.throwIfAborted();
      const sliceSHA256=await sha256(new TextEncoder().encode(JSON.stringify(rows)));
      return validateAlignmentCase({schema:'observability-alignment/1',kind:'alignment',format:'BAM',
        name:bamBlob.name,records:rows,source:{scope:'indexed locus',region:{ref,start,end},
          bamBytes:bamBlob.size,bamSHA256,indexName:indexBlob.name,indexSHA256,sliceSHA256,refLengths}});
    }
  };
}
