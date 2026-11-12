import {alignmentTsv,coverage,inspectQuery,validateAlignment,validateAlignmentCase} from './alignment.js';
import {scaffold,ink} from './molecular-view.js';
import {span} from './science.js';

const descriptions={M:'Aligned · identity unspecified',I:'Insertion',D:'Deletion',N:'Reference skip',S:'Soft clip',
  H:'Hard clip',P:'Padding','=':'Sequence match',X:'Sequence mismatch'};
const color={I:ink.rose,S:ink.violet,X:ink.rose,M:ink.blue,'=':ink.blue};
const stringify=value=>value===null?'unavailable':String(value);
const qTransform=(r,interval,original)=>original&&r.reverse?[r.queryLength-interval[1],r.queryLength-interval[0]]:interval;
const qRange=(r,block,original)=>qTransform(r,block.query,original);

function node(tag,value){
  const el=document.createElement(tag);
  if(value!==undefined)el.textContent=value;
  return el;
}

export function describeQuery(position,original){
  const reference=position.reference!==null?`reference index ${position.reference}`:
    position.insertionAnchor!==null?`insertion at boundary ${position.insertionAnchor}`:'no mapped reference base';
  const index=original?position.originalQuery:position.query;
  return `Query ${index} · ${descriptions[position.op]} · ${reference} · base ${stringify(position.base)} · Phred ${stringify(position.quality)}`;
}

export function makeAlignmentView($,focus){
  let bundle,reads=[],active,settings,depth,focusInfo,currentPage=0;
  function load(value){
    const next=validateAlignmentCase(value);
    const parsed=next.records.map(validateAlignment);
    if(!parsed.some(r=>r.queryLength>0))throw new Error('The selected slice has no records with query-consuming CIGAR operations.');
    bundle=next;reads=parsed;
    $('alignRecord').replaceChildren(...reads.map((r,i)=>
      new Option(`${i+1} · ${r.name} · ${r.ref}:${r.start} · ${r.reverse?'−':'+'}`,String(i))));
    for(const [i,read] of reads.entries())$('alignRecord').options[i].disabled=read.queryLength===0;
    $('alignRecord').value=String(reads.findIndex(r=>r.queryLength>0));
    $('minMapq').value=0;$('nonPrimary').checked=false;$('originalRead').checked=false;
    $('alignment-controls').hidden=false;
    reset();
  }
  function reset(){
    const r=reads[+$('alignRecord').value];
    if(!r)return;
    const bound=bundle.source.region;
    const lo=Math.max(r.start,bound?.start??r.start);
    $('covStart').value=lo;
    $('covEnd').value=Math.min(bound?.end??Infinity,lo+1000,Math.max(lo+1,r.referenceEnd));
    $('queryFocus').value=0;
    currentPage=0;
  }
  function build(){
    if(!bundle)throw new Error('Load a SAM slice or select a BAM and its index.');
    const numeric=id=>$(id).value===''?NaN:Number($(id).value);
    const index=numeric('alignRecord');
    active=reads[index];
    if(!active||active.queryLength<=0)throw new Error('Choose an alignment with a nonempty query.');
    const original=$('originalRead').checked,at=numeric('queryFocus'),lo=numeric('covStart'),hi=numeric('covEnd');
    const options={minMapq:numeric('minMapq'),includeNonPrimary:$('nonPrimary').checked};
    const bound=bundle.source.region;
    if(bound&&(active.ref!==bound.ref||lo<bound.start||hi>bound.end))
      throw new Error(`Coverage for this BAM slice is bounded by the loaded locus ${bound.ref} ${span([bound.start,bound.end])}.`);
    const knownLength=Object.hasOwn(bundle.source.refLengths??{},active.ref)?bundle.source.refLengths[active.ref]:null;
    if(knownLength!==null&&hi>knownLength)throw new Error('The coverage window extends past the reference length.');
    $('queryFocus').max=active.queryLength-1;
    settings={index,lo,hi,original,at,...options};
    depth=coverage(reads,active.ref,lo,hi,options);
    focusInfo=inspectQuery(active,at,original);
    const hiRead=Math.min(active.queryLength,Math.max(100,at+51)),loRead=Math.max(0,hiRead-100);
    const blocks=active.ops.map(b=>[qRange(active,b,original),b.op]);
    const colors=Array(hiRead-loRead).fill(ink.muted);
    for(const [[a,b],op] of blocks)
      for(let i=Math.max(loRead,a);i<Math.min(hiRead,b);i++)colors[i-loRead]=color[op]||ink.muted;
    const geom=scaffold(loRead,hiRead,i=>colors[i-loRead]);
    const soft=active.ops.reduce((n,b)=>n+(b.op==='S'?b.length:0),0);
    const result={read:{name:active.name,ordinal:active.ordinal,flag:active.flag,reverse:active.reverse,
      mapq:active.mapq,ref:active.ref,cigar:active.cigar,start:active.start,end:active.referenceEnd,
      queryLength:active.queryLength,sequenceAvailable:active.sequence!=='*',softClipped:soft,
      hardClipped:active.hardStart+active.hardEnd},position:focusInfo,coverage:depth,operations:active.ops};
    return {geometry:{...geom,cut:null,track:{lo:loRead,hi:hiRead,title:original?'Original query indices':'Stored SEQ query indices',
      rows:[['Aligned',['M','=','X'],'blue'],['Soft clip',['S'],'violet'],['Insertion',['I'],'rose']]
        .map(([label,ops,shade])=>[label,blocks.filter(([,op])=>ops.includes(op)).map(([q])=>q),shade])}},settings,result,
      metrics:[['Reference window',span([lo,hi]),active.ref],
        ['Records counted',`${depth.used} / ${depth.considered}`,'overlapping records under current filters'],
        ['Observed reference bases',`${depth.covered} / ${hi-lo}`,`peak depth ${depth.maxDepth}`]]};
  }
  function paint(){
    if(!depth||$('alignment-evidence').hidden)return;
    const canvas=$('coverage-plot'),w=Math.max(1,canvas.clientWidth);
    const h=78,dpr=Math.min(devicePixelRatio||1,2);
    canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);
    const ctx=canvas.getContext('2d');
    ctx.setTransform(dpr,0,0,dpr,0,0);ctx.fillStyle='#192542';ctx.fillRect(0,0,w,h);
    const [lo,hi]=depth.window, count=Math.min(hi-lo,240);
    const peaks=new Int32Array(count);
    for(const r of depth.runs){
      const a=Math.floor((r.start-lo)*count/(hi-lo));
      const b=Math.ceil((r.end-lo)*count/(hi-lo));
      for(let i=a;i<b;i++)peaks[i]=Math.max(peaks[i],r.depth);
    }
    ctx.fillStyle='#8bcaff';
    for(let i=0;i<count;i++){
      const bar=peaks[i]*60/Math.max(1,depth.maxDepth);
      ctx.fillRect(i*w/count,72-bar,Math.max(1,w/count-1),bar);
    }
    if(focusInfo?.reference!==null&&focusInfo.reference>=lo&&focusInfo.reference<hi){
      ctx.fillStyle='#f5a8bc';ctx.fillRect((focusInfo.reference-lo)*w/(hi-lo),4,2,68);
    }
  }
  function show(){
    $('coverage-label').textContent=`${depth.ref} ${span(depth.window)} · local depth maxima per display bin. Exact depth runs are in TSV and JSON. Excluded: ${depth.excluded.unknownMapq} unknown MAPQ, ${depth.excluded.nonPrimary} non-primary, ${depth.excluded.lowMapq} below threshold.`;
    $('align-query-note').textContent=describeQuery(focusInfo,settings.original);
    const tbl=node('table'),thead=node('thead'),hr=node('tr'),body=node('tbody');
    for(const label of ['CIGAR','Operation','Query interval','Reference interval','Inspect'])hr.append(node('th',label));
    thead.append(hr);tbl.append(thead,body);
    const rows=active.ops.slice(currentPage*120,(currentPage+1)*120);
    for(const b of rows){
      const q=qRange(active,b,settings.original),tr=node('tr');
      for(const cell of [b.length+b.op,descriptions[b.op],span(q),span(b.reference)])tr.append(node('td',cell));
      const last=node('td');
      if(q[0]<q[1]){
        const button=node('button','Select');button.addEventListener('click',()=>focus(q[0]));last.append(button);
      }else last.textContent='No stored base';
      tr.append(last);body.append(tr);
    }
    const nav=node('div');
    if(active.ops.length>120){
      for(const [label,delta,disabled] of [['Previous',-1,currentPage===0],
        ['Next',1,(currentPage+1)*120>=active.ops.length]]){
        const b=node('button',label);b.disabled=disabled;b.addEventListener('click',()=>{currentPage+=delta;show();});nav.append(b);
      }
      nav.append(node('span',` Blocks ${currentPage*120+1}–${currentPage*120+rows.length} of ${active.ops.length}`));
    }
    $('align-table').replaceChildren(tbl,nav);
    paint();
  }
  function inspect(index){
    focusInfo=inspectQuery(active,index,settings.original);settings.at=index;
    $('align-query-note').textContent=describeQuery(focusInfo,settings.original);
    paint();
    return focusInfo;
  }
  window.addEventListener('resize',paint);
  return {load,reset,build,show,inspect,tsv:()=>alignmentTsv(active,depth),get bundle(){return bundle;}};
}
