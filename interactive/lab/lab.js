import {MolecularView, scaffold, atomic, ink} from './molecular-view.js';
import {boundary,compare,validate,span,flip,revcomp,inSpan,complement} from './science.js';
import {readCase, importedGeometry} from './import-case.js';
import {setupGuide} from './guide.js';
import {parseSam,validateAlignmentCase,MAX_CASE_BYTES,MAX_REPLAY_BYTES} from './alignment.js';
import {openIndexedBam,sha256} from './bam-input.js';
import {makeAlignmentView,describeQuery} from './alignment-view.js';

const $=id=>document.getElementById(id);
const data=JSON.parse($('data').textContent);
const {audit,molecule}=data;
const cases=[...audit.cases];
const modes=['edge','compare','atom','import','align'];
const number=id=>Number($(id).value===''?NaN:$(id).value);
const readIds=['start','end','shift','readEnd'];
const pane=(tag,content,cls)=>{
  const node=document.createElement(tag);
  node.textContent=content;
  if(cls)node.className=cls;
  return node;
};
const swatches={blue:'blue',violet:'violet',rose:'rose'};
let mode='edge', view, currentResult, viewSpec, imported, picked, pair, requestedFit=true;
let renderQueued=false, pinned=null;
const guide=setupGuide($,selectPreset);
const alignment=makeAlignmentView($,focusAlignment);
let bamBlob,bam,loadEpoch=0;

function showError(error) {
  $('webgl-error').hidden=false;
  $('webgl-error').textContent='3D requires WebGL 2. The numerical experiments and data export still work.';
  console.warn(error);
}
try{
  view=new MolecularView($('scene'),pick);
  $('scene').addEventListener('webglcontextlost',e=>{e.preventDefault();showError('WebGL context lost');});
}catch(error){showError(error);}

function pick(row) {
  if(mode==='atom'){
    const base=molecule.seq[row.index];
    $('picked').textContent=`Read index ${114+row.index} · strand ${row.strand?'B':'A'} · ${row.strand?complement[base]:base} ${row.atom} · (${row.xyz.map(n=>n.toFixed(3)).join(', ')}) Å`;
  }else if(mode==='compare'){
    const rec=cases[+$('record').value], rev=$('reverse').checked;
    const i=rev?rec.length-1-row.index:row.index;
    const native=rec.segment[i-rec.domain[0]];
    const onShown=(rev!==Boolean(row.strand))?complement[native]??'N':native;
    $('picked').textContent=`Displayed index ${row.index} · source index ${i} · ${onShown}`;
  }else if(mode==='import'){
    const xyz=row.xyz.map(n=>n.toFixed(3)).join(', ');
    const i=viewSpec.points.indexOf(row)+1;
    const a=picked?.xyz;
    const distance=a?Math.hypot(...a.map((n,i)=>n-row.xyz[i])):null;
    $('picked').textContent=`#${i} · ${row.label} · (${xyz}) Å${distance===null?' · select another atom to measure':` · distance ${distance.toFixed(3)} Å`}`;
    $(a?'distanceB':'distanceA').value=i;
    pair=a?{a:{label:picked.label,xyz:a},b:{label:row.label,xyz:row.xyz},distance}:null;
    view?.setMeasurement(picked?.p,row.p);
    reportDistance();
    picked=a?null:row;
  }else if(mode==='align'){
    const position=alignment.inspect(row.index);
    currentResult.result.position=position;
    currentResult.inputs.at=row.index;
    $('queryFocus').value=row.index;
    $('picked').textContent=describeQuery(position,$('originalRead').checked);
  }else{
    $('picked').textContent=`Index ${row.index} · ${row.index<number('readEnd')?'inside the read':'outside the read'}`;
  }
  if(mode!=='import')$('track-position').value=mode==='atom'?114+row.index:row.index;
  view?.select(row.p,Math.max(row.r*1.35,.84));
}

function renderTrack(spec) {
  $('coord-track').hidden=!spec;
  if(!spec)return;
  const {lo,hi,rows,title}=spec;
  $('track-label').textContent=title;
  $('track-domain').textContent=span([lo,hi]);
  const slider=$('track-position');
  slider.min=lo;slider.max=hi-1;slider.value=lo;
  $('track-rows').replaceChildren(...rows.map(([label,interval,color])=>{
    const line=pane('div','','track-row'),bar=pane('div','','track-bar');
    const ranges=Array.isArray(interval[0])?interval:interval.length?[interval]:[];
    for(const [a,b] of ranges){
      const fill=pane('span','',color);
      const left=Math.max(lo,a),right=Math.min(hi,b);
      if(left>=right)continue;
      fill.style.left=100*(left-lo)/(hi-lo)+'%';
      fill.style.width=100*(right-left)/(hi-lo)+'%';bar.append(fill);
    }
    line.append(pane('span',label),bar);
    return line;
  }));
}

function focusAlignment(index) {
  const item=viewSpec?.points.find(point=>point.index===index);
  if(item)pick(item);
  else{$('queryFocus').value=index;redraw();}
}

function renderAlign(){
  if(!alignment.bundle){
    viewSpec={points:[],links:[],frame:15,cut:null};
    $('export').disabled=true;
    present({kind:'alignment_pending'},[['Input','SAM or BAM','Open an example or choose a local file']],[],
      'Open the public example to inspect a reverse-strand read. For BAM, supply its BAI or CSI and choose a locus.',
      'Input files stay on this device.',
      'Indexed alignment input');
    $('experiment').replaceChildren();
    return;
  }
  const next=alignment.build();
  viewSpec=next.geometry;
  const source=alignment.bundle.source;
  present({kind:'alignment_locus_audit',inputs:next.settings,result:next.result},
    next.metrics,[['blue','Aligned read bases'],['violet','Soft clips'],['rose','Insertions or X operations']],
    `Selected ${next.result.read.name}, CIGAR ${next.result.read.cigar}. The reference depth uses ${next.result.coverage.used} records after the displayed filters.`,
    `Evidence universe: ${source.scope??'portable slice'}, ${alignment.bundle.records.length} mapped records in this case. The 3D scaffold is in read coordinates. M records alignment without declaring a match or mismatch.`,
    'Read scaffold · linked to existing alignment');
  $('experiment').replaceChildren(pane('h3','Record identity'),
    pane('p',`FLAG ${next.result.read.flag} · ${next.result.read.reverse?'reverse':'forward'} · MAPQ ${next.result.read.mapq===255?'unavailable':next.result.read.mapq}`),
    pane('p',`${next.result.read.softClipped} stored soft-clipped bases · ${next.result.read.hardClipped} hard-clipped positions absent from SEQ`));
  $('alignment-evidence').hidden=false;
  alignment.show();
}

function reportDistance(){
  if(mode!=='import'||!currentResult)return;
  const threshold=number('distanceCap');
  const ok=$('distanceCap').validity.valid&&Number.isFinite(threshold)&&threshold>=0;
  $('export').disabled=!ok;
  $('input-error').textContent=ok?'':'Use a non-negative distance threshold.';
  const el=$('distance-result');
  currentResult.inputs.threshold=ok?threshold:null;
  currentResult.inputs.atomIndices=pair?['distanceA','distanceB'].map(number):null;
  currentResult.result.measurement=pair&&ok?{...pair,withinThreshold:pair.distance<=threshold}:null;
  el.textContent=pair&&ok?`${pair.a.label} ↔ ${pair.b.label} · ${pair.distance.toFixed(3)} Å · ${pair.distance<=threshold?'within':'beyond'} ${threshold} Å.`:'Select two atoms to evaluate the threshold.';
  guide.update(currentResult,mode);
}

function present(result,metrics,legend,explanation,scope,title) {
  currentResult=result;
  $('scene-title').textContent=title;
  $('metrics').replaceChildren(...metrics.map(([label,value,description])=>{
    const el=pane('div','','card');
    el.append(pane('small',label),pane('strong',value),pane('small',description));
    return el;
  }));
  $('legend').replaceChildren(...legend.map(([color,label])=>pane('span',label,`swatch ${swatches[color]}`)));
  $('explanation').textContent=explanation;
  $('scope').textContent=scope;
  guide.update(result,mode);
}

function button(label,action) {
  const node=pane('button',label);
  node.type='button';
  node.addEventListener('click',action);
  return node;
}

function setBoundary(values) {
  readIds.forEach((id,i)=>$(id).value=values[i]);
  syncReadSlider();
  redraw(true);
}

function syncReadSlider(){
  const L=number('readEnd');
  if(!Number.isInteger(L)||!$('readEnd').validity.valid)return;
  const slider=$('readSlide');
  if(L<+slider.min||L>+slider.max){
    slider.min=Math.max(0,L-55);
    slider.max=Math.min(2000000,L+55);
  }
  slider.value=L;
}

function renderEdge() {
  const [start,end,shift,L]=readIds.map(number);
  for(const id of readIds)if(!$(id).validity.valid)throw new RangeError('Use an integer within the control’s stated range.');
  const r=boundary(start,end,shift,L);
  const baseline=boundary(start,end,0,L);
  const [a,b]=r.anchor;
  const lo=Math.max(0,Math.min(b-26,Math.max(b-91,Math.min(a-4,L-5))));
  const hi=Math.max(lo+36,Math.min(b+91,Math.max(b+7,L+5)));
  const boundedHi=Math.min(hi,lo+120);
  const geometry=scaffold(lo,boundedHi,i=>inSpan(i,r.anchor)?(i<L?ink.blue:ink.rose):ink.muted);
  const cutInView=L>=lo&&L<=boundedHi;
  viewSpec={...geometry,cut:$('cutVisible').checked&&cutInView?geometry.x(L):null,track:{
    lo,hi:boundedHi,title:'Local read coordinates',rows:[
      ['Read coverage',[0,L],'blue'],['Required anchor',r.anchor,'rose']]}};
  const lost=r.outside?(r.observed?` The first missing anchor base is at index ${Math.max(L,a)}.`:' All anchor bases fall outside the read.'):' The full interval is inside the read.';
  present({kind:'coordinate_containment',inputs:{start,end,shift,L},result:r},
    [['Shifted anchor',span(r.anchor),`from ${span(r.unshifted)}`],
     ['Anchor bases observed',`${r.observed} / ${r.width}`,r.contained?'Complete containment':'Incomplete containment'],
     ['First read end that fits',String(b),`current read end ${L}`]],
    [['blue','Observed anchor'],['rose','Missing anchor / read boundary']],
    `The unshifted anchor ${baseline.contained?'fits':'does not fit'}. After a displacement of ${shift} nt, ${r.observed} anchor bases intersect ${span(r.read)}.${lost}`,
    `The 3D scaffold shows local positions ${span([lo,boundedHi])}.${cutInView?'':' The read end falls outside this local view.'} Containment is one admission condition; motif identity, quality and editing outcomes require separate evidence.`,
    'Finite-read coordinate model');
  const el=$('experiment');
  el.replaceChildren(pane('h3','Sweep read end around the threshold'));
  const grid=pane('div','','sweep');
  const sweep=[];
  for(let n=Math.max(0,b-4);n<=b+4;n++){
    const fit=boundary(start,end,shift,n).contained;
    sweep.push({readEnd:n,contained:fit});
    const node=button(String(n),()=>{setBoundary([start,end,shift,n]);});
    node.dataset.fit=String(fit);
    node.setAttribute('aria-label',`Read end ${n}: ${fit?'fits':'does not fit'}`);
    node.setAttribute('aria-pressed',String(L===n));
    node.disabled=n<+$(readIds[3]).min||n>+$(readIds[3]).max;
    grid.append(node);
  }
  currentResult.sweep=sweep;
  el.append(grid,pane('small','Solid underline · fits. Dashed underline · incomplete.'));
  const line=pane('p',pinned?`Pinned read end ${pinned.inputs.L}: ${pinned.result.observed}/${pinned.result.width}. Current: ${r.observed}/${r.width}.`:'Pin a result before changing a parameter.');
  el.append(line,button(pinned?'Replace pinned result':'Pin current result',()=>{pinned=currentResult;redraw();}));
  const contrasts=pane('details','');
  contrasts.append(pane('summary','Frozen CFTR contrast intervals'));
  for(const [i,label] of ['Edge anchor','Inner anchor','Target motifs'].entries()){
    const row=audit.contrasts[i];
    contrasts.append(pane('p',`${label} · [${Number(row.L_pp).toFixed(2)}, ${Number(row.H_pp).toFixed(2)}] pp`));
  }
  contrasts.append(pane('small','Archive results under a common Q30 mask. Slider changes do not recompute treated-minus-mock bounds.'));
  el.append(contrasts);
}

function syncRules(reset=true) {
  const rec=cases[+$('record').value];
  for(const letter of ['A','B']){
    const select=$('rule'+letter);
    select.replaceChildren(...rec.rules.flatMap((rule,i)=>[
      [String(i),rule.title],...rule.otherOptimal.map((interval,j)=>[`${i}:${j}`,`${rule.title} · optimal ${j+1}`])
    ]).map(([v,title])=>new Option(title,v)),new Option('Custom interval','custom'));
    select.value=letter==='A'?(rec.portable?'0':'2'):'1';
    setRule(letter);
  }
  if(reset)redraw(true);
}

function setRule(letter) {
  const select=$('rule'+letter).value;
  if(select==='custom')return;
  const rec=cases[+$('record').value];
  const [i,j]=select.split(':').map(Number);
  const interval=Number.isInteger(j)?rec.rules[i].otherOptimal[j]:rec.rules[i].span;
  $(letter.toLowerCase()+'0').value=interval[0];
  $(letter.toLowerCase()+'1').value=interval[1];
}

function renderCompare() {
  const rec=cases[+$('record').value];
  $('downloadCase').hidden=!rec.portable;
  const A=['a0','a1'].map(number),B=['b0','b1'].map(number);
  const reverse=$('reverse').checked;
  const r=compare(rec,A,B,reverse);
  const [lo,hi]=reverse?flip(rec.domain,rec.length):rec.domain;
  const geometry=scaffold(lo,hi,i=>{
    const a=inSpan(i,r.a),b=inSpan(i,r.b);
    return a&&b?ink.pale:a?ink.blue:b?ink.rose:ink.muted;
  });
  viewSpec={...geometry,cut:null,track:{lo,hi,title:'Displayed record coordinates',
    rows:[['Interval A',r.a,'blue'],['Interval B',r.b,'rose']]}};
  const selectors=[$('ruleA').value,$('ruleB').value];
  const labels=selectors.map(x=>x.includes(':')||x==='custom'?null:rec.rules[+x]?.native??null);
  const agreement=labels.includes(null)?'Native-event labels are available for selected archived rule outputs only.':
    labels[0]===labels[1]?`The archived labels agree (${labels[0]?'native event':'non-event'}).`:'The archived labels disagree.';
  present({kind:'same_record_span_comparison',record:rec.id,inputs:{A,B,reverse,selectors},result:r,archivedLabels:labels},
    [['Interval A',span(r.a),`${r.lengthA} nt`],['Interval B',span(r.b),`${r.lengthB} nt`],
     ['Shared positions',`${r.overlap} nt`,`${r.onlyA} exclusive to A · ${r.onlyB} exclusive to B`]],
    [['blue','A only'],['rose','B only'],['violet','Complementary strand']],
    `The two intervals yield ${r.identicalString?'identical':'different'} strings. ${agreement} ${reverse?'Reverse-complement coordinates use [N − end, N − start).':''}`,
    rec.portable?'The imported FASTA supplies the base string. Candidate spans are user-defined and the complementary strand is computed.':
      'Both rules operate on one record. The schematic shows observed positions and a computed complement. Changing a span does not rerun an alignment or adjudicate the underlying molecule.',
    reverse?'Reverse-complement coordinate view':'Same-record span comparison');
  $('experiment').replaceChildren(pane('h3','Extracted strings'),
    pane('small','A'),pane('pre',r.sequenceA),pane('small','B'),pane('pre',r.sequenceB));
}

function renderAtoms() {
  const start=number('atomStart'),end=number('atomEnd');
  if(![start,end].every(Number.isInteger)||start<114||end>156||start>=end)
    throw new RangeError('Use a non-empty interval within [114, 156).');
  const representation=$('representation').value,strand=$('strands').value;
  const geometry=atomic(molecule,start-114,end-114,strand,representation);
  viewSpec={...geometry,cut:null,track:{lo:start,hi:end,title:'Read positions in the modeled crop',
    rows:[['Modelled span',[start,end],'violet']]}};
  const seq=molecule.seq.slice(start-114,end-114);
  present({kind:'ideal_duplex_inspection',inputs:{start,end,strand,representation},
           result:{visibleAtoms:geometry.points.length,visibleBonds:geometry.links.length,seq}},
    [['Displayed positions',span([start,end]),`${end-start} paired positions in the source model`],
     ['Visible heavy atoms',String(geometry.points.length),'C, N, O and P'],
     ['Displayed covalent links',String(geometry.links.length),representation==='space'?'Hidden in space-filling view':'For atoms within the crop']],
    [['blue','Observed base sequence'],['violet','Computed complement']],
    'This scene uses the 42-base sequence extracted by the programmed-reference rule. Cropping and strand selection change the displayed atoms. Pick an atom to inspect its original model-space coordinates.',
    'The sequence comes from the archived read. Coordinates are idealized B-DNA from AmberTools. The complementary strand is computed.',
    'CFTR · model-space atomic coordinates');
  $('experiment').replaceChildren(pane('h3','Sequence in the selected crop'),
    pane('small','Observed 5′ → 3′'),pane('pre',seq),pane('small','Complement 5′ → 3′'),pane('pre',revcomp(seq)));
}

function renderImport() {
  if(!imported)throw new RangeError('Choose a PDB, XYZ, FASTA or portable case file first.');
  const chain=$('importChain').value;
  const geometry=importedGeometry(imported,chain);
  for(const id of ['distanceA','distanceB']) {
    $(id).max=geometry.points.length;
    if(number(id)>geometry.points.length)$(id).value=geometry.points.length;
  }
  viewSpec={...geometry,cut:null};
  present({kind:'imported_structure_inspection',source:imported.name,format:imported.format,
          inputs:{chain},result:{visibleAtoms:geometry.points.length,explicitBonds:geometry.links.length}},
    [['Format',imported.format,'Processed in this browser'],
     ['Visible atoms',String(geometry.points.length),chain==='*'?'All chains':`Chain ${chain}`],
     ['Explicit bonds',String(geometry.links.length),'From the supplied file']],
    [['blue','Imported coordinates'],['violet','Additional chains']],
    'Pick two atoms in succession to measure their Euclidean separation in ångströms. The chain control changes the visible subset.',
    'Imported positions are read as supplied. The viewer draws explicit connectivity when present. Distance measures coordinate geometry, with scientific interpretation supplied by the case author.',
    'Imported structure · model-space coordinates');
  const distance=pane('p','');distance.id='distance-result';
  $('experiment').replaceChildren(pane('h3','Local file'),pane('p',imported.name),distance,
    button('Download portable case',()=>download(imported,'portable-structure.case.json')));
  reportDistance();
}

function redraw(fit=false) {
  requestedFit ||= fit;
  if(renderQueued)return;
  renderQueued=true;
  requestAnimationFrame(()=>{
    renderQueued=false;
    const refit=requestedFit;
    requestedFit=false;
    $('input-error').textContent='';
    $('picked').textContent='';
    picked=null;
    pair=null;
    $('export').disabled=false;
    try {
      ({edge:renderEdge,compare:renderCompare,atom:renderAtoms,import:renderImport,align:renderAlign})[mode]();
      view?.setModel(viewSpec.points,viewSpec.links,viewSpec.frame,viewSpec.cut,refit);
      renderTrack(viewSpec.track);
      if(mode==='align'&&alignment.bundle)focusAlignment(currentResult.inputs.at);
    } catch(error) {
      $('input-error').textContent=error.message;
      $('export').disabled=true;
      console.warn(error);
    }
  });
}

function switchMode(value) {
  mode=value;
  $('alignment-evidence').hidden=value!=='align'||!alignment.bundle;
  for(const id of modes)$('control-'+id).hidden=id!==mode;
  for(const tab of document.querySelectorAll('button[data-mode]'))
    tab.setAttribute('aria-pressed',String(tab.dataset.mode===mode));
  redraw(true);
}

function selectPreset(name){
  if(name==='alignment'){
    loadInput(data.alignmentExample);return;
  }
  if(name==='edge'){
    switchMode('edge');setBoundary([210,222,3,224]);return;
  }
  if(name==='same-label'||name==='span'){
    updateCases(name==='same-label'?1:0);$('reverse').checked=false;switchMode('compare');return;
  }
  if(name==='atom'){
    $('atomStart').value=114;$('atomEnd').value=156;$('strands').value='both';$('representation').value='bond';
    switchMode('atom');return;
  }
  throw new Error('Unknown preset.');
}

function download(payload,name) {
  downloadText(JSON.stringify(payload,null,2)+'\n',name,'application/json');
}

function downloadText(text,name,type){
  const url=URL.createObjectURL(new Blob([text],{type}));
  const a=document.createElement('a');
  a.href=url;a.download=name;a.click();
  setTimeout(()=>URL.revokeObjectURL(url),1500);
}

function loadInput(input){
  if(input.kind==='alignment'){
    alignment.load(input);
    document.querySelector('[data-mode=align]').hidden=false;
    switchMode('align');
  }else if(input.kind==='sequence'){
    const half=Math.ceil(input.sequence.length/2);
    const record={id:input.name,length:input.sequence.length,domain:[0,input.sequence.length],
      segment:input.sequence,portable:input,rules:[{title:'Full sequence',span:[0,input.sequence.length],otherOptimal:[]},
      {title:'First half',span:[0,half],otherOptimal:[]},{title:'Second half',span:[half,input.sequence.length],otherOptimal:[]}]};
    cases.push(record);updateCases(cases.length-1);switchMode('compare');
  }else{
    imported=input;
    const chains=[...new Set(input.atoms.map(a=>a.chain))];
    $('importChain').replaceChildren(new Option('All chains','*'),...chains.map(c=>new Option(c,c)));
    document.querySelector('[data-mode=import]').hidden=false;
    switchMode('import');
  }
}

function snapshot(){
  if($('export').disabled||!currentResult)throw new Error('Resolve the input error before exporting.');
  const portable=mode==='align'?alignment.bundle:mode==='import'?imported:mode==='compare'?cases[+$('record').value].portable:null;
  const source=portable?{format:portable.format,name:portable.name}:
    {archivedCommit:audit.commit,auditSHA256:data.audit_sha256,cftrModelSHA256:molecule.sha256};
  return {schema:'read-observability-result/1',...currentResult,source,
    ...(portable?{portableInput:portable}:{}),...(mode==='edge'?{showBoundary:$('cutVisible').checked,pinned}:{}),
    ...(view?{camera:view.cameraState()}:{}),appVersion:'0.3.0'};
}

const nextFrame=()=>new Promise(resolve=>requestAnimationFrame(resolve));

async function replay(record){
  if(record?.schema!=='read-observability-result/1'||!record.inputs)throw new Error('Unsupported replay schema.');
  const i=record.inputs;
  const expectedArchive=record.source?.auditSHA256;
  if(expectedArchive&&expectedArchive!==data.audit_sha256)throw new Error('Replay refers to a different archived audit.');
  if(record.source?.cftrModelSHA256&&record.source.cftrModelSHA256!==molecule.sha256)
    throw new Error('Replay refers to a different atomic model.');
  const input=record.portableInput?parseInputJson(record.portableInput):null;
  let measurement;
  if(record.kind==='coordinate_containment'){
    boundary(i.start,i.end,i.shift,i.L);
    const p=record.pinned?.inputs;
    pinned=p?{inputs:p,result:boundary(p.start,p.end,p.shift,p.L)}:null;
    $('cutVisible').checked=record.showBoundary!==false;
    switchMode('edge');setBoundary([i.start,i.end,i.shift,i.L]);
  }else if(record.kind==='same_record_span_comparison'){
    const recIndex=audit.cases.findIndex(r=>r.id===record.record);
    if(input?.kind==='sequence')loadInput(input);
    else if(recIndex>=0)updateCases(recIndex);
    else throw new Error('Replay has no portable sequence input or known archive record.');
    const rec=cases[+$('record').value];
    compare(rec,i.A,i.B,i.reverse);
    for(const [index,letter] of ['A','B'].entries()){
      const selector=i.selectors?.[index]??'custom',candidate=i[letter];
      const element=$('rule'+letter);
      const [rule,alt]=selector.split(':').map(Number);
      const saved=Number.isInteger(alt)?rec.rules[rule]?.otherOptimal[alt]:rec.rules[rule]?.span;
      element.value=saved&&String(saved)===String(candidate)?selector:'custom';
      $(letter.toLowerCase()+'0').value=candidate[0];$(letter.toLowerCase()+'1').value=candidate[1];
    }
    $('reverse').checked=Boolean(i.reverse);switchMode('compare');
  }else if(record.kind==='ideal_duplex_inspection'){
    $('atomStart').value=i.start;$('atomEnd').value=i.end;
    $('strands').value=i.strand;$('representation').value=i.representation;
    if(!$('strands').value||!$('representation').value)throw new Error('Unsupported atomic display options.');
    switchMode('atom');
  }else if(record.kind==='alignment_locus_audit'&&input?.kind==='alignment'){
    if(!Number.isInteger(i.index)||!input.records[i.index])throw new Error('Invalid replay alignment index.');
    loadInput(input);
    $('alignRecord').value=i.index;
    $('covStart').value=i.lo;$('covEnd').value=i.hi;$('minMapq').value=i.minMapq;
    $('nonPrimary').checked=i.includeNonPrimary===true;
    $('originalRead').checked=i.original===true;
    $('queryFocus').value=i.at;
  }else if(record.kind==='imported_structure_inspection'&&input?.kind==='structure'){
    loadInput(input);
    $('importChain').value=i.chain;
    if(!$('importChain').value)throw new Error('Replay chain is unavailable in the input structure.');
    $('distanceCap').value=i.threshold??4;
    measurement=i.atomIndices;
  }else throw new Error('Unsupported experiment kind or missing structure input.');
  await nextFrame();
  if($('export').disabled)throw new Error($('input-error').textContent);
  if(measurement){
    if(!Array.isArray(measurement)||measurement.length!==2||
      !measurement.every(n=>Number.isInteger(n)&&n>=1&&n<=viewSpec.points.length)||measurement[0]===measurement[1])
      throw new Error('Replay references invalid visible atom indices.');
    pick(viewSpec.points[measurement[0]-1]);pick(viewSpec.points[measurement[1]-1]);
  }
  if(record.camera&&view)view.restoreCamera(record.camera);
  return snapshot();
}

function parseInputJson(value){
  if(value?.schema==='observability-alignment/1')return validateAlignmentCase(value);
  return readCase('portable.json',JSON.stringify(value));
}

async function parseLabInput(name,body){
  const lower=name.toLowerCase();
  if(lower.endsWith('.sam'))return parseSam(name,body,{
    ...(globalThis.crypto?.subtle?{samTextSHA256:await sha256(new TextEncoder().encode(body))}:{})});
  if(lower.endsWith('.json'))return parseInputJson(JSON.parse(body));
  return readCase(name,body);
}

function updateCases(selected=0){
  $('record').replaceChildren(...cases.map((c,i)=>new Option(c.id,String(i))));
  $('record').value=String(selected);
  syncRules(false);
}
updateCases();
for(const tab of document.querySelectorAll('button[data-mode]'))
  tab.addEventListener('click',()=>switchMode(tab.dataset.mode));
for(const id of [...readIds,'cutVisible','reverse','atomStart','atomEnd','strands','representation','importChain'])
  $(id).addEventListener('input',()=>{if(id==='readEnd')syncReadSlider();redraw(['atomStart','atomEnd','importChain'].includes(id));});
for(const letter of ['A','B']) {
  $('rule'+letter).addEventListener('change',()=>{setRule(letter);redraw();});
  for(const edge of ['0','1'])$(letter.toLowerCase()+edge).addEventListener('input',()=>{
    $('rule'+letter).value='custom';redraw();
  });
}
$('record').addEventListener('change',()=>syncRules());
$('readSlide').addEventListener('input',()=>{$('readEnd').value=$('readSlide').value;redraw();});
$('witness').addEventListener('click',()=>setBoundary([210,222,3,224]));
$('oneMore').addEventListener('click',()=>{
  const next=number('readEnd')+1;
  if(next<=+$('readEnd').max)setBoundary([number('start'),number('end'),number('shift'),next]);
});
$('inner').addEventListener('click',()=>setBoundary([203,215,3,224]));
$('detail').addEventListener('click',()=>{$('atomStart').value=129;$('atomEnd').value=141;redraw(true);});
$('allAtoms').addEventListener('click',()=>{$('atomStart').value=114;$('atomEnd').value=156;redraw(true);});
$('resetView').addEventListener('click',()=>view?.fit());
$('track-position').addEventListener('input',()=>{
  const index=number('track-position')-(mode==='atom'?114:0);
  const row=viewSpec.points.find(point=>point.index===index);
  if(row)pick(row);
});
$('alignRecord').addEventListener('change',()=>{alignment.reset();redraw(true);});
$('alignmentDemo').addEventListener('click',()=>selectPreset('alignment'));
for(const id of ['covStart','covEnd','minMapq','nonPrimary','originalRead','queryFocus'])
  $(id).addEventListener('input',()=>redraw());
$('exportTsv').addEventListener('click',()=>{
  try{snapshot();downloadText(alignment.tsv(),'alignment-evidence.tsv','text/tab-separated-values');}
  catch(error){$('input-error').textContent=error.message;}
});
$('exportAlignment').addEventListener('click',()=>download(alignment.bundle,'alignment-case.json'));
$('bamIndex').addEventListener('change',async()=>{
  const index=$('bamIndex').files[0],source=bamBlob,epoch=++loadEpoch;
  $('bamLoad').disabled=true;
  bam=null;
  if(!source||!index)return;
  try{
    const opened=await openIndexedBam(source,index);
    if(epoch!==loadEpoch)return;
    bam=opened;
    $('bamRef').replaceChildren(...Object.keys(opened.refLengths).map(name=>new Option(name,name)));
    resetBamWindow();
    $('bamLoad').disabled=false;
    $('import-status').textContent='BAM header and index selected. Choose the locus to read.';
  }catch(error){if(epoch===loadEpoch)$('import-status').textContent=error.message;}
});
function resetBamWindow(){
  const ref=$('bamRef').value;
  $('bamStart').value=0;
  $('bamEnd').value=Math.min(1000,bam.refLengths[ref]);
}
$('bamRef').addEventListener('change',resetBamWindow);
$('bamLoad').addEventListener('click',async()=>{
  const opened=bam,epoch=loadEpoch;
  if(!opened)return;
  $('bamLoad').disabled=true;
  $('import-status').textContent='Reading indexed locus locally…';
  try{
    const input=await opened.query($('bamRef').value,number('bamStart'),number('bamEnd'));
    if(epoch!==loadEpoch)return;
    loadInput(input);
    $('import-status').textContent=`Loaded indexed locus · ${input.records.length} mapped records`;
  }catch(error){if(epoch===loadEpoch)$('import-status').textContent=error.message;}
  finally{if(epoch===loadEpoch)$('bamLoad').disabled=false;}
});
$('distanceCap').addEventListener('input',reportDistance);
$('measureAtoms').addEventListener('click',()=>{
  const [a,b]=['distanceA','distanceB'].map(id=>viewSpec.points[number(id)-1]);
  if(!a||!b||a===b){$('input-error').textContent='Choose two distinct atom numbers from the visible set.';return;}
  picked=null;pick(a);pick(b);
});
$('clearMeasurement').addEventListener('click',()=>{
  pair=null;picked=null;$('picked').textContent='';view?.select(null);view?.setMeasurement();reportDistance();
});
$('runChecks').addEventListener('click',()=>{
  const tests=validate(audit,molecule);
  $('test-report').textContent=tests.map(({test,pass})=>`${pass?'PASS':'FAIL'} · ${test}`).join('\n');
});
$('export').addEventListener('click',()=>download(snapshot(),'read-observability-replay.json'));
$('renderPng').disabled=!view;
$('renderPng').addEventListener('click',()=>{
  try{
    snapshot();
    const png=view.renderPng(number('renderScale'));
    const a=document.createElement('a');
    a.href=png.dataUrl;a.download=`read-observability-${mode}-${png.width}x${png.height}.png`;
    a.click();
  }catch(error){$('input-error').textContent=error.message;}
});
$('loadCase').addEventListener('change',async()=>{
  const file=$('loadCase').files[0];
  if(!file)return;
  const epoch=++loadEpoch;
  $('import-status').textContent='';
  try {
    if(file.name.toLowerCase().endsWith('.bam')){
      bamBlob=file;bam=null;
      $('bamIndex').value='';$('bamRef').replaceChildren();$('bamLoad').disabled=true;
      $('bam-controls').hidden=false;
      document.querySelector('[data-mode=align]').hidden=false;
      switchMode('align');
      $('import-status').textContent='BAM selected. Load its matching index in Read alignments.';
      return;
    }
    const isJson=file.name.toLowerCase().endsWith('.json');
    if(file.size>(isJson?MAX_REPLAY_BYTES:MAX_CASE_BYTES))
      throw new RangeError('Choose a text case up to 6 MB or replay JSON up to 24 MB.');
    $('bam-controls').hidden=true;bam=null;bamBlob=null;
    const body=await file.text();
    if(epoch!==loadEpoch)return;
    const json=isJson?JSON.parse(body):null;
    if(json?.schema==='read-observability-result/1'){
      await replay(json);
      $('import-status').textContent=`Recomputed ${file.name}`;
    }else{
      const input=json?parseInputJson(json):await parseLabInput(file.name,body);
      if(epoch!==loadEpoch)return;
      loadInput(input);
      $('import-status').textContent=`Loaded ${file.name} · ${input.format}`;
    }
  }catch(error){if(epoch===loadEpoch)$('import-status').textContent=error.message;}
  finally{if(epoch===loadEpoch)$('loadCase').value='';}
});
$('downloadCase').addEventListener('click',()=>{
  const input=cases[+$('record').value].portable;
  if(input)download(input,'portable-sequence.case.json');
  else $('import-status').textContent='Load a FASTA or a portable sequence case to export it.';
});
$('provenance').textContent=`Archived audit commit ${audit.commit}. Atomic coordinates are the CFTR ideal B-DNA model. Input hashes are recorded in every archived-case result export.`;
window.readObservabilityLab=Object.freeze({
  snapshot,replay,
  alignmentTsv(){return mode==='align'&&alignment.bundle?alignment.tsv():null;},
  async preset(name){selectPreset(name);await nextFrame();return snapshot();},
  async open(name,text){loadInput(await parseLabInput(name,text));await nextFrame();return snapshot();},
  renderPng(scale){
    snapshot();
    if(!view)throw new Error('WebGL 2 is required for PNG export.');
    return view.renderPng(scale);
  }
});
redraw(true);
