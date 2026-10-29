function renderArgument(){
  const canvas=$('argument-svg');canvas.replaceChildren();
  const id=state.station,scientific=(M.edges||[]).filter(e=>rel(e)!=='reading_order');
  const direct=scientific.filter(e=>endpoint(e,'source')===id||endpoint(e,'target')===id);
  const relevant=new Set([id]);
  for(const e of direct)for(const n of [endpoint(e,'source'),endpoint(e,'target')])if(['formal_statement','formal_derivation'].includes(nodeMap.get(n)?.type))relevant.add(n);
  const conds=(M.conjunctions||[]).filter(c=>relevant.has(c.conclusion||c.target)||direct.some(e=>e.conjunction===c.id));
  const nodeIds=new Set([id]);for(const e of direct){nodeIds.add(endpoint(e,'source'));nodeIds.add(endpoint(e,'target'));}
  for(const c of conds){nodeIds.add(c.id);nodeIds.add(c.conclusion||c.target);for(const p of c.premises||[])nodeIds.add(p);}
  const edges=scientific.filter(e=>nodeIds.has(endpoint(e,'source'))&&nodeIds.has(endpoint(e,'target')));
  const gates=new Set(conds.map(c=>c.id));
  const rightIds=[...nodeIds].filter(n=>n!==id&&!gates.has(n)&&direct.some(e=>endpoint(e,'source')===id&&endpoint(e,'target')===n));
  const leftIds=[...nodeIds].filter(n=>n!==id&&!gates.has(n)&&!rightIds.includes(n));
  const measure=document.createElement('canvas').getContext('2d');
  function linesFor(text,width,size=21,bold=true){
    measure.font=(bold?'700 ':'400 ')+size+'px "DejaVu Sans", sans-serif';
    const lines=[''];for(const word of String(text).split(/\s+/)){
      const i=lines.length-1,candidate=(lines[i]+' '+word).trim();
      if(lines[i]&&measure.measureText(candidate).width>width)lines.push(word);else lines[i]=candidate;
    }return lines;
  }
  const positions=new Map();
  function displayLabel(n){return(nodeMap.get(n)?.type==='refuted_inference'?'Refuted inference: ':'')+label(n);}
  function dimensions(n,w){const gate=gates.has(n),ls=linesFor(displayLabel(n),w-32,21,true);return{w,h:gate?128:Math.max(112,65+ls.length*26.1),lines:ls};}
  let leftY=132,rightY=132;
  for(const n of leftIds){const d=dimensions(n,330);positions.set(n,{x:30,y:leftY,...d});leftY+=d.h+25;}
  for(const n of rightIds){const d=dimensions(n,410);positions.set(n,{x:755,y:rightY,...d});rightY+=d.h+25;}
  const selectedDim=dimensions(id,410),height=Math.max(480,leftY+50,rightY+selectedDim.h+95);
  const centerY=132+(height-182)/2;
  positions.set(id,{x:755,y:Math.max(rightY,centerY-selectedDim.h/2),...selectedDim});
  const gateIds=[...gates];gateIds.forEach((n,i)=>positions.set(n,{x:465,y:Math.max(132,centerY-90+(i-(gateIds.length-1)/2)*205),w:170,h:180,lines:[]}));
  const H=Math.max(height,...[...positions.values()].map(p=>p.y+p.h+45)),W=1200;
  canvas.setAttribute('viewBox',`0 0 ${W} ${H}`);canvas.setAttribute('width',W);canvas.setAttribute('height',H);
  canvas.append(sv('rect',{width:W,height:H,fill:'#fbfcfa'}),sv('title',{},'Scientific argument for '+label(id)),sv('desc',{},'Every edge is an actual typed edge of the source-bound model. A required premise is not sufficient alone. Focus or select a line to read its exact relation and conditions.'));
  const colors={requires:'#765413',derived:'#315c93',support:'#006d72',limit:'#9a3651',defines:'#536c76'};
  function styleFor(r){
    if(r==='requires')return{kind:'requires',dash:'4 4',width:2};
    if(r==='joint_support')return{kind:'requires',dash:'',width:4};
    if(r.includes('derive')||r.includes('computed'))return{kind:'derived',dash:'',width:2.5};
    if(['does_not_establish','qualifies','counterexample_to_interchangeability','distinguishes_scope_from'].includes(r))return{kind:'limit',dash:'10 5',width:2.5};
    if(r==='defines'||r==='motivates')return{kind:'defines',dash:'2 5',width:2};
    return{kind:'support',dash:'',width:2.5};
  }
  const defs=sv('defs');for(const[k,c]of Object.entries(colors)){const mk=sv('marker',{id:'arg-arrow-'+k,markerWidth:8,markerHeight:8,refX:7,refY:4,orient:'auto'});mk.append(sv('path',{d:'M1 1L7 4L1 7Z',fill:c}));defs.append(mk);}canvas.append(defs);
  textSVG(canvas,['SCIENTIFIC ARGUMENT / '+id],30,31,20,'#315c93',true);
  const keys=[['requires (jointly)','requires','4 4',2],['AND joint support','requires','',4],['derive / compute','derived','',2.5],['scope / qualification','limit','10 5',2.5],['define / motivate','defines','2 5',2],['evidence within scope','support','',2.5]];
  keys.forEach((k,i)=>{const col=i%3,row=Math.floor(i/3),x=30+col*390,y=62+row*29;canvas.append(sv('line',{x1:x,x2:x+40,y1:y-4,y2:y-4,stroke:colors[k[1]],'stroke-width':k[3],'stroke-dasharray':k[2]}));textSVG(canvas,[k[0]],x+51,y,17,colors[k[1]]);});
  let edgeLabel=null;
  function inspectEdge(e,g){
    for(const p of canvas.querySelectorAll('.arg-selected-edge'))p.classList.remove('arg-selected-edge');
    g.classList.add('arg-selected-edge');
    const box=$('edge-inspector');box.hidden=false;box.replaceChildren(el('strong',slug(rel(e))),el('p',label(endpoint(e,'source'))+' → '+label(endpoint(e,'target'))),el('p',e.label||e.statement||e.description||M.edgeTypes?.[rel(e)]||''));
    if(e.conditions?.length)box.append(el('p','Conditions: '+e.conditions.map(label).join('; ')));
    if(e.conjunction)box.append(el('p','Conjunction: '+e.conjunction));if(e.scope)box.append(el('p',e.scope));
    if(edgeLabel)edgeLabel.remove();edgeLabel=sv('g',{'class':'argument-edge-label'});edgeLabel.append(sv('rect',{x:385,y:95,width:350,height:32,rx:4,fill:'#ffffff',stroke:'#c8d5d8'}));
    const ls=linesFor(slug(rel(e)),322,16,true);textSVG(edgeLabel,ls,397,116,16,'#183440',true);canvas.append(edgeLabel);
  }
  for(const e of edges){
    const a=positions.get(endpoint(e,'source')),b=positions.get(endpoint(e,'target'));if(!a||!b)continue;
    const r=rel(e),s=styleFor(r),col=colors[s.kind],x1=a.x+a.w,y1=a.y+a.h/2,x2=b.x,y2=b.y+b.h/2;
    let d;if(a.x===b.x)d=`M${a.x+a.w},${y1}C1190,${y1} 1190,${y2} ${b.x+b.w},${y2}`;
    else if(x1>x2)d=`M${a.x},${y1}C${(a.x+b.x+b.w)/2},${y1} ${(a.x+b.x+b.w)/2},${y2} ${b.x+b.w},${y2}`;
    else d=`M${x1},${y1}C${(x1+x2)/2},${y1} ${(x1+x2)/2},${y2} ${x2},${y2}`;
    const g=sv('g',{'data-argument-edge':e.id,'data-relation':r,role:'button',tabindex:0,'aria-label':slug(r)+': '+label(endpoint(e,'source'))+' to '+label(endpoint(e,'target'))});
    g.append(sv('path',{d,fill:'none',stroke:col,'stroke-width':s.width,'stroke-dasharray':s.dash,'marker-end':'url(#arg-arrow-'+s.kind+')'}),sv('path',{d,fill:'none',stroke:'transparent','stroke-width':20}));
    g.addEventListener('click',()=>inspectEdge(e,g));g.addEventListener('focus',()=>inspectEdge(e,g));g.addEventListener('keydown',ev=>{if(ev.key==='Enter'||ev.key===' '){ev.preventDefault();inspectEdge(e,g);}});canvas.append(g);
  }
  for(const[n,p]of positions){
    const gate=gates.has(n),selected=n===id,meta=nodeMap.get(n)||stationMap.get(n)||{},fill=gate?'#faf3e6':['limit','refuted_inference'].includes(meta.type)?'#f9edf0':selected?'#e8f3f0':'#fff',stroke=gate?'#765413':['limit','refuted_inference'].includes(meta.type)?'#9a3651':selected?'#006d72':'#c8d5d8';
    const g=sv('g',{'data-argument-node':n,role:'button',tabindex:0,'aria-label':label(n)});g.append(sv('rect',{x:p.x,y:p.y,width:p.w,height:p.h,rx:9,fill,stroke,'stroke-width':selected?3:1.5}));
    textSVG(g,[n],p.x+15,p.y+26,17,gate?'#765413':'#536c76',true);
    if(gate){const c=conds.find(c=>c.id===n);textSVG(g,['AND'],p.x+15,p.y+61,26,'#765413',true);textSVG(g,[(c.premises||[]).length+' required inputs'],p.x+15,p.y+89,15,'#765413');textSVG(g,linesFor(label(c.conclusion||c.target),p.w-30,15,false),p.x+15,p.y+120,15,'#765413');}
    else textSVG(g,p.lines,p.x+15,p.y+57,21,'#183440',true);
    const inspect=()=>{if(stationMap.has(n)){select(n);return;}const box=$('edge-inspector');box.hidden=false;box.replaceChildren(el('strong',gate?'AND / '+n:label(n)),el('p',meta.statement||meta.title||meta.label||''));if(facts.has(n)){const f=facts.get(n);box.append(el('p',showValue(f)+(f.unit?' '+f.unit:'')));if(f.scope)box.append(el('p',f.scope));for(const ref of f.sourceRefs||[])box.append(el('p',JSON.stringify(ref),'mono'));}if(gate){const c=conds.find(c=>c.id===n);box.append(el('p','Conclusion: '+label(c.conclusion||c.target)),el('p','All required: '+(c.premises||[]).map(label).join(' + ')));}};
    g.addEventListener('click',inspect);g.addEventListener('keydown',ev=>{if(ev.key==='Enter'||ev.key===' '){ev.preventDefault();inspect();}});canvas.append(g);
  }
  if(!edges.length)textSVG(canvas,['This station poses the question. No formal support edge is asserted here.', 'Use the reading route to reach the operational argument.'],70,H-65,22,'#536c76');
  canvas.dataset.renderedEdges=JSON.stringify(edges.map(e=>e.id));canvas.dataset.conjunctions=JSON.stringify([...gates]);
}
