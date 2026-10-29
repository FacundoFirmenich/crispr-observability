'use strict';
// Rendered checks of the public PaperGraph only. Formal article/supplement sources are not required.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const fsp = fs.promises;
const path = require('node:path');
const http = require('node:http');
const crypto = require('node:crypto');
const {chromium} = require(process.env.PLAYWRIGHT_MODULE_PATH || 'playwright-core');
const ROOT = path.resolve(__dirname, '..');
const COMPONENT = path.join(ROOT, 'interactive', 'papergraph');
const OUT = path.resolve(process.env.PAPERGRAPH_QA_OUTPUT || path.join(require('node:os').tmpdir(), 'crispr-papergraph-qa'));
assert.ok(OUT !== ROOT && !OUT.startsWith(ROOT + path.sep), 'Browser QA output must remain outside the public repository');
const CHROME = process.env.CHROME_BINARY || '/usr/bin/google-chrome';
const report = {date:'2026-10-29',status:'running',scope:'Public graphical computational companion: rendered UI, source-model roundtrips and stated interface contracts; not molecular accuracy, independent biological replication or pedagogical efficacy.',browserRoute:'Browser plugin not available. Playwright with local Chrome over loopback HTTP; outgoing application requests blocked. QA evidence is written outside the repository.',viewports:[]};
let browser, server, activePage;
const hash = s => crypto.createHash('sha256').update(s).digest('hex');
const state = p => p.evaluate(()=>window.PAPERGRAPH.getState());
async function bounds(page){
  const expected=page.viewportSize();
  const actual=await page.evaluate(expected=>{
    const shown=n=>!!n.getClientRects().length&&getComputedStyle(n).display!=='none'&&getComputedStyle(n).visibility!=='hidden';
    const bad=[];
    for(const n of document.querySelectorAll('header,.reading-bar,.map-card,.detail,.relation-card,dialog[open]')){
      if(!shown(n))continue;const b=n.getBoundingClientRect();
      if(b.left < -2 || b.right > expected.width+2)bad.push({kind:'element',tag:n.id||n.className,left:b.left,right:b.right});
    }
    const ranges=[];
    for(const n of document.querySelectorAll('.detail,.relation-card,header,.reading-bar,dialog[open]')){
      if(!shown(n))continue;const b=n.getBoundingClientRect(),walker=document.createTreeWalker(n,NodeFilter.SHOW_TEXT);
      while(walker.nextNode()){
        const t=walker.currentNode,p=t.parentElement;
        if(!t.textContent.trim()||!p||!shown(p)||p.closest('svg,select,.live,.skip'))continue;
        const range=document.createRange();range.selectNodeContents(t);
        for(const q of range.getClientRects())if(q.width>0&&(q.left<b.left-2||q.right>b.right+2)){ranges.push({text:t.textContent.trim().slice(0,80),left:q.left,right:q.right,parent:n.id||n.className});break;}
      }
    }
    return {expected,innerWidth,clientWidth:document.documentElement.clientWidth,scrollWidth:document.documentElement.scrollWidth,bad,ranges};
  },expected);
  assert.ok(actual.innerWidth<=expected.width+2&&actual.clientWidth<=expected.width+2&&actual.scrollWidth<=expected.width+2,'Physical viewport expands: '+JSON.stringify(actual));
  assert.deepEqual(actual.bad,[],'A card or control region overflows: '+JSON.stringify(actual));
  assert.deepEqual(actual.ranges,[],'Text is clipped horizontally: '+JSON.stringify(actual));
  return actual;
}
async function capture(page,name){await bounds(page);await page.screenshot({path:path.join(OUT,name+'.png'),fullPage:false});}
async function downloaded(page,id){const promise=page.waitForEvent('download');await page.locator(id).click();const d=await promise;assert.equal(await d.failure(),null);const stream=await d.createReadStream(),chunks=[];for await(const b of stream)chunks.push(b);return{name:d.suggestedFilename(),text:Buffer.concat(chunks).toString('utf8')};}
async function checkStaticFigures(base){
  report.staticFigures=[];
  report.staticFigureExternalRequests=[];
  const context=await browser.newContext({viewport:{width:1280,height:900}});
  await context.route('**/*',route=>{const u=route.request().url();if(u.startsWith('data:')||u.startsWith('blob:')||new URL(u).origin===base)return route.continue();report.staticFigureExternalRequests.push(u);return route.abort();});
  const page=await context.newPage();
  for(const name of ['PaperGraph_Overview.svg','PaperGraph_Article.svg','PaperGraph_Argument.svg']){
    await page.goto(base+'/figures/papergraph/'+name,{waitUntil:'networkidle'});await page.evaluate(()=>document.fonts.ready);
    const failures=await page.evaluate(()=>{
      const bad=[];for(const g of document.querySelectorAll('g.station')){const box=g.querySelector('.station-box').getBBox();for(const t of g.querySelectorAll('text')){const b=t.getBBox();if(b.x<box.x+3||b.x+b.width>box.x+box.width-3||b.y<box.y||b.y+b.height>box.y+box.height){bad.push({station:g.dataset.station,text:t.textContent,bounds:{x:b.x,y:b.y,width:b.width,height:b.height},card:{x:box.x,y:box.y,width:box.width,height:box.height}});}}}
      for(const g of document.querySelectorAll('g[data-node-id]')){const rect=g.querySelector('rect');if(!rect)continue;const box=rect.getBBox();for(const t of g.querySelectorAll('text')){if(!t.textContent.trim())continue;const b=t.getBBox();if(b.x<box.x+2||b.x+b.width>box.x+box.width-2||b.y<box.y||b.y+b.height>box.y+box.height)bad.push({atlasNode:g.dataset.nodeId,text:t.textContent,bounds:{x:b.x,y:b.y,width:b.width,height:b.height},card:{x:box.x,y:box.y,width:box.width,height:box.height}});}}
      const vb=document.documentElement.viewBox.baseVal;for(const t of document.querySelectorAll('text')){const b=t.getBBox();if(b.x< -1||b.x+b.width>vb.width+1||b.y< -1||b.y+b.height>vb.height+1)bad.push({outsideCanvas:true,text:t.textContent,bounds:{x:b.x,y:b.y,width:b.width,height:b.height},canvas:{width:vb.width,height:vb.height}});}
      return bad;
    });
    report.staticFigures.push({name,failures,status:failures.length?'failed':'passed'});
  }
  await context.close();
  assert.deepEqual(report.staticFigures.filter(f=>f.failures.length),[],'A static figure contains text outside its card or SVG canvas');
  assert.deepEqual(report.staticFigureExternalRequests,[],'A static figure requests an external resource');
}
async function runViewport(base,name,viewport,mobile){
  const result={name,viewport,status:'running',checks:[],pageErrors:[],consoleErrors:[],externalRequests:[]};report.viewports.push(result);
  const ctx=await browser.newContext({viewport,isMobile:mobile,hasTouch:mobile,reducedMotion:'reduce',acceptDownloads:true});
  await ctx.route('**/*',route=>{const u=route.request().url();if(u.startsWith('data:')||u.startsWith('blob:')||new URL(u).origin===base)return route.continue();result.externalRequests.push(u);return route.abort();});
  const page=await ctx.newPage();activePage=page;page.on('pageerror',e=>result.pageErrors.push(String(e)));page.on('console',m=>{if(m.type()==='error')result.consoleErrors.push(m.text());});
  await page.goto(base+'/interactive/papergraph/PaperGraph_CRISPR.html',{waitUntil:'networkidle'});
  await page.waitForFunction(()=>window.PAPERGRAPH?.ready===true);
  assert.match(await page.title(),/CRISPR PaperGraph/);
  assert.equal(await page.locator('[role=alert]').count(),0);
  assert.ok((await page.locator('body').innerText()).toLowerCase().includes('papergraph'),'Page does not contain PaperGraph identity');
  assert.equal(await page.locator('#station-picker option').count(),12);
  const raw=await page.locator('#papergraph-data').textContent(),m=JSON.parse(raw);
  assert.equal(hash(raw),await page.evaluate(()=>window.PAPERGRAPH.modelSHA256));
  result.embeddedModelSHA256=hash(raw);result.inputModelSHA256=await page.evaluate(()=>window.PAPERGRAPH.inputModelSHA256);
  assert.equal(await page.evaluate(()=>window.PAPERGRAPH.counts.stations),12);
  assert.equal(m.stations.length,12);
  assert.equal(m.facts.length,49);
  assert.equal(m.nodes.length,92);
  assert.equal(m.edges.length,86);
  assert.equal(m.conjunctions.length,7);
  assert.equal(m.edges.filter(e=>e.relation==='reading_order').length,11);
  assert.equal(m.edges.filter(e=>e.relation!=='reading_order').length,75);
  assert.equal(m.sectionCoverage,undefined,'The private manuscript section map must not be public');
  assert.equal(m.integrationRecommendation,undefined,'Private editorial integration instructions must not be public');
  const inspectPublic = value => {
    if(Array.isArray(value)){value.forEach(inspectPublic);return;}
    if(!value || typeof value!=='object')return;
    for(const key of Object.keys(value))assert.ok(!['literalAnchor','originalPath'].includes(key),'Private reference field found in public model: '+key);
    Object.values(value).forEach(inspectPublic);
  };
  inspectPublic(m);
  assert.ok(!JSON.stringify(m).includes('source_snapshots/'),'Private snapshots must not be embedded in the public model');
  const canonical=JSON.parse(await fsp.readFile(path.join(COMPONENT,'papergraph.json'),'utf8'));
  assert.deepEqual(m,canonical,'The rendered embedded model differs from the public canonical JSON');
  await capture(page,name+'-overview');
  result.checks.push('Nonblank public companion, twelve stations, no error overlay and exact embedded-data SHA');
  const factMap=new Map((Array.isArray(m.facts)?m.facts:Object.entries(m.facts).map(([id,f])=>({id,...f}))).map(f=>[f.id,f]));
  for(const station of m.stations){
    await page.locator('#station-picker').selectOption(station.id);
    assert.equal(await page.locator('#detail-title').innerText(),station.title);
    for(const depth of ['student','technical','audit']){
      await page.locator('[data-depth="'+depth+'"]').click();
      assert.equal((await state(page)).depth,depth);
      assert.equal(await page.locator('[data-depth="'+depth+'"]').getAttribute('aria-pressed'),'true');
      assert.equal((await page.locator('#explanation').innerText()).replace(/\s+/g,' ').trim(),station.levels[depth].replace(/\s+/g,' ').trim());
      assert.equal(await page.locator('#revision').innerText(),station.revisionCondition);
      assert.equal(await page.locator('#consequence').innerText(),station.consequence);
      assert.equal(await page.locator('#facts [data-fact]').count(),station.factIds.length);
      assert.deepEqual(await page.locator('#facts [data-fact]').evaluateAll(ns=>ns.map(n=>n.dataset.fact)),station.factIds);
      assert.deepEqual(await page.evaluate(()=>window.PAPERGRAPH.getSelected()),station);
      for(const id of station.factIds)assert.ok(factMap.has(id));
      await bounds(page);
    }
  }
  result.checks.push('All twelve stations × three depth views preserve source content, fact identity, conditions and consequences');
  const edgeMap=new Map(m.edges.map(e=>[e.id,e]));
  await page.locator('[data-mapmode=argument]').click();
  for(const station of m.stations){
    await page.locator('#station-picker').selectOption(station.id);
    assert.equal((await state(page)).mapmode,'argument');
    assert.equal(await page.locator('#argument-graphic').isVisible(),true);
    const ids=await page.locator('#argument-svg').evaluate(n=>JSON.parse(n.dataset.renderedEdges));
    for(const id of ids){assert.ok(edgeMap.has(id),'Displayed edge missing from source model');assert.notEqual(edgeMap.get(id).relation,'reading_order');}
    if(station.id==='PG05'){
      const refuted=page.locator('#argument-svg [data-argument-node="R_equal_labels_imply_equal_strings"]');
      assert.equal(await refuted.count(),1);
      assert.match(await refuted.textContent(),/Refuted inference:/i,'Refuted status must be explicit text, not colour alone');
    }
    const directGates=m.conjunctions.filter(c=>c.conclusion===station.id),actualGates=await page.locator('#argument-svg').evaluate(n=>JSON.parse(n.dataset.conjunctions));
    for(const c of directGates)assert.ok(actualGates.includes(c.id),'Direct declared conjunction is missing');
    for(const gate of actualGates){const c=m.conjunctions.find(c=>c.id===gate);assert.ok(c,'Rendered conjunction absent from source model');const nodeIds=await page.locator('#argument-svg [data-argument-node]').evaluateAll(ns=>ns.map(n=>n.dataset.argumentNode));assert.ok(nodeIds.includes(c.id));assert.ok(nodeIds.includes(c.conclusion),'Conjunction conclusion is omitted');for(const premise of c.premises)assert.ok(nodeIds.includes(premise),'A conjunctive premise is omitted from the graph');}
    const scientificLayout=await page.locator('#argument-svg').evaluate(svg=>{
      const nodeText=[];for(const g of svg.querySelectorAll('[data-argument-node]')){const box=g.querySelector('rect').getBBox();for(const t of g.querySelectorAll('text')){const b=t.getBBox();if(b.x<box.x+2||b.x+b.width>box.x+box.width-2||b.y<box.y||b.y+b.height>box.y+box.height)nodeText.push({node:g.dataset.argumentNode,text:t.textContent,bounds:{x:b.x,y:b.y,width:b.width,height:b.height},card:{x:box.x,y:box.y,width:box.width,height:box.height}});}}
      const visibleText=t=>{if(!t.getClientRects().length)return false;for(let n=t;n&&n!==svg.parentElement;n=n.parentElement){const s=getComputedStyle(n);if(s.display==='none'||s.visibility==='hidden'||Number(s.opacity)===0)return false;}return true;};
      const edgeText=[...svg.querySelectorAll('[data-argument-edge] text')].filter(visibleText).map(t=>({id:t.closest('[data-argument-edge]').dataset.argumentEdge,text:t.textContent,box:t.getBBox()}));
      const overlaps=[];for(let i=0;i<edgeText.length;i++)for(let j=i+1;j<edgeText.length;j++){const a=edgeText[i],b=edgeText[j];if(a.id===b.id)continue;const w=Math.min(a.box.x+a.box.width,b.box.x+b.box.width)-Math.max(a.box.x,b.box.x),h=Math.min(a.box.y+a.box.height,b.box.y+b.box.height)-Math.max(a.box.y,b.box.y);if(w>3&&h>3)overlaps.push({a:a.id,b:b.id,textA:a.text,textB:b.text,intersection:{width:w,height:h}});}
      return {nodeText,overlaps};
    });
    assert.deepEqual(scientificLayout.nodeText,[],station.id+' dynamic node text exceeds its rectangle');
    assert.deepEqual(scientificLayout.overlaps,[],station.id+' edge labels collide');
    await bounds(page);
  }
  await page.locator('#station-picker').selectOption('PG04');
  const edgeNode=page.locator('#argument-svg [data-argument-edge]').first();await edgeNode.focus();await page.keyboard.press('Enter');assert.equal(await page.locator('#edge-inspector').isVisible(),true);
  assert.ok((await page.locator('#edge-inspector').innerText()).trim());
  await page.locator('#argument-graphic').scrollIntoViewIfNeeded();await capture(page,name+'-scientific-argument');
  const argumentSVG=await downloaded(page,'#export-svg');assert.match(argumentSVG.name,/Argument/);assert.match(argumentSVG.text,/<svg/);
  result.checks.push('Scientific argument mode draws source-model edges only, preserves every declared AND premise, supports keyboard edge inspection and exports its graph');
  await page.locator('[data-mapmode=guide]').click();
  await page.locator('#station-picker').selectOption(m.stations[0].id);
  await page.locator('#next').click();assert.equal((await state(page)).station,m.stations[1].id);
  await page.locator('#prev').click();assert.equal((await state(page)).station,m.stations[0].id);
  const formal=m.stations.find(s=>s.id==='PG04')||m.stations[4];
  if(mobile){await page.locator('#mobile-route [data-station="'+formal.id+'"]').tap();}
  else {const n=page.locator('#graphic [data-station="'+formal.id+'"]:visible').first();await n.focus();await page.keyboard.press('Enter');}
  assert.equal((await state(page)).station,formal.id);
  await page.locator('#focus-path').click();assert.equal((await state(page)).focus,true);
  assert.equal(await page.locator('#focus-path').getAttribute('aria-pressed'),'true');
  await page.locator('#detail').scrollIntoViewIfNeeded();await capture(page,name+'-conditional-bounds');
  const saved=await state(page);await page.reload({waitUntil:'networkidle'});await page.waitForFunction(()=>window.PAPERGRAPH?.ready);assert.deepEqual(await state(page),saved);
  await page.locator('#overview').click();assert.equal((await state(page)).focus,false);
  result.checks.push('Guided previous/next, keyboard or touch selection, focus path, overview and URL state restore');
  await page.locator('#glossary-open').click();assert.equal(await page.locator('#glossary-dialog').evaluate(n=>n.open),true);
  const term=m.glossary[0].term,definition=m.glossary[0].plain||m.glossary[0].definition||m.glossary[0].student||m.glossary[0].text;
  await page.locator('#glossary-search').fill(term);assert.match(await page.locator('#glossary-list').innerText(),new RegExp(term.replace(/[.*+?^${}()|[\]\\]/g,'\\$&'),'i'));
  assert.ok(definition&&await page.locator('#glossary-list').innerText().then(t=>t.includes(definition)),'Glossary definition is absent');
  await capture(page,name+'-glossary');await page.locator('#glossary-close').click();assert.equal(await page.locator('#glossary-dialog').evaluate(n=>n.open),false);
  const svg=await downloaded(page,'#export-svg');assert.match(svg.text,/<svg/);assert.equal(await page.evaluate(s=>new DOMParser().parseFromString(s,'image/svg+xml').querySelectorAll('parsererror').length,svg.text),0);
  const data=await downloaded(page,'#export-data');assert.equal(hash(data.text),hash(raw));
  const view=JSON.parse((await downloaded(page,'#export-view')).text);assert.equal(view.embeddedModelSHA256,hash(raw));assert.deepEqual(view.state,await state(page));
  result.checks.push('Searchable glossary and parseable SVG, exact JSON data and source-bound view exports');
  await page.locator('#station-picker').selectOption(m.stations[8].id);await page.locator('[data-depth=audit]').click();await page.locator('#detail').scrollIntoViewIfNeeded();await capture(page,name+'-sensitivity-audit');
  assert.deepEqual(result.pageErrors,[]);assert.deepEqual(result.consoleErrors,[]);assert.deepEqual(result.externalRequests,[]);
  result.checks.push('No application runtime errors, no external requests, no page/card/text horizontal clipping');
  result.status='passed';await ctx.close();
}
(async()=>{
  await fsp.mkdir(OUT,{recursive:true});
  report.checkerSHA256=hash(await fsp.readFile(__filename));
  report.htmlSHA256=hash(await fsp.readFile(path.join(COMPONENT,'PaperGraph_CRISPR.html')));
  server=http.createServer(async(req,res)=>{try{const pathname=decodeURIComponent(new URL(req.url,'http://localhost').pathname);if(pathname==='/favicon.ico'){res.writeHead(204).end();return;}const p=path.resolve(ROOT,'.'+pathname);if(!p.startsWith(ROOT+path.sep)){res.writeHead(403).end();return;}const b=await fsp.readFile(p);res.writeHead(200,{'Content-Type':p.endsWith('.html')?'text/html; charset=utf-8':p.endsWith('.svg')?'image/svg+xml':'application/json; charset=utf-8','Cache-Control':'no-store'}).end(b);}catch{res.writeHead(404).end();}});
  await new Promise(r=>server.listen(0,'127.0.0.1',r));const base='http://127.0.0.1:'+server.address().port;report.url=base+'/interactive/papergraph/PaperGraph_CRISPR.html';
  browser=await chromium.launch({executablePath:CHROME,headless:true});report.browserVersion=browser.version();
  await checkStaticFigures(base);
  await runViewport(base,'desktop',{width:1280,height:900},false);await runViewport(base,'mobile',{width:390,height:844},true);report.status='passed';
})().catch(async e=>{report.status='failed';report.error=e.stack||String(e);process.exitCode=1;if(activePage&&!activePage.isClosed()){try{await activePage.screenshot({path:path.join(OUT,'failure-state.png'),fullPage:false});}catch{}}}).finally(async()=>{if(browser)await browser.close();if(server)await new Promise(r=>server.close(r));await fsp.writeFile(path.join(OUT,'PAPERGRAPH_BROWSER_QA.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report,null,2));});
