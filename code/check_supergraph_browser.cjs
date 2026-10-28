'use strict';
// Render the public spatial companion in isolated CI over loopback HTTP.
// This checks interface behavior and frozen-data identity, not biological validity.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const fsp = fs.promises;
const http = require('node:http');
const path = require('node:path');
const crypto = require('node:crypto');
const { chromium } = require(process.env.PLAYWRIGHT_MODULE_PATH || 'playwright-core');
const publicDir = path.resolve('interactive');
const file = process.env.SUPERGRAPH_FILE || 'supergraph/crispr_supergraph.html';
const output = path.resolve(process.env.SUPERGRAPH_QA_OUTPUT || '/tmp/crispr-supergraph-qa');
const originalSHA = '29fc7fce6805bb50b4a6b2c8340d745d07552bad4b8fa7d5b359a7d35bad6445';
const referenceURL = 'https://raw.githubusercontent.com/FacundoFirmenich/KwanTube/a63779af2287956ed208f3e36a5c4f5ea1694a18/outputs_data/interactive/epistemic_graph.html';
const referenceSHA = 'c5fa0f9feb407adc5dcb0c1f0a6e14e3fa1506df5dfb36dd12682b12d43eacd3';
const sha = data => crypto.createHash('sha256').update(data).digest('hex');
const report = {
  scope: 'Actual Chrome rendering of the public graph, local data identity and UI interaction; not independent scientific replication.',
  browserPath: 'Browser plugin unavailable; isolated CI Chrome over loopback HTTP. No Cloud Browser file access.',
  flow: 'Open spatial map -> inspect nodes and relations -> filter preparation and lenses -> pan, zoom and export -> preserve the existing explorer.',
  originalExplorerSHA256: originalSHA, viewports: [], status: 'running',
  limitations: ['Chrome only; emulated desktop/mobile, not physical-device validation.', 'Native touch taps plus keyboard and desktop pointer/wheel actions; physical multi-touch pinch is not certified.'],
};
let browser, server, referenceBytes;
async function state(page) { return page.evaluate(() => window.CRISP_SUPERGRAPH.getState()); }
async function visible(page) { return page.evaluate(() => window.CRISP_SUPERGRAPH.getVisible()); }
async function settle(page) { await page.evaluate(() => new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)))); }
async function capture(page, name) { await settle(page); await bounds(page); await page.screenshot({ path: path.join(output, name + '.png'), fullPage: false }); }
async function bounds(page) {
  const expected = page.viewportSize();
  const data = await page.evaluate(expected => {
    const rect = n => { const r = n.getBoundingClientRect(); return { x:r.x,y:r.y,width:r.width,height:r.height,right:r.right,bottom:r.bottom }; };
    const shown = n => getComputedStyle(n).display !== 'none' && getComputedStyle(n).visibility !== 'hidden' && n.getClientRects().length;
    const boxes = [];
    for (const n of document.querySelectorAll('header, .canvas-wrap, #leftPanel, #rightPanel')) {
      if (!shown(n)) continue;
      const r=rect(n);
      if (r.x < -2 || r.right > expected.width+2 || r.y < -2 || r.bottom > expected.height+2) boxes.push({id:n.id||n.className,...r});
    }
    const textOverflow=[];
    for(const panel of document.querySelectorAll('header,#leftPanel,#rightPanel')){
      if(!shown(panel))continue;const b=panel.getBoundingClientRect(),walker=document.createTreeWalker(panel,NodeFilter.SHOW_TEXT);
      while(walker.nextNode()){const n=walker.currentNode,p=n.parentElement;if(!n.textContent.trim()||!p||!shown(p)||p.closest('.sr-only,select,svg'))continue;
        const r=document.createRange();r.selectNodeContents(n);
        for(const q of r.getClientRects())if(q.width>0&&(q.left<b.left-2||q.right>b.right+2)){textOverflow.push({panel:panel.id||'header',text:n.textContent.trim().slice(0,100),left:q.left,right:q.right});break;}
        r.detach();
      }
    }
    return { expected, innerWidth, innerHeight, clientWidth:document.documentElement.clientWidth, scrollWidth:document.documentElement.scrollWidth,
      graph:rect(document.getElementById('graph')), canvas:rect(document.querySelector('.canvas-wrap')), boxes,textOverflow };
  }, expected);
  assert.ok(data.innerWidth <= expected.width+2 && data.clientWidth <= expected.width+2 && data.scrollWidth <= expected.width+2,
    'Physical viewport width expanded: '+JSON.stringify(data));
  assert.deepEqual(data.boxes, [], 'A visible panel or canvas lies outside the physical viewport: '+JSON.stringify(data));
  assert.deepEqual(data.textOverflow,[], 'Visible panel/header text is clipped horizontally: '+JSON.stringify(data));
  assert.ok(data.graph.width >= expected.width*(expected.width < 700 ? .95 : .3), 'Graph is not the primary canvas: '+JSON.stringify(data));
  assert.ok(data.graph.height >= expected.height*.65, 'Graph does not fill the main viewport: '+JSON.stringify(data));
  return data;
}
async function download(page, selector) {
  const promised=page.waitForEvent('download'); await page.locator(selector).click();
  const d=await promised; assert.equal(await d.failure(),null); const stream=await d.createReadStream(),chunks=[];
  for await(const c of stream)chunks.push(c);
  return {name:d.suggestedFilename(),text:Buffer.concat(chunks).toString('utf8')};
}
async function ensurePanel(page, side, open) {
  const key=side==='left'?'leftHidden':'rightHidden';
  if ((await state(page))[key] === open) await page.locator(side==='left'?'#toggleLeft':'#toggleRight').click();
  await settle(page); assert.equal((await state(page))[key],!open);
}
async function selectPicker(page, id) {
  await ensurePanel(page,'left',true); await page.locator('#nodePicker').selectOption(id); await settle(page);
  const s=await page.evaluate(()=>window.CRISP_SUPERGRAPH.getSelected()); assert.equal(s.id,id); assert.equal(s.kind,'node');
  return s;
}
async function runViewport(baseURL,name,viewport,mobile) {
  const result={name,viewport,checks:[],pageErrors:[],consoleErrors:[],consoleWarnings:[],externalRequests:[],status:'running'};
  report.viewports.push(result);
  const context=await browser.newContext({viewport,isMobile:mobile,hasTouch:mobile,reducedMotion:'reduce',acceptDownloads:true});
  const page=await context.newPage();
  page.on('pageerror',e=>result.pageErrors.push(String(e)));
  page.on('console',m=>{if(m.type()==='error')result.consoleErrors.push(m.text());if(m.type()==='warning')result.consoleWarnings.push(m.text());});
  await context.route('**/*',route=>{
    const url=route.request().url();
    if(url.startsWith('data:')||url.startsWith('blob:')||new URL(url).origin===baseURL)return route.continue();
    result.externalRequests.push(url);return route.abort();
  });
  await page.goto(baseURL+'/'+file,{waitUntil:'networkidle'});
  await page.waitForFunction(()=>window.CRISP_SUPERGRAPH?.ready===true);
  assert.match(await page.title(),/CRISPR/i);
  assert.equal(await page.locator('#graph').count(),1);
  assert.ok(await page.locator('#graph g.node').count() >= 20);
  assert.equal(await page.locator('[role=alert]').count(),0);
  const raw=await page.locator('#supergraph-data').textContent(), model=JSON.parse(raw);
  assert.equal(sha(raw),await page.evaluate(()=>window.CRISP_SUPERGRAPH.dataSHA256));
  const counts=await page.evaluate(()=>window.CRISP_SUPERGRAPH.graphCounts);
  assert.equal(counts.nodes,model.nodes.length); assert.equal(counts.edges,model.edges.length); assert.equal(counts.hyperedges,8);
  assert.equal(model.nodes.length,57);
  assert.doesNotMatch(raw,/microsecond_gap|neural_gap|tau_Lindblad|Warm, wet open.system environment/);
  assert.match(await page.locator('#integrityStat').innerText(),/verif|match|SHA|fingerprint/i);
  await bounds(page); await capture(page,name+'-map');
  result.checks.push('Spatial full-viewport graph identity, 57 nodes/8 conjunctions, frozen-data SHA and no unrelated physics metrics');

  // Verify collapsing panels changes the canvas itself on desktop, rather than hiding a card below a report.
  if(!mobile){
    const before=(await bounds(page)).canvas.width;
    await ensurePanel(page,'left',false); const wider=(await bounds(page)).canvas.width; assert.ok(wider>before+150);
    await ensurePanel(page,'right',false); const full=(await bounds(page)).canvas.width; assert.ok(full>=viewport.width-3);
    await capture(page,name+'-canvas-expanded');
    await ensurePanel(page,'left',true); await ensurePanel(page,'right',true);
  }else{
    assert.equal((await state(page)).leftHidden,true); assert.equal((await state(page)).rightHidden,true);
    await ensurePanel(page,'left',true); await capture(page,name+'-filters'); await ensurePanel(page,'left',false);
  }
  result.checks.push('Desktop panel collapse expands the graph; mobile drawers leave the graph visible by default');

  await ensurePanel(page,'left',true);
  const initial=(await visible(page)).nodes.length;
  await page.locator('#showSources').check();await settle(page);
  assert.equal((await visible(page)).nodes.length,57);
  await selectPicker(page,'S_WITNESS_BOUNDS');
  const sourceNode=model.nodes.find(n=>n.id==='S_WITNESS_BOUNDS');
  assert.match(await page.locator('#iMetrics').textContent(),new RegExp(sourceNode.source_identity.sha256));
  assert.ok((await page.locator('#iArtifacts a').getAttribute('href')).includes('/raw_replay/expected/WITNESS_BOUNDS.json'));
  await ensurePanel(page,'left',true);await page.locator('#showSources').uncheck();await settle(page);
  assert.equal((await visible(page)).nodes.length,initial);
  await page.locator('#search').fill('anchor'); await page.waitForFunction(()=>window.CRISP_SUPERGRAPH.getState().search==='anchor');
  assert.ok((await visible(page)).nodes.length<initial);assert.ok((await visible(page)).nodes.length>0);
  await page.locator('#search').fill('no-matching-scientific-node-zzq');
  await page.waitForFunction(()=>window.CRISP_SUPERGRAPH.getState().search==='no-matching-scientific-node-zzq');
  assert.equal((await visible(page)).nodes.length,0);
  await page.locator('#search').fill('');await page.waitForFunction(()=>window.CRISP_SUPERGRAPH.getState().search==='');
  const lenses=page.locator('#lenses input[type=checkbox]'), lensCount=await lenses.count();assert.ok(lensCount>=4);
  await lenses.first().uncheck();await settle(page);assert.equal((await state(page)).lenses.length,lensCount-1);
  await lenses.first().check();await settle(page);assert.equal((await state(page)).lenses.length,lensCount);
  const views=await page.locator('#viewBtns [data-view]').evaluateAll(xs=>xs.map(x=>x.dataset.view));
  assert.ok(views.includes('all')&&views.length>=5);
  for(const v of views){await page.locator(`#viewBtns [data-view="${v}"]`).click();await settle(page);assert.equal((await state(page)).view,v);}
  await page.locator('#viewBtns [data-view="all"]').click();
  const preps=await page.locator('#preparation option').evaluateAll(xs=>xs.map(x=>x.value)); assert.ok(preps.length>=3);
  for(const v of preps){await page.locator('#preparation').selectOption(v);await settle(page);assert.equal((await state(page)).preparation,v);}
  const domains=await page.locator('#domain option').evaluateAll(xs=>xs.map(x=>x.value));assert.ok(domains.length>=2);
  for(const v of domains){await page.locator('#domain').selectOption(v);await settle(page);assert.equal((await state(page)).domain,v);}
  for(const prep of model.preparations)for(const event of model.target_domains){
    await ensurePanel(page,'left',true);
    await page.locator('#preparation').selectOption(prep.id);await page.locator('#domain').selectOption(event.id);await settle(page);
    assert.equal(await page.locator('#domainDefinition').textContent(),event.definition);
    const rows=model.context_metrics.witness_observations.filter(r=>r.preparation===prep.id&&r.event===event.id);
    assert.equal(await page.locator('#contextMetrics .context-card').count(),3);
    for(let i=0;i<rows.length;i++){const r=rows[i],text=await page.locator('#contextMetrics .context-card').nth(i).textContent();
      assert.ok(text.includes(r.run));
      assert.ok(text.includes('['+r.L_pp.toFixed(5)+', '+r.H_pp.toFixed(5)+']%'));
      assert.ok(text.includes('N '+r.N.toLocaleString('en-US')+' · E '+r.E.toLocaleString('en-US')+' · NE '+r.NE.toLocaleString('en-US')+' · U '+r.U.toLocaleString('en-US')));
    }
  }
  await page.locator('#resetView').click();await settle(page);
  result.checks.push('22 optional source nodes and exact source identity, search/empty-state, lenses, every quick view, all 3 preparation × 2 event witnessed partitions and hulls');

  await ensurePanel(page,'left',false);await ensurePanel(page,'right',false);
  const target=page.locator('#graph [data-node-id="C_anchor"]');
  await target.focus();await page.keyboard.press('Enter');await settle(page);
  assert.equal((await page.evaluate(()=>window.CRISP_SUPERGRAPH.getSelected())).id,'C_anchor');
  assert.match(await page.locator('#iTitle').innerText(),/Observation|product|geometry/i);
  await ensurePanel(page,'right',true);await capture(page,name+'-node-inspector');
  if(mobile){
    await ensurePanel(page,'right',false);
    const touch=page.locator('#graph [data-node-id="C_strings"]');await touch.tap();await settle(page);
    assert.equal((await page.evaluate(()=>window.CRISP_SUPERGRAPH.getSelected())).id,'C_strings');
    assert.equal((await state(page)).rightHidden,false);
    await capture(page,name+'-touch-inspector');
  }
  await selectPicker(page,'H_hull');
  assert.match(await page.locator('#iTitle').innerText(),/AND|hull|Unknown|Counts|premise/i);
  await ensurePanel(page,'left',false);await ensurePanel(page,'right',true);
  assert.match(await page.locator('#rightPanel').innerText(),/premis|AND|conjunct/i);
  await capture(page,name+'-conjunction');
  await ensurePanel(page,'right',false);
  const edge=page.locator('#graph [data-edge-id]').first(), edgeId=await edge.getAttribute('data-edge-id');
  await edge.focus();await page.keyboard.press('Enter');await settle(page);
  assert.equal((await page.evaluate(()=>window.CRISP_SUPERGRAPH.getSelected())).id,edgeId);
  assert.equal((await page.evaluate(()=>window.CRISP_SUPERGRAPH.getSelected())).kind,'edge');
  await ensurePanel(page,'right',true);await capture(page,name+'-edge-inspector');
  result.checks.push('Keyboard node selection, native mobile touch tap, explicit conjunction details and selectable relation inspector');

  await ensurePanel(page,'left',false);await ensurePanel(page,'right',false);
  const zoom=(await state(page)).transform.k;
  await page.locator('#zoomIn').click();assert.ok((await state(page)).transform.k>zoom);
  await page.locator('#zoomOut').click();assert.ok((await state(page)).transform.k<zoom*1.2+1e-9);
  await page.locator('#zoomFit').click();await settle(page);
  await page.locator('#graph').focus();const beforePan=(await state(page)).transform;
  await page.keyboard.press('ArrowRight');await settle(page);assert.notEqual((await state(page)).transform.x,beforePan.x);
  await page.keyboard.press('+');await settle(page);assert.ok((await state(page)).transform.k>beforePan.k);
  await page.keyboard.press('0');await settle(page);
  if(!mobile){
    const box=await page.locator('#graph').boundingBox(), before=(await state(page)).transform;
    await page.mouse.move(box.x+10,box.y+box.height-20);await page.mouse.down();
    await page.mouse.move(box.x+75,box.y+box.height-45,{steps:4});await page.mouse.up();await settle(page);
    const after=(await state(page)).transform;assert.ok(after.x!==before.x||after.y!==before.y);
    await page.mouse.wheel(0,-120);await settle(page);assert.ok((await state(page)).transform.k>after.k);
  }
  await page.locator('#zoomFit').click();await settle(page);
  const exportedSVG=await download(page,'#exportSVG');
  assert.match(exportedSVG.text,/<svg/);
  assert.equal(await page.evaluate(t=>new DOMParser().parseFromString(t,'image/svg+xml').querySelectorAll('parsererror').length,exportedSVG.text),0);
  const exportedJSON=await download(page,'#exportJSON'), json=JSON.parse(exportedJSON.text);
  assert.equal(json.embeddedGraphSHA256,sha(raw));
  assert.deepEqual(json.state,await state(page));
  assert.deepEqual(json.context,model.context_metrics.witness_observations.filter(r=>r.preparation===json.state.preparation&&r.event===json.state.domain));
  const savedState=await state(page);assert.ok(new URL(page.url()).hash.startsWith('#view='));
  await page.reload({waitUntil:'networkidle'});await page.waitForFunction(()=>window.CRISP_SUPERGRAPH?.ready===true);await settle(page);
  assert.deepEqual(await state(page),savedState);
  await capture(page,name+'-final-map');
  result.checks.push('Zoom buttons, keyboard pan/zoom/fit, desktop native pointer pan/wheel and valid SVG/JSON downloads');

  assert.deepEqual(result.pageErrors,[]);assert.deepEqual(result.consoleErrors,[]);assert.deepEqual(result.consoleWarnings,[]);assert.deepEqual(result.externalRequests,[]);
  result.checks.push('No application page/console errors or warnings, no external requests, and viewport-bound panels/canvas');
  result.status='passed';await context.close();
}
async function referenceCapture(baseURL){
  try{
    const r=await fetch(referenceURL,{signal:AbortSignal.timeout(15000)});assert.equal(r.status,200);
    const b=Buffer.from(await r.arrayBuffer());assert.equal(sha(b),referenceSHA);referenceBytes=b;
    const context=await browser.newContext({viewport:{width:1280,height:900},reducedMotion:'reduce'}),page=await context.newPage();
    await page.goto(baseURL+'/__reference.html',{waitUntil:'networkidle'});assert.equal(await page.title(),'Interactive Epistemic Map of Claims, Evidence, and Falsification');
    await page.screenshot({path:path.join(output,'reference-published-interface.png'),fullPage:false});await context.close();
    report.reference={status:'captured',url:referenceURL,sha256:referenceSHA,scope:'Pinned repository interface; archived publication differs only in one version metadata string.'};
  }catch(e){report.reference={status:'unavailable',error:String(e),scope:'Optional external reference capture; current CRISPR app checks remain mandatory.'};}
}
(async()=>{
  await fsp.mkdir(output,{recursive:true});
  assert.equal(sha(await fsp.readFile(path.join(publicDir,'epistemic_graph.html'))),originalSHA,'Previously delivered white explorer must remain byte-identical.');
  server=http.createServer(async(req,res)=>{try{
    const name=decodeURIComponent(new URL(req.url,'http://localhost').pathname);
    if(name==='/favicon.ico'){res.writeHead(204).end();return;}
    if(name==='/__reference.html'&&referenceBytes){res.writeHead(200,{'Content-Type':'text/html; charset=utf-8'}).end(referenceBytes);return;}
    const p=path.resolve(publicDir,'.'+name);if(!p.startsWith(publicDir+path.sep)){res.writeHead(403).end();return;}
    const b=await fsp.readFile(p);res.writeHead(200,{'Content-Type':p.endsWith('.html')?'text/html; charset=utf-8':'application/json; charset=utf-8','Cache-Control':'no-store'}).end(b);
  }catch{res.writeHead(404).end();}});
  await new Promise(r=>server.listen(0,'127.0.0.1',r));const baseURL=`http://127.0.0.1:${server.address().port}`;
  browser=await chromium.launch({executablePath:process.env.CHROME_PATH||'/usr/bin/google-chrome',headless:true});report.browserVersion=browser.version();
  await referenceCapture(baseURL);
  await runViewport(baseURL,'desktop',{width:1280,height:900},false);
  await runViewport(baseURL,'mobile',{width:390,height:844},true);
  report.status='passed';
})().catch(e=>{report.status='failed';report.error=e.stack||String(e);process.exitCode=1;}).finally(async()=>{
  if(browser)await browser.close();if(server)await new Promise(r=>server.close(r));
  await fsp.mkdir(output,{recursive:true});await fsp.writeFile(path.join(output,'supergraph-browser-checks.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report,null,2));
});
