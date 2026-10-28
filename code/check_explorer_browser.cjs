'use strict';
// Rendered checks of the public companion only. No paper sources are required.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const fsp = fs.promises;
const http = require('node:http');
const path = require('node:path');
const crypto = require('node:crypto');
const { chromium } = require(process.env.PLAYWRIGHT_MODULE_PATH || 'playwright-core');

const publicDir = path.resolve('interactive');
const output = path.resolve(process.env.EXPLORER_QA_OUTPUT || '/tmp/crispr-explorer-qa');
const report = {
  scope: 'Rendered public numerical companion; not independent scientific replication.',
  browserPath: 'Isolated CI Chrome over a loopback HTTP server.',
  viewports: [], status: 'running',
};
let browser, server;
const nf = new Intl.NumberFormat('en-US');
const fmt = n => nf.format(n);
const pct = n => (100 * n).toFixed(5) + '%';

async function screenshot(page, viewport, view, anchor = null) {
  if (anchor) {
    await page.locator(anchor).evaluate(n => n.scrollIntoView({ block: 'start', behavior: 'instant' }));
  } else {
    await page.evaluate(() => window.scrollTo(0, 0));
  }
  await noOverflow(page);
  await page.screenshot({ path: path.join(output, `${viewport}-${view}.png`), fullPage: false });
}
async function assertPartition(page) {
  const r = await page.evaluate(() => window.CRISP_EXPLORER.getSelected());
  assert.equal(r.E + r.NE + r.U, r.N);
  for (const [key, val] of [['n', r.N], ['e', r.E], ['ne', r.NE], ['u', r.U]]) {
    assert.equal(await page.locator('#count-' + key).innerText(), fmt(val));
  }
  assert.equal(await page.locator('#read-hull').innerText(), `[${pct(r.E / r.N)}, ${pct((r.E + r.U) / r.N)}]`);
  return r;
}
async function downloadText(page, selector) {
  const promise = page.waitForEvent('download');
  await page.locator(selector).click();
  const d = await promise;
  assert.equal(await d.failure(), null);
  const stream = await d.createReadStream();
  const parts = [];
  for await (const chunk of stream) parts.push(chunk);
  return { name: d.suggestedFilename(), text: Buffer.concat(parts).toString('utf8') };
}
async function setRange(page, value) {
  await page.locator('#completion').evaluate((input, q) => {
    input.value = String(q);
    input.dispatchEvent(new Event('input', { bubbles: true }));
  }, value);
}
async function noOverflow(page) {
  // innerWidth can expand with a mobile layout overflow: use the configured
  // physical CSS viewport as the independent width, then inspect actual text.
  const expected = page.viewportSize().width;
  const geometry = await page.evaluate(expectedWidth => {
    const tolerance = 2;
    const panel = document.querySelector('.panel:not([hidden])');
    const excluded = '.tablewrap,.graph-scroller,.tabs'; // deliberate local scrollers
    const rendered = el => {
      const style = getComputedStyle(el);
      return style.display !== 'none' && style.visibility !== 'hidden' && el.getClientRects().length;
    };
    const crosses = r => r.width > 0 && (r.left < -tolerance || r.right > expectedWidth + tolerance);
    const boxes = [];
    for (const el of panel.querySelectorAll('.card,.toolbar,.hypothetical,.number-strip,.legend,.readout-title,button,select,input')) {
      if (!rendered(el) || el.closest(excluded)) continue;
      const rect = el.getBoundingClientRect();
      if (crosses(rect)) boxes.push({ element: el.id || el.className || el.tagName, left: rect.left, right: rect.right, width: rect.width });
    }
    const textOverflow = [];
    const walker = document.createTreeWalker(panel, NodeFilter.SHOW_TEXT);
    while (walker.nextNode()) {
      const text = walker.currentNode, parent = text.parentElement;
      if (!text.textContent.trim() || !parent || !rendered(parent) || parent.closest(excluded)) continue;
      const range = document.createRange(); range.selectNodeContents(text);
      for (const rect of range.getClientRects()) if (crosses(rect)) {
        textOverflow.push({ text: text.textContent.trim().slice(0, 100), left: rect.left, right: rect.right });
        break;
      }
      range.detach();
    }
    return { expectedWidth, innerWidth, clientWidth: document.documentElement.clientWidth,
      scrollWidth: document.documentElement.scrollWidth, visualViewportWidth: window.visualViewport?.width,
      panel: panel.id, boxes, textOverflow };
  }, expected);
  assert.ok(geometry.clientWidth <= expected + 2, `Layout viewport expanded beyond device: ${JSON.stringify(geometry)}`);
  assert.ok(geometry.innerWidth <= expected + 2, `Mobile innerWidth expanded beyond device: ${JSON.stringify(geometry)}`);
  assert.ok(geometry.scrollWidth <= expected + 2, `Page overflows device: ${JSON.stringify(geometry)}`);
  assert.deepEqual(geometry.boxes, [], `Visible panel controls/cards clipped: ${JSON.stringify(geometry)}`);
  assert.deepEqual(geometry.textOverflow, [], `Visible panel text clipped: ${JSON.stringify(geometry)}`);
  return geometry;
}
async function goTab(page, tab) {
  await page.locator('#tab-' + tab).click();
  assert.equal(await page.locator('#tab-' + tab).getAttribute('aria-selected'), 'true');
  assert.equal(await page.locator('#panel-' + tab).isVisible(), true);
  await noOverflow(page);
}

async function runViewport(baseURL, name, viewport, mobile) {
  const result = { name, viewport, checks: [], consoleErrors: [], consoleWarnings: [], pageErrors: [], externalRequests: [] };
  report.viewports.push(result);
  const context = await browser.newContext({ viewport, isMobile: mobile, hasTouch: mobile,
    reducedMotion: 'reduce', acceptDownloads: true });
  const page = await context.newPage();
  page.on('pageerror', err => result.pageErrors.push(String(err)));
  page.on('console', msg => {
    if (msg.type() === 'error') result.consoleErrors.push(msg.text());
    if (msg.type() === 'warning') result.consoleWarnings.push(msg.text());
  });
  await context.route('**/*', route => {
    if (new URL(route.request().url()).origin !== baseURL) {
      result.externalRequests.push(route.request().url());
      return route.abort();
    }
    return route.continue();
  });
  await page.goto(baseURL + '/epistemic_graph.html', { waitUntil: 'networkidle' });
  assert.equal(await page.title(), 'CRISPR observability · Evidence explorer');
  await page.waitForFunction(() => Boolean(window.CRISP_EXPLORER));
  await page.waitForFunction(() => document.getElementById('integrity').textContent === 'Embedded data fingerprint verified');
  assert.equal(await page.locator('[role=alert]').count(), 0);
  assert.match(await page.locator('main').innerText(), /support path is a conjunction/);
  assert.deepEqual(await page.evaluate(() => window.CRISP_EXPLORER.counts),
    { observations: 144, sensitivity: 288, nodes: 21, hyperedges: 8 });
  const raw = await page.locator('#evidence-data').textContent();
  const model = JSON.parse(raw);
  const embeddedSHA = crypto.createHash('sha256').update(raw.trim()).digest('hex');
  assert.equal(await page.evaluate(() => window.CRISP_EXPLORER.modelSHA256), embeddedSHA);
  assert.deepEqual(JSON.parse(await fsp.readFile(path.join(publicDir, 'evidence_model.json'), 'utf8')), model);
  await noOverflow(page);
  await screenshot(page, name, 'map');
  result.checks.push('Identity, nonblank content, no error overlay, frozen-data integrity and source-object identity');

  // Mobile uses the textual map by default; test the same logical path in both forms.
  const edge = model.graph.hyperedges[0];
  const nodeIndex = model.graph.nodes.findIndex(n => n.id === edge.conclusion);
  const node = model.graph.nodes[nodeIndex];
  if (mobile) {
    assert.equal(await page.locator('#textual-map').getAttribute('open'), '');
    await page.locator('#node-list button').nth(nodeIndex).click();
  } else {
    const target = page.locator('#graph-svg g[role=button]').nth(nodeIndex);
    await target.focus();
    await page.keyboard.press('Enter');
  }
  assert.equal(await page.locator('#inspector h2').innerText(), node.label);
  assert.ok(await page.locator('#inspector input[type=checkbox]').count() >= edge.premises.length);
  await page.locator('#inspector input[type=checkbox]').first().uncheck();
  assert.match(await page.locator('#inspector').innerText(), /AND incomplete/);
  assert.ok(new URL(page.url()).hash.includes('withheld='));
  await page.locator('#inspector input[type=checkbox]').first().check();
  assert.equal(await page.locator('#inspector .gate-note.incomplete').count(), 0);
  result.checks.push('Keyboard/touch node selection and AND-premise withholding/restoration');
  if (mobile) await page.locator('#visual-map summary').click();
  await page.locator('#zoom-in').click();
  assert.equal(await page.evaluate(() => window.CRISP_EXPLORER.getState().z), 1.2);
  assert.equal(await page.locator('#graph-svg').evaluate(n => n.style.width), '120%');
  await page.locator('#graph-svg g[role=button]').first().focus();
  await page.keyboard.press('Enter');
  assert.equal(await page.evaluate(() => window.CRISP_EXPLORER.getState().z), 1.2);
  await page.locator('#zoom-fit').click();
  assert.equal(await page.locator('#graph-svg').evaluate(n => n.style.width), '100%');
  await page.locator('#graph-search').fill('unknown');
  const match = await page.locator('#map-count').innerText();
  assert.doesNotMatch(match, /^21 \/ 21/);
  await page.locator('#reset').click();
  assert.equal(new URL(page.url()).hash, '');
  const mapExport = await downloadText(page, '#export-map');
  assert.match(mapExport.text, /<svg/);
  assert.match(mapExport.text, /http:\/\/www.w3.org\/2000\/svg/);
  assert.equal(await page.evaluate(text => new DOMParser().parseFromString(text, 'image/svg+xml').querySelectorAll('parsererror').length, mapExport.text), 0);
  result.checks.push('Zoom persists across selection, fit/reset/search and valid SVG download');

  await page.locator('#tab-map').focus();
  await page.keyboard.press('ArrowRight');
  assert.equal(await page.locator('#tab-readouts').getAttribute('aria-selected'), 'true');
  await screenshot(page, name, 'readouts');
  let selected = await assertPartition(page);
  const runs = new Set();
  for (const prep of model.catalogs.preparations) {
    await page.locator('#read-prep').selectOption(prep.id);
    const options = await page.locator('#read-run option').evaluateAll(nodes => nodes.map(n => n.value));
    for (const run of options) {
      await page.locator('#read-run').selectOption(run);
      selected = await assertPartition(page); runs.add(selected.run);
    }
  }
  assert.equal(runs.size, 9);
  await page.locator('#read-event').selectOption(model.catalogs.events.find(e => e.id !== selected.event).id);
  selected = await assertPartition(page);
  await setRange(page, 37);
  const k = Math.round(selected.U * .37);
  assert.equal(await page.locator('#completion-result').innerText(), pct((selected.E + k) / selected.N));
  assert.match(await page.locator('#completion-detail').innerText(), new RegExp(`${fmt(k)} of ${fmt(selected.U)} unknown`));
  const state = await page.evaluate(() => window.CRISP_EXPLORER.getState());
  await page.reload({ waitUntil: 'networkidle' });
  await page.waitForFunction(() => Boolean(window.CRISP_EXPLORER));
  assert.deepEqual(await page.evaluate(() => window.CRISP_EXPLORER.getState()), state);
  await assertPartition(page);
  const viewExport = await downloadText(page, '#export-view');
  const exported = JSON.parse(viewExport.text);
  assert.equal(exported.modelSHA256, embeddedSHA);
  assert.equal(exported.selectedObservation.id, selected.id);
  assert.equal(exported.state.q, 37);
  const countsExport = await downloadText(page, '#export-counts');
  assert.match(countsExport.text.split('\n')[0], /"N"/);
  assert.match(countsExport.text.split('\n')[0], /"U"/);
  const intervalExport = await downloadText(page, '#export-interval');
  assert.match(intervalExport.text, /<svg/);
  await screenshot(page, name, 'readouts-card', '#read-heading');
  await screenshot(page, name, 'readouts-completion', '.hypothetical');
  await noOverflow(page);
  result.checks.push('Nine raw streams, event change, exact count/hull rendering, attainable integer completion, hash restore, JSON/CSV/SVG exports');

  await goTab(page, 'sensitivity');
  function sensitivityExpected(s) {
    return model.sensitivityContrasts.filter(r => r.preparation === s.sprep && r.event === s.sevent && r.control === s.scontrol && (s.sref === 'both' || r.reference === s.sref));
  }
  let ss = await page.evaluate(() => window.CRISP_EXPLORER.getState());
  let expected = sensitivityExpected(ss);
  assert.equal(expected.length, 16);
  assert.equal(await page.locator('#sens-mobile-cards > div').count(), expected.length);
  assert.match(await page.locator('#sens-caption').innerText(), /16 reference\/configuration hulls/);
  await screenshot(page, name, 'sensitivity');
  await page.locator('#sens-reference').selectOption('mutant');
  await page.locator('#sens-control').selectOption('untreated');
  ss = await page.evaluate(() => window.CRISP_EXPLORER.getState());
  expected = sensitivityExpected(ss);
  assert.equal(expected.length, 8);
  assert.equal(await page.locator('#sens-mobile-cards > div').count(), 8);
  const gridExport = await downloadText(page, '#export-grid');
  assert.equal(gridExport.text.trim().split('\r\n').length, 9);
  const sensExport = await downloadText(page, '#export-sensitivity-svg');
  assert.match(sensExport.text, /<svg/);
  result.checks.push('Sixteen sensitivity configurations, reference/control filters and grid/SVG downloads');

  await goTab(page, 'contract');
  assert.equal(await page.locator('#scope-records').innerText(), '177,006');
  assert.equal(await page.locator('#scope-witness').innerText(), '170,725');
  assert.equal(await page.locator('#scope-unknown').innerText(), '6,281');
  assert.equal(await page.locator('#sources-table tr').count(), model.sources.length);
  const dataExport = await downloadText(page, '#export-model');
  assert.equal(crypto.createHash('sha256').update(dataExport.text).digest('hex'), embeddedSHA);
  await screenshot(page, name, 'contract');
  result.checks.push('Contract arithmetic, authentic aggregate counts, complete source registry and exact embedded-data export');

  assert.deepEqual(result.pageErrors, []);
  assert.deepEqual(result.consoleErrors, []);
  assert.deepEqual(result.consoleWarnings, []);
  assert.deepEqual(result.externalRequests, []);
  result.checks.push('No page/console errors, no network dependency; device-width, card/control and text geometry unclipped');
  result.status = 'passed';
  await context.close();
}

(async () => {
  await fsp.mkdir(output, { recursive: true });
  server = http.createServer(async (req, res) => {
    try {
      const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
      if (pathname === '/favicon.ico') { res.writeHead(204).end(); return; }
      const p = path.resolve(publicDir, '.' + pathname);
      if (!p.startsWith(publicDir + path.sep)) { res.writeHead(403).end(); return; }
      const data = await fsp.readFile(p);
      res.writeHead(200, { 'Content-Type': p.endsWith('.html') ? 'text/html; charset=utf-8' : 'application/json; charset=utf-8', 'Cache-Control': 'no-store' }).end(data);
    } catch { res.writeHead(404).end(); }
  });
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  const baseURL = `http://127.0.0.1:${server.address().port}`;
  report.url = baseURL + '/epistemic_graph.html';
  browser = await chromium.launch({ executablePath: process.env.CHROME_PATH || '/usr/bin/google-chrome', headless: true });
  report.browserVersion = browser.version();
  await runViewport(baseURL, 'desktop', { width: 1280, height: 900 }, false);
  await runViewport(baseURL, 'mobile', { width: 390, height: 844 }, true);
  report.status = 'passed';
})().catch(err => {
  report.status = 'failed'; report.error = err.stack || String(err); process.exitCode = 1;
}).finally(async () => {
  if (browser) await browser.close();
  if (server) await new Promise(resolve => server.close(resolve));
  await fsp.writeFile(path.join(output, 'browser-checks.json'), JSON.stringify(report, null, 2) + '\n');
  console.log(JSON.stringify(report, null, 2));
});
