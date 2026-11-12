'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const { pathToFileURL } = require('node:url');
const { chromium } = require(process.env.PLAYWRIGHT_MODULE_PATH || 'playwright-core');

(async () => {
  const browser = await chromium.launch({
    executablePath: process.env.CHROME_PATH,
    headless: true,
  });
  try {
    const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    const entry = path.resolve('interactive/supergraph/crispr_supergraph.html');
    await page.goto(pathToFileURL(entry).href);
    await page.waitForFunction(() => window.CRISP_SUPERGRAPH?.ready);
    const state = () => page.evaluate(() => window.CRISP_SUPERGRAPH.getState());
    const nav = page.getByRole('navigation', { name: 'Evidence categories' });
    assert.equal(await nav.getByRole('button').count(), 7);
    const navBox = await nav.boundingBox(), graphBox = await page.locator('#graph').boundingBox();
    assert.ok(navBox.y + navBox.height <= graphBox.y + 1);
    await nav.getByRole('button', { name: 'External validation' }).click();
    assert.equal((await state()).activeCategory, 'External validation');
    await nav.getByRole('button', { name: 'Observation mechanisms' }).click();
    const card = page.locator('[data-node-id="C_strings"]');
    const baseline = await card.boundingBox();
    assert.ok(baseline.width > 140);
    const model = await page.locator('#supergraph-data').textContent().then(JSON.parse);
    const expected = async id => {
      const visible = await page.evaluate(() => window.CRISP_SUPERGRAPH.getVisible().edges);
      return model.edges.filter(e => visible.includes(e.id) && [e.source,e.target].includes(id)).map(e => e.id).sort();
    };
    const traced = () => page.locator('[data-trace-edge]').evaluateAll(items => items.map(n => n.dataset.traceEdge).sort());
    await card.click();
    assert.deepEqual(await traced(), await expected('C_strings'));
    assert.equal(await page.locator('.edge-signal').count(), (await expected('C_strings')).length);
    assert.ok(await page.locator('.edge-trace.incoming').count() > 0);
    await page.emulateMedia({ reducedMotion:'reduce' });
    assert.equal(await page.locator('.edge-signal').first().evaluate(e => getComputedStyle(e).display), 'none');
    assert.notEqual(await page.locator('.edge-trace').first().evaluate(e => getComputedStyle(e).display), 'none');
    await page.emulateMedia({ reducedMotion:'no-preference' });
    await page.getByRole('button', { name: 'Fit connections' }).click();
    const framed = await page.locator('g.node[data-focus],g.node[data-context]').evaluateAll(cards => {
      const map = document.getElementById('graph').getBoundingClientRect();
      return cards.every(card => {
        const box = card.getBoundingClientRect();
        return box.left >= map.left - 1 && box.right <= map.right + 1 &&
          box.top >= map.top - 1 && box.bottom <= map.bottom + 1;
      });
    });
    assert.equal(framed, true);
    if (process.env.SUPERGRAPH_QA_OUTPUT) {
      fs.mkdirSync(process.env.SUPERGRAPH_QA_OUTPUT, { recursive: true });
      await page.screenshot({ path: path.join(process.env.SUPERGRAPH_QA_OUTPUT, 'supergraph-connections.png') });
    }
    await nav.getByRole('button', { name: 'Observation mechanisms' }).click();
    const column = page.locator('.column-handle').nth(1);
    const drag = async (locator, dx, dy) => {
      const box = await locator.boundingBox();
      assert.ok(box);
      const x = box.x + box.width / 2, y = box.y + box.height / 2;
      await page.mouse.move(x, y);
      await page.mouse.down();
      await page.mouse.move(x + dx, y + dy, { steps: 8 });
      const id = await locator.getAttribute('data-node-id');
      if (id && !(await state()).locked) {
        await page.waitForFunction(node => document.querySelector('g.node[data-focus]')?.dataset.nodeId === node, id);
        assert.deepEqual(await traced(), await expected(id));
      }
      await page.mouse.up();
      await page.waitForFunction(() => !document.getElementById('graph').classList.contains('dragging'));
    };
    const original = await state();
    const other = page.locator('[data-node-id="C_anchor"]');
    await drag(other, 40, 30);
    assert.deepEqual(await traced(), await expected('C_strings'));
    const moved = await state();
    assert.equal(moved.nodeOffsets.length, 1);
    assert.deepEqual(moved.transform, original.transform);
    await drag(column, 24, 40);
    const shifted = await state();
    assert.equal(shifted.columnOffsets.length, 1);
    assert.deepEqual(shifted.nodeOffsets, moved.nodeOffsets);
    const movedPath = await page.locator('.edge').first().getAttribute('d');
    await page.locator('#lockLayout').click();
    await drag(card, 24, 31);
    const panned = await state();
    assert.notDeepEqual(panned.transform, shifted.transform);
    assert.deepEqual(panned.nodeOffsets, shifted.nodeOffsets);
    assert.deepEqual(panned.columnOffsets, shifted.columnOffsets);
    assert.equal(await page.locator('.edge').first().getAttribute('d'), movedPath);
    await page.locator('#depthNext').click();
    assert.equal((await state()).depth, 0);
    await page.locator('#depthNext').click();
    assert.equal((await state()).depth, 1);
    await page.reload();
    assert.equal((await state()).depth, 1);
    assert.deepEqual((await state()).nodeOffsets, moved.nodeOffsets);
    await page.evaluate(() => {
      const first = document.querySelector('g.node').getAttribute('data-node-id');
      window.CRISP_SUPERGRAPH.selectNode(first);
    });
    assert.match(await page.locator('#iConjunction').innerText(), /linked|entering/i);
    assert.deepEqual(errors, []);
    const css = await page.locator('style').first().textContent();
    assert.doesNotMatch(css, /#(?:4ade80|fbbf24|fb7185|f59e0b)/i);
    if (process.env.SUPERGRAPH_QA_OUTPUT) {
      fs.mkdirSync(process.env.SUPERGRAPH_QA_OUTPUT, { recursive: true });
      await page.screenshot({ path: path.join(process.env.SUPERGRAPH_QA_OUTPUT, 'supergraph-depth.png') });
    }
    const touchPage = await browser.newPage({
      viewport: { width: 390, height: 844 }, isMobile: true, hasTouch: true,
    });
    await touchPage.goto(pathToFileURL(entry).href);
    await touchPage.waitForFunction(() => window.CRISP_SUPERGRAPH?.ready);
    const mobileNav = touchPage.getByRole('navigation', { name: 'Evidence categories' });
    await mobileNav.getByRole('button', { name:'Observation mechanisms' }).tap();
    await touchPage.locator('[data-node-id="C_strings"]').tap();
    assert.equal(await touchPage.evaluate(() => window.CRISP_SUPERGRAPH.getSelected().id), 'C_strings');
    assert.equal(await touchPage.evaluate(() => window.CRISP_SUPERGRAPH.getState().rightHidden), false);
    console.log('PASS: desktop drag and depth, persistent layout, palette, evidence inspector and mobile node tap');
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
