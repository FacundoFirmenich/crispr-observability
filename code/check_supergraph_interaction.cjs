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
    const node = page.locator('g.node').filter({ hasText: 'Conditional' }).first();
    const card = await node.count() ? node : page.locator('g.node').first();
    const column = page.locator('.column-handle').nth(1);
    const drag = async (locator, dx, dy) => {
      const box = await locator.boundingBox();
      assert.ok(box);
      const x = box.x + box.width / 2, y = box.y + box.height / 2;
      await page.mouse.move(x, y);
      await page.mouse.down();
      await page.mouse.move(x + dx, y + dy, { steps: 8 });
      await page.mouse.up();
      await page.waitForFunction(() => !document.getElementById('graph').classList.contains('dragging'));
    };
    const original = await state();
    await drag(card, 40, 30);
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
    console.log('PASS: node and column drag, locked pan, depth, persistence, evidence inspector and palette');
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
