import assert from 'node:assert/strict';
import {readFile,writeFile,access} from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {gzipSync} from 'node:zlib';
import {pathToFileURL} from 'node:url';
import {resolve} from 'node:path';
import {chromium,devices} from 'playwright';

const here=process.cwd();
const url=process.argv[2]??pathToFileURL(resolve(here,'CFTR_OBSERVABILITY_LAB.html')).href;
const local=url.startsWith('file:');
const browser=await chromium.launch({channel:'chrome',headless:true,
  args:['--enable-unsafe-swiftshader','--disable-background-networking','--disk-cache-size=1048576']});
const has=(page,selector,value)=>page.waitForFunction(([s,v])=>document.querySelector(s).textContent.includes(v),[selector,value]);
try{
  const context=await browser.newContext({viewport:{width:1440,height:980},acceptDownloads:true,offline:local});
  const page=await context.newPage();
  const errors=[];
  page.on('pageerror',e=>errors.push(String(e)));
  await page.goto(url);await has(page,'#metrics','11 / 12');
  const network=[];
  page.on('request',r=>{if(/^https?:/.test(r.url()))network.push(r.url());});
  const control=id=>page.locator('#'+id);
  await page.locator('[data-mode=align]').click();
  await control('alignmentDemo').click();await has(page,'#metrics','67 / 67');
  await has(page,'#align-query-note','no mapped reference base');
  await control('originalRead').check();await has(page,'#align-query-note','reference index 10065');
  await control('minMapq').fill('20');await has(page,'#metrics','0 / 67');
  await control('minMapq').fill('0');await has(page,'#metrics','67 / 67');
  await control('loadCase').setInputFiles(resolve(here,'1000genomes_hg00096_chr1.sam'));
  await has(page,'#import-status','SAM');
  assert.equal((await page.evaluate(()=>window.readObservabilityLab.snapshot())).result.read.cigar,'183S67M');
  await control('loadCase').setInputFiles(resolve(here,'1000genomes_hg00096_chr1.bam'));
  await control('bamIndex').setInputFiles(resolve(here,'1000genomes_hg00096_chr1.bam.bai'));
  await has(page,'#import-status','header and index');
  await control('bamRef').selectOption('1');
  await control('bamStart').fill('9999');await control('bamEnd').fill('10066');
  await control('bamLoad').click();await has(page,'#import-status','1 mapped records');
  await control('queryFocus').fill('183');await has(page,'#align-query-note','reference index 9999');
  const indexed=await page.evaluate(()=>window.readObservabilityLab.snapshot());
  assert.ok(indexed.portableInput.source.bamSHA256&&indexed.portableInput.source.indexSHA256);
  const [download]=await Promise.all([page.waitForEvent('download'),control('exportTsv').click()]);
  const tsv=await readFile(await download.path(),'utf8');
  assert.ok(tsv.includes('coverage\t1\t9999\t10066')&&tsv.includes('183S'));
  const replay=await page.evaluate(v=>window.readObservabilityLab.replay(v),indexed);
  assert.deepEqual(replay.result,indexed.result);
  assert.deepEqual(network,[]);
  await page.screenshot({path:'lab-v3-alignment.png',fullPage:true});
  const mobile=await browser.newPage({...devices['Pixel 7'],offline:local});
  await mobile.goto(url);await mobile.locator('[data-mode=align]').click();
  await mobile.locator('#alignmentDemo').click();await has(mobile,'#metrics','67 / 67');
  assert.equal(await mobile.evaluate(()=>document.documentElement.scrollWidth>innerWidth),false);
  const proof={target:local?'offline HTML':'public HTTPS',builtInDemo:true,sam:true,indexedBam:true,
    reverseCoordinate:true,unmappedClip:true,coverageFilter:true,exactTsv:true,replay:true,noUploadRequests:true,
    mobileNoOverflow:true,errors};
  assert.deepEqual(errors,[]);
  await writeFile(local?'BROWSER_V3_CHECKS.json':'BROWSER_V3_HTTPS_CHECKS.json',JSON.stringify(proof,null,2)+'\n');
  console.log(JSON.stringify(proof));
  await context.close();await mobile.close();

  if(local&&await access('V2_BASELINE.html').then(()=>true,()=>false)){
    const paths=['V2_BASELINE.html','CFTR_OBSERVABILITY_LAB.html'];
    const values=paths.map(()=>[]);
    for(let k=0;k<7;k++)for(const i of (k%2?[1,0]:[0,1])){
      const c=await browser.newContext({viewport:{width:1440,height:980},offline:true});
      const p=await c.newPage();
      await p.goto(pathToFileURL(resolve(here,paths[i])).href);
      await has(p,'#metrics','11 / 12');
      values[i].push(Math.round(await p.evaluate(()=>performance.now())));
      await c.close();
    }
    const comparison=await Promise.all(paths.map(async (file,i)=>{
      const bytes=await readFile(file);
      const times=values[i],sorted=[...times].sort((a,b)=>a-b);
      return {version:i===0?'0.2':'0.3',bytes:bytes.length,gzipBytes:gzipSync(bytes).length,
        sha256:createHash('sha256').update(bytes).digest('hex'),timingsMs:times,medianMs:sorted[3]};
    }));
    const result={metric:'Navigation start through first anchor metrics. Same Chromium process, new offline context per trial.',
      scope:'One host, seven interleaved runs per version. No claim about another tool or production hardware.',comparison};
    await writeFile('PORTABILITY_V3_CHECKS.json',JSON.stringify(result,null,2)+'\n');
    console.log(JSON.stringify(result));
  }
}finally{await browser.close();}
