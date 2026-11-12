import {readFile,writeFile} from 'node:fs/promises';
import {build} from 'esbuild';

const bundle=await build({entryPoints:['lab.js'],bundle:true,write:false,format:'iife',minify:true,
  legalComments:'inline',target:['chrome109','firefox115','safari17']});
const html=await readFile('index.template.html','utf8');
const data=await readFile('data.json','utf8');
const code=bundle.outputFiles[0].text;
await writeFile('CFTR_OBSERVABILITY_LAB.html',html.replace('__LAB_DATA__',()=>data.replaceAll('<','\\u003c'))
  .replace('__LAB_SCRIPT__',()=>code.replace(/<\/script/gi,'<\\/script')),'utf8');
