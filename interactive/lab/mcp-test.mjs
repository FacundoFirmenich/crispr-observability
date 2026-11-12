import assert from 'node:assert/strict';
import {Client} from '@modelcontextprotocol/sdk/client/index.js';
import {StdioClientTransport} from '@modelcontextprotocol/sdk/client/stdio.js';
import {createHash} from 'node:crypto';
import {readFile} from 'node:fs/promises';
import {join} from 'node:path';

const workspace=process.cwd();
const client=new Client({name:'lab-integration-test',version:'0.2.0'});
const transport=new StdioClientTransport({command:process.execPath,
  args:[join(workspace,'mcp-server.mjs'),workspace],stderr:'inherit'});
await client.connect(transport);
try{
  const tools=(await client.listTools()).tools.map(t=>t.name).sort();
  assert.deepEqual(tools,['lab_alignment_locus','lab_boundary','lab_compare_spans','lab_normalize_case','lab_render']);
  async function call(name,args){
    const result=await client.callTool({name,arguments:args},undefined,{timeout:95000});
    assert.ok(!result.isError,result.content?.[0]?.text);
    return result.structuredContent;
  }
  const containment=await call('lab_boundary',{start:210,end:222,shift:3,read_end:224});
  assert.equal(containment.result.observed,11);
  const comparison=await call('lab_compare_spans',{archive_read:'SRR29623271.868',a:[114,186],b:[114,153]});
  assert.equal(comparison.result.onlyA,33);
  const normalized=await call('lab_normalize_case',{case_path:'qa-1BNA.cif'});
  assert.equal(normalized.format,'mmCIF');
  assert.ok(normalized.atoms>480);
  const text=(await readFile('qa-1BNA.cif','utf8'));
  const inline=await call('lab_normalize_case',{case_name:'same.cif',case_text:text});
  assert.equal(inline.sha256,normalized.sha256);
  const rendered=await call('lab_render',{scale:4,replay:{schema:'read-observability-result/1',
    kind:'ideal_duplex_inspection',inputs:{start:129,end:141,strand:'both',representation:'bond'}}});
  assert.equal(rendered.result.visibleAtoms,492);
  assert.ok(rendered.width>=3500&&rendered.height>=1500);
  for(const [name,expected] of Object.entries(rendered.files)){
    const bytes=await readFile(join(workspace,rendered.folder,name));
    assert.equal(createHash('sha256').update(bytes).digest('hex'),expected);
    if(name==='scene.png'){
      assert.equal(bytes.readUInt32BE(16),rendered.width);
      assert.equal(bytes.readUInt32BE(20),rendered.height);
    }
  }
  const rejected=await client.callTool({name:'lab_normalize_case',arguments:{case_path:'../CRISPR_ATOMIC_DNA_REVIEW_20261112/1BNA.pdb'}});
  assert.equal(rejected.isError,true);
  assert.match(rejected.content[0].text,/authorized workspace/);
  const audit=await call('lab_alignment_locus',{bam_path:'1000genomes_hg00096_chr1.bam',
    index_path:'1000genomes_hg00096_chr1.bam.bai',reference:'1',locus_start:9999,locus_end:10066,
    query_index:0,original_orientation:true});
  assert.equal(audit.coverage.covered,67);
  assert.equal(audit.position.reference,10065);
  const replay=JSON.parse(await readFile(join(workspace,audit.folder,'alignment-replay.json'),'utf8'));
  const evidence=await call('lab_render',{scale:1,replay});
  assert.equal(evidence.result.covered,67);
  assert.ok(evidence.files['alignment-evidence.tsv']);
  const core=await readFile(join(workspace,audit.folder,'alignment-evidence.tsv'));
  assert.equal(createHash('sha256').update(core).digest('hex'),evidence.files['alignment-evidence.tsv']);
  console.log(JSON.stringify({tools,atomSiteCount:normalized.atoms,inlineInput:true,replay:true,
    workspaceRestriction:true,rendered,alignment:{audit,evidence,sameTsv:true}},null,2));
}finally{await client.close();}
