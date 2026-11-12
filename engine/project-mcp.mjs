import {McpServer} from '@modelcontextprotocol/sdk/server/mcp.js';
import {StdioServerTransport} from '@modelcontextprotocol/sdk/server/stdio.js';
import {realpath,stat} from 'node:fs/promises';
import {dirname,resolve} from 'node:path';
import {z} from 'zod';
import {readManifest,runProject,validateProject} from './project-engine.mjs';
const [requested,...extra]=process.argv.slice(2);
if(!requested||extra.length)throw new Error('Usage: node project-mcp.mjs authorized-workspace');
const workspace=await realpath(requested);
if(!(await stat(workspace)).isDirectory())throw new Error('Workspace must be a directory');
const server=new McpServer({name:'read-observability-projects',version:'0.4.0'});
server.registerTool('lab_project_run',{
  description:'Validate or run a multi-input project manifest using files under the authorized workspace.',
  inputSchema:{manifest_path:z.string(),check_only:z.boolean().default(false)}
},async args=>{
  try{
    const input=await readManifest(workspace,args.manifest_path);
    const result=args.check_only?{valid:true,project:validateProject(input)}:await runProject(input,dirname(resolve(workspace,args.manifest_path)));
    return {content:[{type:'text',text:JSON.stringify(result)}],structuredContent:result};
  }catch(error){return {isError:true,content:[{type:'text',text:error.message}]};}
});
await server.connect(new StdioServerTransport());
