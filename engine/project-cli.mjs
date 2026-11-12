import {dirname, basename, resolve} from 'node:path';
import {z} from 'zod';
import {projectSchema, readManifest, runProject, validateProject} from './project-engine.mjs';
try{
  const [command, requested, ...extra]=process.argv.slice(2);
  if(extra.length)throw new Error('Unexpected arguments');
  if(command==='schema'&&!requested){
    process.stdout.write(JSON.stringify(z.toJSONSchema(projectSchema),null,2)+'\n');
  }else{
    if(!['check','run'].includes(command)||!requested)throw new Error('Usage: node project-cli.mjs check|run manifest.json; or schema');
    const full=resolve(requested),root=dirname(full);
    const input=await readManifest(root,basename(full));
    const out=command==='check'?{valid:true,project:validateProject(input)}:await runProject(input,root);
    process.stdout.write(JSON.stringify(out,null,2)+'\n');
  }
}catch(error){console.error(error.message);process.exitCode=1;}
