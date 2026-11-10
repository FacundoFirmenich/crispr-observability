"use strict";
function selectNodes(graph, query, kind, status) {
 const q=query.trim().toLowerCase();
 return graph.nodes.filter(n=>(!kind||n.kind===kind)&&(!status||n.status===status)&&
  (!q||[n.id,n.label,n.scope,n.details,n.falsifier,...n.artifacts].join(" ").toLowerCase().includes(q)));
}
function neighbors(graph,id) {
 const ids=new Set([id]); for(const e of graph.edges) if(e.source===id||e.target===id){ids.add(e.source);ids.add(e.target);} return [...ids];
}
function summary(graph) {return {nodes:graph.nodes.length,edges:graph.edges.length,claims:graph.nodes.filter(n=>n.kind==="claim").length};}
if(typeof module!=="undefined")module.exports={selectNodes,neighbors,summary};
