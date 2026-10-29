"""Build the public CRISPR PaperGraph from its authored, source-bound model.

SVG coordinates are diagram layout, not measurements or evidential strength.
The builder has no external browser/CDN dependencies. Pillow
is used only to wrap labels using the same DejaVu Sans font used by Inkscape.
"""
from pathlib import Path
from html import escape
import argparse, hashlib, json, re, shutil, subprocess, os
from PIL import ImageFont, Image

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parents[1]
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
BOLD = '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'
COLORS = {'question':'#23495a', 'contract':'#765413', 'mechanism':'#006d72',
          'formal':'#315c93', 'result':'#006d72', 'boundary':'#9a3651',
          'evidence':'#006d72', 'method':'#315c93', 'consequence':'#9a3651'}

def sha(data): return hashlib.sha256(data).hexdigest()
def text_lines(text, width, size, bold=False):
    font = ImageFont.truetype(BOLD if bold else FONT, size)
    lines, current = [], ''
    for word in str(text).split():
        while font.getlength(word) > width:
            if current:lines.append(current);current=''
            chunk=''
            for char in word:
                if chunk and font.getlength(chunk+char)>width:break
                chunk+=char
            lines.append(chunk);word=word[len(chunk):]
        candidate = (current+' '+word).strip()
        if current and font.getlength(candidate) > width:
            lines.append(current); current=word
        else: current=candidate
    if current: lines.append(current)
    return lines

def xml_text(lines,x,y,size,fill='#1b3440',bold=False,lineheight=None):
    lineheight=lineheight or size*1.26
    return '<text font-family="DejaVu Sans" font-size="'+str(size)+'" fill="'+fill+'"'+(' font-weight="700"' if bold else '')+'>'+''.join('<tspan x="'+str(x)+'" y="'+str(round(y+i*lineheight,2))+'">'+escape(str(t))+'</tspan>'for i,t in enumerate(lines))+'</text>'

def role_color(station):
    r=str(station.get('role','')).lower()
    for k,c in COLORS.items():
        if k in r:return c
    if any(x in r for x in ['limit','frontier','scope','negative']): return COLORS['boundary']
    if any(x in r for x in ['derive','bound','proof','math']): return COLORS['formal']
    return COLORS['evidence']

def argument_svg(model):
    """Complete zoomable audit atlas. Layout groups roles, not causal ranks."""
    nodes=model['nodes'];edges=model['edges'];facts={f['id']:f for f in model['facts']}
    stations={s['id']:s for s in model['stations']}
    fact_nodes=[n for n in nodes if n['type']=='evidence_fact']
    groups=[('Evidence facts 1',fact_nodes[0::3],50,370),
            ('Evidence facts 2',fact_nodes[1::3],480,370),
            ('Evidence facts 3',fact_nodes[2::3],910,370),
            ('Declared premises',[n for n in nodes if n['type']=='assumption'],1340,410),
            ('Formal results',[n for n in nodes if n['type']in{'formal_statement','formal_derivation','refuted_inference'}and n['id']not in stations],1810,410),
            ('Joint requirements',[n for n in nodes if n['type']=='conjunction'],2280,360),
            ('Article stations',[n for n in nodes if n['id']in stations],2700,520),
            ('Unestablished claims',[n for n in nodes if n['type']=='limit'],3280,460)]
    included=[n['id']for _,ns,_,_ in groups for n in ns]
    assert set(included)=={n['id']for n in nodes} and len(included)==len(nodes)
    positions={};W=3800
    for title,ns,x,w in groups:
        y=270
        for n in ns:
            name=stations.get(n['id'],{}).get('shortTitle')or n.get('title')or n.get('label')or n['id']
            if n['type']=='refuted_inference':name='REFUTED INFERENCE: '+name
            lines=text_lines(name,w-32,22,True)
            scalar=[]
            if n['id']in facts:
                f=facts[n['id']];value=f['value']
                if not isinstance(value,(dict,list)):
                    value=f'{value:,}'if isinstance(value,int)else str(value)
                    scalar=text_lines(value+(' '+f['unit']if f.get('unit')else''),w-32,20)
                else:scalar=text_lines('Structured data: inspect model',w-32,20)
            header=text_lines(n['id'],w-32,17,True)
            h=53+len(header)*22+len(lines)*28+len(scalar)*25+15
            positions[n['id']]={'x':x,'y':y,'w':w,'h':h,'title':lines,'header':header,'value':scalar,'node':n}
            y+=h+26
    H=max(p['y']+p['h']for p in positions.values())+150
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-labelledby="argument-title argument-description">',
         '<title id="argument-title">PaperGraph: complete typed argument atlas</title>',
         '<desc id="argument-description">Complete source-bound model: all nodes and all typed edges, including a separately marked pedagogical reading route. Role columns do not encode causal rank, truth, confidence or empirical strength. Edge titles contain their declared semantics and conditions. Node titles and embedded metadata retain full objects and facts.</desc>',
         '<metadata id="papergraph-model">'+escape(json.dumps(model,ensure_ascii=False,separators=(',',':')))+'</metadata>',
         f'<rect width="{W}" height="{H}" fill="#fbfcfa"/>',
         xml_text(['PAPERGRAPH / COMPLETE TYPED ARGUMENT ATLAS'],50,58,34,'#183440',True),
         xml_text([f'{len(nodes)} source-bound nodes · {len(edges)} declared links · {len(model["conjunctions"])} AND groups · no added observations'],50,103,25),
         xml_text(['This is a zoomable audit figure, not an A4 summary. Role columns are layout only. Isolated nodes remain isolated when no semantic edge is declared.'],50,145,24,'#536c76')]
    styles={'requires':('#765413','4 5',2),'joint_support':('#765413','',4),
            'reading_order':('#99a9af','2 7',1.5),'does_not_establish':('#9a3651','12 6',2.5),
            'qualifies':('#9a3651','12 6',2.5),'derives_under':('#315c93','',2.5),
            'defines':('#536c76','2 5',2),'motivates':('#536c76','2 5',2),
            'supports_within':('#006d72','',2.5),
            'counterexample_to_interchangeability':('#9a3651','12 6',2.5),
            'distinguishes_scope_from':('#9a3651','12 6',2.5)}
    for i,(rel,(c,dash,sw))in enumerate(styles.items()):
        x=50+(i%6)*615;y=189+(i//6)*35
        out.extend([f'<line x1="{x}" x2="{x+42}" y1="{y-7}" y2="{y-7}" stroke="{c}" stroke-width="{sw}" stroke-dasharray="{dash}"/>',xml_text([rel.replace('_',' ')],x+53,y,18,c)])
    out.append('<defs>')
    for rel,(c,_,_)in styles.items():out.append(f'<marker id="atlas-{rel}" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M1 1L7 4L1 7Z" fill="{c}"/></marker>')
    out.append('</defs>')
    for e in edges:
        a=positions[e['source']];b=positions[e['target']];rel=e['relation'];c,dash,sw=styles[rel]
        x1=a['x']+a['w'];x2=b['x'];y1=a['y']+a['h']/2;y2=b['y']+b['h']/2
        if a['x']==b['x']:
            x1=a['x']+a['w'];x2=b['x']+b['w'];control=min(W-10,x1+45)
            d=f'M{x1},{y1}C{control},{y1} {control},{y2} {x2},{y2}'
        elif x1>x2:
            x1=a['x'];x2=b['x']+b['w'];mid=(x1+x2)/2;d=f'M{x1},{y1}C{mid},{y1} {mid},{y2} {x2},{y2}'
        else:
            mid=(x1+x2)/2;d=f'M{x1},{y1}C{mid},{y1} {mid},{y2} {x2},{y2}'
        out.append(f'<g data-edge-id="{escape(e["id"],quote=True)}" data-relation="{rel}"><title>'+escape(json.dumps(e,ensure_ascii=False))+'</title>'+f'<path d="{d}" fill="none" stroke="{c}" stroke-width="{sw}" stroke-dasharray="{dash}" opacity="{.4 if rel=="reading_order"else .7}" marker-end="url(#atlas-{rel})"/></g>')
    for title,_,x,w in groups:out.append(xml_text(text_lines(title,w,25,True),x,253,25,'#183440',True))
    for n in nodes:
        p=positions[n['id']];typ=n['type'];c='#765413'if typ in{'assumption','conjunction'}else'#9a3651'if typ in{'limit','refuted_inference'}else'#315c93'if typ in{'formal_statement','formal_derivation'}else'#006d72'
        out.extend([f'<g id="node-{escape(n["id"],quote=True)}" data-node-id="{escape(n["id"],quote=True)}"><title>'+escape(json.dumps({'node':n,'fact':facts.get(n['id'])},ensure_ascii=False))+'</title>',
          f'<rect x="{p["x"]}" y="{p["y"]}" width="{p["w"]}" height="{p["h"]}" rx="8" fill="#fff" stroke="{c}" stroke-width="2"/>',
          xml_text(p['header'],p['x']+16,p['y']+28,17,c,True,lineheight=22),
          xml_text(p['title'],p['x']+16,p['y']+48+len(p['header'])*22,22,'#183440',True,lineheight=28),
          xml_text(p['value'],p['x']+16,p['y']+48+len(p['header'])*22+len(p['title'])*28+8,20,c,lineheight=25),'</g>'])
    out.extend([xml_text(['Grey dotted reading-order links carry no inference. Dashed premise links must be taken jointly at their named AND group.'],50,H-85,24,'#536c76'),
         xml_text(['Full source identities, exact quantities, conditions and objects are retained in the embedded JSON metadata and in the companion model.'],50,H-45,24,'#536c76'),'</svg>'])
    return '\n'.join(out)

def station_svg(model, orientation='landscape'):
    portrait=orientation=='portrait'
    W,H=(1600,2660) if portrait else (1760,1320)
    cols=2 if portrait else 4
    margin=50; gap=36 if portrait else 30
    cardW=(W-2*margin-(cols-1)*gap)/cols
    cardH=340 if portrait else 245
    y0=210
    rowgap=12 if portrait else 30
    stations=model['stations']; positions={}
    out=[f'<svg id="papergraph-svg" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-labelledby="pg-title pg-desc">',
         '<title id="pg-title">CRISPR PaperGraph: the argument, evidence and limits</title>',
         '<desc id="pg-desc">Twelve numbered stations give the reading order of the article. Dotted connectors indicate reading order only. The linked model and interactive version retain typed scientific relations, source identities, assumptions and revision conditions. The diagram does not identify molecular accuracy or biological function.</desc>',
         f'<rect width="{W}" height="{H}" fill="#fbfcfa"/>',
         f'<rect x="0" y="0" width="{W}" height="9" fill="#006d72"/>',
         xml_text(['PAPERGRAPH / CRISPR OBSERVABILITY'],margin,53,32 if portrait else 17,'#006d72',True),
         xml_text(['Follow what is observed. Keep every inference conditional.'],margin,105,42 if portrait else 36,bold=True),
         xml_text(['Question → contract → evidence → conditional conclusion → frontier.'],margin,148,32 if portrait else 19),
         xml_text(['Read stations 01–12. Dotted arrows mean “read next”, not “proves”.'],margin,190 if portrait else 184,32 if portrait else 18,'#526873')]
    rows=(len(stations)+cols-1)//cols
    ordered=[]
    for i,st in enumerate(stations):
        row=i//cols; cc=i%cols; col=(cols-1-cc) if row%2 else cc
        x=margin+col*(cardW+gap);y=y0+row*(cardH+rowgap)
        positions[st['id']]={'x':x,'y':y,'w':cardW,'h':cardH}
        ordered.append((x,y))
    out.append('<defs><marker id="reading-arrow" markerWidth="8" markerHeight="8" refX="6" refY="4" orient="auto"><path d="M1 1L6 4L1 7" fill="none" stroke="#9dacb1" stroke-width="1.5"/></marker></defs>')
    for i in range(len(ordered)-1):
        x,y=ordered[i];xx,yy=ordered[i+1]
        if yy==y:
            right=xx>x;a=x+cardW if right else x;b=xx if right else xx+cardW
            d=f'M{a},{y+cardH/2}L{b},{yy+cardH/2}'
        else:
            d=f'M{x+cardW/2},{y+cardH}L{xx+cardW/2},{yy}'
        out.append(f'<path d="{d}" fill="none" stroke="#9dacb1" stroke-width="1.8" stroke-dasharray="4 5" marker-end="url(#reading-arrow)"/>')
    for i,st in enumerate(stations):
        p=positions[st['id']];x,y=p['x'],p['y'];c=role_color(st)
        title=st.get('printTitle',st.get('shortTitle',st['title']));take=st.get('printSummary',st['takeaway']) if portrait else st['takeaway']
        titleSize=40 if portrait else 24; bodySize=33 if portrait else 19
        titleLines=text_lines(title,cardW-54,titleSize,True)
        while len(titleLines)>2 and titleSize>20:
            titleSize-=1;titleLines=text_lines(title,cardW-54,titleSize,True)
        bodyLines=text_lines(take,cardW-54,bodySize)
        maxBody=2 if portrait else 4
        while len(bodyLines)>maxBody and bodySize>(32 if portrait else 16):
            bodySize-=1;bodyLines=text_lines(take,cardW-54,bodySize)
        if len(bodyLines)>maxBody: raise ValueError(f'Takeaway too long for complete print rendering: {st["id"]}')
        role={'adverse_and_control':'INTERPRETIVE CASES'}.get(st.get('role'),str(st.get('role','Station')).replace('_',' ').upper())
        roleLines=text_lines(role,cardW-110,32 if portrait else 14,True)
        evidenceLines=text_lines(st.get('printEvidence',''),cardW-54,32,False)
        sourceLines=text_lines(st.get('printSource','Source anchors in the companion'),cardW-54,32)
        if portrait and (len(evidenceLines)>2 or len(sourceLines)>1): raise ValueError('Evidence/source badge exceeds complete print geometry: '+st['id'])
        out.extend([f'<g class="station" data-station="{escape(st["id"],quote=True)}" role="button" tabindex="0" aria-label="{escape(str(i+1)+". "+st["title"],quote=True)}">',
            f'<rect class="station-box" x="{x}" y="{y}" width="{cardW}" height="{cardH}" rx="11" fill="#ffffff" stroke="#c9d5d8" stroke-width="1.5"/>',
            f'<rect x="{x}" y="{y}" width="6" height="{cardH}" rx="3" fill="{c}"/>',
            xml_text([f'{i+1:02d}'],x+23,y+36 if portrait else y+31,32 if portrait else 16,c,True),xml_text(roleLines[:1],x+89 if portrait else x+63,y+36 if portrait else y+31,32 if portrait else 13,c,True),
            xml_text(titleLines,x+25,y+82 if portrait else y+72,titleSize,bold=True,lineheight=titleSize*1.15),
            xml_text(bodyLines,x+25,y+175 if portrait else y+124,bodySize,lineheight=bodySize*(1.15 if portrait else 1.25)),
            xml_text(evidenceLines,x+25,y+254,32,c,False,lineheight=36.8) if portrait else '',
            xml_text(sourceLines if portrait else [st.get('printSource','Source anchors in the companion')],x+25,y+cardH-10 if portrait else y+cardH-17,32 if portrait else 12,'#526873'),'</g>'])
    baseY=y0+rows*(cardH+rowgap)-rowgap+(32 if portrait else 46)
    if portrait:
        out.extend([f'<line x1="{margin}" y1="{baseY}" x2="{W-margin}" y2="{baseY}" stroke="#c9d5d8"/>',
          xml_text(['A FORMAL DEPENDENCY, NOT A CAUSAL ARROW'],margin,baseY+36,32,'#315c93',True)])
        by=baseY+60
        out.extend([f'<rect x="50" y="{by}" width="650" height="130" rx="8" fill="#faf3e6" stroke="#d5c29e"/>',
          xml_text(['Nonnegative integer partition', 'E + NE + U = N; N > 0', 'Correct labels; unknown labels free'],70,by+35,32,'#765413',lineheight=36),
          f'<path d="M700,{by+65}L850,{by+65}" stroke="#315c93" stroke-width="3" marker-end="url(#reading-arrow)"/>',
          xml_text(['AND'],741,by+47,32,'#315c93',True),
          f'<rect x="865" y="{by}" width="685" height="130" rx="8" fill="#edf1f8" stroke="#b6c8df"/>',
          xml_text(['Attainable: p = (E+k)/N', 'integer k = 0,…,U', 'Hull: [E/N, (E+U)/N]'],890,by+35,32,'#315c93',lineheight=36),
          xml_text(['Conditional theorem: not sampling coverage or molecular truth.'],margin,by+178,32,'#9a3651'),
          xml_text(['Private article-reading guide · Facundo Firmenich'],margin,H-25,32,'#526873'),'</svg>'])
        return '\n'.join(out),positions
    out.extend([f'<line x1="{margin}" y1="{baseY}" x2="{W-margin}" y2="{baseY}" stroke="#c9d5d8"/>',
        xml_text(['SCIENTIFIC RELATIONS LIVE IN THE ARGUMENT MODEL'],margin,baseY+36,16,'#315c93',True),
        xml_text(['Requires jointly: stated assumptions. Derived under: the declared finite-record contract.'],margin,baseY+70,18),
        xml_text(['Supported by: an identified source and operation. Limited by: an unmeasured quantity or changed domain.'],margin,baseY+101,18),
        xml_text(['The overview is a reading map; open a station to inspect these typed links and the evidence they reference.'],margin,baseY+132,17,'#526873'),
        xml_text(['Conditional identification hull: [E/N, (E+U)/N]  ·  Accepted labels must be correct.'],margin,baseY+176,19,'#765413',True),
        xml_text(['No molecular-accuracy, prospective-validation or haplotype-function claim is supplied by this diagram.'],margin,baseY+209,17,'#9a3651'),
        xml_text(['Public diagram · formal text withheld · Facundo Firmenich · source-bound model '+sha(json.dumps(model,sort_keys=True).encode())[:16]],margin,H-25,12,'#526873'),'</svg>'])
    return '\n'.join(out),positions

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--model',type=Path,default=ROOT/'papergraph.json');ap.add_argument('--output',type=Path,default=ROOT);ap.add_argument('--render',action='store_true');ap.add_argument('--figures-output',type=Path,default=REPOSITORY/'figures/papergraph');a=ap.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);fig=a.figures_output;fig.mkdir(parents=True,exist_ok=True)
    modelBytes=a.model.read_bytes();model=json.loads(modelBytes)
    assert len(model['stations'])==12
    assert model['readingPath']==[s['id']for s in model['stations']]
    for st in model['stations']:
        assert all(st['levels'].get(d) for d in ['student','technical','audit'])
    land,positions=station_svg(model,'landscape');port,_=station_svg(model,'portrait');atlas=argument_svg(model)
    (fig/'PaperGraph_Overview.svg').write_text(land);(fig/'PaperGraph_Article.svg').write_text(port);(fig/'PaperGraph_Argument.svg').write_text(atlas)
    raw=json.dumps(model,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c')
    template=(ROOT/'PaperGraph.template.html').read_text()
    template=template.replace('__ARGUMENT_RUNTIME__',(ROOT/'code/argument_runtime.js').read_text())
    html=template.replace('__PAPERGRAPH_DATA__',raw).replace('__MODEL_SHA256__',sha(modelBytes)).replace('__EMBEDDED_SHA256__',sha(raw.encode())).replace('__OVERVIEW_SVG__',land).replace('__PORTRAIT_SVG_JSON__',json.dumps(port).replace('<','\\u003c'))
    assert '__PAPERGRAPH_DATA__'not in html
    assert not re.search(r'<(?:script|link)[^>]+(?:src|href)=["\']https?://',html,re.I)
    target=a.output/'PaperGraph_CRISPR.html';target.write_text(html)
    (a.output/'papergraph.json').write_bytes(modelBytes)
    # Validate declared public sources in place. No source-copy step exists.
    for source in model['sources']:
        if source['access']=='WITHHELD_PRIVATE_SOURCE_IDENTITY':
            assert 'path' not in source
            continue
        rel=Path(source['path']);assert not rel.is_absolute() and '..' not in rel.parts
        sourceFile=REPOSITORY/rel
        assert sourceFile.is_file(),source['id']
        rawSource=sourceFile.read_bytes()
        assert sha(rawSource)==source['sha256'] and len(rawSource)==source['bytes'],source['id']
    receipt={'html_sha256':sha(html.encode()),'model_sha256':sha(modelBytes),'embedded_sha256':sha(raw.encode()),'stations':len(model['stations']),'nodes':len(model.get('nodes',[])),'edges':len(model.get('edges',[])),'svg_overview_sha256':sha(land.encode()),'svg_article_sha256':sha(port.encode()),'svg_argument_sha256':sha(atlas.encode()),'argument_reading_order_links':sum(e['relation']=='reading_order'for e in model['edges']),'argument_semantic_links':sum(e['relation']!='reading_order'for e in model['edges']),'rendered':False,'scope':'Public scientific diagram; full manuscript and supplement withheld. Preserves declared source data and distinguishes reading order from scientific inference.'}
    if a.render:
        exe=shutil.which('inkscape');assert exe,'Inkscape is required for requested vector rendering'
        for name in ['PaperGraph_Overview','PaperGraph_Article','PaperGraph_Argument']:
            png=fig/(name+'.rendering.png');pdf=fig/(name+'.rendering.pdf')
            subprocess.run([exe,str(fig/(name+'.svg')),'--export-type=png','--export-filename='+str(png)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            with Image.open(png)as image:image.load()
            os.replace(png,fig/(name+'.png'))
            subprocess.run([exe,str(fig/(name+'.svg')),'--export-type=pdf','--export-filename='+str(pdf)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            data=pdf.read_bytes();assert data.startswith(b'%PDF-')and b'%%EOF'in data[-100:]
            os.replace(pdf,fig/(name+'.pdf'))
        receipt['rendered']=True
    (a.output/'BUILD_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
