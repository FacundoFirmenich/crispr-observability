# Auditoría del ecosistema para Read Observability Lab

## Alcance y método

Revisión dirigida a funciones transferibles, con 25 referentes de seis familias. Se contrastó documentación oficial o artículos de los autores. Para el laboratorio se usaron el código y las pruebas de ejecución. Las capacidades de los comparadores comerciales y de escritorio se puntuaron documentalmente. Esta cobertura permite decisiones de producto y priorización. No es un censo del mercado, una comparación de rendimiento ni una certificación independiente de cada producto.

Se buscaron visores web y de escritorio, herramientas de alineamiento genómico, analizadores de amplicones, editores comerciales de secuencia, cuadernos reactivos y entornos de visualización y reproducibilidad. La unidad de comparación es una capacidad que un investigador puede ejercer con evidencia trazable. Las hipótesis de investigación originales y la ausencia de verdad molecular independiente conservan el mismo alcance.

## Inventario y decisiones

| Familia y referente | Capacidad documentada | Decisión para el laboratorio | Fuente primaria |
| --- | --- | --- | --- |
| Visor web · Mol* | Representaciones, selección, mediciones, varios formatos y snapshots | Incorporar entrada mmCIF y replay. Evaluar delegación a Mol* para ensamblajes, mapas y superposición. | [Mol* Viewer](https://molstar.org/viewer-docs/query-parameters/) |
| Visor web · RCSB PDB 3D | Sincronización de anotaciones entre estructura y secuencia, comparación y vecindad | Preparar un carril sincronizado que conserve la diferencia entre coordenada del read y posición del modelo. | [RCSB Sequence Annotations in 3D](https://www.rcsb.org/docs/sequence-viewers/sequence-annotations-in-3d) |
| Visor web · iCn3D | Estructura, secuencia, conjuntos y estado enlazable | Adoptar inspección con identidad persistente. Mantener los casos sensibles en replay local. | [NCBI iCn3D Help](https://ncbi.nlm.nih.gov/Structure/icn3d/docs/icn3d_help.html) |
| Biblioteca web · NGL | Formatos estructurales, estilos y trayectorias | Evaluar adaptador de lectura o visualización si aparece un caso que requiere trayectoria. | [NGL File Formats](https://nglviewer.org/ngl/api/manual/file-formats.html) |
| Biblioteca web · 3Dmol.js | Visor embebible, estilos y superficies | Usar como alternativa de integración ligera, según el contrato y el coste de mantenimiento. | [3Dmol.js API](https://www.3dmol.org/doc/) |
| Especificación · MolViewSpec | Descripción declarativa de vistas y anotaciones moleculares | Considerar un exportador de vistas. El replay actual también preserva semántica experimental propia. | [MolViewSpec](https://molstar.org/mol-view-spec-docs/) |
| Escritorio · ChimeraX | Selección molecular, alineamiento estructural, render y comandos | Delegar superposición avanzada y preparar anotaciones exportables cuando existan casos con correspondencia definida. | [UCSF ChimeraX guide](https://plato.cgl.ucsf.edu/chimerax/data/cheat_sheet/ChimeraX%20cheat%20sheet%20by%20Alejandro%20Pereira.pdf) |
| Escritorio · PyMOL | Imágenes, animación, scripting y sesiones. Base abierta y distribución Incentive. | Adoptar exportación de detalle de cámara. Evaluar puente de selección o escena, según el uso editorial. | [PyMOL](https://pymol.org/pymol.html) |
| Escritorio · VMD | Trayectorias, medidas de dinámica y plugins | Mantener como destino para trayectorias. Una escena ideal de CFTR carece de trayectoria experimental. | [VMD plugins](https://www.tcbg.illinois.edu/Research/vmd/plugins/) |
| Cristalografía · CCDC Mercury | Estructuras cristalinas, visualización y análisis | Referencia para medidas y navegación espacial. Empaquetamiento cristalino fuera de los datos CFTR actuales. | [CCDC Mercury](https://www.ccdc.cam.ac.uk/solutions/software/mercury/) |
| Genómica · IGV | Cobertura, bases, indels y contexto por alineamiento | Priorizar un carril de lectura con clipping, orientación y procedencia al admitir alineamientos indexados. | [IGV Alignments](https://igv.org/doc/desktop/UserGuide/tracks/alignments/viewing_alignments_basics/) |
| Genómica · JBrowse 2 | Pistas, clipping y sesiones transportables | Adoptar replay. Futuro adaptador de pistas y loci con contratos explícitos de ensamblaje y coordenadas. | [JBrowse session parameters](https://jbrowse.org/jb2/docs/urlparams/) |
| Genómica · UCSC Genome Browser | Pistas y contexto de ensamblaje | Inspiración para vista multiescala y selección de ventanas. | [UCSC User Guide](https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html) |
| Amplicones · CRISPResso2 | Resultados de edición, clasificación y visualizaciones de amplicones | Adoptar importación de tablas con denominadores y nombres de categorías explícitos. | [CRISPResso2 Help](https://crispresso2.pinellolab.org/help) |
| Amplicones · CRISPR-A | Pipeline y aplicación para análisis de edición y simulación | Mantener identidad de lectura y trazabilidad al integrar salidas. | [Sanvicente-García et al.](https://doi.org/10.1371/journal.pcbi.1011137) |
| Edición de secuencias · Benchling | Alineamientos guardados, vistas de secuencia e importación de archivos | Adoptar coordinación de anotaciones y vista. | [Benchling alignments](https://help.benchling.com/hc/en-us/articles/39922093033869-How-to-create-alignments) |
| Edición de secuencias · SnapGene | Vistas Map y Sequence, controles y anotaciones | Adoptar una navegación directa de región a detalle y exportación de anotaciones. | [SnapGene Views](https://support.snapgene.com/hc/en-us/articles/34768892944916-Changes-to-Map-and-Sequence-Views-in-SnapGene-8-0) |
| Edición de secuencias · Geneious Prime | Gestión de secuencias, alineamientos, análisis y workflows | Usar como referencia para formatos y registro de resultados por caso. | [Geneious manual](https://manual.geneious.com/en/latest/Alignments.html) |
| Cuaderno · Wolfram Mathematica | Controles reactivos enlazados a expresiones y gráficos, estado y distribución interactiva | Adoptar recorrido activo, parámetros visibles y presets comparables. Mantener un contrato declarativo de cálculo. | [Wolfram Manipulate](https://reference.wolfram.com/language/ref/Manipulate.html) |
| Cuaderno · Observable | Celdas reactivas con explicación, código y resultados | Adoptar la secuencia pregunta → intervención → resultado esperado → comprobación. | [Observable Notebooks](https://observablehq.github.io/documentation/notebooks/) |
| Cuaderno · Jupyter/ipywidgets | Controles parametrizados y acceso al código en el cuaderno | Usar como destino para extensiones programáticas y casos repetibles. | [Jupyter Widgets](https://ipywidgets.readthedocs.io/en/stable/examples/Widget%20List.html) |
| App de datos · Streamlit | Estado de sesión y componentes de datos | Adoptar controles de estado explícitos. La distribución actual sigue funcionando con HTML autónomo. | [Streamlit Session State](https://docs.streamlit.io/develop/api-reference/caching-and-state/st.session_state) |
| Visualización general · ParaView | Pipeline de filtros, estado y trazas Python | Adoptar manifest de exportación y verificación. Evaluar presets y lote de escenas. | [ParaView Scripting](https://www.paraview.org/scripting/) |
| Reproducibilidad · Galaxy | Historial de datos, herramientas y salidas | Adoptar separación persistente entre entrada, parámetros, resultado y artefacto. | [Galaxy History](https://training.galaxyproject.org/training-material/topics/galaxy-interface/tutorials/history/tutorial.html) |
| Interoperabilidad · MCP TypeScript SDK | Herramientas tipadas y transporte local por stdio | Implementar ejecución por agente en workspace autorizado y salida de render, replay y hashes. | [MCP TypeScript SDK](https://ts.sdk.modelcontextprotocol.io/server) |

## Scoring comparable por capacidad

Escala ordinal 0–5. 0 corresponde a capacidad ausente de la versión local inspeccionada. 1 indica una aproximación básica. 2, una operación útil con cobertura estrecha. 3, un flujo operativo de punta a punta. 4, cobertura amplia o trazabilidad y pruebas del flujo. 5, un referente maduro en esa dimensión. — indica que el flujo exacto no se puntuó a partir de las fuentes revisadas. Las columnas no son conmensurables y no se promedian.

- **3D** mide profundidad molecular, variedad de representaciones y capacidad estructural.
- **Lectura** mide observación por registro, intervalos y procedencia local de secuencias o alineamientos.
- **Intervención** mide posibilidad de variar entradas y ver cambios calculados en resultados.
- **Replay** mide restituir el estado y/o la procedencia de una ejecución.
- **Didáctica** mide ayuda contextual y guías ejecutables dentro del flujo.
- **IDE** mide automatización documentada en un entorno de desarrollo. La columna distingue el tipo de interfaz en observaciones.

| Referente | 3D | Lectura | Intervención | Replay | Didáctica | IDE | Observación |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Mol* | 5 | 2 | 3 | 5 | 3 | 4 | SDK y estado de visualización. |
| iCn3D | 4 | 2 | 3 | 4 | 3 | 3 | Integración web y secuencia de estructura. |
| ChimeraX | 5 | 1 | 4 | 4 | 3 | 4 | Comandos y extensiones de escritorio. |
| IGV / JBrowse 2 | — | 5 | 2 | 4 | 3 | 3 | Pistas, evidencia de alineamiento y sesiones. |
| CRISPResso2 / CRISPR-A | — | 4 | 3 | 3 | 2 | 4 | Pipeline y salidas de cuantificación. |
| Benchling / Geneious Prime | — | 4 | 3 | 4 | 3 | — | Alineamientos y gestión de datos. |
| Wolfram Mathematica | 4 | — | 5 | 4 | 5 | 4 | Autoría programable y cuadernos computables. |
| Observable / Jupyter | 3 | — | 5 | 4 | 4 | 4 | Entorno reactivo y código editable. |
| ParaView | 5 | — | 5 | 5 | 3 | 5 | Filtros, estado y automatización Python. |
| Lab anterior | 2 | 3 | 3 | 1 | 2 | 0 | Cálculo local y exportación de resultado. |
| **Lab actualizado** | **3** | **3** | **4** | **4** | **4** | **3** | Replay y PNG ampliado. Cuatro herramientas MCP stdio. |

Los grupos son lecturas comparativas de los flujos documentados, con diferencias entre productos individuales. El 3 en lectura del laboratorio refleja un panel de registros y secuencias y la ausencia actual de pistas BAM/CRAM. La ventaja defendible es más estrecha: intervención local explícita sobre contención, correspondencia de strings separada de acuerdo de etiquetas, resultado reproducible y explicación vinculada a esos datos. La revisión no establece prioridad mundial ni valida la exactitud de edición sobre la molécula.

## Integración y secuencia de trabajo

1. **Integrado**. Recorrido con comprobaciones vivas, ayuda por control, importación del primer modelo mmCIF, línea de medición 3D, replay local y PNG de alta resolución. El puente MCP usa la aplicación HTML, guarda el PNG y su replay y emite hashes verificables.
2. **Siguiente brecha para la investigación**. Vista de secuencia y coordenadas vinculada al modelo, con un importador acotado de alineamientos y adaptación comprobada del denominador de cada herramienta. El primer control es distinguir coordenada del read, referencia y modelo sin equivalencias implícitas.
3. **Interoperabilidad especializada**. Evaluar Mol* o un exportador MolViewSpec para los casos estructurales exigentes. Evaluar IGV.js/JBrowse para lectura genómica. Reservar el motor propio para observabilidad, pruebas y procedencia.
4. **Publicación reproducible**. Mantener dos salidas separadas, el artefacto renderizado y el estado ejecutable. Añadir lote de renders y una guía de importadores solo tras verificar semánticas sobre casos externos.

Esta auditoría inspecciona funciones de visualización, análisis, reproducibilidad e integración. Las piezas propietarias inspiran requisitos e interfaces. Su código, recursos y formatos no documentados conservan sus licencias y límites.
