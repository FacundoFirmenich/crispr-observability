# Compañero computacional público

Versión **0.4.0 / MP04**, 29 de octubre de 2026.

El repositorio contiene código, datos numéricos, trazabilidad y figuras editables.
El manuscrito completo y su suplemento se conservan de forma privada hasta el
depósito que decida el autor.

Ejecutar `python3 code/run_review_checks.py` desde la raíz para reproducir las
comprobaciones declaradas. Las figuras están en `figures/`; los diecinueve CSV
en `tables/`. La reejecución completa de alineamientos se explica por separado
en `engine_replay/README.md`.

Descargar `interactive/epistemic_graph.html` y abrirlo en el navegador permite
explorar el mapa de evidencias, las dependencias conjuntivas, los conteos y cotas,
y la sensibilidad histórica. Funciona sin conexión y sin dependencias externas.
Sus controles modifican la vista o completan desconocidos hipotéticamente; no
añaden datos biológicos. La reconstrucción y sus límites están documentados en
`interactive/README.md`.

Los agradecimientos vigentes están en `ACKNOWLEDGMENTS.md`. Las afiliaciones y
los registros ORCID mantienen sus categorías originales en
`metadata/AUTHOR_PROFILE_VERIFIED.json`.

La segunda herramienta, `interactive/supergraph/crispr_supergraph.html`, conserva
la interfaz espacial oscura del compañero publicado de *Beyond Redfield*: filtros
laterales, mapa navegable e inspector de evidencias. Es una adaptación separada
bajo GPLv3. El explorador blanco anterior permanece idéntico y disponible.
El código, la procedencia exacta y las instrucciones de reconstrucción del nuevo
mapa están en `interactive/supergraph/README.md`; los datos científicos conservan
sus condiciones de origen.

La tercera herramienta es el [PaperGraph](interactive/papergraph/PaperGraph_CRISPR.html):
una guía del argumento científico con doce estaciones y tres profundidades de
lectura sobre los mismos hechos, condiciones y límites. El recorrido de lectura
se distingue expresamente de las relaciones científicas. Los grupos AND indican
premisas conjuntas; una arista aislada no sustituye al conjunto de requisitos.
El atlas completo, la figura vertical y la vista horizontal están en
[figures/papergraph/](figures/papergraph/), con SVG editables y versiones PNG/PDF.
Descargar el HTML y abrirlo localmente permite usarlo sin conexión; el visor de
código de GitHub no ejecuta una aplicación HTML.

La edición pública del PaperGraph conserva resúmenes científicos y vínculos a
los datos ya públicos. No incluye el manuscrito, el suplemento, sus instantáneas
ni citas de su texto privado. Sus identidades se declaran sin confundirlas con
fuentes públicamente accesibles. Las instrucciones propias están en
[interactive/papergraph/README.md](interactive/papergraph/README.md).

Los grafos de argumentos tienen antecedentes directos, como
[Micropublications (2014)](https://doi.org/10.1186/2041-1480-5-28).
Aquí la aportación es la integración concreta de lectura, evidencia y límites
en este proyecto. La utilidad pedagógica es una hipótesis por evaluar con
lectores; las comprobaciones del programa no la demuestran. Ninguna de las tres
herramientas establece precisión molecular ni validación biológica independiente.
