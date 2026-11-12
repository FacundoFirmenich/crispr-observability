# Guía para investigadores y comunidad

[English](RESEARCHER_GUIDE.md) · [Caso guía](LEADING_CASE.es.md)

## Elegir la tarea científica

| Tarea | Entrada | Resultado |
| --- | --- | --- |
| Entender o cuestionar una afirmación | Supergraph y PaperGraph | Premisas, fuentes, alcance y mapa exportable |
| Examinar recuentos históricos exactos | Explorador epistémico numérico | Recuentos originales e intervalos condicionados |
| Manipular secuencias/estructuras o examinar alineamientos | Read Observability Lab | Escena 3D, coordenadas, JSON/replay, PNG y TSV de alineamiento cuando corresponda |
| Usar varias entradas o conectar un IDE | Motor CLI/API/MCP | Manifiesto versionado e informe con procedencia |
| Recalcular resultados publicados | Autoejecutable, Python, Streamlit o cuaderno | Cálculo delimitado y recibo de ejecución |

Cada interfaz responde una pregunta. La escena no reemplaza al cálculo del
alineamiento, y el hash identifica los datos, no autentica la etiqueta biológica.
Empiece por el caso guía y declare después el contrato de su propio estudio.

## Abrir y navegar

La portada enlaza las herramientas; descargue los HTML para trabajar sin red.
El laboratorio necesita WebGL 2. Arrastre para orbitar, utilice la rueda para
acercar y los controles de cámara para recuperar la orientación. El Supergraph
permite mover tarjetas, columnas y el mapa, resaltar relaciones y seguir un
recorrido conceptual con reproducción y pausa. PaperGraph ofrece tres niveles
de explicación: Student, Technical y Audit.

El laboratorio tiene cinco modos: contención de anclas, tramos competidores,
ADN ideal, estructuras importadas y alineamientos. Cambie una variable, examine
la explicación y el resultado, y guarde el replay. Las vistas de datos
conservados y los cambios contrafactuales mantienen alcances distintos. Los
controles internos del visor molecular conservan sus nombres científicos en
inglés; esta guía y el dashboard proporcionan orientación en castellano.

## Incorporar datos

| Entrada | Contrato |
| --- | --- |
| FASTA | Una secuencia, 2–360 bases ACGTN, inspección local |
| PDB / mmCIF textual / XYZ | Coordenadas normalizadas; hasta 30.000 átomos y 6 MB de texto |
| SAM | Registro, CIGAR, referencia y ventana explícitos |
| BAM con BAI/CSI | Lectura indexada local; BAI contrastado con fixtures públicos; cobertura de pruebas CSI limitada |
| JSON portátil | Esquemas versionados de secuencia, estructura, alineamiento o replay |
| Manifiesto | Hasta 32 entradas y 64 operaciones; comprobación SHA-256 opcional |

Las coordenadas empiezan en cero y usan intervalos semiabiertos. Las distancias
atómicas se expresan en angstroms. La ventana admite hasta 50.000 posiciones y
2.500 registros. El visor mantiene un caso activo; el motor admite entradas
múltiples. CRAM y procesamiento genómico completo no están implementados.
Consulte la [guía de formatos](../interactive/lab/README.md) y el
[contrato de uniones](../engine/GUIA_EVIDENCIA_SPLICING_v0.5.md).

La importación normaliza a JSON; no preserva necesariamente todos los campos
de cada formato original. Conserve fuente y hash junto al caso. «No observado»
en una ventana acotada no equivale a ausencia en toda la biblioteca.

## Trabajo programático y MCP

Instale con `npm ci --prefix engine`. Valide con
`node engine/project-cli.mjs check <manifiesto>` y ejecute con `run`.
El mismo núcleo funciona desde JavaScript y desde el servidor stdio:
`node engine/project-mcp.mjs <directorio-autorizado>`. Expone `lab_project_run`
y restringe las rutas al espacio autorizado. El puente MCP del visor permite
normalizar, calcular anclas y tramos, inspeccionar alineamientos y exportar
renders. Consulte [MCP](../interactive/lab/MCP.md).

El agente debe conservar manifiesto, parámetros y fuentes de cada resultado.
Pedir una operación no significa haberla ejecutado. Un nombre de archivo
coincidente no demuestra identidad de muestra ni molécula.

## Reproducción, Streamlit y Docker

Los paquetes nativos v0.9.0 incluyen runtime Python y entradas congeladas:
Windows x86-64 `START_HERE.bat`, Linux x86-64 `START_HERE.sh`, macOS Intel/Apple
silicon `START_HERE.command`. Cada binario corresponde a su plataforma. Las
mejoras actuales de fuente son posteriores a esa release: consulte el commit
antes de atribuirle capacidades nuevas a un paquete antiguo.

El dashboard tiene controles en castellano e inglés, evidencia, Supergraph,
laboratorio, reproducción por alcance, ledger y guías:

```sh
python -m pip install -r reproducibility/requirements-reviewer.txt
python -m streamlit run reproducibility/reviewer_app.py
```

O utilice el contenedor:

```sh
docker build -f Dockerfile.reviewer -t crispr-reviewer .
docker run --rm -p 127.0.0.1:8501:8501 -v crispr-review-results:/results crispr-reviewer
```

Abra `http://localhost:8501`. El volumen conserva resultados. Docker incluye
visores y reproducción Python; el motor Node tiene instalación separada y
versionada. Quick/counts/full tienen alcances diferentes que constan en el
recibo. Elija almacenamiento persistente adecuado. Ningún lanzador inicia
Colab ni Kaggle.

## Hardware y accesibilidad real

El desarrollo local y la revisión editorial se realizaron en un portátil
convencional de 2016, Intel Core i5-6200U a 2,30 GHz, según informa el autor.
Las campañas grandes utilizaron recursos alojados y los binarios se construyen
en CI nativa. Esto describe el contexto de desarrollo; no prueba que todo
análisis completo se ejecute en cualquier ordenador de diez años. El HTML es
pequeño, pero el render depende del navegador y WebGL. La reproducción Python
congelada no requiere una GPU dedicada; recursos y tiempo dependen del alcance.

## Participar en la investigación

La infraestructura se usa o adapta en [tres líneas](RESEARCH_LINES.es.md).
Para informar un problema aporte commit, comando/manifiesto, hash, coordenadas,
resultado esperado y observado y unidad de evidencia. Publique sólo datos que
pueda compartir. Use issues para fallos reproducibles, discussions para
preguntas y propuestas, y pull requests para cambios contrastados. La
[hoja de ruta](COMMUNITY_ROADMAP.es.md) describe vías concretas.

Son aportes útiles un caso independiente documentado, una prueba adversa de
coordenadas, un adaptador verificado, una explicación accesible o un recibo de
reproducción externa. Los resultados desfavorables también cuentan. Citación,
uso independiente y reproducción son resultados distintos; no se sustituyen
con insignias, promoción o recuentos de issues.

## Publicar y citar

Conserve entradas originales, commit, manifiesto, informe, unidad de evidencia,
resultados irresueltos y versión del visor. Cite componente y versión exactos.
El DOI 10.5281/zenodo.24121127 corresponde a v0.8.1, no al laboratorio posterior.
El código actual original es GPL-3.0-only; la extensión de modelos generativos
retiene MIT y cada dependencia conserva su licencia. El manuscrito y el
suplemento privados no forman parte de las releases públicas de código.
