# Read Observability Lab v0.4 · Guía de usuario del motor multientrada

El motor ejecuta análisis de secuencias, coordenadas estructurales y alineamientos a partir de un manifiesto por estudio. Un mismo manifiesto puede usar archivos de varios formatos, comparar segmentos de distintas entradas y producir un informe con resultados y procedencia por operación. Cada estudio conserva su propio espacio de nombres. Podés reutilizar el mismo identificador en dos proyectos sin compartir archivos ni resultados entre ellos.

Esta guía sirve para conectar un paper nuevo con sus insumos existentes. El motor se ejecuta por terminal, desde JavaScript o mediante MCP en un IDE. El visor web existente conserva su flujo interactivo de un caso activo.

## 1. Preparar el estudio

En un entorno remoto con acceso a tus insumos y Node.js 22, obtené el código de la rama `lab-engine-multientrada-v04` de `FacundoFirmenich/crispr-observability` y ejecutá, desde la raíz del repositorio:

~~~sh
npm ci --prefix engine --no-fund --no-audit
~~~

En la carpeta del estudio, guardá `project.json` junto a sus entradas, organizadas en subcarpetas si te conviene. En el manifiesto usá rutas relativas a esa carpeta, con `/` incluso en Windows. Los IDs de entrada y operación deben ser únicos en su lista y empezar por una letra. Admiten letras, números, puntos, guiones y guiones bajos, hasta 64 caracteres. Podés declarar hasta 32 entradas y ordenar de 1 a 64 operaciones.

| Tipo | Archivos | Qué identifica |
| --- | --- | --- |
| `sequence` | FASTA de un registro o JSON `observability-case/1` de secuencia | Cadena de ADN de 2 a 360 bases, alfabeto ACGTN |
| `structure` | PDB, mmCIF/CIF, XYZ de un fotograma o JSON de estructura compatible | Lista normalizada de átomos y coordenadas en ångströms |
| `alignment` | SAM o JSON de alineamiento compatible con el laboratorio | Registros de alineamiento, con referencia, coordenadas, CIGAR y banderas |
| `indexed-bam` | BAM con índice BAI o CSI, declarados como dos archivos | Acceso por intervalo a una referencia del encabezado BAM |

Un archivo textual tiene un límite de 6 MB. Para loci, el máximo es 50.000 bases y el informe admite hasta 2.500 registros por caso. La lectura indexada aplica además límites de 24 MiB para el rango de lectura y 32 MB para el índice. Podés añadir `sha256` en minúsculas al objeto de cada archivo para exigir una coincidencia exacta. Para BAM, esta opción lee el archivo completo en streaming, con un costo adicional de E/S. Las consultas ordinarias usan acceso indexado.

## 2. Escribir y ejecutar el manifiesto

Ejemplo de configuración. Ajustá los nombres de archivo, referencias e intervalos a los insumos del otro estudio.

~~~json
{
  "schema": "observability-project/1",
  "project_id": "mathematical-biology",
  "datasets": [
    {"id":"reference","type":"sequence","file":{"path":"inputs/reference.fa"}},
    {"id":"observed","type":"sequence","file":{"path":"inputs/observed.fa"}},
    {"id":"structure","type":"structure","file":{"path":"inputs/model.cif"}},
    {"id":"reads","type":"indexed-bam",
     "bam":{"path":"inputs/reads.bam"},"index":{"path":"inputs/reads.bam.bai"}}
  ],
  "jobs": [
    {"id":"inspect-structure","op":"inspect","input":"structure"},
    {"id":"compare-strings","op":"compare-sequence-spans",
     "a":{"input":"reference","span":[2,14]},
     "b":{"input":"observed","span":[4,16],"reverse_complement":true}},
    {"id":"distance","op":"atom-distance","input":"structure","atoms":[0,1]},
    {"id":"locus","op":"alignment-locus","input":"reads",
     "locus":{"ref":"chr1","span":[10000,10080]},"alignment_index":0,
     "query_index":0,"min_mapq":20,"original_orientation":true},
    {"id":"finite-window","op":"anchor-containment",
     "anchor":[20,28],"shift":3,"read_end":30}
  ]
}
~~~

~~~sh
node engine/project-cli.mjs check /workspace/study/project.json
node engine/project-cli.mjs run /workspace/study/project.json > /workspace/study/report.json
node engine/project-cli.mjs schema > /workspace/study/project.schema.json
~~~

`check` valida el manifiesto y sus referencias internas. `run` abre las entradas usadas y verifica su formato, contenido y límites. Conservá la salida en la carpeta remota del estudio. Una ruta inexistente, un checksum distinto o una operación incompatible detiene la ejecución con el ID correspondiente y un código de salida distinto de cero. Las entradas declaradas sin uso figuran como `not-read`.

## 3. Elegir la operación e interpretar la respuesta

Los intervalos son de base cero y semiabiertos `[inicio, fin)`. El inicio se incluye y el fin se excluye. Un intervalo `[2,14)` contiene 12 posiciones. El nombre de una referencia debe coincidir con el encabezado del alineamiento, incluyendo el prefijo si lo hay.

| Operación | Resultado e interpretación |
| --- | --- |
| `inspect` | Devuelve un resumen de la entrada cargada, como longitud de la secuencia, número de átomos, registros SAM o referencias en el encabezado BAM. |
| `compare-sequence-spans` | Extrae ambos intervalos, aplica el reverso complementario únicamente si se solicita y comprueba igualdad de cadenas. Devuelve longitud de intersección si ambos intervalos pertenecen a la misma entrada. Para entradas distintas, devuelve `overlap: null`. La igualdad de cadenas carece de autoridad para establecer correspondencia molecular entre archivos. |
| `anchor-containment` | Desplaza el intervalo de anclaje según `shift` y cuantifica qué parte queda en la ventana `[0,read_end)`. Sirve para probar observabilidad geométrica. |
| `atom-distance` | Calcula distancia euclídea en Å entre dos índices de la lista normalizada. El resultado identifica ambos átomos. Elegí los índices según esa lista. El lector estructural toma el primer modelo cuando el formato lo permite y aplica sus reglas de selección de posiciones alternativas. |
| `alignment-locus` | Devuelve el registro seleccionado, cobertura por ventana, proyección opcional de una posición de lectura y evidencia tabulada. `alignment_index` elige una posición en la lista devuelta por el lector. `min_mapq` e `include_non_primary` gobiernan el cálculo de cobertura. `original_orientation` cambia la interpretación de `query_index`; la referencia sigue su propio sistema de coordenadas. |

Si contrastás dos observaciones, registrá en el manifiesto dos operaciones con iguales coordenadas y criterios cuando el diseño lo permita. Identificá explícitamente la entrada en cada una y verificá el alcance antes de comparar. En SAM y BAM, un índice numérico igual puede seleccionar registros distintos. Comprobá el identificador, referencia, posición y CIGAR devueltos. La cobertura y la coincidencia de cadenas son resultados operacionales de los datos ingresados. La independencia de muestras, el enlace preamplificación y la exactitud biológica requieren evidencia propia.

## 4. Integrar con un IDE agéntico o código

Servidor MCP por stdio:

~~~sh
node /ruta/al/repositorio/engine/project-mcp.mjs /workspace/estudios
~~~

Configurá el cliente con `node` como comando y las dos rutas anteriores como argumentos. El servidor expone `lab_project_run`. Para validar, llamalo con `{"manifest_path":"estudio-a/project.json","check_only":true}`. Para ejecutar, usá `check_only:false` o dejá ese campo fuera. Las rutas se resuelven dentro del workspace autorizado del servidor. Cada manifiesto resuelve sus entradas respecto de su propia carpeta, igual que la CLI. Podés abrir varios estudios en un servidor con una raíz común y conservar rutas de datos independientes.

Desde JavaScript, importá `readManifest` y `runProject` de `engine/project-engine.mjs` y ejecutá `await runProject(await readManifest(root, "project.json"), root)`. El motor devuelve el informe como objeto. La aplicación que lo invoca decide dónde persistirlo o cómo visualizarlo.

## 5. Custodia y revisión

El informe `observability-project-report/1` conserva la versión del motor, el ID del estudio, el SHA-256 del manifiesto normalizado, las entradas usadas, los parámetros efectivos y las salidas en orden. El digest del manifiesto describe la configuración normalizada. Los hashes de fuente documentan el contenido. En BAM indexado, el informe del lector conserva la procedencia del índice y de los bloques consultados. El hash integral del BAM se verifica cuando se declara en la entrada. Archivá juntos manifiesto, insumos autorizados, informe, commit y condiciones de ejecución.

Para entregar el estudio a otro chat, pasale esta guía y las rutas remotas autorizadas. Pedile que identifique formatos y coordenadas, arme un manifiesto independiente, corra `check` y `run`, y devuelva el informe junto con la interpretación de cada operación. Si los insumos requieren un formato distinto, acordá primero una conversión trazable y preservá el original.
