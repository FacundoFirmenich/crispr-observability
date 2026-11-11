# Cierre científico MP34–MP35: validación real y frontera discriminante

**Estamos mejor posicionados para sostener y delimitar la investigación, pero A+B+C no queda cerrado como puerta de envío.** Se completaron la campaña controlada MP34 y la evaluación real MP35 en Colab CPU, con checkpoints continuables conservados por lote. MP34 usó 128 loci de calibración y 256 de evaluación, con 27.648 vistas simuladas emparejadas. MP35 utilizó HG005: 64 loci de calibración y otros 64 de evaluación, 2.539 registros lectura–locus en total. Sus lecturas fueron adquiridas en septiembre de 2024 y los archivos del benchmark GIAB v4.2.1 tienen fecha de septiembre de 2021. Esta separación aporta verdad externa para el genotipo de los objetivos seleccionados; no convierte el benchmark en una observación independiente de cada molécula ni en validación funcional de edición génica.

**El resultado real es una ausencia de mejora por ampliar el contexto bajo la regla congelada.** En las 1.266 observaciones de evaluación, tanto con 150 como con 300 nucleótidos hubo 508 decisiones unívocas concordantes y 758 ambiguas; no hubo decisiones unívocas discordantes ni conjuntos vacíos. Las transiciones emparejadas fueron todas idénticas. Con el consenso por locus fijado antes de los resultados, se resolvieron 25/64 objetivos: 10/32 variantes presentes y 15/32 ausentes; quedaron 39/64 sin adjudicar. No corresponde vender el cero de discordancias como exactitud universal ni retirar las abstenciones del denominador. En MP34 sí hubo una mejora secundaria considerable de resolución, pero su criterio principal no pasó uniformemente: la regla absoluta de contexto largo tuvo cobertura basal239/256, inferior al95% descriptivo requerido. La ganancia de aquella simulación no se transportó a esta evaluación real.

**La nueva explicación sustantiva es una limitación de la regla calibrada, y además apareció una limitación de independencia.** El margen relativo congelado fue4 en ambas longitudes. En las732 observaciones de indels de2–4 nucleótidos, la diferencia entre puntuaciones quedó dentro de ese margen y ambas alternativas permanecieron admitidas; el contexto extra no cambió esa decisión. En indels de5–30 nucleótidos se resolvieron508/534 observaciones y quedaron26 ambiguas. Es un diagnóstico posterior a los resultados, sin modificar el análisis principal. Separadamente, la auditoría prefijada encontró152 identificadores de lectura compartidos entre calibración y evaluación, que afectaban177 registros de evaluación. La separación de loci no había separado todas las moléculas largas. La sensibilidad secundaria excluyó esos177 registros, mantuvo los umbrales y los64 objetivos originales y dejó55 loci con al menos10 lecturas: ambas longitudes siguieron idénticas, con399 decisiones unívocas concordantes y657 ambiguas entre1.056 registros evaluables;20/64 loci se resolvieron y44/64 quedaron no adjudicados. Los otros33 registros retenidos pertenecían a nueve loci con menos de10 lecturas y no se contaron como evaluación por lectura. Esta sensibilidad no transforma retrospectivamente la partición original en un ensayo independiente.

**Las comprobaciones técnicas respaldan que estos resultados pertenecen al método ejecutado y no a un desacuerdo entre dos implementaciones de distancia.** Edlib recalculó65.452 distancias de ambas campañas sin diferencias; una segunda extracción verificó32 lecturas en16 loci prefijados por hash. Todos los contextos analizados quedaron dentro del intervalo de referencia de las plantillas según el CIGAR almacenado; no se observó una salida de ese intervalo que explicara este resultado. Antes del análisis se corrigieron dos obstáculos de acceso: referencia CRAM distinta de la de UCSC y desfase de la cabecera respecto del índice publicado. La referencia exacta fue verificada, el índice original se preservó y la reparación se contrastó con cabeceras y decodificación efectiva. No se atribuye la causa histórica del desfase al proveedor sin confirmación. Los lotes MP34 sumaron aproximadamente93 segundos; los128 loci MP35 sumaron256 segundos de ejecución incluyendo acceso remoto. Son tiempos de esas operaciones, no la duración total de preparación, custodia y auditoría.

**La consecuencia para el preprint es incorporar evidencia favorable y adversa con un alcance más preciso, no declarar superado el cuello de botella.** Ya existe una evaluación externa real de variantes, pero la utilidad discriminante es parcial y el análisis original presenta dependencia entre particiones. El siguiente diseño decisivo debe separar calibración y evaluación por regiones y por identificadores de lectura, y distinguir la decisión a nivel de variante de la exigencia de acertar todas las lecturas de un locus. Una agregación de evidencia por variante podría evaluarse prospectivamente, manteniendo las abstenciones por lectura; todavía es una propuesta, no un resultado. No se retocaron los umbrales ni se recalibró sobre el conjunto evaluado. Los manuscritos siguen privados y no se ha regenerado en este cierre un PDF del paper ni del supplementary. El código, resultados y recibos se conservan en el repositorio privado de desarrollo y en Drive; la publicación de un manuscrito actualizado exige integrar esta adjudicación y comprobar nuevamente sus afirmaciones.

## Resultados reales: denominadores explícitos

| Resultado | Primario150nt | Primario300nt | Sensibilidad sin lecturas compartidas150/300nt |
|---|---:|---:|---:|
| Registros evaluables lectura–locus |1266|1266|1056 en ambas|
| Unívocos concordantes |508|508|399 en ambas|
| Unívocos discordantes |0|0|0 en ambas|
| Ambiguos |758|758|657 en ambas|
| Vacíos |0|0|0 en ambas|
| Loci evaluables / planificados |64/64|64/64|55/64|
| Loci resueltos / planificados |25/64|25/64|20/64|
| Positivos resueltos / positivos planificados |10/32|10/32|8/32|
| Negativos resueltos / negativos planificados |15/32|15/32|12/32|
| Ganancia por ampliar el contexto |—|0|0|

La resolución por lectura del análisis primario fue40,13%; por locus39,06%. En la sensibilidad, la resolución por locus sobre los64 planificados fue31,25%. Los porcentajes positivos y negativos son fracciones correctamente adjudicadas sobre objetivos planificados con abstención explícita, no estimaciones de rendimiento poblacional. El panel estaba equilibrado y enriquecido para indels simples; no establece precisión universal ni sensibilidad de descubrimiento de variantes.

## Scoring de preparación de esta validación

Escala de evidencia explícita: **0=no demostrado/no cumplido;1=parcial o limitado;2=cumplido en el alcance declarado**. No es un ranking de revistas ni una nueva puntuación editorial del manuscrito, que no se ha revisado íntegramente en este cierre.

| Criterio | Antes de MP34–35 | Ahora | Interpretación |
|---|---:|---:|---|
| Verificación independiente completa de las distancias de estas campañas |0|2|65.452 comparaciones con segunda implementación, cero diferencias|
| Evaluación real contra benchmark de otro individuo |0|1|HG005 ejecutado; alcance genotípico, un individuo y una librería|
| Discriminación útil en el panel real |0|1|25/64 loci en primario;20/64 en sensibilidad|
| Mejora por contexto adicional en lecturas reales |0|0|Ninguna decisión cambió|
| Separación estricta de lecturas entre calibración y evaluación primaria |0|0|152 IDs compartidos; sensibilidad secundaria conservada|
| Puerta A+B+C plenamente satisfecha para envío |0|0|Falta una evaluación prospectiva con partición adecuada y cierre editorial|

## Reproducibilidad y custodia

El resultado primario MP35 permanece inalterado. El protocolo, selección, código y corrección del índice están enlazados por hashes en los checkpoints. Los archivos se descargaron desde Colab y sus manifiestos se comprobaron completamente. Los recibos de subida a Drive acreditan la carga, no una lectura completa posterior de los objetos en Drive. No se borraron originales ni PDF. La licencia del código original de esta línea es MIT. Las fuentes primarias externas conservan sus condiciones propias.

Fuentes de interpretación técnica: especificación CRAM https://samtools.github.io/hts-specs/CRAMv3.pdf ; implementación independiente de distancia https://github.com/Martinsos/edlib . Datos: https://epi2me.nanoporetech.com/giab-2025.01/ y los archivos GIAB HG005 v4.2.1 identificados en `MP35_SOURCE_MANIFEST.json` y `MP35_BENCHMARK_SOURCE_DATES.json`.
