# MP34: contexto observable, discriminación y límite de cobertura

La investigación está mejor posicionada para explicar el mecanismo de discriminación, pero no ha cerrado todavía la validación biológica independiente. Se completó en Colab CPU una campaña nueva en chr22, con 128 loci de calibración y 256 de evaluación separados antes de observar resultados. Los 13.824 registros fuente generaron 27.648 vistas cortas/largas, emparejadas mediante el mismo historial de edición. Se verificaron los manifiestos, se recontaron las doce condiciones desde las distancias almacenadas y se reconstruyeron las 27.648 vistas a partir de sus operaciones literales. Esta reconstrucción no equivale a recalcular independientemente todas las distancias de alineamiento. Los datos son simulaciones controladas sobre secuencias genómicas reales; no son nuevas lecturas experimentales.

El resultado favorable es concreto: la regla combinada de distancia absoluta y relativa, con contexto de 300 bases fuente, aumenta mucho las decisiones unívocas correctas frente a 150. En el régimen basal pasa de 1.575 a 3.058 sobre 3.072, y de una decisión incorrecta a cero. Con variantes adicionales omitidas del catálogo pasa de 1.290 a 2.779 correctas, sin decisiones incorrectas en ambas condiciones. La explicación compatible con este experimento es que añadir secuencia observada aporta información discriminante cuando la vista corta no distingue suficientemente los candidatos. No demuestra que el mismo rendimiento se transporte a secuenciación real ni que cualquier ampliación del contexto mejore todos los indicadores.

La evidencia adversa impide declarar éxito global: el contraste principal preespecificado usaba la regla absoluta. En basal, su cobertura empírica por locus baja de 249/256 a 239/256 (93,36%), por debajo del umbral descriptivo del 95%; en el régimen de más ruido sí cumple el criterio. La condición combinada es secundaria y no puede reemplazar retrospectivamente al análisis principal. Además, con más ruido, sus errores unívocos aumentan de uno a dos al ampliar el contexto. El fallo de un umbral empírico en una muestra finita no refuta por sí solo una garantía conformal marginal bajo intercambiabilidad; sí impide presentar esta campaña como cumplimiento uniforme del criterio registrado. Cualquier ajuste posterior necesita una nueva evaluación reservada.

El siguiente paso decisivo sigue siendo B con lecturas reales y después C a nivel de variantes. Se priorizó HG005, cuyas lecturas candidatas fueron adquiridas en septiembre de 2024, frente al benchmark GIAB v4.2.1 anterior: esta separación temporal es más defendible que inferir independencia de la fecha de publicación de un conjunto HG002. La admisión detectó una referencia CRAM distinta de la secuencia UCSC inicialmente utilizada; no se interpretó ninguna lectura con esa discrepancia. El acceso y la identidad exacta de la referencia deben resolverse antes de cerrar B/C. Tampoco se equiparará concordancia con un genotipo de benchmark a verdad independiente de cada molécula. Los manuscritos permanecen privados y este cierre no constituye autorización científica de envío.

## Resultados completos

Cada fila contiene 3.072 vistas de evaluación; cobertura por locus significa que las doce vistas del régimen retienen su etiqueta fuente. Las longitudes 150/300 se refieren a bases fuente antes del ruido, no a longitudes observadas necesariamente exactas. «Vacío» significa que no se admite ningún candidato; «ambiguo», que se admiten ambos.

| Régimen | Contexto | Regla | Correctos unívocos | Incorrectos unívocos | Ambiguos | Vacíos | Loci cubiertos /256 |
|---|---:|---|---:|---:|---:|---:|---:|
| Basal | 150 | Absoluta | 679 | 0 | 2386 | 7 | 249 |
| Basal | 300 | Absoluta | 1202 | 0 | 1852 | 18 | 239 |
| Basal | 150 | Combinada | 1575 | 1 | 1489 | 7 | 248 |
| Basal | 300 | Combinada | 3058 | 0 | 11 | 3 | 253 |
| Más ruido | 150 | Absoluta | 252 | 0 | 2817 | 3 | 253 |
| Más ruido | 300 | Absoluta | 412 | 0 | 2656 | 4 | 252 |
| Más ruido | 150 | Combinada | 1272 | 1 | 1799 | 0 | 255 |
| Más ruido | 300 | Combinada | 2776 | 2 | 291 | 3 | 251 |
| Catálogo incompleto | 150 | Absoluta | 553 | 0 | 2511 | 8 | 249 |
| Catálogo incompleto | 300 | Absoluta | 1089 | 0 | 1963 | 20 | 237 |
| Catálogo incompleto | 150 | Combinada | 1290 | 0 | 1774 | 8 | 249 |
| Catálogo incompleto | 300 | Combinada | 2779 | 0 | 292 | 1 | 255 |

## Scoring de seguimiento

Escala editorial orientativa de 0–10, no una medición ni probabilidad de aceptación. Se conserva la valoración previa para no convertir una mejora de alcance controlado en validación biológica o utilidad general ya establecida.

| Dimensión | Antes | Ahora | Cambio sustentado |
|---|---:|---:|---|
| Calidad científica global | 8,5 | 8,5 | Más evidencia y mejor delimitación; el criterio principal no pasa uniformemente |
| Rigor metodológico | 8,7 | 8,7 | Separación prospectiva calibración/evaluación; resultado adverso preservado |
| Validación computacional | 8,7 | 8,7 | Nuevo chr22, emparejamiento exacto y recuento completo; alcance todavía controlado |
| Utilidad discriminante general | 6,5 | 6,5 | Gran mejora secundaria en simulación, pendiente transporte real |
| Validación biológica independiente | 3,0 | 3,0 | Fuente candidata mejor admitida; C no ejecutada todavía |

## Custodia y reproducción

Código y protocolo congelados: SHA-256 del protocolo `3584ab1354beee1ec4eaeb10ffe7a2245e7aa2e2a9820242633acecc2377dadc`; trabajador `fe4e63b2271795df0a6163c5f0541f851035f4af6c1ea4a99ee4cd488642fb8c`. Archivo final `MP34_test_7.zip`: 1,782 MB; SHA-256 `ced659bfea9b5ff95991bea949090b817e3ec09805aa19678870feea6e7741b0`. Los doce lotes se conservaron en Drive; los recibos de carga no son una lectura completa verificada de los objetos en Drive. La descarga desde Colab y los manifiestos internos sí fueron verificados. No se borraron originales ni PDF. El tiempo de cálculo de los doce lotes fue aproximadamente 93 segundos; no incluye preparación, transferencias ni validación de fuentes reales.
