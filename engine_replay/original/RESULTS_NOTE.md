# Sensibilidad real de alineamiento: cierre verificable

Se completaron y verificaron 144 ejecuciones del motor oficial CRISPResso2 2.3.4, commit `1fc4b629e5400a48d6569b0ab4b8d1766f33636f`, sin modificar su fuente. Se utilizaron las nueve corrientes originales, con 177.006 registros raw y el mismo filtro externo de media Phred+33 ≥30. No se añadieron observaciones biológicas.

Los 18 baselines frescos, dos referencias por nueve corridas, coinciden exactamente con los resultados congelados: ambos eventos moleculares proyectados y las tablas nativas completas de cuantificación. Cada ejecución fue comprobada mediante identidad individual raw, dos implementaciones de proyección, recuento desde FASTQ anotado y recuento desde tabla de alelos. Las cantidades nativas se cotejaron con las clases anotadas del motor. `BASELINE_VALIDATION.json` y `VERIFICATION.json` documentan esos cierres.

La base personalizada es apertura −20, extensión −2 y puntuación mínima 60. Las seis perturbaciones OAT varían sólo apertura −10/−30, extensión −1/−4 o score 50/70, manteniendo el resto de la configuración y la ventana congelados. Son 108 perturbaciones más 18 baselines. El perfil stock sin guía añade 18 trabajos separados y no entra en los rangos OAT: cambia simultáneamente opciones omitidas de alineamiento y ventana. Las decisiones se fijaron antes de obtener estos resultados, sobre datos retrospectivos ya conocidos; no constituyen un holdout prospectivo.

Los contrastes PE−mock conservan cota inferior positiva en las 14 combinaciones de referencia y configuración OAT. El menor extremo inferior observado por preparación es:

| Preparación | Cadena programada exacta (pp) | Péptido local nativo (pp) |
|---|---:|---:|
| 1 | 15,90186372 | 16,01852148 |
| 2 | 28,70347550 | 29,05793755 |
| 3 | 24,92238411 | 25,16093373 |

Son extremos de envolventes de identificación sobre N raw, no intervalos de confianza. Los seis mínimos ocurren con referencia mutante y score 70. No demuestran exactitud molecular de las etiquetas.

La estabilidad de esos contrastes agregados coexiste con sensibilidad individual. Dentro de las siete configuraciones OAT de referencia mutante, 48 IDs raw presentan al menos una etiqueta E y una NE, tanto para cadena exacta como para péptido. Con referencia programada, los conflictos resueltos son cero para exacto y uno para péptido. Entre los dos baselines de distinta referencia ya existían 135 IDs conflictivos para exacto y 137 para péptido; al incluir las 14 configuraciones se observan 161 y 163, respectivamente: 26 IDs adicionales por evento respecto de esos baselines incluidos. Estos recuentos son por ID distinto dentro del ámbito declarado. La suma de 63 conflictos péptido frente al baseline en tablas de transición cuenta comparaciones solapadas y no equivale a 63 IDs independientes.

Los cambios de puntuación mínima producen cambios de resolución frente al baseline, sin conflictos E↔NE en esta grilla. Las penalizaciones de gaps sí pueden producir conflictos resueltos. `TRANSITIONS.json` conserva su dirección por corrida; `JOINT_STABILITY.json` conserva los conteos conjuntos y `JOINT_LABELS.jsonl.gz` las firmas por los 177.006 IDs originales, con leyenda separada. `SENSITIVITY_RANGES.json` contiene mínimos y máximos E/NE/U por corrida y referencia; sus extremos son marginales y no deben combinarse como si ocurrieran simultáneamente en una configuración.

En todos los 18 trabajos stock, los eventos proyectados conservan exactamente los conteos del baseline correspondiente. Las ventanas nativas efectivas, sin guía y con exclusiones upstream de 15 bases a cada extremo, son posiciones 15–208 en referencia mutante y 15–211 en programada. El incentivo configurado es 1, pero su vector efectivo es cero en estos trabajos sin guía. El porcentaje nativo MODIFIED cambia: en SRR29623272/mutante pasa de 28,45459952% a 42,61498810%; en programada, de 81,64155432% a 85,86439334%. Su denominador son las lecturas alineadas y su predicado es modificación respecto de la referencia dentro de su ventana; no es el evento de péptido CFTR nativo.

Se preservaron dos ejecuciones fallidas por una entrada Q30 local truncada después de su preparación y una salida marcada completada pero detectada truncada. No se interpretaron. La entrada fue regenerada desde raw autenticado hasta igualar exactamente su SHA preparado; las ejecuciones afectadas se repitieron. La causa de esas truncaciones no está establecida. Las ejecuciones posteriores usaron staging en `/tmp`, cierre del motor, auditoría gzip/zip e identidad raw, fsync y comprobación de hashes antes y después de publicar su carpeta final. `failed_attempts/`, `INPUT_REPAIR.json` y los manifiestos conservan la historia.

La evidencia nueva establece sensibilidad computacional finita de una proyección local: cambiar alineamientos puede alterar etiquetas individuales aun cuando el contraste siga positivo. El evento, el denominador y la configuración deben viajar juntos al comparar resultados. No adjudica cuál alineamiento es verdad molecular ni extiende la conclusión a función, eficacia, población o datos futuros. P1/M8/v26 y S1 originales permanecen intactos. La siguiente frontera publicable sigue siendo la utilidad externa y una adjudicación independiente cuando se reclame exactitud.
