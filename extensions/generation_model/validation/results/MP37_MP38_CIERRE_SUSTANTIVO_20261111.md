# Validacion externa y adjudicacion de los dos limites

Checkpoint de 2026-11-11. No equivale a cierre de toda la investigacion.

La investigacion queda mejor posicionada: se completo un transporte congelado a otra adquisicion de HG005 y a otro individuo, HG001. El problema de carecer de un segundo individuo tiene ahora evidencia experimental directa dentro del dominio evaluado. La validacion de edicion con verdad ortogonal enlazada al mismo fragmento sigue NO ESTIMADA; localizar estudios candidatos no modifica ese estado.

| Panel, ventanas de 300 nt | Consenso previo, tolerancia 4 | Consenso 80%, tolerancia 0 | Mediana, tolerancia 0 |
|---|---|---|---|
| MP37, HG005 PAW88001 | 17 concordantes / 0 discordantes / 47 U | 61 / 0 / 3 | 64 / 0 / 0 |
| MP38, HG001 PAW79146 | 17 / 0 / 47 | 59 / 0 / 5 | 63 / 0 / 1 |

Cada panel tiene 64 loci planificados. Los umbrales se fijaron antes de calcular los resultados. Los dos criterios primarios pasan en MP37; en MP38 pasa mediana y falla consenso, porque 59/64 no alcanza el minimo del 95%. No se intercambian sus veredictos. Las ventanas de 150 nt conservan los resultados de mediana y consenso 0; consenso 4 en HG001 da 16/0/48. Estos son recuentos descriptivos en paneles dirigidos de genotipos homocigotos y no una precision universal o poblacional.

El avance mantiene la evidencia local. HG001 chr20:16520126 GTT>G sigue sin resolverse: 8 preferencias concordantes, 8 contrarias y 4 empates en 20 lecturas. En MP37, tres loci resueltos solo por mediana tienen 13 preferencias concordantes y 7 contrarias. Una decision local concordante no autentica cada lectura. La identidad exacta demostrada para esta regla es que, con arrepentimientos no negativos y umbral cero, mediana cero requiere una mayoria estricta de valores cero; los empates pueden apoyar ambos alelos. No es equivalente a consenso singleton del 80%.

Se verificaron las 10.204 distancias con un motor independiente, sin discrepancias, y las 1.024 decisiones con una caracterizacion independiente por recuentos. Una segunda extraccion de la fuente, con un procedimiento distinto, verifico 64 lecturas en 32 loci. MP37 conserva 1.271 identidades unicas; MP38 contiene 1.280 registros lectura-locus pero 1.274 lecturas unicas. Seis lecturas aparecen en dos loci cada una. Al retirar sus 12 registros en una sensibilidad posterior, sin reemplazo, no cambia ninguno de los ocho resultados agregados de HG001. Esa sensibilidad no transforma los datos originales en observaciones independientes.

## Fuentes candidatas y autoridad permitida

| Fuente | Evidencia localizada | Estado y limite |
|---|---|---|
| ConSeqUMI, GSE288938, DOI 10.1093/nar/gkaf1304 | 42 muestras, consensos publicos, ejemplo IRF3; recuperacion por UMI y Sanger descrita para plasmidos/barcodes | No se ha demostrado que los ejemplos CRISPR sean los mismos ejemplos validados por Sanger. El TAR de consensos inspeccionado no contiene trazas AB1. Falta enlace maquina-legible UMI-traza-evento. |
| Park et al., PRJNA780655, DOI 10.1126/sciadv.abo7676 | 60 runs; nueve constructos de ADN con barcode y Sanger descritos en metodos | Los titulos inspeccionados no identifican expresamente los runs de los nueve estandares. Zenodo 6805011/6805013 contiene codigo, no esos datos crudos identificados. |
| Simpson et al., PRJNA916868, DOI 10.1016/j.ymthe.2023.01.004 | 42 runs; ATXN2, nCATS sin PCR, amplicones y ddPCR | Candidato de validacion ortogonal de cuantificacion por muestra, no de etiqueta independiente por fragmento. Transgen humano en cerebro de raton. |
| IDMseq, PRJNA606194, DOI 10.1186/s13059-020-02143-8 | UMI y secuenciacion multiplataforma | El articulo indica subconjuntos distintos de moleculas entre Illumina y Nanopore. No adjudicar como verdad multiplataforma del mismo fragmento. |
| SAFE, Dryad 10.5061/dryad.905qfttp6 | Clones humanos, SNP flanqueante y ddPCR | Potencial verdad a nivel de clon y numero de copias; no etiqueta de cada fragmento. |
| Garvan HG001, gtgseq | Otra preparacion/fuente, BAM e indice publicos y cabecera consultada | No ejecutado. La cabecera no tiene M5 de referencia; resolver referencia o congelar remapeo explicito antes de admitir. |

Fuentes primarias: https://pmc.ncbi.nlm.nih.gov/articles/PMC12700093/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC9586483/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC10014281/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC7444080/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC10018342/ ; https://registry.opendata.aws/gtgseq/ . La inspeccion dirigida de metadatos y metodos no equivale a lectura exhaustiva de todos los suplementos.

## Evaluacion de madurez (juicio, no metrica de precision)

| Dimension, escala 0-10 | Anterior MP36 | Ahora MP38 | Motivo |
|---|---:|---:|---|
| Validacion externa de genotipo en el dominio estudiado | 6,5 | 7,5 | Nuevo individuo y otra adquisicion sin recalibracion; dominio limitado a panel dirigido. |
| Utilidad discriminante | 7,5 | 8,0 | Mayor resolucion transferida; consenso falla en HG001 y una jurisdiccion sigue U. |
| Validacion biologica de edicion / mismo fragmento | No estimada | No estimada | No existe aun enlace ortogonal admitido para ese endpoint. No corresponde inventar un score de precision. |
| Calidad cientifica de conjunto | 8,7 | 8,8 | Se cierra un limite de transporte; no se compensa aritmeticamente el limite molecular. |

No se ha demostrado imposibilidad universal de cerrar la validacion molecular. Si una biblioteca historica carece de un identificador que enlace mediciones independientes y no conserva material recuperable, mas computo no reconstruye esa identidad. Un conjunto nuevo con el enlace adecuado puede validar el metodo prospectiva o retrospectivamente en su propia jurisdiccion, sin reescribir los datos antiguos. La siguiente accion decisiva es admitir un enlace real barcode/UMI-Sanger-evento o una comparacion ddPCR por muestra con endpoints y denominadores armonizados; los dos objetivos tienen distinta autoridad.

Los resultados MP37/MP38, sus fallos y sensibilidades quedan registrados. El manuscrito y suplemento ingleses y la extension castellana se actualizaron en privado; la traduccion castellana completa del nucleo historico sigue pendiente. El paquete privado no es una version final de envio. No se publico ningun manuscrito en GitHub. Los resultados experimentales nuevos no hacen desaparecer los resultados adversos MP34/MP35 ni autentican los fragmentos originales de HEK3/CFTR.
