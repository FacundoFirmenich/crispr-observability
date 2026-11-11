# Cierre científico MP33: catálogo, calibración y ambigüedad

MP33 se ejecutó realmente por la CLI instalada de Google Colab, en la sesión CPU `crispr-mp33-cpu`. Queda mejor localizado el cuello de botella científico: completar el catálogo repara parte de los errores, mientras que cambiar solamente los umbrales no los repara. Las 24 comparaciones cruzaron dos catálogos de puntuación con dos catálogos de calibración, tres regímenes y dos niveles nominales. Se verificaron íntegramente 36.864 decisiones a partir de las puntuaciones congeladas, se recalcularon 16 umbrales y se reprodujeron exactamente los seis brazos anteriores de MP32. Es un análisis retrospectivo emparejado sobre MP31: no añade observaciones independientes, nuevas moléculas ni evidencia biológica ortogonal.

Resultado principal con nivel nominal del 95%, 128 bloques y 1.536 lecturas del régimen con variantes omitidas. Cada bloque exige incluir el alelo generador en sus doce lecturas; el porcentaje de bloques observado no es una garantía de cobertura por locus.

| Catálogo para puntuar | Catálogo para calibrar | Bloques cubiertos | Únicas correctas | Únicas erróneas | Ambiguas | Vacías |
|---|---|---:|---:|---:|---:|---:|
| Canónico | Canónico | 126/128 | 794 | 5 | 737 | 0 |
| Canónico | Ampliado | 25/128 | 633 | 5 | 718 | 180 |
| Ampliado | Canónico | 127/128 | 802 | 1 | 733 | 0 |
| Ampliado | Ampliado | 127/128 | 802 | 1 | 733 | 0 |

La comparación cruzada distingue el mecanismo: la mejora aparece al incorporar las secuencias omitidas incluso conservando los umbrales canónicos. Recalibrar con el catálogo ampliado y aplicar ese umbral más estricto al catálogo incompleto conserva los cinco errores y provoca 180 abstenciones vacías. Entre las cinco respuestas erróneas originales, dos pasan a correctas, dos pasan a ambiguas y una sigue siendo errónea. Además, siete respuestas antes ambiguas pasan a correctas y una antes correcta pasa a ambigua. No equivale a haber descubierto haplotipos: el catálogo ampliado conoce por construcción las variantes utilizadas para generar estas lecturas.

Un diagnóstico posterior a los resultados recalculó con una implementación independiente de distancia los cinco casos adversos, tanto observados como sin corrupción simulada. En dos lecturas centradas, el catálogo incompleto prefiere el alelo incorrecto incluso sin ruido; el ampliado recupera el verdadero. En otras dos, la secuencia limpia del fragmento es compatible con ambos alelos y completar el catálogo restaura la ambigüedad. En el caso persistente `chr21_31794476`, la lectura limpia puntúa 0 frente a ambos alelos del catálogo ampliado; una inserción y una sustitución simuladas producen distancias 2 al alelo generador y 1 al alternativo. Por tanto, conocer el catálogo completo no basta para convertir una preferencia de distancia en identificación de procedencia. Este análisis de cinco casos fue seleccionado después de observar los errores y no constituye una muestra independiente de validación.

El panel con ruido alto conserva los resultados de MP32 en los cuatro cruces de catálogo al 95%: 125/128 bloques, 596 únicas correctas, cero únicas erróneas y tres vacías. Los ocho brazos de controles externos rechazan las 1.536 lecturas respectivas. Son controles externos lejanos, no un panel de alta homología. La siguiente decisión científica es evaluar una regla explícita de abstención ante fragmentos compatibles con varios alelos y contrastarla en loci separados del desarrollo, junto con negativos genómicos cercanos seleccionados sin mirar los resultados del clasificador. No se ha resuelto aún esa frontera ni la obstrucción de procedencia de bases de MP27.

Scoring de conjunto orientativo sobre diez, manteniendo la rúbrica anterior: juicio comparativo, no probabilidad de publicación. Este diagnóstico aumenta la explicación del mecanismo; no justifica inflar la valoración global por reutilizar las mismas observaciones.

| Dimensión | Antes: MP31–MP32 | Después: MP33 | Interpretación |
|---|---:|---:|---|
| Calidad científica global | 8,5 | 8,5 | Mecanismo mejor delimitado; falta confirmación independiente |
| Solidez metodológica | 8,7 | 8,7 | Descomposición factorial y recuento independiente |
| Evidencia computacional | 8,7 | 8,7 | Reutilización explícita del mismo panel |
| Utilidad discriminante demostrada | 6,5 | 6,5 | Reparación parcial con catálogo oracular; queda un error |
| Validación biológica independiente | 3,0 | 3,0 | Sin nuevas mediciones moleculares |

Los tres bloques se descargaron por CLI, se verificaron y se guardaron en Drive antes de continuar; los archivos pesan aproximadamente 0,58, 1,14 y 1,70 MB. La lectura completa de las copias de Drive no se ha verificado: no se afirma identidad criptográfica de esas copias ni se autoriza limpieza local por su sola presencia. Los archivos descargados desde Colab sí coinciden con los hashes emitidos por la ejecución. Se conserva el fallo inicial de lectura por prefijo del ZIP y su corrección anterior a las decisiones. La sesión propia se detuvo al terminar. Los manuscritos permanecen privados y no se regeneraron en esta ejecución.
