# Comparación HiFi y Nanopore en los mismos 113 loci de HG003 y HG004

Se puntuaron los mismos 113 loci fijados para Nanopore con los BAM HiFi de chr20 distribuidos en el caso DeepTrio de PrecisionFDA v2. Los grupos de lectura declaraban HG003 y HG004. El diccionario chr20 contiene nombre y longitud, sin MD5. Se aplicaron las plantillas GRCh38 verificadas de la ejecución anterior y las mismas reglas de selección, puntuación y abstención.

| Donante y tramo | Nanopore C/D/U | HiFi C/D/U | HiFi tras mayoría exacta C/D/U |
|---|---:|---:|---:|
| HG003 largo | 14/7/2 | 23/0/0 | 19/0/4 |
| HG003 corto | 32/0/0 | 32/0/0 | 32/0/0 |
| HG004 largo | 18/8/0 | 26/0/0 | 24/0/2 |
| HG004 corto | 31/0/1 | 32/0/0 | 32/0/0 |

HiFi concordó con GIAB en los 113 loci, incluidos los 15 discordantes y los tres sin resolver en Nanopore. En el locus `chr20:52989164` de HG004, 16 de 20 lecturas HiFi mostraron las 18 copias esperadas y cuatro mostraron 19. En Nanopore, 11 mostraron 17 y ninguna mostró 18. La mayoría exacta conservó las dos llamadas centrales, que eran discordantes entre bibliotecas. El contraste comprueba genotipos en las mismas muestras y loci. La tecnología, el llamado de bases y el alineamiento previo varían conjuntamente.

La ejecución está vinculada al mismo hash de selección y de código que Nanopore. Se comprobaron los 113 loci y se reconstruyeron ambas decisiones de cada locus a partir de los puntajes y recuentos guardados. Fuente primaria de los archivos: https://github.com/IntelLabs/open-omics-deepvariant/blob/r1.5/docs/deeptrio-pacbio-case-study.md#download-hg002-hg003-and-hg004-bam-files .
