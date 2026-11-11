# Réplica en HG003 y HG004

Se ejecutaron 113 loci fijados antes de leer los resultados de secuenciación. Se excluyeron los 128 loci de la campaña anterior y los loci de desarrollo ya registrados. Los controles de referencia se buscaron primero en HG002, cuyo catálogo no aportó controles aislados para HG003. Se usó entonces el catálogo GIAB de HG005, de otro trío. Antes de puntuar se fijó un máximo de 16 por estrato y se conservaron todos los controles elegibles si eran menos. Los controles largos fueron siete en HG003 y diez en HG004. No se sustituyeron loci tras ver lecturas o puntajes.

| Donante | Longitud mayor | Concordantes / discordantes / sin resolver | Discordancias entre llamadas emitidas | Aciertos retirados por mayoría exacta | Errores retirados |
|---|---|---:|---:|---:|---:|
| HG003 | ≥15 | 14 / 7 / 2 | 7/21 | 13/14 | 7/7 |
| HG003 | <15 | 32 / 0 / 0 | 0/32 | 13/32 | No estimable |
| HG004 | ≥15 | 18 / 8 / 0 | 8/26 | 18/18 | 7/8 |
| HG004 | <15 | 31 / 0 / 1 | 0/31 | 4/31 | No estimable |

En la réplica hubo 15 discordancias entre 47 llamadas en tramos largos y cero entre 63 en cortos. Al juntar los paneles disjuntos de los cuatro donantes como recuento descriptivo, resultan 23/107 y 3/126, respectivamente. Los donantes corresponden a dos tríos y los paneles tienen cupos deliberados por longitud y genotipo. El agrupamiento no estima una tasa poblacional.

En los cuatro paneles la mayoría exacta retiró 25 de 26 errores iniciales y 113 de 207 aciertos. En los 113 loci de tramo largo dejó una llamada concordante, una discordante y 111 abstenciones. La discordancia conservada fue `chr20:52989164 C>CA` en HG004. El genotipo esperado es el alternativo de 18 copias. Entre las 20 lecturas seleccionadas, 11 mostraron exactamente las 17 copias del alelo de referencia y nueve tuvieron otra longitud. La mayoría exacta no evita una medición de longitud discordante que se repite en la mayoría de las lecturas.

Se verificaron las identidades de archivo y referencia, los 113 loci de esta réplica y los 128 previos, sus puntajes y la reconstrucción de cada decisión. La réplica aportó 1.100 pares lectura–locus de HG003 y 1.153 de HG004. Los datos por lectura y el agregado están en `HOMOPOLYMER_HG003_HG004_RESULT.json` y `HOMOPOLYMER_REPLICATION_COMPARISON.json`. La ejecución anterior conserva su código fijado en el commit `69601a1`.
