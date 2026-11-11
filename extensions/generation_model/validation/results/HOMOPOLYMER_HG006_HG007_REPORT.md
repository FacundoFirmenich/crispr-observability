# Transferencia de HG002 a HG006 y HG007

Se decodificaron y puntuaron los 128 loci congelados en `HOMOPOLYMER_FROZEN_PANEL.json`. El análisis recuperó 1.263 pares lectura–locus de HG006 y 1.280 de HG007. Se comprobaron los hashes de la selección, el inventario, el código y la referencia, la correspondencia de los 128 loci y la reconstrucción de cada decisión desde los puntajes y los recuentos guardados. No hubo identificadores de lectura compartidos entre los dos donantes.

| Donante | Tramo del alelo mayor | Correctas / discordantes / sin resolver con la regla central | Errores retirados por mayoría exacta | Aciertos retirados por mayoría exacta |
|---|---|---:|---:|---:|
| HG006 | ≥15 | 25 / 5 / 2 | 5/5 | 25/25 |
| HG006 | <15 | 30 / 2 / 0 | 2/2 | 11/30 |
| HG007 | ≥15 | 27 / 3 / 2 | 3/3 | 27/27 |
| HG007 | <15 | 30 / 1 / 1 | 1/1 | 2/30 |

El estrato largo tuvo una proporción de discordancia mayor en ambos donantes. Entre llamadas emitidas fueron 8/60 en los tramos largos y 3/63 en los cortos. La asociación descriptiva observada en HG002 se repitió, con tres errores adicionales en tramos cortos. El umbral de 15 bases sirve como señal de riesgo en estos paneles y carece de alcance universal. Se usaron dos progenitores del mismo trío y paneles de loci distintos, seleccionados por longitud y genotipo.

La mayoría estricta de lecturas con longitud exacta del alelo llamado retiró 11/11 discordancias y 65/112 aciertos. El coste fue especialmente claro en tramos largos: quedaron 64/64 sin resolver, incluidos los 52 aciertos iniciales. En cortos quedaron 47 aciertos y 17 abstenciones. Entre las dos clases del estrato largo hubo seis discordancias en alternativos confirmados y dos en loci con alelo de referencia. El filtro pierde utilidad como veto general en estas bibliotecas.

Los alias de HG003 y HG004 quedaron enlazados con los nombres GIAB mediante la hoja de muestra, el resumen de adquisición y el grupo de lectura del CRAM. Coinciden identificador de corrida, flow cell y hora. La prueba figura en `HOMOPOLYMER_ALIAS_READBACK.json`. Se mantiene el inventario anterior porque su hash forma parte del panel congelado. La comprobación de alias habilita preparar una selección adicional con esos donantes. Los resultados presentes corresponden a HG006 y HG007.
