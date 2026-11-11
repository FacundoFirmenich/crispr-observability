# MP36: resolución local, mecanismo y anidaciones conservadas

El proyecto queda mejor posicionado en resolución de genotipos del panel y, sobre todo, en explicación del mecanismo. La ejecución completa de MP36 terminó con 64 loci de evaluación en 26 regiones separadas de las 30 regiones de calibración, sin identificadores de lectura compartidos entre particiones. La regla primaria resolvió 63/64 loci concordantemente con el benchmark, frente a 18/64 del consenso transportado, sin llamadas discordantes por locus. El criterio descriptivo congelado pasó. Esto se refiere a variantes simples homocigotas de un individuo y una biblioteca; no establece precisión genómica general ni verdad molecular de HEK3/CFTR.

La adjudicación local impide atribuir indebidamente todo el salto a las medianas. Los 45 loci recuperados contienen indels de 2–4 nucleótidos. A 300 nt, cambiar el umbral relativo de cuatro a cero, manteniendo el consenso del 80%, ya resuelve 62/64; la mediana añade un locus. En ese locus, chr20:48196173 GGA>G, catorce lecturas favorecen el genotipo esperado y seis la alternativa. La mediana concuerda con el benchmark, pero no elimina esa discrepancia. En otros 20 loci concordantes también hay preferencias contrarias: son 35 preferencias contrarias dentro de 21 loci concordantes. No se presentan como errores moleculares autenticados. El único locus no evaluable, chr20:25434443 T>TTTTTC, tiene seis lecturas frente al mínimo fijo de diez; su estado abierto se conserva.

La granularidad aporta una observación adicional que el resumen primario no muestra: en chr20:40147874 G>GGA, ampliar el contexto de 150 a 300 nt aumenta el soporte concordante de 15/20 a 16/20 y permite alcanzar el consenso del 80%, aunque la mediana no cambie. La auditoría de anidaciones conserva 1.254 registros como universo común, sus 114 conjuntos nombrados y las intersecciones especificadas. Las ventanas anidadas no se cuentan como réplicas. En doce grupos de resúmenes coincidentes, contados por longitud, iguales conteos, etiquetas y medianas coexisten con diferentes cantidades de preferencias contrarias o empates. Es evidencia observada de pérdida de detalle al resumir, no una demostración de superioridad universal de una arquitectura ni una estimación de información de Fisher. El diagnóstico de conjuntos se añadió después de conocer los resultados y queda marcado como retrospectivo.

Las 10.064 distancias almacenadas coinciden con una implementación independiente; la segunda extracción desde la fuente reprodujo 32 lecturas de 16 loci seleccionados por hash. Se verificaron los cuatro brazos primarios, las particiones, la aritmética local, las 2.508 filas de ventanas y los 153 archivos del paquete científico completo. Los resultados adversos de MP34 y MP35 permanecen: este resultado no transforma retrospectivamente sus criterios fallidos ni su dependencia en éxitos. Las garantías condicionales por jurisdicción, la replicación entre individuos o bibliotecas y la verdad enlazada al mismo fragmento siguen abiertas. La cobertura empírica 63/64 no autoriza una garantía individual para cada locus.

La consecuencia para el paper es una afirmación más fuerte y más precisa: podemos localizar una recuperación real de resolución, separar sus causas operacionales y conservar la evidencia minoritaria que la conclusión agregada no reconstruye. Manuscrito y suplemento ingleses incorporan las nuevas secciones; las nuevas secciones de resultados y métodos también están completas en castellano. La traducción íntegra del núcleo histórico al castellano y la adaptación de todas las figuras históricas siguen pendientes y no se presentan como terminadas. Los manuscritos permanecen fuera del repositorio. El siguiente paso científico decisivo es validar en datos reservados cualquier regla local nueva y obtener soporte independiente para las jurisdicciones actualmente no evaluables, sin relajar retrospectivamente sus umbrales.

## Comparación pareada sobre los mismos 64 loci

| Regla | Contexto | Concordantes | Discordantes | No resueltos |
|---|---:|---:|---:|---:|
| Consenso transportado, umbral relativo 4 | 150 nt | 18 | 0 | 46 |
| Consenso diagnóstico, umbral relativo 0 | 150 nt | 61 | 0 | 3 |
| Mediana primaria, umbral relativo 0 | 150 nt | 63 | 0 | 1 |
| Consenso transportado, umbral relativo 4 | 300 nt | 18 | 0 | 46 |
| Consenso diagnóstico, umbral relativo 0 | 300 nt | 62 | 0 | 2 |
| Mediana primaria, umbral relativo 0 | 300 nt | 63 | 0 | 1 |

La calibración de medianas no otorga calibración independiente al consenso diagnóstico de umbral cero. Cada brazo conserva su contrato.

## Scoring del avance, con alcance explícito

Rúbrica de madurez de evidencia: 0 = no establecido; 1 = parcial o diagnóstico; 2 = verificado en el dominio declarado. No son probabilidades de aceptación editorial ni una puntuación global de calidad.

| Dimensión | Antes de cerrar MP36 | Ahora | Jurisdicción y límite |
|---|---:|---:|---|
| Separación de lecturas calibración/evaluación | 0 en MP35 primario | 2 en MP36 | No equivale a independencia entre muestras biológicas |
| Resolución local de indels pequeños | 1 | 2 | 45 loci de 2–4 nt de este panel; no sensibilidad genómica |
| Atribución del salto a umbral frente a mediana | 0 | 1 | Descomposición secundaria sobre las mismas observaciones |
| Conservación de discordancias e intersecciones | 1 | 2 | Registro exacto y reconstruible; sin nueva replicación |
| Verificación computacional independiente | 2 | 2 | Distancias completas; segunda extracción de un subconjunto fijo |
| Cobertura condicional por locus | 0 | 0 | No establecida |
| Verdad independiente del mismo fragmento HEK3/CFTR | 0 | 0 | No aportada por HG005 |

## Artefactos y custodia

- `results/MP36_RESULT.json`: resultado primario intacto.
- `results/MP36_MECHANISM_AUDIT.json`: descomposición secundaria.
- `results/MP36_LOCAL_JURISDICTION_AUDIT.json` y `MP36_LOCAL_READ_SEGMENTS.csv`: evidencia por locus y ventana.
- `results/MP36_SET_INCLUSIONS_AND_PROJECTIONS.json`: pertenencias, inclusiones e intersecciones exactas; diagnóstico retrospectivo.
- `docs/MP36_CARTOGRAFIA_LOCAL.html`: explorador local completo.
- `docs/MP36_MAPA_EPISTEMICO.html`: grafo de afirmaciones, premisas y límites.
- `MP36_COMPLETE_AUDITED.zip`: 29,8 MB; lectura íntegra local verificada; Drive confirmó la carga. La lectura íntegra remota desde Drive no se realizó. No se borró ningún archivo.
