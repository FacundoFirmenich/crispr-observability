# Caso guía: una base ausente cambia lo que podemos observar

[English](LEADING_CASE.md) · [Guía para investigadores](RESEARCHER_GUIDE.es.md)

El recorrido conecta el límite físico de una lectura, un fallo de correspondencia
de secuencias y una intervención computacional auditable. Utiliza dos lecturas
CFTR auténticas y las coordenadas de anclaje del análisis conservado. Las fuentes,
transformaciones e identidades están en
[`PROVENANCE.json`](../examples/cftr-leading-case/PROVENANCE.json).

## 1. Observar el mecanismo

Abra [Read Observability Lab](../interactive/lab/CFTR_OBSERVABILITY_LAB.html) en un
navegador con WebGL 2 y seleccione **Read boundary**. El ancla [210,222),
desplazada +3 posiciones dentro de una lectura que termina en 224, pasa a
[213,225): solamente se observan 11 de sus 12 posiciones. Con desplazamiento
cero vuelve a estar contenida. Restaure +3 y aumente el final a 225: también se
recupera la contención.

El último cambio es contrafactual y geométrico: no crea una base medida ni
demuestra su identidad o calidad. En el análisis CFTR completo, cambiar la regla
de observación bajo una máscara de calidad común cambia un extremo inferior de
contraste de −10,09177 a +13,23196 puntos porcentuales. Ese contraste empírico
conservado es distinto del cálculo geométrico y no es una eficiencia predicha por
la escena tridimensional.

## 2. Examinar las bases auténticas

En **Competing spans**, elija `SRR29623270.2091`: [114,156) frente a [114,201)
produce cadenas de 42 y 87 bases; 45 posiciones pertenecen exclusivamente al
tramo largo. Las etiquetas locales también difieren en este caso.

Para observar acuerdo de etiquetas con desacuerdo de cadenas, elija
`SRR29623271.868`: [114,186) frente a [114,153), 72 frente a 39 bases. Las dos
etiquetas nativas conservadas son falsas, pero las cadenas son distintas. La
proyección óptima alternativa [114,186) recupera la igualdad de cadenas con la
primera regla. Elegir un traceback y certificar todos los candidatos de un
dominio completo son operaciones diferentes.

**Atomic model** permite girar, seleccionar y medir una geometría ideal asociada
a la secuencia. Las coordenadas experimentales importadas conservan su propia
procedencia. Ninguna de esas representaciones demuestra que la conformación
visualizada estuviera presente en el material secuenciado.

## 3. Ejecutar y conservar

Desde la raíz del repositorio, con Node.js 22 o posterior:

```sh
npm ci --prefix engine --no-fund --no-audit
node engine/project-cli.mjs check examples/cftr-leading-case/project.json
node engine/project-cli.mjs run examples/cftr-leading-case/project.json > leading-case-result.json
node examples/cftr-leading-case/verify.mjs
```

Elija un destino persistente apropiado para su equipo. La instalación necesita
conexión; el ejemplo se ejecuta después sin red. `check` valida el contrato;
`run` lee los archivos, comprueba sus SHA-256 y calcula siete operaciones. El
último comando verifica los resultados esperados.

| Operación | Resultado esperado | Significado |
| --- | --- | --- |
| Ancla sin desplazar | 12/12; contenida | Condición geométrica original |
| +3; final 224 | 11/12; una posición fuera | Pérdida dependiente del producto |
| +3; final 225 | 12/12; contenida | Contrafactual declarado |
| Primera lectura auténtica | 42 frente a 87; distintas | Desacuerdo de secuencia observada |
| Segunda lectura auténtica | 72 frente a 39; distintas | La etiqueta no identifica la cadena |
| Óptimo alternativo | Cadenas de 72 bases idénticas | Importa el candidato; falta certificar el dominio completo |

## 4. Relacionar el resultado con la afirmación

En el [Supergraph](../interactive/supergraph/crispr_supergraph.html), siga el
recorrido conceptual: geometría, correspondencia y cotas condicionadas. Examine
las premisas y fuentes. El explorador epistémico numérico permite consultar los
recuentos. Las aristas declaran dependencias; la evidencia procede del cálculo
y de los registros enlazados.

## 5. Probar una idea nueva

Copie el manifiesto en un estudio nuevo, conserve el original y cambie una
operación identificada. Declare si modifica una medición, una regla de análisis
o un contrafactual. Mantenga coordenadas comparables, resultados irresueltos y
denominador original. Exporte JSON/replay junto al PNG: una imagen aislada no
reproduce la intervención. Este caso demuestra geometría, correspondencia local,
procedencia y repetibilidad computacional; la exactitud de edición de la molécula
original exige una observación independiente enlazada a esa unidad.
