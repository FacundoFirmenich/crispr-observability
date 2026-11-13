# Del contrato de observación a las decisiones locales

La ampliación comparativa del 13 de noviembre de 2026 incorpora ejecuciones nativas nuevas, transiciones de admisión y una cartografía completa del coste de decisión en un panel finito. Conserva las versiones anteriores y todos los regímenes adversos. Son experimentos computacionales retrospectivos sobre datos públicos autenticados.

## Qué cambia para el investigador

**Un porcentaje estable puede ocultar una población distinta.** Sobre los mismos 175.010 registros Q30 de nueve bibliotecas CFTR, ampliCan 1.32.1 admite 112.822 con el anclaje computacional [0,20), 174.439 con [4,24) y 174.510 con [8,28). Estos anclajes no se presentan como los cebadores experimentales originales. Su desplazamiento también cambia la frontera de análisis que excluye el cebador; por tanto, modifica conjuntamente admisión y dominio de observación.

**Un cambio neto pequeño puede ocultar transiciones opuestas.** El segundo desplazamiento incorpora 209 registros y excluye 138: ganancia neta de 71. La correspondencia utiliza familias de secuencias observadas completas y sus multiplicidades invariables, dentro de cada biblioteca. No demuestra identidad molecular original ni una ampliación anidada de la evidencia.

**Una política local puede ayudar en un régimen de costes y perjudicar en otro.** Cada uno de cuatro donantes se evalúa con decisiones aprendidas únicamente de los otros tres. Con coste seis por discordancia y uno por abstención, la política local conserva 124 concordancias y cuatro discordancias y deja 113 observaciones sin resolver. Su pérdida es 137, frente a 164 de la retención inicial, 152 de la abstención fija y 165 del enrutamiento agregado. La enumeración racional completa delimita la ventaja frente a esos tres comparadores en 17/3 < coste <= 20/3, exclusivamente para este panel. Los costes adversos se conservan.

La conexión científica une geometría de observación, identidad de la evidencia y decisión selectiva. No constituye un ranking universal de exactitud ni una política óptima poblacional. La exactitud biológica sobre el fragmento original conserva su límite y queda fuera de los requisitos de cierre de esta versión.

## Reproducir e inspeccionar

Desde la raíz del repositorio:

```sh
python extensions/comparative_execution/check_results.py
python extensions/comparative_execution/local_policy.py --output /destino/politica-nueva.json
python extensions/comparative_execution/cost_frontier.py --output /destino/frontera-nueva.json
python extensions/comparative_execution/collect_native.py --output /destino/anclajes-nuevos.json
```

Los destinos deben ser nuevos. Los resultados pequeños están en `extensions/comparative_execution/results/`; los archivos nativos completos permanecen en la release `comparison-evidence-43493655364`. El colector lee los archivos de ampliCan en memoria. Los contenedores se ejecutaron en CPU de GitHub Actions. Los tiempos pertenecen a etapas concretas y no comparan velocidades de herramientas con tareas distintas.

CRISPResso2 completó nueve ejecuciones nativas y la reconstrucción de sus registros. CRISPR-A completó nueve llamadas de variantes tras una reparación limitada a serialización; la reparación posterior de la interfaz se verificó en nueve reejecuciones, conservando idénticas las seis tablas de control. Se conserva la auditoría en `NATIVE_EXECUTION_ADJUDICATION.md`.

## Usar archivos propios

En Streamlit, desplegá **Analice su FASTQ**, cargá tratamiento y control, revisá la configuración editable de anclajes y generá el informe HTML/JSON. El procesamiento se realiza en memoria del equipo que aloja la aplicación; un alojamiento remoto recibe los archivos. La interfaz rechaza archivos mayores de 20 MB o 250.000 registros en vez de truncarlos. Para conjuntos mayores está la CLI en streaming. El observador por anclajes y la certificación del dominio completo de alineamientos mantienen contratos distintos.

El Supergraph incorpora las jurisdicciones favorables y adversas y las transiciones de admisión. El grafo numérico y PaperGraph preservan sus identidades históricas; sus documentos acompañantes explican la ampliación actual. La interpretación comienza por definir qué evento se desea medir.
