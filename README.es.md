# CRISPR Observability

[English](README.md) · [Abrir el entorno de investigación](https://facundofirmenich.github.io/crispr-observability/)

**Seguir la afirmación. Examinar la evidencia. Cambiar parámetros. Reproducir el resultado.**

Una edición puede desplazar fuera de una lectura finita el ancla necesaria para
detectarla. Dos reglas pueden coincidir en una etiqueta y extraer bases distintas.
CRISPR Observability hace visibles esos mecanismos y certifica objetivos dentro
de dominios de interpretación explícitos. El método se acompaña de un entorno
de evidencia y un laboratorio computacional reutilizable.

| Entrada | Capacidad | Utilidad científica |
| --- | --- | --- |
| [Epistemic Supergraph](interactive/supergraph/crispr_supergraph.html) | Recorrido conceptual, premisas, fuentes, tarjetas/columnas móviles, relaciones iluminadas | Entender y contrastar de qué depende una conclusión |
| [Read Observability Lab](interactive/lab/CFTR_OBSERVABILITY_LAB.html) | Secuencias y estructuras 3D, anclas, tramos, SAM/BAM, importación, JSON/replay/PNG | Ensayar intervenciones computacionales y conservar resultados |
| [Reproducción](reproducibility/LEEME_ES.md) | Autoejecutables Windows/macOS/Linux, quick/counts/full, Streamlit, Docker y cuaderno | Recalcular con alcance y ledger explícitos |

Descargue los HTML mediante **Download raw file** de GitHub o use la portada.
El laboratorio requiere WebGL 2; no necesita enviar sus archivos a un servidor.
El explorador epistémico numérico y PaperGraph ofrecen vistas complementarias,
con sus respectivos alcances. Las geometrías ideales y las coordenadas importadas
conservan su procedencia; el visor no predice dinámica molecular ni eficiencia.

## Empezar con evidencia real

El [caso guía CFTR](docs/LEADING_CASE.es.md) conecta dos lecturas auténticas,
una intervención geométrica, siete operaciones y resultados exactos. La
[guía de investigadores](docs/RESEARCHER_GUIDE.es.md) explica importación, MCP,
exportación, reproducción y contribuciones. El [motor multientrada](engine/README.md)
ejecuta manifiestos con hasta 32 entradas y 64 operaciones por CLI, API y MCP;
el visor mantiene un caso activo.

La infraestructura se utiliza o adapta en **CRISPR Observability**, **Biogen /
BioGen S1 —observabilidad local del splicing—** y **ChronoWitness —biología
matemática y dinámica temporal local—**. [Sus estados](docs/RESEARCH_LINES.es.md)
son diferentes: laboratorio nativo, motor de uniones aplicado y adaptador temporal
separado. La integración del último con el visor común sigue pendiente. Son
tres líneas del programa, no tres adopciones externas independientes.

## Acceso y reproducción

Con Python 3.11+ y Git, en un destino persistente adecuado:

```sh
git clone https://github.com/FacundoFirmenich/crispr-observability.git
cd crispr-observability
python -B reproducibility/reproduce.py --mode quick --output replay-quick
```

Para abrir sin instalar Python, use la [release nativa v0.9.0](https://github.com/FacundoFirmenich/crispr-observability/releases/tag/v0.9.0):
`START_HERE.bat` en Windows x86-64, `START_HERE.sh` en Linux x86-64 o
`START_HERE.command` en macOS Intel/Apple silicon. El recibo `REPRODUCTION.json`
declara qué se ejecutó. El dashboard Streamlit y el contenedor Docker están
documentados en la guía. No se inicia Colab ni Kaggle.

El desarrollo local utilizó un **portátil convencional de 2016, Intel Core
i5-6200U a 2,30 GHz**, según informa el autor. Las campañas grandes usaron
cómputo alojado. Esto describe el hardware, no establece que toda campaña
completa funcione en cualquier equipo antiguo.

## Resultados y límites

| Hallazgo | Resultado | Alcance |
| --- | --- | --- |
| Censura por geometría | Ancla [213,225) en lectura de 224; extremo inferior CFTR −10,09177 → +13,23196 pp | Mismo universo y máscara; el intervalo anterior contiene cero |
| Etiqueta frente a cadena | 83 cadenas desiguales entre 170.808 registros conjuntos; 18 con etiqueta nativa coincidente | Contraejemplo de correspondencia; no 83 errores moleculares adjudicados |
| Certificación local | Un objetivo se resuelve si todos los candidatos del dominio completo dan el mismo valor determinado | Certificado relativo al modelo y abstención explícita |
| Coste de retirar llamadas | Regla germinal de cuatro donantes: retira 25/26 discordancias y 113/207 concordancias en 241 loci | Concordancia por locus bajo selección y regla congeladas |

Las discrepancias locales y resultados adversos se conservan. Consulte el mapa
de evidencia y las instrucciones de reproducción antes de interpretar una cifra.

## Versiones y contribuciones

El DOI [10.5281/zenodo.24121127](https://doi.org/10.5281/zenodo.24121127) identifica
v0.8.1. Los paquetes v0.9.0 y el desarrollo posterior tienen sus propios commits.
Las capacidades nuevas no se atribuyen retrospectivamente a esos archivos.
Manuscrito y suplemento permanecen privados para el depósito del autor.

La [hoja de ruta comunitaria](docs/COMMUNITY_ROADMAP.es.md) invita a casos nuevos,
reproducciones independientes, adaptadores verificados y correcciones científicas.
Use Issues para fallos reproducibles, Discussions para preguntas y pull requests
para mejoras contrastadas. Aporte commit, fuente, coordenadas, comando y resultado.
Se informará adopción independiente cuando exista evidencia de ella.

El software original actual usa **GPL-3.0-only**; la extensión independiente
`extensions/generation_model/` conserva **MIT**. Se mantienen licencias previas,
avisos de terceros y derechos de los datos y publicaciones originales.
