# Reejecución y alcance de la sensibilidad de alineamiento

Este directorio contiene código y resultados propios. El motor CRISPResso2, su fuente, su entorno virtual y sus binarios se obtuvieron privadamente del upstream oficial y se excluyen de la entrega. Consultar su licencia antes de cualquier uso. No se modificó el motor.

Identidad requerida: etiqueta `v2.3.4`, commit `1fc4b629e5400a48d6569b0ab4b8d1766f33636f`. La licencia leída tiene SHA256 `183a65307ac1e05357c3129275c34ab45c29034e29fe2fe80fb3c9ab353a2fa3`. `ENGINE_ENVIRONMENT.json` registra versiones y coincidencia de los archivos fuente con el entorno histórico. La construcción usó el C generado por upstream, `CC=gcc`, y un entorno virtual separado con las dependencias científicas compartidas; no instaló Cython.

Desde la raíz del proyecto, preparar el motor en el entorno autorizado:

```bash
git clone --depth 1 --branch v2.3.4 https://github.com/pinellolab/CRISPResso2.git _mp02_engine/CRISPResso2
git -C _mp02_engine/CRISPResso2 rev-parse HEAD
python3 -m venv --system-site-packages _mp02_engine/env
_mp02_engine/env/bin/python -m pip install jinja2==3.1.6
CC=gcc _mp02_engine/env/bin/python -m pip install --no-deps --no-build-isolation _mp02_engine/CRISPResso2
_mp02_engine/env/bin/CRISPResso --version
git -C _mp02_engine/CRISPResso2 diff --exit-code
```

Las entradas son las nueve corrientes FASTQ originales de `publication_verification/P1/advance_m7/replay/project/evidence/v22/v21/sources/fastq/`. Las referencias, diseño histórico y resultados congelados se leen en `restored/m2/advance_m2/`. No se cambian esos archivos.

`GRID.json` y `GRID_FREEZE.json` fijaron antes de los nuevos resultados la base personalizada `(-20,-2,60)` y seis variantes de un parámetro por vez: apertura `-10/-30`, extensión `-1/-4`, puntuación mínima `50/70`. Las dos referencias por nueve corridas por siete configuraciones producen 126 trabajos completos del motor; las ventanas permanecen fijas en esta grilla. `STOCK_PROFILE.json` y su congelación agregaron, por separado y antes de su ejecución, 18 trabajos del perfil sin guía con opciones de alineamiento/ventana omitidas para usar sus defaults upstream. Este perfil cambia varias opciones y no se incluye en los rangos OAT. La base personalizada no se denomina default estándar.

La preparación original se realizó con:

```bash
python3 publication_02/sensitivity/run_grid.py prepare
python3 publication_02/sensitivity/run_profiles.py prepare
```

Esos pasos conservan las primeras congelaciones y rechazan sobrescribirlas. En una reproducción nueva, usar una copia sin outputs, conservar las congelaciones originales como antecedentes y registrar una nueva ejecución sin borrar resultados previos. Los comandos exactos por trabajo quedan en `executions/`; sus rutas absolutas se deben adaptar a la ubicación del proyecto. El entorno original reconstruye las entradas Q30 deterministas y exige identidad del flujo secuencia/calidad raw, recuentos raw y Q30 coincidentes con los registros congelados.

La secuencia de ejecución y comprobación es:

```bash
python3 publication_02/sensitivity/run_grid.py pilot --workers 1
python3 publication_02/sensitivity/checked_run.py baseline --workers 3
python3 publication_02/sensitivity/summarize_grid.py baseline
python3 publication_02/sensitivity/staged_run.py stock --workers 3
python3 publication_02/sensitivity/staged_run.py oat --workers 3
python3 publication_02/sensitivity/summarize_grid.py all
```

El método endurecido para nuevas ejecuciones es el de `staged_run.py`: escribe cada salida en `/tmp/mp02_engine_jobs`, espera el cierre del motor, verifica gzip/zip e identidad de todas las lecturas Q30, sincroniza archivos, copia a una carpeta temporal cerrada y verifica hashes antes y después del renombrado final. Los registros de trabajo incluyen los hashes de cada salida. `checked_run.py` autentica antes y después el SHA binario y el flujo secuencia/calidad Q30. El perfil baseline original tuvo dos fallos por una entrada local truncada y una salida completada pero truncada; todos se conservaron en `failed_attempts/` y se reemplazaron para el análisis mediante reejecución verificable. `INPUT_REPAIR.json` exige que la regeneración de la entrada tenga exactamente el SHA de su primer manifiesto. La causa de las truncaciones no está establecida.

`summarize_grid.py` compara los 18 baselines con los recuentos y las tablas nativas congeladas. Para cada trabajo, reconstruye los eventos desde las bases observadas alineadas, mantiene fallos de calidad/alineamiento/cobertura como U dentro de N raw y coteja dos rutas de recuento: FASTQ anotado y tabla de alelos detallada. Guarda configuraciones efectivas, `include_idxs`, `exclude_idxs` y vector efectivo de incentivo de gaps. También conserva transiciones individuales contra baseline, firmas conjuntas por lectura y ejemplos de conflictos. Los intervalos son envolventes de identificación `E/N` a `(E+U)/N`; no son intervalos de confianza.

La frecuencia nativa `MODIFIED` del software queda separada de los eventos de cadena programada exacta y péptido local nativo. La prueba mide sensibilidad computacional finita, retrospectiva y en la misma genealogía experimental. No aporta verdad molecular independiente, datos biológicos nuevos, eficacia ni validación poblacional.
