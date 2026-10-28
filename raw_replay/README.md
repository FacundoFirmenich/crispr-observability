# P1: replay mínimo desde los nueve FASTQ auténticos

Requiere Python 3.11 o posterior y su biblioteca estándar. No requiere instalación,
red, CRISPResso, artículos, archivos históricos ni paquetes de entrenamiento.

```sh
python run_raw_replay.py
```

El destino predeterminado es `results/generated/`. Se puede elegir otro con
`--output /ruta/nueva`. El runner conserva los FASTQ y las fuentes congeladas.
No lee los resultados esperados hasta terminar la reconstrucción científica.

Se reconstruyen y contrastan con resultados originales tres productos:

1. `M1_TARGET_COUNTS.json`: los dos eventos del observador target de M1 en los
   nueve streams primarios; trims 0/8 y cuatro políticas estructurales/de calidad
   (144 filas). Las políticas se conservan separadas.
2. `RAW_UNIVERSE_IDENTITIES.json`: nueve tamaños, hashes de los archivos FASTQ
   auténticos y hashes lógicos ordenados de secuencia+calidad. Además se verifica
   que los identificadores de lectura sean únicos dentro de cada corrida.
3. `WITNESS_BOUNDS.json`: 18 filas de recuentos, 12 contrastes con aritmética
   racional exacta, nueve dominios de correspondencia y las abundancias de ambas
   tablas congeladas. Se mantienen 177.006 registros, 170.725 testigos y 6.281 U.

Cada registro raw se observa en ambas orientaciones. El target exige un único par
de motivos observados; la proyección por coordenadas utiliza los alineamientos
preexistentes de las dos referencias, con flancos inmediatos observados. Sólo
se acepta como testigo la misma cadena codificante en las tres representaciones.
Las abundancias de todos los alelos deben corresponder exactamente a lecturas raw.
U conserva tanto fallos de calidad como ausencia o discordancia de correspondencia.

`INPUT_MANIFEST.json` identifica cada entrada, su función, procedencia y alcance
de derechos. Una alteración de secuencia/calidad, ZIP, referencia o resultado
esperado detiene el runner por identidad; el recuento raw, los hashes lógicos,
la conservación de abundancias y los resultados originales son comprobaciones
adicionales. Estos hashes no son firmas criptográficas ni evidencia independiente
de verdad molecular. Cambiar deliberadamente código, manifiesto y goldens juntos
queda fuera de esa garantía de integridad.

El runner sí recalcula observadores desde raw y proyecciones sobre alineamientos
congelados. No reejecuta el motor de alineamiento, no reconstruye entrenamientos
históricos, no aporta otro humano ni otra cohorte y no calibra verdad molecular,
función o eficacia terapéutica. Esta implementación mínima adapta los mismos
contratos históricos y los compara con sus resultados originales; no representa
una validación experimental independiente ni la intervención de otro analista
humano.

Los FASTQ y tablas son evidencia de fuente y productos derivados; el aviso de
datos separado prevalece sobre cualquier licencia de código del repositorio
integrador. Este compañero permanece separado hasta su integración expresa.

## Comprobación de sensibilidad de la verificación

```sh
python check_tamper_detection.py
```

Esta comprobación usa copias temporales. Cambia una base de un FASTQ auténtico,
un recuento de una tabla congelada y un valor esperado, y exige rechazo. También
actualiza el hash de transporte de la tabla alterada en la copia temporal para
comprobar que la conservación de abundancias detecta un fallo semántico que ya
no es un simple desacuerdo de SHA. Nunca modifica los insumos del paquete.
