# Reproducción científica sin conexión

Extraé el paquete correspondiente a tu sistema. Conservá las carpetas de
evidencia y del entorno junto al ejecutable. Desde una terminal en esa carpeta:

- Windows (PowerShell): `.\crispr-reproduce.exe --output .\replay-quick`.
- Linux y macOS: `./crispr-reproduce --output ./replay-quick`.

El paquete ya incluye Python. El directorio de salida debe ser nuevo.
Comprobá el archivo descargado con la suma SHA-256 de su release. Consultá
`REPRODUCTION.json` para ver qué se ejecutó, con qué entradas y hasta qué
nivel de evidencia.

- `--mode quick` comprueba integridad y grafo, vuelve a resolver el caso
  auténtico de Clone3 y reconstruye los intervalos exactos de CFTR y las
  decisiones germinales a partir de las observaciones archivadas.
- `--mode counts` recompone los recuentos de los 439.589 pares COL17A1.
- `--mode full` añade el recálculo de 20.260 casos referencia–segmento,
  81.112 testigos y la comprobación con el oráculo independiente acotado.

La comparación germinal mantiene los denominadores, discordancias y
abstenciones por donante y estrato. Sus puntajes y llamadas de flanco son
entradas archivadas. La exactitud de edición por fragmento original
requiere observaciones enlazadas independientes.
