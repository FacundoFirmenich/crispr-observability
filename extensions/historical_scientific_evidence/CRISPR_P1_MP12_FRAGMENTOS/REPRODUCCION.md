# Reproducción privada de MP12

Descomprimir el ZIP de evidencia en un directorio vacío y ejecutar desde su raíz con Python 3.11 o posterior. Los datos adjuntos proceden de los depósitos originales y conservan sus términos. El árbol incluye los FASTQ de los 21 accesos, los dos FASTQ del control técnico y las tablas históricas estrictamente necesarias para clasificar los casos adversos. Las órdenes vuelven a escribir las salidas en el árbol extraído. Conservar el ZIP para cotejar con la ejecución archivada.

```sh
python -B CRISPR_P1_MP12_FRAGMENTOS/audit_fragments.py
python -B CRISPR_P1_MP12_FRAGMENTOS/verify_evidence.py
```

`FRAGMENT_AUDIT.json` contiene la criba de catálogos, la evidencia por acceso y el conteo de todos los cruces. `CLUSTER_LINKS.csv.gz` conserva las 54.429 coincidencias físicas de una misma observación con los dos identificadores de archivo. `READBACK_VERIFICATION.json` verifica cada cruce a partir de los FASTQ, comprueba el control pareado y resume los archivos publicados que se revisaron. Los archivos `*RECEIPT.json` de `sources/` documentan las fuentes recuperadas. Las consultas de adquisición adicionales se detallan en `SUPPLEMENTAL_SOURCE_RECEIPT.json`.
