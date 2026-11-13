# Empiece por la pregunta y obtenga un informe de decisión

**Pregunta:** ¿estos registros resuelven un contraste positivo entre tratado y mock para el evento de secuencia que he definido?

En Streamlit, elija **Empezar aquí**. El ejemplo CFTR auténtico ya está cargado. Lea los cuatro pasos y descargue el informe HTML, el contrato de entrada y los resultados JSON. Los grafos y el laboratorio 3D quedan disponibles para profundizar; no son requisitos para entender el primer resultado.

Desde la terminal basta Python 3.10+ y el repositorio, sin instalar paquetes:

```sh
python code/observability_report.py demo --output /su/ubicacion/persistente/revision-cftr
```

1. **Defina** evento, referencia, coordenadas, calidad y dominio de candidatos. ADN programado, péptido local y modificación genérica son preguntas distintas.
2. **Conserve** el universo original: E + NE + U = N. Filtrar no permite cambiar silenciosamente el denominador.
3. **Interprete** las cotas: un extremo inferior positivo conserva el signo para todas las completaciones binarias de U bajo el supuesto de etiquetas aceptadas correctas. No es un intervalo de confianza.
4. **Elija** el siguiente contraste: flancos ausentes, cadenas ambiguas, sensibilidad a errores de etiqueta o evidencia independiente enlazada.

Para sus datos, edite el contrato descargado y ejecute:

```sh
python code/observability_report.py analyse su-estudio.json --output /su/ubicacion/persistente/revision
```

Esta entrada recibe recuentos ya clasificados: no transforma cualquier FASTQ en etiquetas ni certifica automáticamente todos los alineamientos. Valida estructura y aritmética; la procedencia de sus declaraciones requiere evidencia. No sobrescribe informes anteriores.

El ejemplo muestra además ocho jurisdicciones donante–longitud del contraste germinal existente. Hay que informar discordancias retiradas y concordancias sacrificadas. Sin llamadas retenidas, la discordancia retenida no es estimable. Para la pérdida explícita `L = c_D D + c_U U`, el agregado favorece la regla mayoritaria solo si `c_D/c_U > 5,52`, con coste positivo de abstención. Los costes los define la aplicación; no se ajustan para favorecer el método. Las diferencias locales se conservan.

El puente funcional CFTR recupera medidas de canal auténticas. La correlación con secuencias exige autenticar donante, cultivo/pool, condición, tiempo y biblioteca/ensayo; el rótulo «réplica 1» no demuestra por sí solo esa correspondencia.
