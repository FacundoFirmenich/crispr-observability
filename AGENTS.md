# Conservación íntegra y recuperación de investigación

Instrucción expresa de Facundo Firmenich, ratificada el 1 de noviembre de 2026.

Conservar los archivos reales de investigación: fuentes, datos, configuraciones,
resultados, intermedios costosos, registros, intentos fallidos y entregas. Un hash,
un índice, un registro SQL, una cita conversacional o un enlace temporal no
sustituyen los bytes del archivo. Los hashes sirven para comprobar identidades
de archivos que efectivamente se han recuperado.

Antes de una ejecución costosa debe existir un destino persistente autorizado.
Guardar por etapas y antes del cierre. Los archivos propios de un repositorio
se sincronizan con su remoto; los materiales privados y los grandes resultados
que Git no incluya se guardan mediante el destino persistente autorizado. El
manuscrito y el suplemento CRISPR completos permanecen privados para el depósito
del autor. No publicar un resultado privado para resolver una falta de respaldo.

Al retomar, recuperar primero los originales desde sus ubicaciones persistentes.
Verificar recepción, apertura, tamaño e identidad. Para entregas científicas
costosas, comprobar también que los archivos enumerados dentro de los paquetes
están presentes y corresponden a sus inventarios. Distinguir copia temporal,
guardado remoto confirmado y descarga posterior realmente verificada.

Si falta un archivo, buscarlo en las versiones y respaldos autorizados. Si no se
recupera, declararlo FALTANTE. No presentar su hash como recuperación, no inventar
su contenido y no lanzar una recomputación para sustituir silenciosamente el
original. Una futura ejecución expresamente autorizada sería un nuevo resultado
con su fecha y procedencia; no se convertiría en el archivo histórico perdido.
La instrucción actual no autoriza recomputar la investigación como reparación.

No declarar una entrega completa sin confirmar cada guardado. Registrar por
archivo su ubicación persistente, identidad, estado y último acceso verificado.
Una fila de base de datos es un índice salvo que contenga también el archivo.
Si se conservan BLOBs en SQL, verificar que se pueden extraer sus bytes reales.
Conservar los fallos como fallos; la recuperación de un archivo defectuoso no
lo convierte en una ejecución científica exitosa.

No cambiar resultados, estados aprobatorios o fechas para aparentar continuidad.
La integridad de los archivos no constituye por sí sola una validación del método,
una medida de exactitud o una puntuación científica. Las limitaciones deben
revisarse sobre la evidencia íntegra, sin añadir barreras editoriales inventadas.

No incluir credenciales en respaldos. No prometer control sobre la limpieza de
proveedores ni activación automática en todas las conversaciones. Identificar
el alcance efectivamente instalado y los respaldos recuperables.
