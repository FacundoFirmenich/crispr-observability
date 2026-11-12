# Evidencia local de splicing con Read Observability Lab v0.5

El motor multientrada añade `junction-evidence` para comprobar, en un SAM o BAM indexado ya generado, si dos o más uniones requeridas aparecen juntas en el alineamiento de un mate identificado. Cada resultado corresponde a la ventana, referencia y archivo declarados. La operación se ejecuta desde la CLI, la API JavaScript y el servidor MCP del motor de proyectos.

## Manifiesto mínimo

En este ejemplo, el SAM contiene los registros inspeccionados. Sustituí los valores por los del manifiesto congelado del estudio. El locus y las uniones usan coordenadas de base cero y fin excluido.

~~~json
{
  "schema": "observability-project/1",
  "project_id": "splicing-audit",
  "datasets": [
    {"id":"condition-a","type":"alignment","file":{"path":"inputs/condition-a.sam"}}
  ],
  "jobs": [
    {
      "id":"local-evidence-a",
      "op":"junction-evidence",
      "input":"condition-a",
      "locus":{"ref":"chr1","span":[90,160]},
      "read_name":"fragment-id",
      "read_group":"read-group-id",
      "mate":"first",
      "junctions":[[104,114],[118,138]]
    }
  ]
}
~~~

`mate` admite `first`, `last` o `unpaired`, según los flags SAM. Si omitís `read_group` y el nombre aparece en varios grupos del locus, el motor exige desambiguación. Para BAM, declará `type:"indexed-bam"` y los objetos `bam` e `index` con sus rutas. Podés fijar SHA-256 por entrada. Repetí el job sobre otra entrada únicamente con una correspondencia explícita de bibliotecas e identidades, conservada fuera de este operador.

~~~sh
npm ci --prefix engine --no-fund --no-audit
node engine/project-cli.mjs check /workspace/study/project.json
node engine/project-cli.mjs run /workspace/study/project.json > /workspace/study/report.json
~~~

## Lectura del informe

- `exactAnyReported` es `observed` si un único registro del mate seleccionado contiene todos los intervalos requeridos como operaciones CIGAR `N`. Cada alineamiento se comprueba por separado. Una deleción `D` conserva su identidad y queda fuera del recuento de uniones.
- `exactUniquePrimary` exige un registro testigo primario y `NH=1`. Devuelve `indeterminate-nh` ante un testigo primario con `NH` ausente o no informativo. Devuelve `indeterminate-primary` si hay primarios competidores para el mate en la ventana y al menos uno cumple la conjunción.
- `not-observed-in-input` describe el archivo o la ventana inspeccionada. La disponibilidad de registros, los filtros previos del archivo, el denominador de biblioteca y el significado de `NH` en el alineador son parte del contrato de cada estudio.
- El informe conserva flags, CIGAR, `NH`, `HI`, `RG`, referencia y posición del mate, `TLEN` y el MAPQ original. Para `MAPQ=255`, la calidad de mapeo evaluable se representa con `null`. El Phred de base se extrae por separado en cada flanco inmediatamente alineado de una unión y es `null` cuando falta. Cada coincidencia queda vinculada al registro que la contiene.
- Las reglas de cobertura por posición siguen en `alignment-locus`. Su filtro de MAPQ carece de autoridad sobre los endpoints de `junction-evidence`.

## Traspaso para un replay BIO

1. Fijá identidad compuesta de biblioteca, referencia, condición, fragmento, grupo y mate. Conservá la versión y configuración del alineador.
2. Compará las observaciones por identidad contra las decisiones congeladas. Registrá testigos persistentes, recuperados y perdidos. Cuando falte un insumo o una identidad verificable, registrá la abstención por fragmento.
3. Aplicá las reglas propias de BIO para contexto, calidad de flancos y denominadores. Conciliá cada desacuerdo con los registros originales antes de agregar. Los recuentos históricos BIO14–16 conservan su propia jurisdicción.
4. Versioná por separado la comprobación de integridad, la reconstrucción de tablas, el replay por fragmento y la ejecución desde fuentes crudas. La paridad con el adjudicador BIO requiere insumos de ese estudio y sigue siendo una prueba aparte.

Referencias de formato: [SAM/BAM v1.6](https://samtools.github.io/hts-specs/SAMv1.pdf) y [etiquetas SAM](https://samtools.github.io/hts-specs/SAMtags.pdf).
