# CRISPR-A template export compatibility repair

The repaired staging run 43494837211 reaches variant calling but fails in `write.csv(t_info[3:length(t_info)], ...)`: its template result is a nonrectangular list (observed row lengths 1, 2, 895). The failure is a CSV serialization error, not a numerical accuracy result. The unmodified native run and its log remain evidence.

A separate compatibility execution replaces only this export statement. It saves the complete `t_info` object as RDS and writes its exported components as a long table with component/index/value columns. It does not change `t_info`, template counts, labels, alignments, filters or later computations. The original and patched script hashes and exact replacement are archived. Any completed results are attributed to the pinned upstream version plus this explicit export-only patch, never to an unmodified upstream success. The pipeline's downstream process declarations do not consume `_template-reads.csv`.
