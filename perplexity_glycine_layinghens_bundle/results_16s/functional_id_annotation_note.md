# Functional ID Annotation Note

- `PWY-*` IDs in the PICRUSt2 pathway output are MetaCyc/BioCyc pathway identifiers, not KEGG pathway IDs.
- `EC:*` IDs are Enzyme Commission identifiers. These were checked against KEGG enzyme naming when possible, and cross-checked against IUBMB/ExPASy.
- Only IDs with externally verified names were annotated. Uncertain entries were intentionally left blank rather than guessed.
- `EC:3.6.3.54` was flagged as a transferred EC entry; the current accepted entry is `EC 7.2.2.8` (P-type Cu(+) transporter).
- `EC:3.4.24.11` is exactly mapped to neprilysin, but that signal is biologically suspicious in bacterial functional prediction and should be interpreted cautiously.

## Generated annotated files
- Pathway_kw_annotated.csv
- Pathway_dose_trend_annotated.csv
- Pathway_pc_mimicry_annotated.csv
- EC_kw_annotated.csv
- EC_dose_trend_annotated.csv
- EC_pc_mimicry_annotated.csv