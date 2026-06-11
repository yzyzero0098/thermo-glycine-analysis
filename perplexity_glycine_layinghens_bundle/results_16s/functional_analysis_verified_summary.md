# Functional Analysis Verified Summary

- `PWY-*` IDs are `MetaCyc/BioCyc pathway identifiers`, not KEGG pathway IDs.
- `EC:*` IDs are `Enzyme Commission` identifiers. These were re-checked against KEGG enzyme pages when possible and cross-checked against ExPASy/IUBMB.
- Uncertain or outdated identifiers were not force-labeled. Transferred EC entries were flagged explicitly.

## Top verified pathway-level KW hits

- `PWY-7219`: adenosine ribonucleotides de novo biosynthesis
- `ANAGLYCOLYSIS-PWY`: glycolysis III (from glucose)
- `PWY-621`: sucrose degradation III (sucrose invertase)
- `P122-PWY`: heterolactic fermentation
- `P42-PWY`: incomplete reductive TCA cycle
- `UDPNAGSYN-PWY`: UDP-N-acetyl-D-glucosamine biosynthesis I
- `PWY4FS-7`: phosphatidylglycerol biosynthesis I (plastidic)
- `PWY4FS-8`: phosphatidylglycerol biosynthesis II (non-plastidic)
- `PWY-5178`: toluene degradation IV (aerobic) (via catechol)
- `PWY-7208`: superpathway of pyrimidine nucleobases salvage

## Top verified pathway-level glycine dose trends

- `PWY-5384`: sucrose degradation IV (sucrose phosphorylase)
- `PWY0-1061`: superpathway of L-alanine biosynthesis
- `GLUCOSE1PMETAB-PWY`: glucose and glucose-1-phosphate degradation
- `GLYCOLYSIS-E-D`: Entner-Doudoroff pathway
- `P461-PWY`: hexitol fermentation to lactate, formate, ethanol and acetate
- `PWY-7371`: 1,4-dihydroxy-6-naphthoate biosynthesis II
- `PWY-7374`: 1,4-dihydroxy-6-naphthoate biosynthesis I
- `PWY-6728`: methylaspartate cycle
- `AST-PWY`: L-arginine degradation II (AST pathway)
- `PWY-7456`: beta-(1,4)-mannan degradation

## Top verified EC-level KW hits

- `EC:4.1.1.47`: tartronate-semialdehyde synthase
- `EC:3.6.3.54`: transferred entry; current accepted entry is `EC 7.2.2.8`, P-type Cu(+) transporter
- `EC:2.7.7.61`: citrate lyase holo-[acyl-carrier protein] synthase
- `EC:1.5.1.3`: dihydrofolate reductase
- `EC:2.1.1.45`: thymidylate synthase
- `EC:4.2.1.11`: phosphopyruvate hydratase (enolase)
- `EC:1.3.98.1`: dihydroorotate oxidase (fumarate)
- `EC:3.4.24.11`: neprilysin; exact EC match but biologically suspicious in bacterial PICRUSt2 prediction
- `EC:1.11.1.15`: transferred entry; do not treat as a stable current enzyme label

## Top verified EC-level glycine dose trends

- `EC:2.4.1.10`: levansucrase
- `EC:3.1.1.11`: pectinesterase
- `EC:3.1.3.4`: phosphatidate phosphatase
- `EC:3.1.3.81`: transferred entry to `EC 3.6.1.75`
- `EC:3.2.1.139`: alpha-glucuronidase
- `EC:4.2.2.10`: pectin lyase
- `EC:1.1.1.202`: 1,3-propanediol dehydrogenase
- `EC:1.1.1.205`: IMP dehydrogenase
- `EC:1.11.1.1`: NADH peroxidase
- `EC:2.4.1.7`: sucrose phosphorylase
- `EC:2.4.2.52`: triphosphoribosyl-dephospho-CoA synthase

## Output files

- `tables/Pathway_kw_annotated.csv`
- `tables/Pathway_dose_trend_annotated.csv`
- `tables/Pathway_pc_mimicry_annotated.csv`
- `tables/EC_kw_annotated.csv`
- `tables/EC_dose_trend_annotated.csv`
- `tables/EC_pc_mimicry_annotated.csv`
