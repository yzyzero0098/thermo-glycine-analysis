# EC enzyme follow-up summary

This summary reuses the saved PICRUSt2-derived functional tables and prioritizes features that support one or more of the following:
- omnibus group differences,
- monotonic dose trend,
- partial PC-like directionality/convergence.

Important limitation: these are **predicted functions from 16S-based PICRUSt2**, not directly measured microbial functional fluxes.

## Top omnibus-difference features
```csv
feature,display_name,kw_pvalue
EC:4.1.1.47,tartronate-semialdehyde synthase,0.0074754825567466
EC:2.7.7.61,citrate lyase holo-[acyl-carrier protein] synthase,0.0154063414916086
EC:1.5.1.3,dihydrofolate reductase,0.0167695133412348
EC:2.1.1.45,thymidylate synthase,0.0173331840150143
EC:4.2.1.11,phosphopyruvate hydratase (enolase),0.0173457333524049
```

## Top dose-trend features
```csv
feature,display_name,spearman_rho,spearman_pvalue
EC:3.1.3.4,phosphatidate phosphatase,-0.894427190999916,0.0027136820350938
EC:3.2.1.139,alpha-glucuronidase,0.894427190999916,0.0027136820350938
EC:2.7.7.61,citrate lyase holo-[acyl-carrier protein] synthase,-0.87287156094397,0.0046592149439939
EC:2.4.2.52,triphosphoribosyl-dephospho-CoA synthase,-0.87287156094397,0.0046592149439939
EC:2.8.3.10,citrate CoA-transferase,-0.87287156094397,0.0046592149439939
```

## Top PC-like features
```csv
feature,display_name,best_dose,min_delta
EC:3.1.3.4,phosphatidate phosphatase,G50,0.0321774292049148
EC:3.2.1.139,alpha-glucuronidase,G25,0.0964795199081951
EC:1.1.1.205,IMP dehydrogenase,G25,0.097819480501118
EC:2.8.3.10,citrate CoA-transferase,G50,0.105750603431453
EC:6.3.5.5,EC:6.3.5.5,G25,0.106773894985567
```
