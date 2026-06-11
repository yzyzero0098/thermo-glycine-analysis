# MetaCyc pathway follow-up summary

This summary reuses the saved PICRUSt2-derived functional tables and prioritizes features that support one or more of the following:
- omnibus group differences,
- monotonic dose trend,
- partial PC-like directionality/convergence.

Important limitation: these are **predicted functions from 16S-based PICRUSt2**, not directly measured microbial functional fluxes.

## Top omnibus-difference features
```csv
feature,display_name,kw_pvalue
PWY-7219,adenosine ribonucleotides de novo biosynthesis,0.0183791641764655
ANAGLYCOLYSIS-PWY,glycolysis III (from glucose),0.0188955677404606
UDPNAGSYN-PWY,UDP-N-acetyl-D-glucosamine biosynthesis I,0.0247731348150176
PWY4FS-7,phosphatidylglycerol biosynthesis I (plastidic),0.0249698610363008
PWY4FS-8,phosphatidylglycerol biosynthesis II (non-plastidic),0.0249698610363008
```

## Top dose-trend features
```csv
feature,display_name,spearman_rho,spearman_pvalue
GLUCOSE1PMETAB-PWY,glucose and glucose-1-phosphate degradation,-0.763762615825973,0.027396043909112
PWY-5860,PWY-5860,0.377964473009227,0.355917683749582
PWY-5862,PWY-5862,0.377964473009227,0.355917683749582
PWY-6123,PWY-6123,0.327326835353989,0.428691507235313
PWY-7219,adenosine ribonucleotides de novo biosynthesis,-0.218217890235992,0.603645056510136
```

## Top PC-like features
```csv
feature,display_name,best_dose,min_delta
ANAGLYCOLYSIS-PWY,glycolysis III (from glucose),G50,0.383482811887984
PHOSLIPSYN-PWY,PHOSLIPSYN-PWY,G50,0.38643347874878
PWY-5686,PWY-5686,G50,0.391517144207378
PWY4FS-7,phosphatidylglycerol biosynthesis I (plastidic),G50,0.394851614058567
PWY4FS-8,phosphatidylglycerol biosynthesis II (non-plastidic),G50,0.394851614058567
```
