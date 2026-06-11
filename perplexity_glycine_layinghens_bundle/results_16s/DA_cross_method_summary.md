# Differential Abundance Cross-Method Summary

## Current status

Additional DA analyses were run using:

- `ANCOMBC` genus-level compositional analysis
- `MaAsLin2` 4-group treatment model
- `MaAsLin2` glycine-dose model (`NC`, `G25`, `G50`)

The main point is that the results are directionally informative, but the two methods behave differently:

- `ANCOMBC` in the currently available environment is relatively aggressive and returns many zero-adjusted significant signals.
- `MaAsLin2` is much more conservative in this dataset and does not retain taxa at `q <= 0.25`, although several biologically plausible taxa rank near the top by raw p-value.

## Recurrently highlighted taxa

Across the current analyses, the following genera repeatedly appeared near the top:

- `Romboutsia`
- `Turicibacter`
- `Clostridia UCG-014`
- `Muribaculaceae`
- `Bacteroides`
- `Lachnospiraceae NK4A136 group`
- `Enterococcus`
- `Gallibacterium`

These taxa also overlap meaningfully with the earlier nonparametric analyses from the existing `Functional analysis` folder, especially the dose-response and glycine-mimicry outputs.

## ANCOMBC interpretation

### NC-reference model

The strongest genus-level signals versus `NC` included:

- `Romboutsia` increased in both `PC` and `G50`
- `Muribaculaceae` increased in `G25`
- `Turicibacter` increased in `PC` and `G50`
- `Clostridia UCG-014` increased in `PC` and `G50`
- `Bacteroides` increased in `G25`
- `Lachnospiraceae NK4A136 group` increased in `G25`
- `Gallibacterium` increased in `G50`

This pattern is broadly consistent with the previous interpretation that glycine does not produce a uniform diversity restoration effect, but instead shifts a subset of treatment-responsive taxa, some of which move toward a `PC`-like state.

### PC-reference model

When `PC` was used as the reference, the strongest contrasts included:

- `Romboutsia`: `NC` lower than `PC`
- `Turicibacter`: `NC` lower than `PC`
- `Clostridia UCG-014`: `NC` lower than `PC`
- `Lachnospiraceae NK4A136 group`: `G25` higher than `PC`
- `Muribaculaceae`: `G25` higher than `PC`

This suggests that some glycine-responsive taxa do not simply “restore toward PC” in a uniform way. Instead, glycine appears to create a partially overlapping but not identical ecological configuration.

## Glycine dose model

The clearest linear glycine-dose signal in `ANCOMBC` was:

- `Romboutsia`: positive glycine-dose coefficient, `q = 0.000276`

Several other taxa showed nominal dose-related behavior but did not survive multiple testing:

- `Turicibacter`
- `Clostridia UCG-014`
- `Gallibacterium`
- `Propioniciclava`

This is important because it supports the earlier Spearman-based trend result where `Romboutsia`, `Enterococcus`, and `Clostridia UCG-014` showed dose-related movement.

## MaAsLin2 interpretation

`MaAsLin2` did not produce significant taxa after FDR correction in either the 4-group model or the glycine-dose model. However, the top raw p-value taxa are biologically coherent:

### Treatment model top raw signals

- `Lactobacillus` lower in `PC`
- `Turicibacter` higher in `PC`
- `Romboutsia` higher in `PC`
- `Clostridia UCG-014` higher in `PC`
- `Muribaculaceae` higher in `G25`
- `Aeriscardovia` lower in `G25`

### Glycine-dose model top raw signals

- `Romboutsia`
- `Clostridia UCG-014`
- `Enterococcus`
- `Turicibacter`
- `Gallibacterium`

Even though these are not FDR-significant in `MaAsLin2`, the overlap with the earlier KW trend analysis and with `ANCOMBC` makes them useful for biological interpretation and candidate prioritization.

## Practical interpretation for the manuscript

The strongest defensible statement at this stage is:

- `Romboutsia` is the most robust glycine dose-responsive candidate across methods.
- `Turicibacter`, `Clostridia UCG-014`, and `Enterococcus` remain biologically plausible glycine-responsive taxa, but their statistical support is method-dependent.
- `Muribaculaceae`, `Bacteroides`, and `Lachnospiraceae NK4A136 group` are relevant to the glycine-mimicry framing, especially for `G25`.

## Caution

These DA results should be written carefully in the paper:

- `ANCOMBC` signals are stronger, but they come from an older `ANCOM-BC` implementation available in the current environment.
- `MaAsLin2` is more conservative and does not confirm strong FDR-significant associations.
- Therefore, the safest wording is:
  - some taxa were **consistently highlighted across multiple analytical approaches**
  - `Romboutsia` showed the clearest glycine dose-related pattern
  - other taxa such as `Turicibacter`, `Clostridia UCG-014`, and `Enterococcus` should be framed as **supportive or trend-level candidates**

## Files

- [ANCOMBC_genus_ref_NC.csv](/home/yzyzero/glycine_16s_analysis/da_functional_phenotype_20260323/tables/ANCOMBC_genus_ref_NC.csv)
- [ANCOMBC_genus_ref_PC.csv](/home/yzyzero/glycine_16s_analysis/da_functional_phenotype_20260323/tables/ANCOMBC_genus_ref_PC.csv)
- [ANCOMBC_genus_glydose_linear.csv](/home/yzyzero/glycine_16s_analysis/da_functional_phenotype_20260323/tables/ANCOMBC_genus_glydose_linear.csv)
- [all_results.tsv](/home/yzyzero/glycine_16s_analysis/da_functional_phenotype_20260323/maaslin2_treat/all_results.tsv)
- [all_results.tsv](/home/yzyzero/glycine_16s_analysis/da_functional_phenotype_20260323/maaslin2_glydose/all_results.tsv)
