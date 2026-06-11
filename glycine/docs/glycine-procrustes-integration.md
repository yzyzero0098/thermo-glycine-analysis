# Glycine heat-stress 프로젝트의 publication-defensible host–microbiome integration 전략 재구성

## Executive Summary
이번 재구성의 핵심 결론은 간단하다. **이 데이터는 sample-level latent integration을 주분석으로 밀기에는 구조가 약하고, treatment-centroid geometry를 중심으로 읽는 것이 가장 방어적**이다.[1][2][4][5][6][33][34][35] 로컬 파일을 다시 교차검증한 결과, 기존에 더 넓게 언급되던 host-linked matched set 17개와 달리, **현재 보존된 Bray-Curtis / weighted UniFrac ordination에 실제로 남아 있고 RNA-seq와 exact pairing되는 host–microbiome sample은 14쌍**이었다 (`experiments/glycine-procrustes-integration/crosswalk/matched_host_microbiome_retained16.csv`, `experiments/glycine-procrustes-integration/crosswalk/matched_host_microbiome_retained16_counts.csv`; [4][11]). 따라서 Procrustes는 **treatment-centroid level을 primary**, exact matched sample-level Procrustes는 **secondary / exploratory**로 두는 것이 타당하다.[5][6][33][34][35]

실제 계산 결과도 이 방향을 지지했다. **centroid-level Procrustes는 exact 4-centroid permutation test에서 통계적으로 지지되는 concordance를 보여주지 않았고, Bray-Curtis에서만 방향성 수준의 해석 여지가 있었으며 weighted UniFrac에서는 그 정렬이 더 약했다.** summary 값은 `experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_summary.csv`에 기록되어 있으며, host full25 vs Bray-Curtis는 disparity 0.4306, exact permutation p=0.2917, host full25 vs weighted UniFrac는 disparity 0.7845, p=0.8333이었다.[6][22] sample-level Procrustes 역시 exact retained 14쌍 기준으로 주결론을 뒷받침할 만큼 강하지 않았다 (`experiments/glycine-procrustes-integration/procrustes/sample_procrustes_summary.csv`; Bray-Curtis disparity 0.8173, p=0.1731; weighted UniFrac disparity 0.9913, p=0.9830).[6][27] 따라서 최종 통합 서사는 “강한 sample-level host–microbiome coupling”이 아니라, **treatment에 따른 방향성 있는 비교는 가능하지만 host–microbiome concordance를 강하게 주장하지 않는 구조**로 재설계하는 것이 안전하다.[2][5][6]

이 요약은 **non-significant small-centroid permutation 결과를 바탕으로 한 보수적 방향성 해석**이며, 확인적 결론이라기보다 reviewer-safe한 integration framing으로 읽어야 한다.[5][6][22][23][24][25][26]

또한 직접 **G50 vs PC** 관점을 넣어보면, 이전의 단순한 “G50이 PC를 mimic한다”는 해석은 그대로 유지되기 어렵다. **Transcriptome에서는 G50이 PC와 완전히 겹치지 않고 더 강하고 부분적으로 다른 host state**로 보인다. 이는 `experiments/glycine-procrustes-integration/procrustes/host_full25_centroid_pairwise_distances.csv`와 `.../host_matched16_centroid_pairwise_distances.csv`, 그리고 `.../transcriptome_direct_contrast/G50_vs_PC_PyDESeq2_summary.csv`에서 직접 확인된다.[7][8][14][15][16] 반면 **microbiome centroid 수준에서는 direct G50–PC distance가 Bray-Curtis와 weighted UniFrac 모두에서 NC–PC distance보다 작다** (`experiments/glycine-procrustes-integration/procrustes/bray_centroid_pairwise_distances.csv`, `.../weighted_unifrac_centroid_pairwise_distances.csv`).[3][7][20][21] 다만 **global host–microbiome geometry concordance**는 weighted UniFrac에서 더 약했다.[6][22][24][26] **Phenotype의 현재 보존 자산은 treatment mean 반복 구조**이므로 sample-level 해석은 불가하지만, treatment-profile distance 기준으로는 **G25와 G50 모두 PC 쪽으로 이동하며, 오히려 G25가 G50보다 PC에 더 가까운 패턴**이 나타난다 (`experiments/glycine-procrustes-integration/phenotype/phenotype_pairwise_standardized_distances.csv`).[3][10][31] 따라서 출판용 서사는 “전체 층위에서의 단순한 PC mimicry”보다, **G50의 distinct high-dose host state + microbiome의 partial PC closeness + phenotype의 partial PC-approach**가 더 방어적이다. 여기서 “distinct high-dose”는 threshold-like behavior를 직접 증명한 결론이 아니라, 현재 구조와 direct contrast가 시사하는 working interpretation이다.[7][10]

---

## 1. 연구 목적과 설계 원칙

이번 작업의 목적은 다음 질문에 답하는 것이었다.[1][2][4][5][6][7]

1. 현재 보존된 glycine heat-stress 프로젝트 자산에서 **가장 publication-defensible한 host–microbiome integration 전략**은 무엇인가?
2. **Procrustes**를 주 통합 프레임으로 쓸 때, sample-level이 아니라 **treatment-centroid level**이 더 타당한가?
3. microbiome ordination을 **Bray-Curtis primary**, **weighted UniFrac sensitivity**로 재구성했을 때 treatment geometry가 일관적인가?
4. **G50 vs PC direct comparison**을 transcriptome, microbiome, phenotype 각 층위에서 어떻게 읽어야 하는가?

이번 재구성에서 유지한 분석 규칙은 세 가지였다.[1][5][6][33][34][35]

- **sample-level correlation 또는 latent integration을 강제로 만들지 않는다.**[1][5][33][34][35]
- **centroid-level Procrustes를 primary로 둔다.**[1][5][6][33]
- **Evidence / Inference / Limitation을 분리해서 쓴다.**[1][2][5]

---

## 2. 데이터 재구성과 실제 사용 가능한 분석 단위

### 2.1 Modality별 보존 자산
- **16S**: `~/glycine_16s_analysis/16s_preprocessing` 아래에 raw/preprocessing asset과 Bray-Curtis / UniFrac QIIME2 artifact가 실제로 존재했다.[2]
- **RNA-seq**: `~/glycine_rnaseq_analysis_20260324/tables/` 아래에 `vst_matrix.csv`, count matrix, sample metadata, NC-referenced DE table이 보존되어 있었다.[3]
- **Phenotype**: `~/glycine_phenotype/integrated_metadata_clean.csv`와 workbook/PDF가 있었지만, phenotype-integrated CSV는 treatment mean을 sample row에 반복한 구조였다.[3][10]

### 2.2 Exact matched set은 17이 아니라 14 retained pairs
기존 continuity artifact에서는 host-linked matched sample이 17개로 요약되어 있었지만, 이번에 **실제 preserved microbiome ordination에 남아 있는 sample**과 RNA-seq metadata를 `Treat + Replicate` 기준으로 다시 맞춘 결과, **exact retained host–microbiome pair는 14개**였다 (`experiments/glycine-procrustes-integration/crosswalk/matched_host_microbiome_retained16_counts.csv`).[4][11] 여기서 파일명에 남아 있는 `matched16` 표기는 **retained microbiome ordination basis가 16 sample**이라는 뜻이지, exact host–microbiome pair 수가 16이라는 뜻은 아니다.[4][11]

- NC = 3
- PC = 4
- G25 = 3
- G50 = 4

즉,
- broader host-linked integration set = 17 (`experiments/glycine-matched-samples/matched_group_counts.csv`를 반영한 과거 요약; [3][10])
- current Bray/weighted UniFrac ordination-supported exact host–microbiome set = 14 (`experiments/glycine-procrustes-integration/crosswalk/matched_host_microbiome_retained16_counts.csv`; [4][11])

으로 구분해야 한다.

### Evidence
- `notes/glycine-procrustes-integration-crosswalk.md`[4]
- `experiments/glycine-procrustes-integration/crosswalk/matched_host_microbiome_retained16.csv`[11]
- `experiments/glycine-procrustes-integration/crosswalk/matched_host_microbiome_retained16_counts.csv`[11]

### Inference
이 차이는 단순 bookkeeping 문제가 아니라, **sample-level integration claim의 방어력 자체를 낮춘다**. 따라서 sample-level Procrustes나 sCCA/DIABLO를 본문 주분석으로 두기 어렵다.[1][4][5][33][34][35]

### Limitation
이 14쌍은 **현재 보존된 16S ordination artifact 기준**의 retained exact pairing이다. raw 16S 전체를 다른 rarefaction/depth 정의로 다시 재실행하면 retained set이 달라질 가능성은 있다. 다만 현재 제출 가능한 증거는 보존 artifact 기준이어야 한다.[2][4]

---

## 3. Transcriptome ordination과 G50 vs PC direct comparison

### 3.1 Host PCA 재구성
`vst_matrix.csv`로부터 host PCA를 다시 계산했다.[3][7]

- Axis1 explained variance = **31.65%** (`experiments/glycine-procrustes-integration/transcriptome_ordination/host_pca_full25_variance.csv`; [13])
- Axis2 explained variance = **11.51%** (`experiments/glycine-procrustes-integration/transcriptome_ordination/host_pca_full25_variance.csv`; [13])
- Axis3 explained variance = **7.23%** (`experiments/glycine-procrustes-integration/transcriptome_ordination/host_pca_full25_variance.csv`; [13])

Primary host geometry는 full host set 25개를 기준으로 잡고, retained host–microbiome exact pair 14개 subset은 sensitivity로 병행했다.[3][7][12][13] 다만 이 primary host centroid는 microbiome centroid와 **동일한 sample frame**에서 추정된 것이 아니므로, full25-based centroid Procrustes는 strict matched-cohort integration이 아니라 **treatment-level geometric comparison**으로 해석해야 한다. exact retained 14-pair에 대응하는 host centroid는 sensitivity로 함께 제시했다.[4][5][6]

### 3.2 Host centroid geometry
#### Full host set 기준 centroid distance
- NC–PC = **28.69** (`experiments/glycine-procrustes-integration/procrustes/host_full25_centroid_pairwise_distances.csv`; [14])
- G25–PC = **4.82** (`experiments/glycine-procrustes-integration/procrustes/host_full25_centroid_pairwise_distances.csv`; [14])
- G50–PC = **22.66** (`experiments/glycine-procrustes-integration/procrustes/host_full25_centroid_pairwise_distances.csv`; [14])
- G50–NC = **44.86** (`experiments/glycine-procrustes-integration/procrustes/host_full25_centroid_pairwise_distances.csv`; [14])

#### Retained exact-pair host subset 기준 centroid distance
- NC–PC = **21.51** (`experiments/glycine-procrustes-integration/procrustes/host_matched16_centroid_pairwise_distances.csv`; [15])
- G25–PC = **3.22** (`experiments/glycine-procrustes-integration/procrustes/host_matched16_centroid_pairwise_distances.csv`; [15])
- G50–PC = **30.78** (`experiments/glycine-procrustes-integration/procrustes/host_matched16_centroid_pairwise_distances.csv`; [15])
- G50–NC = **46.57** (`experiments/glycine-procrustes-integration/procrustes/host_matched16_centroid_pairwise_distances.csv`; [15])

두 host 정의 모두에서 **G25가 PC에 가장 가까운 dose**였고, **G50는 PC와 완전히 합쳐지지 않는 별도의 더 큰 이동**을 보였다.[7][14][15]

### 3.3 Direct G50 vs PC differential comparison
보존된 count matrix로 **PyDESeq2** direct contrast를 새로 실행했다.[7]

- comparison: **G50 vs PC**
- samples: G50=5, PC=7 (`experiments/glycine-procrustes-integration/transcriptome_direct_contrast/G50_vs_PC_PyDESeq2_summary.csv`; [16])
- genes tested: **21,392** (`experiments/glycine-procrustes-integration/transcriptome_direct_contrast/G50_vs_PC_PyDESeq2_summary.csv`; [16])
- significant genes at FDR<0.05 and |log2FC|>=1: **12** (`experiments/glycine-procrustes-integration/transcriptome_direct_contrast/G50_vs_PC_PyDESeq2_summary.csv`; [16])
  - up in G50: 11 (`.../G50_vs_PC_PyDESeq2_summary.csv`; [16])
  - up in PC: 1 (`.../G50_vs_PC_PyDESeq2_summary.csv`; [16])

이 contrast는 null은 아니지만, 기존 보존 자산에서 가장 강했던 `NC vs G50` 급의 대규모 차이는 아니다.[3][7][16][17] 즉 **G50와 PC는 가까운 면이 있지만 동일하지는 않다**.[7][10]

### Evidence
- `experiments/glycine-procrustes-integration/transcriptome_ordination/host_pca_full25_coordinates.csv`[12]
- `experiments/glycine-procrustes-integration/transcriptome_ordination/host_pca_full25_variance.csv`[13]
- `experiments/glycine-procrustes-integration/procrustes/host_full25_centroid_pairwise_distances.csv`[14]
- `experiments/glycine-procrustes-integration/procrustes/host_matched16_centroid_pairwise_distances.csv`[15]
- `experiments/glycine-procrustes-integration/transcriptome_direct_contrast/G50_vs_PC_PyDESeq2_summary.csv`[16]
- `experiments/glycine-procrustes-integration/transcriptome_direct_contrast/G50_vs_PC_PyDESeq2_full.csv`[17]

### Inference
Transcriptome 층위에서는 **G50를 단순한 PC mimic**로 쓰기보다, **distinct high-dose host state**로 쓰는 것이 훨씬 안전하다. 현재 결과는 threshold-like behavior와 양립 가능하지만, 그것을 직접 검정한 것은 아니다. G25는 host geometry 기준으로 오히려 PC convergence가 더 강하다.[7][10][14][15]

### Limitation
- direct contrast는 **PyDESeq2**로 새로 수행한 결과이며, 기존 archived DESeq2 output과 같은 실행환경/R package stack에서 재실행한 것은 아니다.[7][16][17]
- pathway-level sample score matrix는 현재 보존되어 있지 않아, 이번 round의 host ordination primary는 PCA다.[3][7][10]

---

## 4. Microbiome ordination: Bray-Curtis primary, weighted UniFrac sensitivity

### 4.1 재구성 가능한 microbiome space
현재 workspace에는 Bray-Curtis와 weighted UniFrac의 QIIME2 artifact가 모두 남아 있어서, **기존 ordination geometry를 artifact에서 직접 추출**할 수 있었다.[2][8] 단, visible shell script의 rarefaction depth와 saved artifact provenance 사이에는 불일치가 있었다.[2]

### 4.2 Bray-Curtis centroid distance
- NC–PC = **0.2390** (`experiments/glycine-procrustes-integration/procrustes/bray_centroid_pairwise_distances.csv`; [20])
- G25–PC = **0.2982** (`experiments/glycine-procrustes-integration/procrustes/bray_centroid_pairwise_distances.csv`; [20])
- G50–PC = **0.1888** (`experiments/glycine-procrustes-integration/procrustes/bray_centroid_pairwise_distances.csv`; [20])

Bray-Curtis에서는 **G50가 NC보다 PC에 더 가까웠다**.[7][20]

### 4.3 weighted UniFrac centroid distance
- NC–PC = **0.0674** (`experiments/glycine-procrustes-integration/procrustes/weighted_unifrac_centroid_pairwise_distances.csv`; [21])
- G25–PC = **0.0741** (`experiments/glycine-procrustes-integration/procrustes/weighted_unifrac_centroid_pairwise_distances.csv`; [21])
- G50–PC = **0.0402** (`experiments/glycine-procrustes-integration/procrustes/weighted_unifrac_centroid_pairwise_distances.csv`; [21])

weighted UniFrac에서도 direct distance 질문만 보면 **G50가 NC보다 PC에 더 가까웠다**.[7][21]

### 4.4 두 거리계의 의미 차이
그러나 이 둘은 같은 정도의 story를 주지 않는다.[2][7]

- **Direct G50→PC closeness 기준**: Bray-Curtis와 weighted UniFrac 모두에서 G50–PC 거리는 NC–PC 거리보다 작다.[7][20][21]
- **Global host–microbiome geometry concordance 기준**: weighted UniFrac는 Bray-Curtis보다 더 약한 Procrustes 정렬을 보인다.[6][22][24][26]

즉,
- “G50는 PC에 microbiome적으로 일부 가까워진다”는 문장은 두 거리계 모두에서 가능하지만,[7][20][21]
- “microbiome 전체 geometry가 host와 강하게 concordant하다”는 문장은 두 거리계 모두에서 방어가 어렵고, 특히 weighted UniFrac에서 더 약하다.[2][6][22][24][26]

### Evidence
- `experiments/glycine-procrustes-integration/microbiome_ordination/bray_pcoa_coordinates.csv`[18]
- `experiments/glycine-procrustes-integration/microbiome_ordination/weighted_unifrac_pcoa_coordinates.csv`[19]
- `experiments/glycine-procrustes-integration/procrustes/bray_centroid_pairwise_distances.csv`[20]
- `experiments/glycine-procrustes-integration/procrustes/weighted_unifrac_centroid_pairwise_distances.csv`[21]

### Inference
Microbiome 층위에서 가장 안전한 해석은 **G50의 partial PC convergence**다. 다만 이 convergence는 **distance metric-dependent**하고, 전체 multi-omic concordance를 강하게 주장할 만큼 robust하지는 않다.[2][6][7]

### Limitation
- 현재 사용한 microbiome geometry는 **보존된 ordination artifact**를 재사용한 것이다.[2][8]
- script provenance와 saved artifact provenance 사이 rarefaction depth mismatch가 있어, 완전한 end-to-end rerun provenance는 아직 미완이다.[2]

---

## 5. Primary integration result: centroid-level Procrustes

### 5.1 왜 centroid-level이 primary인가
이 프로젝트에서 exact retained host–microbiome pair는 14개뿐이고, phenotype은 sample-level로 독립적이지 않다.[3][4][10][11] 또한 문헌/공식 method source에 따르면:
- **Procrustes**는 matched rows를 전제하고,[1][33]
- **sCCA**는 paired matrices를 전제하며,[1]
- **DIABLO**는 multiple omics blocks measured on the same samples를 전제한다.[1][34][35]

현재 데이터 구조에서 가장 정직한 matched row는 **individual sample이 아니라 treatment centroid**다.[1][4][5][6]

### 5.2 Exact-permutation centroid-level Procrustes results
#### Host full25 vs Bray-Curtis
- disparity = **0.4306** (`experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_summary.csv`; [22])
- exact permutation p = **0.2917** (`experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_host_full25_vs_bray_permutations.csv`; [23])

#### Host full25 vs weighted UniFrac
- disparity = **0.7845** (`experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_summary.csv`; [22])
- exact permutation p = **0.8333** (`experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_host_full25_vs_weighted_unifrac_permutations.csv`; [24])

#### Host matched14-equivalent subset vs Bray-Curtis
- disparity = **0.5116** (`experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_summary.csv`; [22])
- exact permutation p = **0.4167** (`experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_host_matched16_vs_bray_permutations.csv`; [25])

#### Host matched14-equivalent subset vs weighted UniFrac
- disparity = **0.8386** (`experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_summary.csv`; [22])
- exact permutation p = **0.9167** (`experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_host_matched16_vs_weighted_unifrac_permutations.csv`; [26])

### 5.3 해석
이 결과는 “strongly significant multi-omic concordance”를 지지하지 않는다.[5][6][22][23][24][25][26] 현재 evidence로는:
- **Bray-Curtis에서 방향성 있는 treatment-geometry 정렬을 시사하는 정도**[6][22][23][25]
- **weighted UniFrac에서는 그 정렬이 더 약하거나 거의 지지되지 않는 정도**[6][22][24][26]

수준이 가장 안전하다.

### Evidence
- `experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_summary.csv`[22]
- `experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_host_full25_vs_bray_permutations.csv`[23]
- `experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_host_full25_vs_weighted_unifrac_permutations.csv`[24]
- `experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_host_matched16_vs_bray_permutations.csv`[25]
- `experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_host_matched16_vs_weighted_unifrac_permutations.csv`[26]

### Inference
Main text의 integration conclusion은 **“host와 microbiome의 treatment geometry 사이에 방향성 있게 비교 가능한 패턴은 보이지만, exact permutation 기준으로 통계적으로 지지되는 concordance는 확인되지 않았고 weighted UniFrac에서는 그 정렬이 더 약하다”** 정도로 제한해야 한다.[5][6]

### Limitation
4 treatment centroid만을 쓰는 exact permutation은 가능한 permutation 수 자체가 24개뿐이다. 즉 통계 검정력은 본질적으로 낮다. 하지만 이는 오히려 **과장된 significance claim을 막아주는 보수적 장치**이기도 하다.[5][6][22][23][24][25][26]

---

## 6. Secondary / exploratory result: sample-level Procrustes

Exact retained 14 host–microbiome pair에서 sample-level Procrustes를 exploratory로만 실행했다.[4][5][6]

- **Bray-Curtis:** disparity = 0.8173, permutation p = 0.1731 (`experiments/glycine-procrustes-integration/procrustes/sample_procrustes_summary.csv`; [27])
- **weighted UniFrac:** disparity = 0.9913, permutation p = 0.9830 (`experiments/glycine-procrustes-integration/procrustes/sample_procrustes_summary.csv`; [27])

이 결과는 sample-level integration을 main claim으로 올릴 만큼 강하지 않다.[5][6][27]

### Evidence
- `experiments/glycine-procrustes-integration/procrustes/sample_procrustes_summary.csv`[27]
- `experiments/glycine-procrustes-integration/procrustes/sample_procrustes_matched16_bray_paired_coordinates.csv`[28]
- `experiments/glycine-procrustes-integration/procrustes/sample_procrustes_matched16_weighted_unifrac_paired_coordinates.csv`[29]

### Inference
sample-level 결과는 **주장을 강화하기보다 restraint를 걸어주는 sensitivity check**로 쓰는 편이 낫다.[5][6]

### Limitation
pair 수가 14로 작고 group imbalance도 있어, 이 결과는 exploratory beyond primary 수준을 넘기기 어렵다.[4][5][6]

---

## 7. Dose-response와 PC mimicry의 재분류

## 7.1 Transcriptome
- **G25**: PC에 가장 가까운 geometry (`experiments/glycine-procrustes-integration/procrustes/host_full25_centroid_pairwise_distances.csv`, `.../host_matched16_centroid_pairwise_distances.csv`; [14][15])
- **G50**: 더 큰 host-state transition (`.../host_full25_centroid_pairwise_distances.csv`; [14])
- direct `G50 vs PC`는 non-null이지만 소규모 (`.../transcriptome_direct_contrast/G50_vs_PC_PyDESeq2_summary.csv`; [16])

→ **distinct high-dose host response**. threshold-like 해석은 현재 결과와 양립 가능하지만 working hypothesis 수준에 머문다.[7][10]

## 7.2 Microbiome
- Bray-Curtis와 weighted UniFrac 모두 direct G50–PC distance는 NC–PC보다 작음 (`.../bray_centroid_pairwise_distances.csv`, `.../weighted_unifrac_centroid_pairwise_distances.csv`; [20][21])
- 하지만 weighted UniFrac는 전체 concordance를 약화시킴 (`.../centroid_procrustes_summary.csv`; [22])

→ **partial G50-to-PC convergence**, but **not robust global concordance**.[2][6][7]

## 7.3 Phenotype
현재 available phenotype은 treatment mean 반복 구조이므로 group-profile standardized distance로만 평가했다.[3][10]

- NC–PC = **8.66** (`experiments/glycine-procrustes-integration/phenotype/phenotype_pairwise_standardized_distances.csv`; [31])
- G25–PC = **5.76** (`experiments/glycine-procrustes-integration/phenotype/phenotype_pairwise_standardized_distances.csv`; [31])
- G50–PC = **6.61** (`experiments/glycine-procrustes-integration/phenotype/phenotype_pairwise_standardized_distances.csv`; [31])

즉, 둘 다 NC보다 PC에 가깝지만, **G25가 G50보다 더 PC-like**하다.[7][10][31]

### Figure: G25/G50의 PC distance를 NC–PC baseline 대비 비교

![PC distance ratio chart](../glycine-procrustes-integration-pc-distance-ratios.png)

**Caption.** 각 층위에서 G25 또는 G50의 PC까지의 거리(`distance to PC`)를 동일 층위의 `NC–PC` 거리로 나눈 비율. 1보다 작으면 NC보다 PC에 더 가깝고, 1보다 크면 NC보다 더 멀다. Transcriptome은 full host PCA centroid 기준, microbiome은 Bray-Curtis 및 weighted UniFrac centroid 기준, phenotype은 group-level standardized profile distance 기준이다. 기반 데이터는 `experiments/glycine-procrustes-integration/procrustes/direct_pc_closeness_summary.csv`와 `experiments/glycine-procrustes-integration/phenotype/phenotype_pairwise_standardized_distances.csv`에 있다.[30][31]

### Evidence
- `notes/glycine-procrustes-integration-dose-response.md`[7]
- `experiments/glycine-procrustes-integration/procrustes/direct_pc_closeness_summary.csv`[30]
- `experiments/glycine-procrustes-integration/phenotype/phenotype_pairwise_standardized_distances.csv`[31]

### Inference
전체 층위를 하나의 monotonic dose-response로 묶기 어렵다. 현재 evidence는 더 잘게 나눠서 써야 한다.[7][10]

### Limitation
Phenotype은 group-level summary distance이고, microbiome은 distance metric에 따라 story 강도가 달라진다. 따라서 “cross-layer monotonicity”는 과장 표현이 된다.[2][3][7][10]

---

## 8. 가장 방어적인 biological narrative

가장 방어적인 서사는 다음과 같다.[3][7][10][36][37]

1. **G50는 host transcriptome에서 가장 강한 remodeling을 유도한다.**[3][7][14][16]
2. 이 host response는 **PC의 단순 재현이라기보다 더 강하고 부분적으로 다른 state**다.[7][14][15][16]
3. **Microbiome은 G50에서 PC 쪽으로 부분적 convergence**를 보이지만, 이는 **global concordance 기준에서는** Bray-Curtis에서 더 해석 가능하고 weighted UniFrac에서는 더 제한적이다.[2][7][20][21][22]
4. **Phenotype은 두 glycine dose 모두 PC 쪽으로 이동**하나, 현재 보존 자산에서는 G25가 더 PC-like하게 보인다.[3][7][10][31]
5. 따라서 integrative conclusion은 **“strong host–microbiome mimicry”**가 아니라,
   **“host 쪽에서 더 분명하게 관찰되는 distinct high-dose remodeling, microbiome의 partial PC closeness, phenotype-level partial PC approach”** 정도가 적절하다. 여기서 “host 쪽에서 더 분명하다”는 표현은 공통 척도로 정량 비교한 결과가 아니라, 현재 분석 프레임에서의 질적 요약이다.[5][6][7]

### Evidence
위 3–7절의 계산 결과와 prior phenotype / RNA-seq continuity artifact.[3][7][10][14][15][16][20][21][22][27][31]

### Inference
주논문 Results에서 microbiome은 **supportive layer**로 두고, phenotype + RNA-seq를 중심축으로 두는 전략이 가장 안전하다.[5][36][37]

### Limitation
microbiome과 phenotype 모두 host transcriptome만큼 강한 독립적 주장을 지지하지는 않는다. 특히 phenotype은 raw per-bird values 부재, microbiome은 metric sensitivity가 문제다.[2][3][10]

---

## 9. 추천 Figure 구조

## Main text

### Figure 1. Study design + sample overlap + modality availability
- treatment design
- RNA-seq / 16S / phenotype sample counts
- exact retained host–microbiome pair count = **14** 명시[4][11]
- integration analysis hierarchy 명시[5][6]

### Figure 2. Host transcriptome headline
- PCA[12][13]
- 대표 pathway/module 결과 또는 preserved GSEA summary[9][38][40]
- `G50 vs PC` direct contrast 요약[16][17]

### Figure 3. Phenotype anchor
- reviewer에게 중요한 phenotype endpoint만 정리[3][10]
- PC, G25, G50의 treatment-level proximity를 보조로 제시[31]

### Figure 4. Microbiome ordination and centroid geometry
- Bray-Curtis PCoA primary[18][20]
- weighted UniFrac inset 또는 side panel sensitivity[19][21]
- centroid 표시
- G50–PC direct distance annotation[20][21]

### Figure 5. Primary integration figure
- centroid-level Procrustes (host PCA centroid vs Bray-Curtis centroid)[22][23]
- weighted UniFrac sensitivity 결과를 함께 요약[22][24][26]
- sample-level Procrustes는 main text가 아니라 supplement 또는 작은 inset[27]

### Figure 6. Working model
- observed / associated / inferred 화살표를 구분한 summary cartoon[5][36][37]

## Supplement
- exact sample crosswalk full table[4][11]
- weighted UniFrac full sensitivity panel[19][21][24][26]
- sample-level Procrustes full result[27][28][29]
- all permutation distributions[23][24][25][26]
- direct `G50 vs PC` PyDESeq2 full table summary[16][17]
- phenotype standardized distance derivation[31]
- any exploratory taxa/pathway bridge heatmaps[9][10]

---

## 10. Main text vs supplement recommendation

## Main text에 두어야 할 것
- host transcriptome main story[3][7]
- phenotype anchor[3][10][31]
- Bray-Curtis centroid geometry[20]
- centroid-level Procrustes primary result[22][23]
- explicit limitation sentence on imperfect matching[4][5][6]

## Supplement로 보내야 할 것
- sample-level Procrustes[27][28][29]
- weighted UniFrac full details[19][21][24][26]
- exhaustive centroid permutations[23][24][25][26]
- full DE table / full taxa detail[17]
- phenotype distance derivation table[31]

이 분리는 reviewer에게 “integration을 과장하지 않았다”는 인상을 준다.[5][36][37]

---

## 11. Exact executable scripts and artifacts

## Primary analysis script
- `experiments/glycine-procrustes-integration/run_analysis.py`[32]

## Reproducible execution command
```bash
uv run --with pandas --with numpy --with scipy --with scikit-learn --with pydeseq2 \
  experiments/glycine-procrustes-integration/run_analysis.py
```

## Key output directories
- `experiments/glycine-procrustes-integration/crosswalk/`[11]
- `experiments/glycine-procrustes-integration/transcriptome_ordination/`[12][13]
- `experiments/glycine-procrustes-integration/microbiome_ordination/`[18][19]
- `experiments/glycine-procrustes-integration/procrustes/`[14][15][20][21][22][23][24][25][26][27][28][29][30]
- `experiments/glycine-procrustes-integration/phenotype/`[31]
- `experiments/glycine-procrustes-integration/transcriptome_direct_contrast/`[16][17]

---

## 12. Reviewer-safe limitation statements

다음 문장들이 본문/Discussion에서 가장 안전하다.[1][5][6][33][34][35][36][37]

1. **On pairing**  
   Host transcriptome and microbiome profiles were not available for every animal, and the currently retained exact host–microbiome overlap for the preserved ordinations was limited to 14 pairs; therefore, cross-omics integration was interpreted primarily at the treatment-centroid level.[4][5][11]

2. **On phenotype**  
   The currently archived phenotype-linked table contains treatment-level repeated values rather than preserved per-bird measurements, so phenotype was used as a treatment-profile layer rather than an independent sample-level integration block.[3][10]

3. **On concordance strength**  
   Procrustes analysis suggested at most partial host–microbiome treatment-geometry concordance, and this signal weakened under weighted UniFrac sensitivity analysis.[6][22][24][26]

4. **On causality**  
   These cross-layer associations do not establish causal or directional coupling between microbial shifts and host transcriptional responses.[5][36][37]

5. **On interpretation**  
   Accordingly, we interpret the integration as supportive evidence of coordinated treatment-associated remodeling, rather than proof of strong sample-level host–microbiome coupling.[5][6]

---

## Open Questions
1. 16S를 현재 artifact provenance와 일치하는 depth/filtered table lineage로 완전 재실행하면 retained set이 달라지는가?[2]
2. raw per-bird phenotype values를 복구할 수 있는가?[3][10]
3. pathway-level sample score matrix를 별도 구축하면 host centroid narrative가 더 안정적으로 정리되는가?[3][9][10][38][40]
4. Bray-Curtis에서 보이는 G50–PC convergence가 taxa shortlist 수준에서도 일관적인가?[2][9][10]

---

## Bottom line
현재 preserved evidence 기준으로 가장 강한 publication-defensible integration 전략은 다음이다.[1][5][6]

- **Primary:** treatment-centroid-level Procrustes[5][6][22]
- **Primary microbiome distance:** Bray-Curtis[2][20]
- **Sensitivity:** weighted UniFrac[2][21][22]
- **Secondary/exploratory only:** sample-level Procrustes on exact retained pairs[4][27]
- **G50 vs PC conclusion:** transcriptome에서는 단순 mimicry가 아니라 stronger/distinct host state, microbiome에서는 partial convergence, phenotype에서는 group-profile 기준 partial PC approach[7][14][15][16][20][21][31]

즉, 이 프로젝트의 최적 출판 전략은 **phenotype + transcriptome을 중심축으로 두고, microbiome integration은 centroid-level geometry와 direct G50–PC distance로 제한적으로 제시하는 것**이다.[5][6][7][36][37] 다만 이 결론은 **non-significant small-centroid permutation 결과와 imperfect pairing을 전제로 한 보수적 방향성 해석**이며, strong host–microbiome concordance를 확인한 결과로 쓰면 안 된다.[5][6][22][23][24][25][26]

---

## Sources
1. `notes/glycine-procrustes-integration-research-methods.md` — file:///home/yzyzero/feynman/notes/glycine-procrustes-integration-research-methods.md
2. `notes/glycine-procrustes-integration-research-16s.md` — file:///home/yzyzero/feynman/notes/glycine-procrustes-integration-research-16s.md
3. `notes/glycine-procrustes-integration-research-host.md` — file:///home/yzyzero/feynman/notes/glycine-procrustes-integration-research-host.md
4. `notes/glycine-procrustes-integration-crosswalk.md` — file:///home/yzyzero/feynman/notes/glycine-procrustes-integration-crosswalk.md
5. `notes/glycine-procrustes-integration-design-decision.md` — file:///home/yzyzero/feynman/notes/glycine-procrustes-integration-design-decision.md
6. `notes/glycine-procrustes-integration-research-framing.md` — file:///home/yzyzero/feynman/notes/glycine-procrustes-integration-research-framing.md
7. `notes/glycine-procrustes-integration-dose-response.md` — file:///home/yzyzero/feynman/notes/glycine-procrustes-integration-dose-response.md
8. `outputs/glycine-omics-reanalysis.md` — file:///home/yzyzero/feynman/outputs/glycine-omics-reanalysis.md
9. `outputs/glycine-gsea-followup.md` — file:///home/yzyzero/feynman/outputs/glycine-gsea-followup.md
10. `outputs/glycine-log-interpretation.md` — file:///home/yzyzero/feynman/outputs/glycine-log-interpretation.md
11. `matched_host_microbiome_retained16` artifacts — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/crosswalk/matched_host_microbiome_retained16.csv ; file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/crosswalk/matched_host_microbiome_retained16_counts.csv
12. `host_pca_full25_coordinates.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/transcriptome_ordination/host_pca_full25_coordinates.csv
13. `host_pca_full25_variance.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/transcriptome_ordination/host_pca_full25_variance.csv
14. `host_full25_centroid_pairwise_distances.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/host_full25_centroid_pairwise_distances.csv
15. `host_matched16_centroid_pairwise_distances.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/host_matched16_centroid_pairwise_distances.csv
16. `G50_vs_PC_PyDESeq2_summary.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/transcriptome_direct_contrast/G50_vs_PC_PyDESeq2_summary.csv
17. `G50_vs_PC_PyDESeq2_full.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/transcriptome_direct_contrast/G50_vs_PC_PyDESeq2_full.csv
18. `bray_pcoa_coordinates.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/microbiome_ordination/bray_pcoa_coordinates.csv
19. `weighted_unifrac_pcoa_coordinates.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/microbiome_ordination/weighted_unifrac_pcoa_coordinates.csv
20. `bray_centroid_pairwise_distances.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/bray_centroid_pairwise_distances.csv
21. `weighted_unifrac_centroid_pairwise_distances.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/weighted_unifrac_centroid_pairwise_distances.csv
22. `centroid_procrustes_summary.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_summary.csv
23. `centroid_procrustes_host_full25_vs_bray_permutations.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_host_full25_vs_bray_permutations.csv
24. `centroid_procrustes_host_full25_vs_weighted_unifrac_permutations.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_host_full25_vs_weighted_unifrac_permutations.csv
25. `centroid_procrustes_host_matched16_vs_bray_permutations.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_host_matched16_vs_bray_permutations.csv
26. `centroid_procrustes_host_matched16_vs_weighted_unifrac_permutations.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/centroid_procrustes_host_matched16_vs_weighted_unifrac_permutations.csv
27. `sample_procrustes_summary.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/sample_procrustes_summary.csv
28. `sample_procrustes_matched16_bray_paired_coordinates.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/sample_procrustes_matched16_bray_paired_coordinates.csv
29. `sample_procrustes_matched16_weighted_unifrac_paired_coordinates.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/sample_procrustes_matched16_weighted_unifrac_paired_coordinates.csv
30. `direct_pc_closeness_summary.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/procrustes/direct_pc_closeness_summary.csv
31. `phenotype_pairwise_standardized_distances.csv` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/phenotype/phenotype_pairwise_standardized_distances.csv
32. `run_analysis.py` — file:///home/yzyzero/feynman/experiments/glycine-procrustes-integration/run_analysis.py
33. ade4 `procuste` documentation — https://rdrr.io/cran/ade4/man/procuste.html
34. mixOmics DIABLO documentation — https://mixomics.org/mixdiablo/
35. mixOmics N-integration methods documentation — https://mixomics.org/methods/n-integration-methods/
36. Li et al. (2022), *Integrated multi-omics of the gastrointestinal microbiome and ruminant host reveals metabolic adaptation underlying early life development* — https://link.springer.com/article/10.1186/s40168-022-01396-8
37. Fernandes et al. (2026), *Microbiome by transcriptome interactions triggered by a switch to an alternative diet in Nellore cattle* — https://www.nature.com/articles/s41598-025-29588-w
38. fgsea package page — https://bioconductor.org/packages/fgsea/
39. rocker/r-ver Docker image documentation — https://hub.docker.com/r/rocker/r-ver
40. msigdbr package site — https://igordot.github.io/msigdbr/
