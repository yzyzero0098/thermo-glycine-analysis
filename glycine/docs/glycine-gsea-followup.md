# Glycine omics follow-up: file check, matched-sample crosswalk, and preranked GSEA

작성일: 2026-04-03

## Summary
이번 턴에서는 세 가지를 실제로 진행했다.

1. **이전 턴에서 저장했다고 한 파일들이 실제로 존재하는지 확인했다.**
2. **RNA-seq–16S exact matched 17 samples crosswalk를 재생성하고, phenotype 메타데이터의 repeated group-mean 패턴을 다시 확인했다.**
3. **Docker 기반 R 4.3.3 환경에서 preranked fgsea(GO:BP)를 실제 실행했다.**

핵심적으로는:

- 이전에 보고된 주요 파일들은 **모두 존재**했다.
- matched sample은 다시 확인해도 **총 17개(NC=4, PC=5, G25=3, G50=5)**였다.
- 16S metadata의 phenotype 관련 26개 변수는 이번 재점검에서도 **모두 treatment group 내부에서 상수**였다. 따라서 현재 테이블만으로는 individual-level phenotype correlation을 강하게 주장하기 어렵다.
- 우선순위였던 **NC vs G25 preranked GSEA는 실제로 신호가 나온다.** 즉, G25 contrast는 hard DEG가 적더라도 pathway-level로는 빈 contrast가 아니다.

---

## 1. 파일 존재 확인
아래 파일들은 모두 workspace에서 확인됐다.

### 세션/보고서
- `notes/session-logs/2026-04-02-glycine-omics-reanalysis.md`
- `outputs/glycine-omics-reanalysis.md`
- `outputs/glycine-omics-reanalysis.provenance.md`

### RNA-seq 요약 산출물
- `experiments/glycine-rnaseq-summary/README.md`
- `experiments/glycine-rnaseq-summary/deseq2_threshold_sensitivity_counts.csv`
- `experiments/glycine-rnaseq-summary/mds_group_dispersion_summary.csv`
- `experiments/glycine-rnaseq-summary/outlier_key_scenarios.csv`
- `experiments/glycine-rnaseq-summary/NC_vs_G25_top30_by_padj.csv`

### lab notebook
- `CHANGELOG.md`

따라서 이전 턴의 로그/보고 체인은 현재 workspace에서 그대로 이어서 사용할 수 있다.

---

## 2. 새로 만든 matched-sample crosswalk
이번 턴에서 아래 산출물을 추가 생성했다.

- `experiments/glycine-matched-samples/README.md`
- `experiments/glycine-matched-samples/matched_17_samples_full_metadata.csv`
- `experiments/glycine-matched-samples/matched_17_samples_integration_ready.csv`
- `experiments/glycine-matched-samples/matched_group_counts.csv`
- `experiments/glycine-matched-samples/phenotype_repeat_pattern_summary.csv`
- `experiments/glycine-matched-samples/build_matched_crosswalk.py`

### matching rule
- RNA-seq sample metadata: `sample_metadata.csv`
- 16S integrated metadata: `integrated_metadata_clean.csv`
- matching key: **`Treat` + `Replicate`**
- 16S는 **`Tissue == Jejunum_Content`**로 제한

### exact matched counts
| Group | n |
|---|---:|
| NC | 4 |
| PC | 5 |
| G25 | 3 |
| G50 | 5 |
| Total | 17 |

이 값은 이전 보고와 동일하다.

### integration-ready sample pairs
대표적으로 아래와 같이 정규화되었다.

- `P1-2-JM-G25-R` ↔ `P1.2.JC.G25`
- `P1-7-JM-G50-R` ↔ `P1.7.JC.G50`
- `P1-7-JM-NC-R` ↔ `P1.7.JC.NC`
- `P1-5-JM-PC-R` ↔ `P1.5.JC.PC`

즉, 다음 단계의 pathway–taxa integration table은 이 crosswalk를 기준으로 바로 작성할 수 있다.

---

## 3. phenotype metadata 재확인
`phenotype_repeat_pattern_summary.csv`를 새로 계산한 결과:

- audit 대상 변수: **26개**
- treatment group 내부에서 상수인 변수: **26개 전부**

즉, 이번 재확인에서도 현재 16S 통합 metadata 안의 phenotype 값은 **실측 개인값이라기보다 group-level repeated value**로 읽히는 구조였다.

### 해석
이 결과는 이전 판단을 강화한다.

- 지금 상태에서 phenotype을 넣은 상관/예측/매개효과 모델을 강하게 쓰면 위험하다.
- phenotype은 현재로서는 **ecological / group-level annotation** 정도로만 제한하는 것이 안전하다.
- 따라서 다음 통합 단계의 우선순위는 phenotype 예측형 모델이 아니라,
  - **RNA pathway score ↔ taxa abundance**
  - **group-level directional consistency**
  - **dose-pattern concordance**
  중심이 더 적절하다.

---

## 4. Docker 기반 preranked GSEA 실제 실행
이전 턴에서 못 했던 GSEA를 이번에는 **Docker 기반 R 환경**으로 실제 실행했다.

### 실행 환경
- Docker image: `rocker/r-ver:4.3.3`
- 추가 설치: `libcurl4-openssl-dev`, `libssl-dev`, `libxml2-dev`
- R packages: `fgsea`, `msigdbr`, `data.table`, `dplyr`, `curl`, `BiocManager`

### 실행 스크립트
- `experiments/glycine-gsea-docker/run_preranked_fgsea.R`

### 결과 폴더
- `experiments/glycine-gsea-docker/results/README.md`
- `experiments/glycine-gsea-docker/results/fgsea_go_bp_summary.csv`
- `experiments/glycine-gsea-docker/results/NC_vs_G25_fgsea_go_bp.csv`
- `experiments/glycine-gsea-docker/results/NC_vs_G25_fgsea_go_bp_top30_sig.csv`
- `experiments/glycine-gsea-docker/results/NC_vs_G50_fgsea_go_bp.csv`
- `experiments/glycine-gsea-docker/results/NC_vs_G50_fgsea_go_bp_top30_sig.csv`
- `experiments/glycine-gsea-docker/results/NC_vs_PC_fgsea_go_bp.csv`
- `experiments/glycine-gsea-docker/results/NC_vs_PC_fgsea_go_bp_top30_sig.csv`

### 설정
- ranking statistic: 기존 DESeq2 full result table의 **Wald `stat`**
- gene sets: `msigdbr` **Gallus gallus GO:BP**
- pathway size filter: **10–500 genes**

중요하게, 이것은 **fresh DESeq2 rerun이 아니라 기존 full table을 이용한 ranked GSEA**다. 따라서
- DEG cutoff 의존성을 줄이는 목적에는 적합하지만,
- 디자인식 자체를 새로 돌린 것은 아니다.

---

## 5. GSEA 핵심 결과
`fgsea_go_bp_summary.csv` 기준 요약은 아래와 같다.

| Comparison | Ranked genes | FDR < 0.05 pathways | FDR < 0.10 pathways |
|---|---:|---:|---:|
| NC vs G25 | 17,363 | 182 | 355 |
| NC vs G50 | 17,363 | 813 | 1,115 |
| NC vs PC | 17,363 | 271 | 391 |

### 핵심 해석 1: NC vs G25는 pathway-level로는 비어 있지 않다
이건 이번 턴에서 가장 중요한 추가 확인이다.

기존 gene-level cutoff에서는 `NC vs G25`가 매우 약했지만,
preranked GSEA에서는 **유의 pathway가 다수 검출**됐다.

따라서 `G25는 효과가 없다`가 아니라,

> **G25는 gene-level hard DEG에는 불리하지만 ranked pathway signal은 존재하는 contrast**

로 읽는 쪽이 더 타당하다.

### 핵심 해석 2: 방향성은 기존 functional narrative와 상당히 맞는다
NC vs G25 / NC vs G50 / NC vs PC 공통으로 상위 positive hit에 다음 축이 반복된다.

- `GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION`
- `GOBP_COLLAGEN_FIBRIL_ORGANIZATION`
- `GOBP_BASEMENT_MEMBRANE_ORGANIZATION`
- `GOBP_CELL_SUBSTRATE_ADHESION`
- branching / morphogenesis 관련 항목

반대로 negative 쪽에는 반복적으로

- `GOBP_GOLGI_VESICLE_TRANSPORT`
- `GOBP_ENDOPLASMIC_RETICULUM_TO_GOLGI_VESICLE_MEDIATED_TRANSPORT`
- `GOBP_VESICLE_ORGANIZATION`
- `GOBP_MACROAUTOPHAGY`
- `GOBP_CELLULAR_RESPIRATION`, `GOBP_OXIDATIVE_PHOSPHORYLATION`

축이 나온다.

즉, 기존에 정리했던
- extracellular matrix / adhesion / remodeling
- local jejunal structural transition
- trafficking / organelle / metabolism 계열 재배치

서술과 대체로 합치된다.

### 핵심 해석 3: G50가 가장 강하다는 기존 결론도 유지된다
`NC vs G50`는 gene-level에서도 강했고, GSEA에서도 pathway 수가 가장 많다.
이 점은 기존 결론과 일관적이다.

---

## 6. NC vs G25에서 특히 주목할 포인트
`NC_vs_G25_fgsea_go_bp_top30_sig.csv` 기준으로, 가장 먼저 정리할 항목은 아래 두 축이다.

### positive enrichment
- extracellular encapsulating structure organization
- collagen fibril organization
- basement membrane organization
- cell-substrate adhesion
- bone / branching structure morphogenesis 관련 일부 항목

### negative enrichment
- Golgi vesicle transport
- ER-to-Golgi transport
- vesicle organization
- macroautophagy / regulation of autophagy
- cellular respiration / ER stress 관련 항목

### 의미
이 패턴은 `NC vs G25`가 “아무것도 안 변했다”가 아니라,

- **구조적/ECM-adhesion 축은 위로 움직이고**
- **vesicle trafficking / 일부 proteostasis·metabolic 축은 아래로 움직이는**

재배치가 존재함을 시사한다.

다만 GO term 중복이 매우 크기 때문에, 논문 본문에는 그대로 20–30개를 나열하기보다
**중복 term cluster를 4–6개 narrative module로 묶는 작업**이 필요하다.

---

## 7. 이번 결과로 업데이트되는 결론
이번 턴까지 합치면, RNA-seq 쪽 결론은 더 명확해진다.

1. **main analysis는 full sample 유지가 맞다.**
2. **NC vs G25는 threshold-sensitive contrast다.**
3. **하지만 pathway-level로는 실제 신호가 있다.**
4. 따라서 다음 단계의 우선순위는 여전히
   - outlier 제거
   - 더 복잡한 supervised integration
   가 아니라,
   - **GSEA 기반 pathway condensation**
   - **matched 17 sample 기준 pathway–taxa integration**
   이다.

---

## 8. 바로 다음 추천 실행 순서
### 우선순위 1
`NC_vs_G25_fgsea_go_bp.csv`와 `NC_vs_G50_fgsea_go_bp.csv`에서
중복 GO term을 묶어 **condensed pathway module table**을 만든다.

예시 모듈:
- ECM / collagen / basement membrane
- adhesion / migration / vasculature-like remodeling
- vesicle trafficking / Golgi
- autophagy / proteostasis
- oxidative phosphorylation / energy metabolism

### 우선순위 2
이번에 만든 `matched_17_samples_integration_ready.csv`를 기반으로
**pathway–taxa integration table**을 만든다.

이때는
- RNA: contrast-level pathway module direction 또는 sample-level pathway score
- microbiome: cross-method-supported taxa only
  - 우선 `Romboutsia`
  - supportive: `Turicibacter`, `Clostridia UCG-014`, `Enterococcus`

구조로 가는 것이 가장 방어적이다.

### 우선순위 3
그 다음에만 exploratory로
- GSVA/ssGSEA sample score
- module-level correlation
- WGCNA
- sCCA / DIABLO
를 검토한다.

---

## Conclusion
이번 턴에서 확인된 가장 중요한 추가 사실은 이것이다.

> **NC vs G25는 gene-level hard DEG는 적지만, preranked GSEA에서는 분명한 pathway-level 신호가 존재한다.**

따라서 지금의 병목은 outlier 제거가 아니라,
**GSEA 결과를 narrative module로 정리하고 matched 17 samples 기준의 보수적 pathway–taxa integration으로 넘어가는 것**이다.

---

## Sources
- Local RNA-seq results directory: `file:///home/yzyzero/glycine_rnaseq_analysis_20260324`
- Local 16S bundle metadata: `file:///home/yzyzero/perplexity_glycine_layinghens_bundle/inputs_16s/integrated_metadata_clean.csv`
- Local matched-sample outputs: `file:///home/yzyzero/feynman/experiments/glycine-matched-samples`
- Local GSEA outputs: `file:///home/yzyzero/feynman/experiments/glycine-gsea-docker/results`
- Docker image documentation: https://hub.docker.com/r/rocker/r-ver
- fgsea package: https://bioconductor.org/packages/fgsea/
- msigdbr package: https://igordot.github.io/msigdbr/
