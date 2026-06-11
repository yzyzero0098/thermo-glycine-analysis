# Glycine omics reanalysis update

작성일: 2026-04-02

## 요약
오늘 요청한 내용 중 **지금 이 환경에서 실제로 확인·수행 가능한 부분**을 우선 정리했다.

핵심 결론은 다음과 같다.

1. **오늘 질문 내용과 결과는 로그로 남겼다.**
2. RNA-seq은 이미 존재하는 결과와 outlier sensitivity sweep를 재검토했고, **주 분석에서는 전 샘플 유지(full)가 가장 방어 가능**했다.
3. `NC vs G25`에서 DEG가 거의 안 나오는 현상은 단순히 “테이블이 비어서”가 아니라,
   - 현재 기준인 `FDR < 0.05` 및 `|log2FC| >= 1`에서는 매우 적게 남지만,
   - **완화 기준에서는 near-threshold 유전자들이 꽤 존재**하는 구조였다.
4. 따라서 **G25에서 DEG가 적은 이유를 outlier 1–2개 제거로 해결하려는 접근은 주 분석으로는 방어력이 약하다.**
5. 기존 결과 기준으로는 `NC vs G50`에서 GO/KEGG 신호가 강하고, `NC vs G25`는 보수적 기준에서 약하다.
6. microbiome 쪽은 이미 상당한 downstream 결과가 존재하며, **DA는 방법 간 일치 여부를 중심으로 해석해야 한다.**
7. phenotype은 현재 확인된 메타데이터 기준으로 **반복된 group mean 패턴**이 많아서, phenotype 통합 분석은 강하게 주장하면 위험하다.
8. 통합분석은 가능하되, 현재 표본 구조에서는 **단순 상관/거리 기반 + pathway 단위 요약 + 보수적 모듈 접근**이 우선이고, **sCCA/DIABLO/TkNA/MintTea류는 과적합 방지 장치 없이는 본문 주 분석으로 바로 올리기 어렵다.**

---

## 1. 오늘 요청 내용 로그
오늘 사용자가 요청한 핵심은 아래였다.

- 오늘 질문한 내용과 결과를 로그로 남기기
- RNA-seq에서 기존 기준(`|logFC| >= 1`, `FDR < 0.05`) 재검토
- MDS 기준 샘플별/군집별 분포 확인
- outlier 제거를 시도해 batch effect 가능성을 점검
- 특정 비교에서 DEG가 잘 안 나오는 현상의 원인 점검
- 이어서 DEG, GO, KEGG, pathway 중심 functional analysis 정리
- 가능하면 GSEA도 포함
- microbiome 분석은 기존 flow와 기준을 최대한 유지하되, ALDEx2 / ANCOM-BC2 / MaAsLin2 같은 compositional DA 도구를 함께 보강
- phenotype과의 통합 가능성, host–microbe interaction, local jejunal 관점의 통합 전략 정리
- 과적합을 피하는 통합 모델 제안

이 요청 내용과 현재까지의 결과는 아래 파일들에 반영했다.

### 오늘 생성/갱신된 작업 로그
- `notes/session-logs/2026-04-02-glycine-omics-reanalysis.md`
- `experiments/glycine-rnaseq-summary/README.md`
- `outputs/glycine-omics-reanalysis.md` (이 문서)

---

## 2. 데이터 상태와 즉시 가능한 분석 범위

### 확인된 로컬 자산
- RNA-seq: `/home/yzyzero/glycine_rnaseq_analysis_20260324`
- 16S: `/home/yzyzero/glycine_16s_analysis`
- 통합 번들: `/home/yzyzero/perplexity_glycine_layinghens_bundle`
- phenotype 보조 폴더: `/home/yzyzero/glycine_phenotype`

### 현재 바로 가능한 것
- RNA-seq 결과/QC/DEG/pathway 재점검
- 16S diversity/DA/functional 결과 재점검
- RNA-seq–16S 샘플 crosswalk 작성
- host–microbe integration 전략 수립
- phenotype 제한점 평가

### 현재 바로 안 되는 것
- raw FASTQ부터의 full reprocessing
- DESeq2/edgeR를 현재 세션에서 처음부터 재실행하는 작업
  - 이유: 현재 기본 환경에 **R / Bioconductor가 설치되어 있지 않음**
  - 따라서 이번 업데이트는 **기존 결과 + 새 요약 계산 + 기존 sweep 검증** 중심이다.

이 점은 중요하다. 아래에서 말하는 RNA-seq 핵심 결론은 **기존에 이미 계산된 DESeq2/edgeR 결과 파일과 outlier sensitivity 결과를 재검토한 것**이며, 이 세션에서 de novo DESeq2를 다시 돌렸다고 말할 수는 없다.

---

## 3. RNA-seq 재점검

## 3.1 현재 기준과 baseline 확인
기존 보고 기준은 다음과 같다.

- **DEG 기준:** `FDR < 0.05`, `|log2FC| >= 1`
- 방법 비교: DESeq2 + edgeR

기존 결과 파일에서 baseline DEG 수는 아래와 같았다.

| 비교 | DESeq2 | edgeR |
|---|---:|---:|
| NC vs PC | 86 | 44 |
| NC vs G25 | 4 | 0 |
| NC vs G50 | 884 | 944 |

이 수치는 `DEG_method_comparison.csv`와 일치하는지 다시 확인했고, **일치했다.**

---

## 3.2 MDS 기준 샘플/군집 분포 확인
MDS 주좌표 파일(`MDS_main_figure_coordinates.csv`)을 이용해 군집 중심과 군내 분산을 재요약했다.

### 군집별 중심거리 요약
| 그룹 | n | 평균 centroid 거리 | 최대 centroid 거리 | nearest-centroid 정확도 | 기존 outlier flag 수 |
|---|---:|---:|---:|---:|---:|
| G25 | 6 | 0.243 | 0.487 | 0.667 | 0 |
| G50 | 5 | 0.619 | 1.092 | 0.600 | 0 |
| NC | 7 | 0.656 | 0.986 | 0.571 | 1 |
| PC | 7 | 0.317 | 0.659 | 0.429 | 2 |

### 군집 간 centroid 거리
| 그룹쌍 | centroid distance |
|---|---:|
| G50 vs NC | 0.973 |
| G25 vs NC | 0.705 |
| NC vs PC | 0.660 |
| G25 vs G50 | 0.538 |
| G50 vs PC | 0.464 |
| G25 vs PC | 0.106 |

### 해석
- **G25와 PC의 centroid 거리가 매우 짧다 (`0.106`)**. 즉, MDS 상에서 G25는 PC와 상당히 가까운 방향성을 보인다.
- 반면 **G50는 NC와 가장 멀다 (`0.973`)**. 즉, transcriptome 전체 구조상 G50 변화량이 가장 크다.
- PC와 G25가 MDS 상에서 가깝다는 점은, **“G25에서 DEG가 안 나온다 = 변화가 없다”**로 읽으면 안 된다는 뜻이다. 오히려 **NC 대비 변화량은 작지만, 방향은 PC-like일 가능성**이 있다.
- G50는 구조적으로 큰 이동이 보이므로, 현재 DEG 수가 많은 것도 MDS 구조와 일치한다.

즉, 현재 RNA-seq은 단순한 `NC -> G25 -> G50` 선형 증폭이라기보다,
- **G25 = PC 근접 but 작은 변화량**
- **G50 = 강한 전사체 재배치**
라는 해석이 더 맞다.

---

## 3.3 Outlier 제거 / batch effect 의심 검토
기존 outlier screening과 exhaustive removal sweep를 다시 확인했다.

### MDS 기반 candidate
기존 outlier check에서 주로 거론된 샘플은 아래였다.
- `P1-2-JM-NC-R`
- `P1-2-JM-PC-R`
- `P1-5-JM-PC-R`

확장 sweep에서는 `P1-7-JM-G50-R` 제거 시나리오도 함께 강하게 검토되었다.

### 기존 QC 요약
후보 샘플들은 trimming / unique mapping / overall alignment가 대체로 모두 양호했다.
즉, **명백한 technical failure sample로 자동 제외할 근거는 약했다.**

### 핵심 제거 시나리오 결과
| 시나리오 | 권고 | 핵심 영향 |
|---|---|---|
| full | **best retained** | 기준선. 가장 방어 가능 |
| drop `P1-5-JM-PC-R` | sensitivity-only | MDS relief는 있으나 주분석 제외 근거는 약함 |
| drop `P1-2-JM-PC-R` | not defensible | `NC vs PC` edgeR DEG가 사실상 붕괴 |
| drop `P1-2-JM-NC-R` | not defensible | `NC vs PC`와 `NC vs G50` 둘 다 약화 |
| drop `P1-7-JM-G50-R` | not defensible | `NC vs G50` 신호가 급격히 붕괴 |

### 특히 중요한 점
- `P1-7-JM-G50-R` 제거 시 `NC vs G50` edgeR DEG가 **944 -> 82**로 급감했다.
- `P1-2-JM-PC-R` 제거 시 `NC vs PC` edgeR DEG가 **44 -> 0**이 되었다.
- `P1-2-JM-NC-R` 제거 시 `NC vs PC` edgeR DEG가 **44 -> 5**로 줄었다.

### 결론
**이번 데이터에서 outlier 제거는 주 분석 전략으로 채택하기 어렵다.**

왜냐하면:
1. QC 실패 샘플이라고 단정할 근거가 약하고,
2. 제거 시 biological contrast 자체가 붕괴하거나 크게 왜곡되며,
3. “batch correction”이 아니라 **원래 있던 생물학적 신호를 잘라내는 결과**로 해석될 가능성이 크기 때문이다.

따라서 가장 정직한 결론은 다음이다.

> **주 분석은 full sample 유지로 두고, outlier removal은 sensitivity analysis로만 제시하는 것이 적절하다.**

---

## 3.4 왜 `NC vs G25`에서 DEG가 적은가?
이 질문이 이번 RNA-seq 파트의 핵심이었다.

기존 기준(`FDR < 0.05`, `|log2FC| >= 1`)에서는:
- DESeq2: 4개
- edgeR: 0개

하지만 threshold sensitivity를 다시 계산해보면:

| 비교 | FDR<0.05 & |log2FC|>=1 | FDR<0.10 & |log2FC|>=1 | FDR<0.05 & |log2FC|>=0.5 | FDR<0.10 & |log2FC|>=0.5 |
|---|---:|---:|---:|---:|
| NC vs G25 | 4 | 40 | 7 | 96 |
| NC vs PC | 86 | 175 | 429 | 669 |
| NC vs G50 | 884 | 1188 | 1796 | 2395 |

### 해석
이건 매우 중요하다.

`NC vs G25`는
- 완전히 신호가 없는 비교가 아니라,
- **보수적 다중검정 기준에서 살아남는 유전자가 적은 구조**다.

즉,
- 효과크기가 작거나,
- 표본 수 대비 분산이 크거나,
- 여러 유전자가 **near-threshold**에 몰려 있어
- 최종 FDR 0.05를 통과하지 못하는 형태다.

### G25의 near-threshold 예시
상위 유전자들을 보면,
- `IL15`는 기준 통과
- 여러 유전자(`RSPH14`, `STEAP1`, `SYAP1`, `MPHOSPH6`, `ZBTB14` 등)는 **padj 0.05~0.07대**에 위치
- 일부는 `|log2FC|`가 작아서 현재 cutoff에서 제외

즉, G25 문제는 “데이터 이상”이라기보다,

> **현재 cutoff와 표본 구조에서 신호가 얕고 분산에 민감한 contrast**

로 보는 쪽이 더 타당하다.

### 이번 단계에서의 권고
논문용 main DEG는 기존 기준을 유지하되,
추가로 아래를 같이 제시하는 것이 좋다.

1. **main table:** `FDR < 0.05`, `|log2FC| >= 1`
2. **supplementary sensitivity:** `FDR < 0.10`, `|log2FC| >= 1`
3. **ranked pathway 분석:** threshold 의존성을 줄이기 위해 preranked 방식(GSEA/fgsea 계열)을 병행

즉, G25는 gene-level hard cutoff보다 **pathway-level signal**이 더 중요할 가능성이 높다.

---

## 3.5 GO / KEGG existing result 해석
현재 확인된 RNA-seq enrichment 파일은 `NC vs G50`에 대해 존재했다.

### GO BP 상위 신호 (`NC vs G50`)
상위 항목은 대체로 다음 축에 몰린다.
- anatomical structure morphogenesis
- cell adhesion
- circulatory system / vasculature development
- cell migration / cell motility
- angiogenesis
- axon guidance / neuron projection guidance

### KEGG 상위 신호 (`NC vs G50`)
상위 항목은 대체로 다음과 같다.
- Cytoskeleton in muscle cells
- Cell adhesion molecule interaction
- Hormone signaling
- TGF-beta signaling pathway
- Integrin signaling
- Efferocytosis

### 해석
`NC vs G50`는 단순한 “stress marker 몇 개 변화”보다는,
**morphogenesis / adhesion / vascular-development-like remodeling / signaling reprogramming** 쪽의 전사체 재배치가 강하다.

이건 현재 glycine 프로젝트를
- 단순 DEG 나열형보다는
- **pathway-centered / tissue remodeling / local host-state transition**
으로 서술해야 한다는 기존 방향과 잘 맞는다.

---

## 3.6 GSEA 관련 판단
이번 세션에서 **새로운 GSEA 결과 파일은 확인되지 않았다.**

따라서 지금 말할 수 있는 건 두 가지다.

1. **GSEA 결과가 이미 저장돼 있다고 확인된 것은 아니다.**
2. 그러나 현재 구조상 G25에서는 오히려 **preranked GSEA가 꼭 필요**하다.
   - 이유: hard threshold DEG가 너무 적어서 ORA만으로는 정보 손실이 큼

### 다만 현재 blocker
- 현재 세션에는 R/Bioconductor 기반 재실행 환경이 없음
- 따라서 `fgsea`, `clusterProfiler GSEA`, `camera`, `roast` 등은 이번 턴에서 직접 실행하지 못했다

### 따라서 추천
다음 실행 우선순위에서
- **DESeq2 statistic / Wald stat 기반 preranked GSEA**
를 최우선 추가 항목으로 두는 것이 맞다.

---

## 4. 16S microbiome 분석 현황과 판단
사용자가 제시한 분석 flow는 이미 현재 로컬 결과와 상당 부분 일치한다.

### 이미 확인된 기존 16S 분석 축
- alpha diversity
- beta diversity (PCoA / PERMANOVA)
- taxonomy composition
- genus-level differential abundance (KW)
- ANCOM-BC 계열
- MaAsLin2
- PC mimicry / glycine-dose trend
- PICRUSt2 기반 기능 예측
- phenotype correlation / Mantel

즉, microbiome 쪽은 “새로 설계만 해야 하는 상태”가 아니라,
**이미 여러 차례 downstream analysis가 돌아간 상태**다.

---

## 4.1 Diversity와 community-level 결과
현재 읽은 표 기준으로 전역 Bray-Curtis PERMANOVA는:
- `p = 0.157`

즉, community 전체 구조 차이는 **강하게 유의하다고 보기 어렵다**.

반면 genus-level 비모수 검정에서는 상위 후보가 존재한다.
예:
- `Lactobacillus`
- `Enterococcus`
- `Rothia`
- `Romboutsia`

이 구조는 흔하다.
즉,
- 전체 community shift는 약하거나 분산에 묻히지만,
- 일부 주요 taxa는 treatment-responsive일 수 있다.

따라서 논문에서는
- “global beta diversity restored” 같은 강한 문장보다는
- **specific taxa restructuring** 중심 서술이 더 안전하다.

---

## 4.2 DA cross-method 요약
기존 `DA_cross_method_summary.md`를 확인한 결과:

### 반복적으로 상위에 등장한 taxa
- `Romboutsia`
- `Turicibacter`
- `Clostridia UCG-014`
- `Muribaculaceae`
- `Bacteroides`
- `Lachnospiraceae NK4A136 group`
- `Enterococcus`
- `Gallibacterium`

### 방법별 해석
- **ANCOM-BC 계열**: 비교적 공격적으로 많은 신호를 반환
- **MaAsLin2**: 훨씬 보수적이며, 현재 데이터에서는 FDR 기준에서 강한 통과 taxa가 많지 않음

### 가장 강한 메시지
현재 로컬 결과만 놓고 가장 방어적으로 쓸 수 있는 문장은:

> **`Romboutsia`는 glycine dose-related pattern이 가장 견고한 후보다.**

그 외:
- `Turicibacter`
- `Clostridia UCG-014`
- `Enterococcus`

은 **방법 의존적인 supportive candidate**로 쓰는 것이 적절하다.

---

## 4.3 phenotype 통합의 가장 큰 제한점
현재 `phenotype_integrity_check.csv`와 통합 메타데이터를 보면, 다수 phenotype이
- 각 treatment group에서 unique value가 1개
- 즉 **반복된 group mean 패턴**
을 보인다.

이건 매우 중요하다.

### 의미
지금 phenotype–microbiome 상관이나 Mantel 결과는
- **individual-level association**이 아니라
- 사실상 **group-level ecological association**에 가깝다.

### 따라서 무엇이 위험한가?
- WGCNA module–phenotype correlation을 강하게 주장하는 것
- sample-level host–microbe–phenotype mediation을 단정하는 것
- phenotype 기반 supervised multiblock model을 “예측 모델”처럼 내세우는 것

### 무엇은 가능한가?
- exploratory integration
- treatment-level trend 해석
- stress/liver/serum domain 단위의 보조 해석

즉, phenotype 파트는 **강한 인과 주장의 기반이 아니라 context layer**로 쓰는 게 안전하다.

---

## 5. transcriptome–microbiome 샘플 매칭 상태
RNA-seq 25개와 16S 19개의 sample ID를 정규화해서 비교한 결과:

- RNA-seq 총 25
- 16S 총 19
- **정확히 매칭되는 공통 샘플 17개**

### matched sample 수
| 그룹 | matched n |
|---|---:|
| NC | 4 |
| PC | 5 |
| G25 | 3 |
| G50 | 5 |

### 해석
이건 통합분석 설계에서 결정적이다.

- G25 matched sample이 **3개**뿐이라서,
- 복잡한 supervised integration은 과적합 위험이 매우 높다.

즉,
- `sCCA`, `DIABLO`, `MintTea`, `TkNA` 같은 방법을 바로 메인 분석으로 넣기보다,
- **feature 수를 강하게 줄인 pathway/taxa 단위 요약 + permutation 검증**이 선행돼야 한다.

---

## 6. 통합분석 전략: 지금 “바로 가도 되는 것” vs “조심할 것”

## 6.1 바로 가도 되는 분석 (Go)
### A. Pathway-level transcriptome 요약
가장 먼저 추천한다.

이유:
- gene-level보다 차원이 낮아짐
- G25 같은 약한 contrast에서도 신호를 살리기 좋음
- microbiome과 연결할 때 해석이 쉬움

권장 예:
- GO/KEGG/Reactome gene set score
- GSVA/ssGSEA/fgsea leading-edge summary
- top pathway eigengene 또는 module score

### B. Taxa–pathway correlation
- 입력: 소수의 상위 pathway score + 상위 treatment-responsive taxa
- 방법: Spearman + permutation / leave-one-out robustness
- 논문 톤: **exploratory host–microbe coordination**

### C. Mantel / Procrustes / RV-like distance association
- transcriptome: pathway-level 또는 VST/PCA space
- microbiome: Bray-Curtis / Weighted UniFrac
- treatment stratified 또는 all matched samples

이건 현재 표본 구조에서 비교적 안전하다.

### D. PC mimicry / convergence-divergence score
이 프로젝트에 가장 잘 맞는다.

예:
- RNA-seq pathway별로 `PC vs NC`와 `G25/G50 vs NC` 방향 비교
- microbiome taxa별로 동일 비교
- “PC mimicry score”와 “divergent compensation score”를 분리

이건 이 논문의 고유한 framing과도 가장 잘 맞다.

---

## 6.2 조건부로 가능한 분석 (Caution)
### A. WGCNA
**Transcriptome 단독 WGCNA**는 가능하다.
하지만 주의점이 있다.

- 25개 샘플이면 아주 넉넉하진 않다
- trait 연결은 phenotype raw individual 값이 아니라면 약해진다
- 따라서 추천은:
  - **RNA-seq module eigengene 생성**
  - module–taxa association
  - module–treatment / glycine dose / PC proximity association

즉,
- **WGCNA on transcriptome**: 가능
- **WGCNA on phenotype-integrated multiblock**: 현재는 비권장

### B. sCCA / DIABLO / block.sPLS
이론적으로는 매우 매력적이다.
하지만 현재 matched sample이 17개이고 G25는 3개다.

따라서 하려면 반드시:
- feature prefiltering 강하게 수행
- pathway/taxon 단위로 축소
- repeated CV / permutation 기반 성능 검증
- tuning 결과의 stability 확인
- main finding이 아니라 **sensitivity / exploratory figure**로 배치

즉,
> 할 수는 있지만 지금 표본 구조에서 메인 근거로 쓰면 위험하다.

### C. TkNA / MintTea
최근 방법론적으로는 흥미롭다.
하지만 현재 데이터셋에는 다음 문제가 있다.
- matched n이 작음
- phenotype raw sample-level이 아님
- external validation cohort 없음

따라서 이 방법들은 **지금 당장 본문 메인 분석으로 넣기보다**,
- 후속 확장 분석이나
- supplementary exploratory network
로 제한하는 게 좋다.

---

## 6.3 지금은 피하는 것이 좋은 것 (Stop)
### A. phenotype를 강한 supervised label로 쓴 과적합 모델
반복된 group mean이 많기 때문에,
- random forest / PLS-DA / sparse classifier 등으로 phenotype 예측을 강하게 주장하면 위험하다.

### B. mediation / causal chain 강한 주장
현재 구조로는
- glycine -> microbiome -> transcriptome -> phenotype
같은 경로를 인과적으로 확정할 수 없다.

가능한 표현은:
- associated with
- consistent with
- suggests coordinated host–microbe remodeling
정도다.

---

## 7. 최근 방법 논문을 반영한 추천
이번 세션에서 최근 방법론/응용 논문도 확인했다.

### 실질적으로 참고할 만한 축
1. **MintTea (Nature Communications 2024)**
   - sGCCA 기반 intermediate integration + 반복 subsampling + consensus module
   - 장점: 다중 오믹 coherent module 발굴
   - 단점: 작은 표본에서 tuning과 stability 검증이 필수
   - 현재 프로젝트 적용: **직접 본문 주 분석보다는 설계 참고용**

2. **TkNA (Nature Protocols 2024)**
   - host–microbiota network 기반 causal candidate prioritization 프레임
   - 장점: transkingdom network 관점이 명확
   - 단점: 표본 수, 재현성, 네트워크 안정성 요구가 높음
   - 현재 적용: **후속 네트워크 확장용**

3. **Chicken duodenal host–microbiota integration paper (Microbiome 2025)**
   - duodenal mucosal genes와 microbiota의 coordinated role을 chicken에서 직접 다룸
   - 장점: poultry + mucosal + host–microbe라는 점에서 컨셉이 매우 가깝다
   - 현재 적용: **문장 framing과 local interaction justification에 매우 유용**

### 이 프로젝트에 맞춘 현실적 결론
최근 유수 저널 방법을 그대로 끌어오더라도,
현재 데이터에서는 아래 순서가 가장 낫다.

1. **pathway-level host state 정의**
2. **microbiome treatment-responsive taxa 정의**
3. **matched 17 samples 기준 simple robust association**
4. **module/network는 exploratory로만 확장**

즉,
- 먼저 단순하고 방어적인 분석으로 골격을 세우고,
- 그 다음 복잡한 통합 모델은 부가 분석으로 넣는 것이 최선이다.

---

## 8. 내가 권하는 다음 실제 실행 순서

### Phase 1 — RNA-seq main analysis 정리
1. full sample 유지
2. main DEG: `FDR < 0.05`, `|log2FC| >= 1`
3. supplementary threshold sensitivity 표 추가
4. `NC vs G25`는 **hard DEG 부족 contrast**로 규정
5. `NC vs G50` GO/KEGG는 pathway 중심으로 서술
6. **preranked GSEA 추가 실행**

### Phase 2 — transcriptome pathway layer 만들기
1. DESeq2 stat 기반 preranked GSEA
2. GSVA/ssGSEA 또는 pathway score
3. pathway별 PC mimicry / dose-response 요약

### Phase 3 — microbiome main analysis 정리
1. 기존 alpha/beta/taxonomy 결과 재정리
2. DA는 KW + ANCOM-BC + MaAsLin2의 **intersection / directional concordance** 중심
3. `Romboutsia`를 가장 강한 dose-responsive candidate로 우선 배치
4. `Turicibacter`, `Clostridia UCG-014`, `Enterococcus`는 supportive candidate

### Phase 4 — integration
1. matched 17 samples만 사용
2. pathway score × taxa correlation
3. Mantel / Procrustes
4. PC mimicry / convergence-divergence matrix
5. 이후 필요시 transcriptome WGCNA eigengene와 taxa 연계

### Phase 5 — only if reviewers need it
1. sparse CCA / DIABLO exploratory
2. TkNA 또는 module network
3. permutation / CV / stability plot을 반드시 동반

---

## 9. 이번 업데이트의 한계
- 이번 세션에서는 **기존 결과 검증 + 새 요약 계산** 중심이었다.
- 현재 환경에는 R/Bioconductor가 없어 **DESeq2/edgeR/GSEA를 처음부터 재실행하지 못했다.**
- phenotype는 현재 확인된 파일 기준으로 반복된 group mean 패턴이 많아, 강한 sample-level 통합 해석은 제한된다.
- 따라서 지금 문서의 핵심 가치는:
  - 현재 결과를 어떤 톤으로 써야 하는지,
  - 어떤 분석은 main으로 가고 어떤 분석은 supplementary로 내려야 하는지,
  - 어디서 과적합과 과해석이 생기는지,
  를 정리했다는 데 있다.

---

## 권장 다음 단계
가장 먼저 바로 이어서 할 작업은 다음 두 가지다.

1. **RNA-seq preranked GSEA 실행 환경 확보 후 실제 실행**
2. **matched 17 samples 기준 pathway–taxa integration 테이블 작성**

이 두 단계가 끝나면,
- transcriptome–microbiome–phenotype을 단순 병렬 결과가 아니라
- **local jejunal host–microbe coordination**으로 더 강하게 묶을 수 있다.

---

## Sources
### Local files
- `/home/yzyzero/glycine_rnaseq_analysis_20260324/tables/sample_metadata.csv`
- `/home/yzyzero/glycine_rnaseq_analysis_20260324/mds_main_figure/MDS_main_figure_coordinates.csv`
- `/home/yzyzero/glycine_rnaseq_analysis_20260324/outlier_sensitivity_sweep/outlier_sensitivity_report.md`
- `/home/yzyzero/glycine_rnaseq_analysis_20260324/outlier_sensitivity_sweep/scenario_ranking.csv`
- `/home/yzyzero/glycine_rnaseq_analysis_20260324/tables/DEG_method_comparison.csv`
- `/home/yzyzero/glycine_rnaseq_analysis_20260324/tables/NC_vs_G25_DESeq2_full.csv`
- `/home/yzyzero/glycine_rnaseq_analysis_20260324/tables/NC_vs_G50_GO_BP.csv`
- `/home/yzyzero/glycine_rnaseq_analysis_20260324/tables/NC_vs_G50_KEGG.csv`
- `/home/yzyzero/perplexity_glycine_layinghens_bundle/inputs_16s/integrated_metadata_clean.csv`
- `/home/yzyzero/glycine_16s_analysis/phyloseq_results_20260323/tables/beta_diversity_permanova.csv`
- `/home/yzyzero/glycine_16s_analysis/phyloseq_results_20260323/tables/genus_differential_kruskal.csv`
- `/home/yzyzero/glycine_16s_analysis/da_functional_phenotype_20260323/DA_cross_method_summary.md`
- `/home/yzyzero/glycine_16s_analysis/da_functional_phenotype_20260323/integration_status_report.md`
- `/home/yzyzero/glycine_16s_analysis/da_functional_phenotype_20260323/tables/phenotype_integrity_check.csv`
- `experiments/glycine-rnaseq-summary/README.md`
- `experiments/glycine-rnaseq-summary/deseq2_threshold_sensitivity_counts.csv`
- `experiments/glycine-rnaseq-summary/mds_group_dispersion_summary.csv`
- `experiments/glycine-rnaseq-summary/outlier_key_scenarios.csv`
- `experiments/glycine-rnaseq-summary/NC_vs_G25_top30_by_padj.csv`

### Web / paper sources
- MintTea / multi-omic module integration: https://www.nature.com/articles/s41467-024-46888-3
- TkNA protocol: https://www.nature.com/articles/s41596-024-00960-w
- Review on microbiome multi-omics methods: https://pmc.ncbi.nlm.nih.gov/articles/PMC11739165/
- Chicken duodenal mucosal gene–microbiota integration paper: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-025-02054-5
