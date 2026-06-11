# Glycine 작업 로그 기반 해석: 현재까지 결과가 시사하는 jejunal biology와 재분석 우선순위

## Executive Summary
현재 workspace에 남아 있는 glycine 관련 로그와 산출물을 종합하면, 이 프로젝트의 핵심 해석은 비교적 분명하다. 현재 프로젝트의 기존 로그와 분석 문맥에서는 모든 treatment가 shared heat-stress 배경 위에 놓인 것으로 해석되어 왔고, 이 전제 아래에서는 glycine의 효과를 “heat stress를 없앴다”기보다는 **heat-stressed jejunum이 취하는 mucosal state를 얼마나, 어떤 방향으로 재배치했는가**로 읽는 편이 더 타당하다.[9][10] 내부 결과만 놓고 보면, **G50은 가장 강한 transcriptomic remodeling 상태**이고, **G25는 hard DEG 기준에서는 약하지만 pathway 수준에서는 분명한 신호가 있는 저진폭 상태**다.[3][4][6][9][10]

현재까지 가장 방어적으로 말할 수 있는 생물학적 해석은 다음과 같다. 이 데이터셋은 직접적인 barrier 기능이나 redox flux를 측정한 것은 아니지만, 외부 heat-stress poultry 및 intestinal epithelial 문헌과 나란히 읽으면 glycine-associated host state는 **redox buffering / cytoprotection / epithelial survival / 부분적인 barrier-support와 정합적인 방향**으로 해석할 수 있다.[11][12][13][17][18][19][20][21][22][23] 여기에 더해 transcriptome에는 **ECM·collagen·basement membrane·adhesion-like remodeling** 축이 반복적으로 보이고, 반대 방향으로는 **Golgi/vesicle trafficking, ER–Golgi transport, 일부 autophagy/proteostasis/energy metabolism** 축의 상대적 down-enrichment가 보인다.[6][7][8] 다만 이 두 번째 층의 해석, 특히 ECM-remodeling과 autophagy/vesicle trafficking은 아직 직접 관측이라기보다 **문헌으로 보강된 해석 가설**에 가깝다.[14][15][21][24]

따라서 지금 단계의 최적 전략은 outlier 제거 같은 구조 흔들기보다, **(1) pathway term condensation, (2) matched 17 samples 기반 pathway score 구축, (3) 짧은 taxa shortlist와의 보수적 pathway↔taxa association 분석**으로 넘어가는 것이다.[1][2][5][25][26][27][28][29][30] 반대로 phenotype 기반 개체 수준 통합 주장이나 DIABLO/sCCA 같은 고차원 잠재모형을 본문 주 분석으로 올리는 것은 현재 데이터 구조상 방어력이 약하다.[2][27][30]

## 1. 내부 로그가 이미 말해주고 있는 것

### 1.1 설계 해석의 중심은 heat stress 유무가 아니라 treatment 간 mucosal state 차이
현재 프로젝트의 기존 로그와 재분석 문맥에서는 이 데이터를 **shared heat-stress 배경 위의 treatment 비교**로 해석해 왔다.[9][10] 따라서 가장 적절한 질문은 “glycine이 heat stress를 완화했는가”보다, **NC, PC, G25, G50이 같은 stress 배경에서 서로 얼마나 유사하거나 다른 jejunal state에 도달했는가**이다.[9][10]

### 1.2 현재 내부 evidence에서 가장 일관된 구조
내부 산출물은 아래 세 가지를 일관되게 지지한다.[3][4][6][9][10]

1. **G50는 가장 큰 host transcriptome shift를 보인다.**  
   DEG 수, MDS 구조, preranked GSEA pathway 수 모두 이 결론과 맞는다.[3][4][6][9][10]
2. **G25는 null contrast가 아니라 threshold-sensitive contrast다.**  
   hard cutoff에서는 거의 안 남지만, near-threshold gene이 많고 pathway 수준에서는 유의한 신호가 존재한다.[3][6][7][10]
3. **이 데이터는 gene list보다 pathway/module 수준에서 읽어야 한다.**  
   특히 G25는 그 점이 더 분명하다.[3][6][10][26][28][29]

## 2. 현재까지의 결과에서 직접 관측된 것

### 2.1 샘플 구조와 해석 한계
- RNA-seq과 16S의 exact matched sample은 **총 17개**다.[1]
- 그룹별로는 **NC=4, PC=5, G25=3, G50=5**다.[1]
- 통합 메타데이터의 phenotype 관련 26개 변수는 모두 treatment group 내에서 상수였다.[2]

이 말은 곧, 현재 phenotype은 개체 수준 독립 측정값이 아니라 **treatment annotation에 가까운 값**이라는 뜻이다.[2] 따라서 phenotype correlation이나 phenotype-mediated mechanism을 강하게 주장하는 순간 pseudoreplication 문제가 생긴다.[2]

### 2.2 G25는 gene-level로 약하지만 signal-free는 아니다
내부 threshold sensitivity 결과는 `NC vs G25`가 아래 구조임을 보여준다.[3]

- `FDR < 0.05 & |log2FC| >= 1`: 4 genes
- `FDR < 0.10 & |log2FC| >= 1`: 40 genes
- `FDR < 0.10 & |log2FC| >= 0.5`: 96 genes
- `rawP < 0.05 & |log2FC| >= 1`: 485 genes

이 패턴은 G25가 “아예 반응이 없다”라기보다, **효과 크기가 얕고 분산/다중검정에 민감한 contrast**라는 뜻이다.[3][10] 즉, hard DEG만으로 G25를 해석하면 biology를 놓칠 가능성이 크다.[3][6][10]

### 2.3 GSEA는 G25의 pathway-level signal을 실제로 보여준다
실제로 preranked fgsea를 돌린 결과, FDR<0.05 기준 pathway 수는 다음과 같다.[6][10]

- `NC vs G25`: 182
- `NC vs PC`: 271
- `NC vs G50`: 813

이 결과는 두 가지를 동시에 말한다.[6][10][28]

1. **G50가 가장 큰 remodeling 상태**라는 기존 해석은 유지된다.[6][10]
2. **G25도 pathway-level로는 분명히 non-null**이다.[6][7][10]

### 2.4 반복되는 pathway 축
내부 GSEA 결과에서 contrast 전반에 걸쳐 가장 반복적으로 보이는 축은 아래와 같다.[6][7][8]

**양의 방향**
- external encapsulating structure organization
- collagen fibril organization
- basement membrane organization
- cell-substrate / cell-matrix adhesion

**음의 방향**
- Golgi vesicle transport
- ER-to-Golgi vesicle transport
- vesicle organization
- macroautophagy / regulation of autophagy
- cellular respiration / oxidative phosphorylation

이건 현재 host response가 단순히 “몇 개 stress gene이 올랐다/내렸다”가 아니라, **조직 구조·접착·matrix 관련 상태 전이**와 **secretory/proteostasis/energy-handling 상태 전이**를 함께 포함할 가능성을 시사한다.[6][7][8][14][15]

### 2.5 outlier 제거는 해결책이 아니다
outlier sensitivity 결과에서 full-sample 유지가 가장 방어적이었다.[5][9] 특정 샘플 제거는 MDS를 일부 정리해 보일 수 있지만, 주요 contrast 자체를 붕괴시키기도 했다.[5] 따라서 현재 병목은 outlier가 아니라 **해석 프레임과 차원 축소 방식**이다.[5][9][10]

## 3. 외부 문헌을 붙였을 때 가능한 생물학적 해석

### 3.1 가장 정합적인 해석 가설: glycine은 heat-stressed jejunum에서 cytoprotective / redox-buffering / partial barrier-support 쪽으로 읽힌다
이 해석은 현재 가장 강하다.[11][12][13][17][18][19][20][21][22][23]

직접 poultry heat-stress + glycine 문헌에서는 glycine이
- villus height,
- VH:CD,
- goblet cell,
- 일부 tight-junction/barrier-associated 지표,
- antioxidant status
를 개선하는 쪽으로 반복적으로 보고되었다.[11][12][13] 비가금 intestinal epithelial 모델과 piglet/LPS 모델에서도 glycine은
- glutathione 보존,
- ROS 감소,
- apoptosis 감소,
- mTOR-linked protein synthesis 유지,
- ER-stress-related barrier injury 완화
와 연결된다.[17][18][19][20][21][22][23]

따라서 현재 glycine host signal을 해석할 때 가장 먼저 놓을 수 있는 biological frame은 다음이다.

> glycine-associated host state는 **상피 세포 생존, redox buffering, 부분적 barrier maintenance, 점막 안정화와 정합적인 방향**으로 해석될 수 있다. 다만 이는 현재 코호트에서 직접 기능 검증된 결론이 아니라 외부 문헌을 붙인 해석이다.[11][12][13][17][18][19][20][21][22][23]

이 해석은 내부 결과의 G25/G50 pathway structure와 외부 장 상피·heat stress 문헌이 가장 잘 만나는 지점이다.[6][7][8][11][12][13][14][15][17][18][19][20][21]

### 3.2 그다음으로 강한 해석: proteostasis / ER-stress support
외부 문헌에서 heat-stressed chicken jejunum은 HSP, protein processing in ER, glutathione metabolism, inflammatory signaling 변화가 반복적으로 보고된다.[14][15][16] glycine 쪽 문헌도 intestinal cell 수준에서 mTORC1 의존적 ER-stress 완화, apoptosis 감소, barrier 유지와 연결된다.[17][20][21]

따라서 현재 내부 결과에서 보이는
- ER-to-Golgi transport,
- vesicle organization,
- response to ER stress,
- topologically incorrect protein response,
- macroautophagy-related term
같은 축은, 최소한 **proteostasis / secretory-state 관련 상대적 전사체 재배치**로는 읽을 수 있다.[6][7][8][14][15][21]

다만 여기서 주의할 점은, 이것이 곧바로 “autophagy가 활성/억제되었다” 또는 “vesicle trafficking이 좋아졌다/나빠졌다”는 직접 결론은 아니라는 점이다.[21] 세포 조성 변화, 성장 상태 변화, 번역 부담 차이 같은 대안 설명도 가능하다. 현 시점에서 더 안전한 문장은 아래다.

> glycine-associated host state는 heat-stressed jejunal mucosa의 **proteostasis and secretory-state signature가 상대적으로 다르게 배치된 상태**를 시사할 수 있다.[14][15][21]

### 3.3 조건부 해석: ECM / basement membrane / adhesion-like remodeling
내부 GSEA에서 가장 반복되는 positive theme은 ECM/collagen/basement membrane/adhesion이다.[6][7][8] 이건 무시하기 어려운 signal이다. 외부 문헌상 glycine은 collagen의 주요 구성 아미노산이며, mTORC1-ATF4-serine/glycine axis는 collagen biosynthesis 수요와 연결될 수 있다.[24] 또 glycine은 epithelial junction distribution 변화와 barrier remodeling과도 연결된다.[18]

그렇지만 여기에는 불확실성이 남아 있다. 현재 observed pathway가 의미하는 것이
- epithelial repair,
- matrix remodeling,
- vascular-like remodeling,
- barrier basement membrane reinforcement,
- 또는 broader tissue restructuring
중 무엇인지는 아직 분리되지 않았다.[6][7][8][24]

따라서 이 축은 다음처럼 쓰는 것이 가장 안전하다.

> 현재 glycine-responsive host transcriptome은 **ECM/adhesion/basement-membrane-like remodeling**을 강하게 시사하지만, 이것이 정확히 어떤 조직학적 과정인지까지는 아직 직접 확정할 수 없다.[6][7][8][18][24]

### 3.4 G25와 G50의 생물학적 차이: same direction 가능성(미검정), amplitude 차이는 관측된다
내부 결과를 보면 G25는 PC 쪽과 transcriptomic 위치가 가깝고, G50는 훨씬 큰 이동을 보인다.[4][9] 동시에 G25도 GSEA에서는 non-null이다.[6][7][10] 이 구조는 단순 선형 dose effect보다 아래 그림과 **양립 가능**하다.

- **G25**: 방향성은 존재하지만 진폭이 작고 threshold-sensitive한 “subtle coordinated adjustment” 가능성
- **G50**: 훨씬 더 넓은 host-state remodeling을 동반하는 “high-amplitude transition” 가능성

다만 현재는 G25와 G50의 방향 일치성을 정량 검정한 결과가 제시된 것은 아니므로, **둘 다 반응하되 G50가 더 강한 상태**라는 수준까지가 직접적으로 더 안전하다.[3][4][6][9][10]

## 4. observation 과 inference를 분리하면 어디까지 말할 수 있나

### 4.1 observation으로 말할 수 있는 것
1. G50는 G25보다 훨씬 큰 transcriptomic remodeling을 보인다.[3][4][6][9][10]
2. G25는 hard DEG는 적지만 pathway-level signal은 존재한다.[3][6][7][10]
3. positive pathway는 ECM/collagen/basement membrane/adhesion-like theme에 집중된다.[6][7][8]
4. negative pathway는 Golgi/vesicle trafficking, ER–Golgi transport, autophagy/proteostasis-like, energy metabolism theme에 반복된다.[6][7][8]
5. matched host–microbiome sample은 17개뿐이다.[1]
6. 현재 phenotype metadata는 개체 수준 통계에 직접 쓰기 어렵다.[2]

### 4.2 inference로만 말해야 하는 것
1. ECM/adhesion 축이 실제 barrier repair인지, stromal remodeling인지, vascular remodeling인지.[6][7][8][24]
2. vesicle/autophagy 축이 glycine에 의한 “개선”인지, 단순 state shift인지.[6][7][8][21]
3. G25가 실제로 PC mimicry인지, 아니면 부분적 convergence인지.[4][9]
4. 특정 taxa가 특정 host pathway를 매개한다는 인과적 그림.[1][26][27][30]
5. phenotype 개선이 host–microbe coordination의 결과라는 매개 모델.[2][27][30]

## 5. 현재 데이터가 시사하는 가장 타당한 biological model

### 5.1 제안 모델
heat-stressed jejunum에서 glycine은 우선 **상피 세포의 stress-handling capacity와 연관된 host state**를 조정했을 가능성이 있다. 그 결과는
- redox buffering,
- epithelial survival,
- 부분적 barrier/mucosal support,
- proteostasis / secretory-state signature의 변화
와 정합적인 방향으로 읽을 수 있다.[11][12][13][17][18][19][20][21][22][23] 이 변화가 충분히 커지면 transcriptome에서는 **ECM/adhesion/basement membrane-like remodeling**이 보이고, 동시에 **Golgi/vesicle/autophagy/energy metabolism 축의 상대적 재배치**가 동반될 수 있다.[6][7][8][14][15][21][24]

이 모델에서 G25는 같은 방향의 저강도 적응 상태일 가능성이 있고, G50는 더 넓고 깊은 고강도 적응 상태일 가능성이 있다. 다만 이는 현 단계에서 관측과 문헌을 합친 작동 가설이며, 방향 일치성 자체를 직접 검정한 결과는 아직 없다.[3][4][6][9][10]

### 5.2 왜 이 모델이 현재 가장 타당한가
- 내부 RNA-seq/GSEA 구조와 모순되지 않는다.[3][4][6][7][8][9][10]
- heat-stressed poultry jejunum literature와 맞는다.[11][12][13][14][15][16]
- glycine intestinal biology literature와도 대체로 맞는다.[17][18][19][20][21][22][23]
- phenotype 및 matched-n 제약을 넘어서 과도한 host–microbe 인과 해석을 하지 않는다.[1][2][27][30]

## 6. 다시 분석해야 하는 것

### 6.1 최우선 재분석
#### A. GO redundancy condensation
지금 GSEA term은 너무 중복적이다. raw term count를 그대로 서술하면 biology가 부풀려진다. 우선 해야 할 일은 `rrvgo` 또는 유사 도구로 term을 묶어 **4–6개의 narrative module**로 정리하는 것이다.[25]

권장 모듈 초안:
1. ECM / collagen / basement membrane
2. adhesion / migration / epithelial remodeling
3. Golgi / ER-to-Golgi / vesicle trafficking
4. proteostasis / ER stress / autophagy-adjacent response
5. oxidative phosphorylation / cellular respiration / energy metabolism

#### B. matched 17 samples 기반 pathway score 만들기
각 sample에 대해 위 모듈 점수를 만들어야 한다. GSVA/ssGSEA/mean z-score 등 몇 가지 정의를 비교하되, **짧은 pathway set**에서 안정적인지 먼저 확인하는 것이 좋다.[1][26]

#### C. pathway↔taxa integration table
그 다음에만 matched 17개에서
- host: condensed pathway score
- microbiome: prevalence/DA evidence가 있는 짧은 taxa shortlist
를 연결해야 한다. 여기서 중요한 것은 **low-dimensional planned analysis**여야 한다는 점이다.[1][26][27][30]

### 6.2 중간 우선순위 재분석
#### A. microbiome taxa shortlist 재정의
현재는 `Romboutsia`가 가장 견고한 candidate이고, `Turicibacter`, `Clostridia UCG-014`, `Enterococcus` 등은 supportive candidate로 보인다. 재분석에서는 abundance/prevalence/DA concordance를 기준으로 shortlist를 더 엄격히 고정해야 한다.[9][10][30]

#### B. pathway score 정의 민감도 점검
GSVA, ssGSEA, leading-edge aggregate, mean z-score 중 어떤 방식에서도 주요 결론이 유지되는지 확인해야 한다. 이건 reviewer 방어에 도움이 된다.[26]

#### C. correlation robustness
association는 permutation 또는 bootstrap으로 흔들어봐야 한다. n=17이기 때문에 부드럽고 예쁜 heatmap보다 **불안정성 자체를 같이 보여주는 것**이 더 정직하다.[26][27][30]

### 6.3 supplementary로만 고려할 것
#### A. DIABLO / sPLS / sCCA
이 방법들은 할 수는 있지만, 본문 주 분석 근거로 쓰기엔 표본 수가 너무 작다. pathway↔taxa association가 안정적일 때만 exploratory figure로 추가하는 것이 좋다.[27]

#### B. WGCNA
RNA-seq 단독 module 탐색은 가능하지만, phenotype와 직접 연결하거나 host–microbiome 전체의 강한 모듈 네트워크 주장을 하는 순간 과해석 위험이 커진다.[2][27][30]

### 6.4 지금은 피해야 하는 것
1. phenotype repeated values를 개체 수준 trait처럼 넣는 분석.[2][30]
2. gene-by-taxon full matrix correlation.[26][27][30]
3. microbiome ecological network를 강하게 해석하는 분석.[27][30]
4. host–microbe causal mediation narrative.[2][27][30]
5. outlier 제거를 main result로 재구성하는 시도.[5][9]

## 7. Open Questions
1. ECM/adhesion positive axis는 실제 조직학적으로 무엇을 의미하는가?[6][7][8][24]
2. vesicle/autophagy/ER-stress negative axis는 adaptive relief인가, 아니면 다른 형태의 burden shift인가?[6][7][8][21]
3. G25의 PC-like 위치는 실제 mimicry인가, 부분적 convergence인가?[4][9]
4. host pathway 변화와 microbiome taxa 변화 사이에 reproducible association가 존재하는가?[1][26][27][30]
5. phenotype 원자료가 복구되면 현재의 phenotype limitation이 해소되는가?[2]

## 8. Recommended Next Steps
1. GSEA 결과를 semantic-similarity 기반으로 **4–6 modules**로 condense한다.[25]
2. matched 17 samples에 대해 module별 **sample-level pathway score**를 만든다.[1][26]
3. microbiome에서 **짧은 taxa shortlist**를 고정한다.[27][30]
4. pathway score × taxa abundance 간 **보수적 association table**을 만든다.[26][27][30]
5. 그 결과가 안정적일 때만 DIABLO/sPLS류를 supplementary로 검토한다.[27]
6. phenotype 원자료가 있다면 repeated-value issue를 먼저 해결한 뒤 별도 통합을 다시 평가한다.[2][30]

## Sources
1. `matched_group_counts.csv` (local artifact URL): file:///home/yzyzero/feynman/experiments/glycine-matched-samples/matched_group_counts.csv
2. `phenotype_repeat_pattern_summary.csv` (local artifact URL): file:///home/yzyzero/feynman/experiments/glycine-matched-samples/phenotype_repeat_pattern_summary.csv
3. `deseq2_threshold_sensitivity_counts.csv` (local artifact URL): file:///home/yzyzero/feynman/experiments/glycine-rnaseq-summary/deseq2_threshold_sensitivity_counts.csv
4. `mds_group_dispersion_summary.csv` (local artifact URL): file:///home/yzyzero/feynman/experiments/glycine-rnaseq-summary/mds_group_dispersion_summary.csv
5. `outlier_key_scenarios.csv` (local artifact URL): file:///home/yzyzero/feynman/experiments/glycine-rnaseq-summary/outlier_key_scenarios.csv
6. `fgsea_go_bp_summary.csv` (local artifact URL): file:///home/yzyzero/feynman/experiments/glycine-gsea-docker/results/fgsea_go_bp_summary.csv
7. `NC_vs_G25_fgsea_go_bp_top30_sig.csv` (local artifact URL): file:///home/yzyzero/feynman/experiments/glycine-gsea-docker/results/NC_vs_G25_fgsea_go_bp_top30_sig.csv
8. `NC_vs_G50_fgsea_go_bp_top30_sig.csv` (local artifact URL): file:///home/yzyzero/feynman/experiments/glycine-gsea-docker/results/NC_vs_G50_fgsea_go_bp_top30_sig.csv
9. `outputs/glycine-omics-reanalysis.md` (local artifact URL): file:///home/yzyzero/feynman/outputs/glycine-omics-reanalysis.md
10. `outputs/glycine-gsea-followup.md` (local artifact URL): file:///home/yzyzero/feynman/outputs/glycine-gsea-followup.md
11. Deng et al., 2023, *Dietary glycine supplementation prevents heat stress-induced impairment of antioxidant status and intestinal barrier function in broilers*: https://pmc.ncbi.nlm.nih.gov/articles/PMC9827071/
12. Kwon et al., 2024, *Effect of increasing supplementation of dietary glycine on growth performance, meat quality, liver characteristics, and intestinal health in broiler chickens raised under heat stress conditions*: https://pmc.ncbi.nlm.nih.gov/articles/PMC11490916/
13. Won et al., 2023, *Effect of individual or combination of dietary betaine and glycine on productive performance, stress response, liver health, and intestinal barrier function in broiler chickens raised under heat stress conditions*: https://pmc.ncbi.nlm.nih.gov/articles/PMC10232888/
14. Zhu et al., 2024, *Transcriptome analysis of jejunal mucosal tissue in breeder hens exposed to acute heat stress*: https://pmc.ncbi.nlm.nih.gov/articles/PMC11617225/
15. Kim et al., 2022, *Integrated transcriptome analysis for the hepatic and jejunal mucosa tissues of broiler chickens raised under heat stress conditions*: https://jasbsci.biomedcentral.com/articles/10.1186/s40104-022-00734-y
16. Nanto-Hara et al., 2020, *Heat Stress Directly Affects Intestinal Integrity in Broiler Chickens*: https://pubmed.ncbi.nlm.nih.gov/33132728/
17. Wang et al., 2014, *Glycine stimulates protein synthesis and inhibits oxidative stress in pig small intestinal epithelial cells*: https://pubmed.ncbi.nlm.nih.gov/25122646/
18. Li et al., 2016, *Glycine Regulates Expression and Distribution of Claudin-7 and ZO-3 Proteins in Intestinal Porcine Epithelial Cells*: https://pubmed.ncbi.nlm.nih.gov/27029941/
19. Howard et al., 2010, *Glycine transporter GLYT1 is essential for glycine-mediated protection of human intestinal epithelial cells against oxidative damage*: https://pmc.ncbi.nlm.nih.gov/articles/PMC2849964/
20. Xu et al., 2018, *Glycine Relieves Intestinal Injury by Maintaining mTOR Signaling and Suppressing AMPK, TLR4, and NOD Signaling in Weaned Piglets after Lipopolysaccharide Challenge*: https://pmc.ncbi.nlm.nih.gov/articles/PMC6073676/
21. Yang et al., 2022, *Glycine represses endoplasmic reticulum stress-related apoptosis and improves intestinal barrier by activating mammalian target of rapamycin complex 1 signaling*: https://pmc.ncbi.nlm.nih.gov/articles/PMC8669258/
22. Zhang et al., 2022, *Protective effects of glycine against lipopolysaccharide-induced intestinal apoptosis and inflammation*: https://doi.org/10.1007/s00726-021-03011-w
23. Ji et al., 2022, *Glycine regulates mucosal immunity and the intestinal microbial composition in weaned piglets*: https://doi.org/10.1007/s00726-021-02976-y
24. Selvarajah et al., 2019, *mTORC1 amplifies the ATF4-dependent de novo serine-glycine pathway to supply glycine during TGF-β1-induced collagen biosynthesis*: https://pmc.ncbi.nlm.nih.gov/articles/PMC6584619/
25. `rrvgo` vignette: https://ssayols.github.io/rrvgo/articles/rrvgo.html
26. `GSVA` Bioconductor vignette: https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html
27. mixOmics DIABLO method page: https://mixomics.org/mixdiablo/
28. Wu & Smyth, 2012, *Camera: a competitive gene set test accounting for inter-gene correlation*: https://pmc.ncbi.nlm.nih.gov/articles/PMC3458527/
29. Wu et al., 2010, *ROAST: rotation gene set tests for complex microarray experiments*: https://pmc.ncbi.nlm.nih.gov/articles/PMC2922896/
30. ANCOM-BC2 vignette: https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html

## Downgraded or uncited claims
1. **ECM/remodeling 해석은 유지하되 강도를 낮췄다.** 내부 GSEA에서는 반복적으로 관찰되지만, 이것이 barrier repair인지 stromal/vascular remodeling인지까지는 직접 입증되지 않아 “강하게 시사” 수준으로 정리했다.[6][7][8][24]
2. **autophagy/vesicle trafficking 해석은 ‘개선’ 진술을 피했다.** 현재 자료와 문헌은 proteostasis/ER-stress burden 재조정까지는 지지하지만, autophagic flux나 trafficking efficiency의 직접 개선까지는 지지하지 않는다.[14][15][21]
3. **phenotype-mediated mechanism 관련 문구는 방어적으로 낮췄다.** 현재 phenotype 값은 group 내 반복값이라 개체 수준 독립 증거로 쓰기 어렵다.[2][30]
4. **G25의 ‘PC mimicry’ 표현은 확정하지 않았다.** 현재는 MDS상 proximity와 pathway-level non-null signal까지만 직접 지지되므로, mimicry 대신 partial convergence 가능성으로 남겼다.[4][9][10]
5. **직접 근거를 찾지 못한 문장은 남기지 않았다.** 본 brief에는 연구 노트에서 URL 또는 명시적 내부 artifact 근거를 확인하지 못한 주장을 별도로 남기지 않았다.
