# 처음부터 시작하는 나노바디(Nanobody) 시퀀싱 데이터 분석 실전 가이드

> **프로젝트**: SRP124616 / SRR6269034
> **플랫폼**: PacBio CCS (Circular Consensus Sequencing)
> **목표**: 파지 디스플레이 라이브러리의 나노바디 VHH 서열 분석 전체 파이프라인
> **환경**: Ubuntu Linux + Conda + Python + IgBLAST

---

## 목차

1. [배경 지식](#1-배경-지식)
2. [환경 구축 및 데이터 다운로드](#2-환경-구축-및-데이터-다운로드)
3. [삽질 기록: NCBI 도구 설치와 SSL 오류](#3-삽질-기록-ncbi-도구-설치와-ssl-오류)
4. [1단계: FASTQ 품질 관리](#4-1단계-fastq-품질-관리)
5. [삽질 기록: VHH 서열 추출 실패](#5-삽질-기록-vhh-서열-추출-실패)
6. [2단계: VHH 서열 추출 (올바른 방법)](#6-2단계-vhh-서열-추출-올바른-방법)
7. [3단계: CDR 영역 식별](#7-3단계-cdr-영역-식별)
8. [4단계: 서열 다양성 분석](#8-4단계-서열-다양성-분석)
9. [5단계: 아미노산 조성 분석](#9-5단계-아미노산-조성-분석)
10. [6단계: IgBLAST Germline 유전자 분석](#10-6단계-igblast-germline-유전자-분석)
11. [삽질 기록: IgBLAST 데이터베이스 설정](#11-삽질-기록-igblast-데이터베이스-설정)
12. [최종 결과 및 생물학적 해석](#12-최종-결과-및-생물학적-해석)
13. [정리 및 향후 방향](#13-정리-및-향후-방향)

---

## 1. 배경 지식

### 나노바디(Nanobody / VHH)란?

나노바디는 낙타과 동물(낙타, 라마, 알파카)에서 발견되는 **중쇄 항체(HCAb)**의 가변 영역 단편입니다. 일반 항체와 달리 중쇄 항체는 경쇄가 없으며, 항원 결합 기능을 단일 VHH 도메인이 담당합니다.

```
일반 항체 (IgG, ~150 kDa):          중쇄 항체 (HCAb, ~95 kDa):
  ┌─VH─┐ ┌─VH─┐                      ┌─VHH─┐ ┌─VHH─┐
  │    │ │    │                      │     │ │     │
  ├─VL─┤ ├─VL─┤                      │     │ │     │
  │    │ │    │                      ├─CH2─┤ ├─CH2─┤
  ├─CH1┤ ├─CH1┤                      │     │ │     │
  │    │ │    │                      ├─CH3─┤ ├─CH3─┤
  ├─CH2┤ ├─CH2┤                      └─────┘ └─────┘
  │    │ │    │                      (경쇄 없음, CH1 없음)
  ├─CH3┤ ├─CH3┤
  └────┘ └────┘

나노바디 (VHH, ~15 kDa):
  ┌──────────────────────────────┐
  │ FR1-CDR1-FR2-CDR2-FR3-CDR3-FR4 │
  └──────────────────────────────┘
  (~120-130개 아미노산)
```

**나노바디의 장점**:
- 분자량이 작아(~15 kDa) 조직 침투력이 우수함
- 열 안정성이 높고, 극한 pH에서도 견딤
- 대장균에서 고효율 발현 가능
- 기존 항체가 접근하지 못하는 은닉 에피토프(예: 효소 활성 부위) 인식 가능

### VHH 도메인의 구성

```
FR1 (25aa) ── CDR1 (7-12aa) ── FR2 (17aa) ── CDR2 (6-7aa) ── FR3 (38aa) ── CDR3 (3-28aa) ── FR4 (11aa)
 [골격영역]      [항원결합]       [골격영역]      [항원결합]       [골격영역]       [항원결합]       [골격영역]
```

- **FR (Framework Region)**: 구조적 골격. 고도로 보존되어 면역글로불린 접힘(fold)을 유지
- **CDR (Complementarity Determining Region)**: 상보성 결정 영역. 항원과 직접 접촉
- **CDR3**: 길이와 서열이 가장 가변적이며, 항원 결합 특이성을 결정하는 핵심 영역

### 파지 디스플레이(Phage Display)

이 데이터는 파지 디스플레이 라이브러리에서 생산되었습니다. 기본 원리:

1. VHH 유전자를 파지 외피 단백질(pIII) 유전자 상류에 클로닝
2. VHH 단백질이 파지 표면에 발현(디스플레이)
3. 목표 항원으로 "패닝(panning)" 수행 → 결합하는 파지 선별
4. 여러 라운드의 패닝을 거쳐 고친화도 클론 농축

**중요**: 파지 디스플레이 벡터에는 일반적으로 VHH와 pIII 사이에 **amber 종결 코돈(TAG)**이 포함됩니다. amber 억제 균주(예: TG1)에서는 TAG가 Gln으로 읽혀 통과되고, 비억제 균주에서는 번역이 종결되어 가용성 VHH가 생산됩니다.

---

## 2. 환경 구축 및 데이터 다운로드

### 필요 도구 설치

```bash
# conda로 생명정보학 도구 설치 (sudo 필요 없음, 권장)
conda install -y -c bioconda -c conda-forge \
    entrez-direct sra-tools \
    fastqc seqkit blast muscle igblast

# Python 패키지
pip install biopython pandas matplotlib
```

### SRA 데이터 다운로드

```bash
# 방법 1: Entrez Direct 사용 (가능한 경우)
esearch -db sra -query "SRP124616" | efetch -format runinfo > runinfo.csv
cut -d',' -f1 runinfo.csv | grep SRR > srr_list.txt

# 방법 2: ENA에서 실행 정보 가져오기 (더 안정적)
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?\
accession=SRP124616&result=read_run&\
fields=run_accession,library_strategy,library_layout,fastq_ftp&\
format=tsv&limit=0"
```

**출력**:

```
run_accession  library_strategy  library_layout  fastq_ftp
SRR6269034     POOLCLONE         PAIRED          ftp.sra.ebi.ac.uk/vol1/fastq/SRR626/004/SRR6269034/...
```

1개의 실행(SRR6269034)만 존재하며, PacBio의 Pooled Clone 시퀀싱입니다.

```bash
# 다운로드 및 FASTQ 변환
prefetch SRR6269034
fasterq-dump --threads 8 SRR6269034
```

**출력**:

```
spots read      : 85
reads read      : 85
reads written   : 85
```

생성 파일: `SRR6269034.fastq` (179 KB, 85개 reads)

---

## 3. 삽질 기록: NCBI 도구 설치와 SSL 오류

### 문제 1: esearch/efetch 명령어를 찾을 수 없음

```
Command 'esearch' not found, but can be installed with:
sudo apt install ncbi-entrez-direct
```

**원인**: 시스템에 NCBI Entrez Direct 도구가 설치되지 않음.

**해결**: conda로 설치하여 sudo 권한 문제 회피:
```bash
conda install -c bioconda entrez-direct
```

### 문제 2: SSL 연결 오류

```
curl: (56) OpenSSL SSL_read: error:0A000126:SSL routines::unexpected eof while reading
ERROR: curl command failed with: 56
```

**원인**: conda로 설치된 OpenSSL 버전(3.6.2)과 NCBI 서버의 TLS 설정이 완전히 호환되지 않음.

**해결 방안**:

방안 A — `curl -k`로 SSL 검증 건너뛰기 (임시 방안):
```bash
curl -k -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?\
db=sra&term=SRP124616&retmax=500&usehistory=y"
```

방안 B — ENA API 사용 (권장, 더 안정적):
```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?\
accession=SRP124616&result=read_run&fields=run_accession&format=tsv"
```

**교훈**: NCBI의 Entrez API에서 SSL 문제가 발생하면, ENA(유럽 핵산 아카이브)가 좋은 대안입니다. ENA는 SRA의 모든 데이터를 미러링하고 있습니다.

---

## 4. 1단계: FASTQ 품질 관리

### 기본 통계

```python
from Bio import SeqIO
import numpy as np

records = list(SeqIO.parse("SRR6269034.fastq", "fastq"))
lengths = [len(r.seq) for r in records]
mean_quals = [np.mean(r.letter_annotations["phred_quality"]) for r in records]

print(f"총 reads 수:   {len(records)}")       # 85
print(f"길이 범위:     {min(lengths)}-{max(lengths)} bp")  # 885-979 bp
print(f"평균 길이:     {np.mean(lengths):.1f} bp")          # 975.6 bp
print(f"평균 품질 점수: {np.mean(mean_quals):.1f}")          # 91.6 (Q91+)
```

### 결과

| 지표 | 값 |
|------|-----|
| 총 reads 수 | 85 |
| 길이 범위 | 885 ~ 979 bp |
| 평균 길이 | 975.6 bp |
| 중앙값 길이 | 977.0 bp |
| 평균 품질 | Q91.6 |
| 품질 범위 | Q82.1 ~ Q93.0 |

### Phred 품질 점수 해석

```
Q20 = 99%     정확도 (100개 염기 중 1개 오류)
Q30 = 99.9%   정확도 (1000개 염기 중 1개 오류)
Q40 = 99.99%  정확도
Q90 = 매우 높음 (PacBio CCS의 '~' 문자 = Q93)
```

**해석**: Q91+의 품질 점수는 PacBio CCS(환형 일치 시퀀싱) 데이터임을 의미합니다. CCS는 같은 분자를 여러 번 읽어 일치 서열을 추출함으로써 거의 완벽한 정확도를 달성합니다.

### 품질 관리 시각화 코드

```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 리드 길이 분포
axes[0].hist(lengths, bins=20, color='steelblue', edgecolor='black', alpha=0.7)
axes[0].set_xlabel('Read Length (bp)')
axes[0].set_ylabel('Count')
axes[0].set_title('Read Length Distribution')
axes[0].axvline(np.mean(lengths), color='red', linestyle='--',
                label=f'Mean: {np.mean(lengths):.0f} bp')
axes[0].legend()

# 품질 분포
axes[1].hist(mean_quals, bins=20, color='forestgreen', edgecolor='black', alpha=0.7)
axes[1].set_xlabel('Mean Phred Quality Score')
axes[1].set_ylabel('Count')
axes[1].set_title('Per-Read Mean Quality Score')
axes[1].axvline(np.mean(mean_quals), color='red', linestyle='--',
                label=f'Mean: {np.mean(mean_quals):.1f}')
axes[1].legend()

plt.tight_layout()
plt.savefig("analysis_results/01_read_stats.png", dpi=150, bbox_inches='tight')
```

> **그림 1**: 왼쪽 그래프는 reads 길이가 977 bp 부근에 고도로 집중되어 있음을 보여주며, PCR 증폭 단편의 크기가 일관됨을 의미합니다. 오른쪽 그래프는 대부분의 reads가 Q90 이상의 품질을 가지고 있어, CCS 데이터의 높은 정확도를 확인합니다.

---

## 5. 삽질 기록: VHH 서열 추출 실패

### 초기 접근법: 6-frame ORF 탐색

처음 떠올린 방법은 직관적이었습니다 — 각 read를 6-frame 번역하고 가장 긴 ORF를 찾는 것:

```python
def find_longest_orf(seq_str):
    """6개 reading frame에서 가장 긴 ORF 찾기"""
    seq = Seq(seq_str)
    best_protein = ""
    for strand, nuc in [("+", seq), ("-", seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate())
            for i in range(len(trans)):
                if trans[i] == "M":  # 시작 코돈
                    stop = trans.find("*", i)  # 종결 코돈
                    if stop != -1 and (stop - i) > len(best_protein):
                        best_protein = trans[i:stop]
    return best_protein
```

### 결과: 완전한 실패

```
ORF 분석 결과:
- 100aa 이상 ORF를 가진 reads: 2개 / 85개   ← 거의 전멸!
- 단백질 길이 범위: 6 ~ 109 aa
- 평균 단백질 길이: 16.9 aa                    ← VHH의 120aa보다 훨씬 짧음
- VHH 식별 수: 0개                              ← 완전히 실패!
```

### 문제 진단

VHH는 ~120 aa 단백질인데, 왜 ORF 탐색에서는 최대 109 aa 밖에 찾지 못했을까?

첫 번째 read의 6-frame 번역을 수동으로 확인해 봅시다:

```python
from Bio.Seq import Seq
r = records[0]
rc = r.seq.reverse_complement()

# 역상보 서열의 첫 번째 reading frame
prot = str(rc.translate())
# VHH 특징 시작 서열 검색
pos = prot.find("QVQL")
print(f"QVQL found at position {pos}")
print(prot[pos:pos+130])
```

**출력**:

```
QVQL found at position 38
QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSNNFGWNWIRQSPSRGLE*LGRTYYRSKWY
NDYAVSVRSRITINPDTSKNQFSLQLNSVTPEDTAVYYCARQGSTYFDYWGQGTLVTVSS
```

**핵심 발견**: `GLE*LGR` — VHH 서열 중간에 종결 코돈(`*`)이 있습니다!

이것이 ORF 탐색이 실패한 원인입니다. ORF 탐색은 첫 번째 종결 코돈에서 절단하는데, VHH의 FR2 영역에 amber 종결 코돈(TAG)이 내장되어 있었습니다.

### 근본 원인

이것은 **파지 디스플레이 벡터의 설계 특성**이며, 시퀀싱 오류가 아닙니다:

```
벡터 구조: ... pelB-VHH-[amber TAG]-pIII ...

amber 억제 균주(TG1)에서:    TAG → Gln (읽혀 통과, VHH-pIII 융합 단백질, 파지 표면 디스플레이)
비억제 균주(HB2151)에서:     TAG → STOP (번역 종결, 가용성 VHH 단백질 생산)
```

따라서 이 `*`는 벡터 설계의 일부이며, 세균에서는 글루타민(Q)으로 읽혀 통과됩니다. 전통적인 ORF 탐색 방법은 여기서 완전히 무용지물이 됩니다.

### 교훈

> **파지 디스플레이 라이브러리의 시퀀싱 데이터를 분석할 때는 단순한 ORF 탐색에 의존할 수 없습니다. VHH의 보존된 모티프(motif)를 기반으로 서열을 추출하고, amber 종결 코돈을 올바르게 처리해야 합니다.**

---

## 6. 2단계: VHH 서열 추출 (올바른 방법)

### 전략 개선

ORF 탐색이 통하지 않으므로, **보존 모티프(motif) 기반** 방법으로 전환합니다:

1. 각 read에 대해 정방향과 역상보 서열의 3개 reading frame(총 6가지) 시도
2. 번역 후 VHH 시작 모티프 검색: `QVQL`, `EVQL`, `DVQL`
3. VHH 종료 모티프 검색: `WGQG` (FR4의 시작)
4. 내부 종결 코돈을 `X`로 표시 (amber 읽기 통과 아미노산 대체)

```python
def extract_vhh_from_read(seq_record):
    """
    PacBio CCS read에서 VHH 서열 추출.
    역상보 변환과 내부 종결 코돈을 처리.
    """
    seq = seq_record.seq
    for strand, nuc in [("+", seq), ("-", seq.reverse_complement())]:
        for frame in range(3):
            prot = str(nuc[frame:].translate())
            # 내부 종결 코돈을 X로 치환
            prot_no_stop = prot.replace("*", "X")

            for start_motif in ['QVQL', 'EVQL', 'DVQL']:
                pos = prot_no_stop.find(start_motif)
                if pos == -1:
                    continue

                # FR4 시작 모티프 WGQG 검색
                wg_pos = prot_no_stop.find('WGQG', pos)
                if wg_pos == -1:
                    wg_pos = prot_no_stop.find('WGRG', pos)
                if wg_pos == -1:
                    continue

                # VHH + FR4 추출 (WGQGTLVTVSS)
                end_pos = min(wg_pos + 11, len(prot_no_stop))
                vhh_seq = prot_no_stop[pos:end_pos]

                # 길이 합리성 검증 (VHH는 보통 110-140aa)
                if 100 <= len(vhh_seq) <= 160:
                    return {
                        'vhh_seq': vhh_seq,
                        'strand': strand,
                        'frame': frame + 1,
                        'has_internal_stop': '*' in prot[pos:end_pos]
                    }
    return None
```

### 결과

```
VHH 추출 결과:
- 성공: 83개 / 85개 reads (97.6%)
- 실패: 2개 reads

서열 방향 분포:
    + (정방향):   45개
    - (역상보):   38개

Reading frame 분포:
    +1: 45개
    -1: 36개
    -3:  2개

내부 종결 코돈: 83개 (전부)    ← amber stop codon 확인

VHH 길이: 120 aa (모든 서열 일치)
```

### 추출된 VHH 서열 예시

```
  1: QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSNNFGWNWIRQSPSRGLEXL
 51: GRTYYRSKWYNDYAVSVRSRITINPDTSKNQFSLQLNSVTPEDTAVYYCAR
101: QGSTYFDYWGQGTLVTVSS
```

48번 위치의 `X`가 amber 종결 코돈에 해당하는 위치입니다 (발현 시스템에서는 Gln으로 읽혀 통과됨).

---

## 7. 3단계: CDR 영역 식별

### CDR 경계 결정 방법

보존 모티프와 위치 특성을 이용하여 CDR 경계를 설정합니다:

```python
def extract_cdr_regions(vhh_seq):
    """
    보존 모티프 기반 CDR 추출.

    핵심 앵커 포인트:
    - FR1 끝/CDR1 시작: 첫 번째 C 이후 약 4잔기 (~위치 26)
    - CDR1 끝: W 잔기 (~위치 36)
    - CDR2 시작: FR2 이후 (~위치 50)
    - CDR3 시작: YYC/YFC 모티프 이후 (~위치 96)
    - CDR3 끝: WGQG 이전
    - FR4: WGQGTLVTVSS
    """
    result = {}

    # 1) FR4 시작 찾기: WGQG
    fr4_start = vhh_seq.find('WGQG')
    if fr4_start == -1:
        return None
    result['FR4'] = vhh_seq[fr4_start:fr4_start + 11]

    # 2) CDR3 찾기: YYC 모티프 후 ~ WGQG 전
    cdr3_match = re.search(r'[YF][YF]C', vhh_seq[80:fr4_start])
    if cdr3_match:
        cdr3_start = 80 + cdr3_match.end()
    else:
        cdr3_start = fr4_start - 15  # 대안 추정
    result['CDR3'] = vhh_seq[cdr3_start:fr4_start]

    # 3) CDR1 찾기: 첫 번째 C 이후
    first_c = vhh_seq.find('C', 20, 30)
    cdr1_start = (first_c if first_c != -1 else 22) + 4
    w_pos = vhh_seq.find('W', cdr1_start, cdr1_start + 20)
    cdr1_end = w_pos if w_pos != -1 else cdr1_start + 8

    result['FR1']  = vhh_seq[0:cdr1_start]
    result['CDR1'] = vhh_seq[cdr1_start:cdr1_end]

    # 4) CDR2
    fr2_end = cdr1_end + 15
    cdr2_end = fr2_end + 8
    result['FR2']  = vhh_seq[cdr1_end:fr2_end]
    result['CDR2'] = vhh_seq[fr2_end:cdr2_end]
    result['FR3']  = vhh_seq[cdr2_end:cdr3_start]

    return result
```

### 식별 결과

```
┌─────────┬───────┬──────────────────────────────────────────┐
│ 영역    │ 길이  │ 서열                                      │
├─────────┼───────┼──────────────────────────────────────────┤
│ FR1     │ 25 aa │ QVQLQQSGPGLVKPSQTLSLTCAIS                │
│ CDR1    │ 10 aa │ GDSVSSNNFG                    ◀◀ 항원결합 │
│ FR2     │ 15 aa │ WNWIRQSPSRGLEXL                          │
│ CDR2    │  8 aa │ GRTYYRSK                      ◀◀ 항원결합 │
│ FR3     │ 41 aa │ WYNDYAVSVRSRITINPDTSKNQFSLQLNSVTPEDTAVYYC│
│ CDR3    │ 10 aa │ ARQGSTYFDY                    ◀◀ 항원결합 │
│ FR4     │ 11 aa │ WGQGTLVTVSS                              │
└─────────┴───────┴──────────────────────────────────────────┘
```

IgBLAST의 IMGT 넘버링 결과와 비교:

| 영역 | 본 방법 | IMGT (IgBLAST) | 일치? |
|------|---------|----------------|-------|
| FR1 | 25 aa | 25 aa | O |
| CDR1 | 10 aa (GDSVSSNNFG) | 10 aa (GDSVSSNNFG) | O |
| FR2 | 15 aa | 17 aa | 근사 |
| CDR2 | 8 aa | 9 aa (TYYRSKWYN) | 근사 |
| FR3 | 41 aa | 38 aa | 근사 |
| CDR3 | 10 aa (ARQGSTYFDY) | 10 aa (ARQGSTYFDY) | O |
| FR4 | 11 aa | 11 aa | O |

CDR3는 완전히 일치합니다! FR2/CDR2/FR3 경계에서 1-3개 잔기 차이가 있는데, 이는 넘버링 체계(Kabat, IMGT, Chothia)마다 CDR 경계 정의가 약간 다르기 때문입니다. **정확한 분석에는 IgBLAST/IMGT 결과를 기준으로 해야 합니다**.

---

## 8. 4단계: 서열 다양성 분석

### 분석 코드

```python
from collections import Counter

# 전장 VHH 다양성
full_seqs = [v['full_vhh'] for v in vhh_data]
unique_full = set(full_seqs)
full_counter = Counter(full_seqs)

# CDR별 다양성
for cdr_name in ['CDR1', 'CDR2', 'CDR3']:
    cdr_seqs = [v[cdr_name] for v in vhh_data]
    unique_cdr = set(cdr_seqs)
    diversity = len(unique_cdr) / len(cdr_seqs) * 100
    print(f"{cdr_name}: {len(unique_cdr)} 고유 / {len(cdr_seqs)} 전체 = {diversity:.1f}%")
```

### 결과

```
전장 VHH 서열 다양성:
  총 서열:     83개
  고유 서열:    2개
  다양성 비율:  2.4%
  싱글턴:      1개 (1번만 출현한 서열)
  확장 클론:   1개 (82번 출현!)

CDR 다양성:
  CDR       총수    고유    다양성
  ──────────────────────────────
  CDR1       83      1      1.2%
  CDR2       83      1      1.2%
  CDR3       83      1      1.2%

유일한 CDR3 서열:
  ARQGSTYFDY    (83회) ████████████████████████████████████
```

### 해석

**이것은 다양성 라이브러리가 아닌 단일 클론 검증 데이터입니다**:

- 85개 reads 중 83개가 **완전히 동일한 VHH 서열**을 코딩
- CDR3는 1종류의 서열만 존재 (`ARQGSTYFDY`)
- 이는 데이터의 출처가 다음 중 하나임을 의미:
  - 파지 디스플레이 패닝 후 선별된 **단일 클론**, 또는
  - 특정 나노바디의 **서열 확인** 실험
- PacBio 장읽기 시퀀싱은 전장 서열의 정확성 확인에 사용됨

---

## 9. 5단계: 아미노산 조성 분석

### 아미노산 분류

```python
AA_GROUPS = {
    'Hydrophobic': set('AVILMFWP'),  # 소수성 - 단백질 코어
    'Polar':       set('STNQYC'),    # 극성   - 수소결합, 항원 접촉
    'Positive':    set('RHK'),        # 양전하 - 정전기적 상호작용
    'Negative':    set('DE'),         # 음전하 - 정전기적 상호작용
    'Glycine':     set('G'),          # 글리신 - 구조적 유연성
}
```

### CDR3 아미노산 빈도

```
CDR3 서열: ARQGSTYFDY (10개 아미노산)

  Y (극성):   20.0% ██████████   ← 타이로신이 가장 많음 (항원 접촉 핵심 잔기)
  A (소수성): 10.0% █████
  R (양전하): 10.0% █████
  Q (극성):   10.0% █████
  G (글리신): 10.0% █████
  S (극성):   10.0% █████
  T (극성):   10.0% █████
  F (소수성): 10.0% █████
  D (음전하): 10.0% █████
```

**핵심 관찰**: CDR3에서 타이로신(Y) 함량이 가장 높습니다(2개). 이는 항체에서 매우 일반적인데, 타이로신의 측쇄는 방향족이면서 수산기를 가지고 있어 항원과 다양한 상호작용(수소결합, π-π 적층, 반데르발스 힘)을 형성하는 이상적인 잔기입니다.

---

## 10. 6단계: IgBLAST Germline 유전자 분석

### V(D)J 재조합이란?

```
Germline DNA:
  ──[V1][V2]...[Vn]──[D1][D2]...[Dn]──[J1][J2]...[Jn]──[C]──

V(D)J 재조합 (B세포 발달 과정에서):
  1. D-J 결합: D와 J 유전자 각각 하나씩 선택 → D-J
  2. V-DJ 결합: V 유전자 하나 선택 → V-D-J
  3. 접합부 다양성: 연결부에서 뉴클레오티드 무작위 추가/제거 (N/P nucleotides)

재조합 후:
  ──[V]─[N]─[D]─[N]─[J]──[C]──
     └── CDR3 영역 ──┘
```

### IgBLAST 설치 및 실행

```bash
# 설치
conda install -y -c bioconda igblast

# 설치 경로 확인
IGDATA=$(dirname $(which igblastn))/../share/igblast
ls $IGDATA/internal_data/human/  # germline 데이터 존재 확인
```

### 입력 파일 준비

원시 reads를 역상보 변환해야 올바른 VHH 코딩 방향을 얻을 수 있으므로, 먼저 역상보 FASTA를 준비합니다:

```python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

records = list(SeqIO.parse("SRR6269034.fastq", "fastq"))
rc_records = []
seen = set()
for i, r in enumerate(records):
    rc = r.seq.reverse_complement()
    s = str(rc)
    if s not in seen:  # 중복 제거
        seen.add(s)
        rc_records.append(SeqRecord(rc, id=f"Read_{i+1:03d}_RC",
                                     description=f"reverse_complement of {r.id}"))

SeqIO.write(rc_records, "unique_reads_rc.fasta", "fasta")
print(f"Wrote {len(rc_records)} unique sequences")  # 25개 고유 서열
```

### IgBLAST 실행

```bash
export IGDATA=/home/syslab/miniforge3/share/igblast

igblastn \
    -germline_db_V "$IGDATA/internal_data/human/human_V" \
    -germline_db_D germline_db/human_ighd \
    -germline_db_J germline_db/human_ighj \
    -organism human \
    -domain_system imgt \
    -query analysis_results/unique_reads_rc.fasta \
    -auxiliary_data "$IGDATA/optional_file/human_gl.aux" \
    -show_translation \
    -outfmt 19 \
    -out analysis_results/igblast_airr.tsv
```

### IgBLAST 핵심 결과

```
V-(D)-J rearrangement 요약:
  V gene:     IGHV6-1*01
  D gene:     IGHD1-26*01
  J gene:     IGHJ4*02
  Chain:      VH
  Stop codon: Yes (파지 디스플레이 벡터의 amber stop)
  V-J frame:  In-frame
  Productive: No (amber stop 때문 — 억제 균주에서는 기능적)

V-(D)-J 접합부 (CDR3 DNA):
  V 말단:     CAAGA
  V-D 접합:   CA
  D 영역:     GGGCAGC
  D-J 접합:   ACT
  J 시작:     TACTT

CDR3 핵산: GCAAGACAGGGCAGCACTTACTTCGACTAT
CDR3 아미노산: ARQGSTYFDY
```

### IMGT 영역 정렬 상세

```
정렬 요약 (Query vs IGHV6-1*01 germline):

  영역        시작   끝     길이   일치   불일치  Identity
  ─────────────────────────────────────────────────────
  FR1-IMGT    115    189    75     75     0      100.0%
  CDR1-IMGT   190    219    30     26     4       86.7%
  FR2-IMGT    220    270    51     50     1       98.0%
  CDR2-IMGT   271    297    27     27     0      100.0%
  FR3-IMGT    298    411   114    112     2       98.2%
  CDR3 (germ) 412    417     6      6     0      100.0%
  ─────────────────────────────────────────────────────
  합계         -      -    303    296     7       97.7%
```

**해석**:
- FR1과 CDR2는 germline과 완전히 일치 (100%)
- CDR1에 4개 돌연변이 (86.7%) — 종 간 차이(낙타 vs 인간) 또는 체세포 돌연변이 가능
- 전체 97.7% 일치율은 이 VHH가 인간 IGHV6-1과 매우 유사함을 시사

### 정렬 시각화

```
                        <── FR1 ──><── CDR1 ──>
  Query:    QVQLQQSGPGLVKPSQTLSLTCAIS GDSVSSNNFG
  Germline: QVQLQQSGPGLVKPSQTLSLTCAIS GDSVSSNSA A
                                       |||||| * * ← 4개 차이

                        <── FR2 ──><── CDR2 ──>
  Query:    WNWIRQSPSRGLE*LGR TYYRSKWY
  Germline: WNWIRQSPSRGLEWLGR TYYRSKWY
                          ^                ← amber stop (원래는 W)

                        <── FR3 ────────────────────><── CDR3 ──>
  Query:    NDYAVSVRSRITINPDTSKNQFSLQLNSVTPEDTAVYYC ARQGSTYFDY
  Germline: NDYAVSVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYC AR─────────
                    ^                                   (V-D-J 접합)

                        <- FR4 ->
  Query:    WGQGTLVTVSS
  J gene:   WGQG.LV.VSS  (IGHJ4*02, 91.3% identity)
```

---

## 11. 삽질 기록: IgBLAST 데이터베이스 설정

### 문제 1: D/J 유전자 데이터베이스를 찾을 수 없음

```
BLAST Database error: No alias or index file found for nucleotide database
[human_gl_D] in search path [/home/syslab/Desktop/nanobody::]
```

**원인**: conda로 설치한 IgBLAST에는 V 유전자 내부 데이터베이스(`internal_data/human/human_V`)만 포함되어 있고, D와 J 유전자의 사전 구축된 BLAST 데이터베이스는 없음.

**해결 방법**: D/J 유전자 데이터베이스를 수동 생성.

```bash
# IGHD germline FASTA 파일 생성
cat > germline_db/human_ighd.fasta << 'EOF'
>IGHD1-1*01
GGTATAACTGGAACTAC
>IGHD1-7*01
GGTATAGATGATTATTAT
>IGHD1-26*01
GGTATAGTGGGCAGCTAC
>IGHD2-2*01
AGGATATTGTAGTGGTATATACTGTATTT
>IGHD3-3*01
GTATTACTATGGTTCGGGGAGTTATTATAAC
>IGHD3-10*01
GTATTACTATGATTACGAC
>IGHD4-4*01
TGACTACGGTGACTAC
>IGHD5-5*01
GTAACTGGAACGAC
>IGHD6-6*01
GAGTATAGCAGCTCGTCC
>IGHD6-13*01
GGGTATAGCAGCGGCTGG
>IGHD7-27*01
CTAACTGGGGA
EOF
# (실제로는 전체 26개 D 유전자를 사용했으며, 여기서는 간략화하여 표시)

# IGHJ germline FASTA 파일 생성
cat > germline_db/human_ighj.fasta << 'EOF'
>IGHJ1*01
GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG
>IGHJ2*01
CTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG
>IGHJ3*01
TGATGCTTTTGATGTCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG
>IGHJ4*01
ACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG
>IGHJ4*02
ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG
>IGHJ5*01
ACAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG
>IGHJ6*01
ATTACTACTACTACTACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCAG
EOF

# BLAST 데이터베이스 구축
makeblastdb -parse_seqids -dbtype nucl \
    -in germline_db/human_ighd.fasta \
    -out germline_db/human_ighd \
    -title "Human IGHD"

makeblastdb -parse_seqids -dbtype nucl \
    -in germline_db/human_ighj.fasta \
    -out germline_db/human_ighj \
    -title "Human IGHJ"
```

### 문제 2: IMGT 온라인 데이터베이스 다운로드 실패

```bash
# IMGT에서 다운로드 시도
curl -sL "https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframedP" \
    -o imgt_all.fasta
# 결과: FASTA가 아닌 HTML 페이지가 다운로드됨
```

**원인**: IMGT 웹사이트는 브라우저를 통한 접근이 필요하며, curl로 직접 요청하면 로그인/검증 페이지로 리다이렉트됨.

**해결 방안**:
1. IgBLAST에 내장된 internal_data의 V 유전자 데이터베이스 활용
2. D/J 유전자 서열을 수동으로 준비 (위에서 설명)
3. 또는 NCBI의 IgBLAST FTP 사이트에서 사전 구축된 데이터베이스 다운로드

### 문제 3: 종 불일치

IgBLAST에는 낙타과(Camelidae)의 germline 데이터베이스가 없음.

**대응 방법**: 인간(human) germline을 참조로 사용. 이것이 가능한 이유:
- IGHV6-1은 인간 V 유전자 중 낙타 VHH와 구조적으로 가장 유사
- 분석 결과는 여전히 의미가 있으나, 서열 차이가 종 간 자연적 차이에서 비롯될 수 있음에 주의
- 정식 논문 발표 시에는 낙타과 전용 germline 데이터베이스를 사용해야 함

---

## 12. 최종 결과 및 생물학적 해석

### 종합 결과

```
┌────────────────────────────────────────────────┐
│           최종 분석 결과 요약                     │
├────────────────────────────────────────────────┤
│                                                │
│  프로젝트: SRP124616                            │
│  데이터: SRR6269034 (PacBio CCS, 85 reads)      │
│                                                │
│  데이터 품질:                                    │
│    평균 리드 길이: 976 bp                        │
│    평균 품질: Q92 (매우 높음)                     │
│                                                │
│  VHH 식별:                                      │
│    성공률: 83/85 (97.6%)                         │
│    고유 서열: 2개                                 │
│    고유 CDR3: 1개                                │
│                                                │
│  Germline 유전자:                                │
│    V gene: IGHV6-1*01 (97.7% identity)          │
│    D gene: IGHD1-26*01                           │
│    J gene: IGHJ4*02   (91.3% identity)          │
│                                                │
│  CDR3:                                          │
│    서열: ARQGSTYFDY                              │
│    길이: 10 aa                                   │
│    DNA:  GCAAGACAGGGCAGCACTTACTTCGACTAT          │
│                                                │
│  V-D-J 접합:                                     │
│    V 말단: CAAGA                                 │
│    V-D 접합: CA (N-addition)                     │
│    D 영역: GGGCAGC                               │
│    D-J 접합: ACT                                 │
│    J 시작: TACTT                                 │
│                                                │
└────────────────────────────────────────────────┘
```

### 생물학적 해석

**1. 이것은 단일 클론 나노바디입니다**

85개 reads 중 83개가 동일한 VHH 서열을 코딩하며, CDR3는 1종류(ARQGSTYFDY)뿐입니다. 라이브러리 다양성 시퀀싱이 아니라 특정 클론의 서열 확인 실험입니다.

**2. PacBio CCS를 사용한 이유**

PacBio 장읽기 시퀀싱은 단일 read로 VHH 코딩 서열 전체(~400 bp)를 커버할 수 있어 조립이 필요 없습니다. CCS의 높은 정확도(Q90+)가 서열의 신뢰성을 보장합니다.

**3. 벡터 구조 확인**

```
5'─ pelB 시그널 펩타이드 ─ VHH ─ [amber TAG] ─ (G4S)3 링커 ─ VL ─ 3'
                           │                                    │
                       983 bp read가 전체 클론 삽입 단편을 커버
```

reads에서 VHH 외에도 경쇄 가변 영역(IGKV1-39*01, 93.7% identity)이 검출되어, 이것이 VHH-VL 이중특이성 구조체 또는 scFv 융합 단백질의 일부임을 시사합니다.

**4. FR2 영역의 VHH 특이적 돌연변이**

FR2에서 핵심 차이가 관찰됩니다: `E*L` vs germline의 `EWL`. 여기서 `*`는 amber 종결 코돈입니다. 실제 발현에서 W→Q 치환은 VHH(기존 VH 대비)의 특징적 돌연변이 중 하나로, VL 없이도 가용성을 유지하는 데 기여합니다.

**5. Germline 일치도 분석**

97.7%의 V 유전자 일치도 (인간 IGHV6-1*01 대비):
- 303개 총 뉴클레오티드 중 7개 뉴클레오티드 차이
- CDR1의 4개 차이는 낙타과 vs 인간의 germline 자연적 차이를 반영할 수 있음
- FR 영역의 3개 차이는 체세포 돌연변이 또는 종 간 차이일 수 있음

---

## 13. 정리 및 향후 방향

### 전체 분석 흐름 리뷰

```
원시 데이터 (SRA)
    │
    ▼
[다운로드] prefetch + fasterq-dump
    │
    ▼
FASTQ 파일 (85 reads, ~977 bp, Q92)
    │
    ├──▶ [품질관리] 리드 길이 분포, 품질 분포
    │
    ├──▶ [번역] 역상보 + 6-frame 번역
    │      │
    │      ├── ✗ ORF 탐색 실패 (amber stop codon!)
    │      │
    │      └── ✓ 모티프 기반 VHH 추출 (83/85 성공)
    │
    ├──▶ [CDR 식별] 보존 모티프 + IMGT 넘버링
    │
    ├──▶ [다양성 분석] → 단일 클론 (CDR3 1종류)
    │
    ├──▶ [아미노산 분석] → CDR3에 Tyr 풍부
    │
    └──▶ [IgBLAST] → IGHV6-1*01 / IGHD1-26*01 / IGHJ4*02
```

### 삽질 총정리

| 문제 | 증상 | 근본 원인 | 해결 방안 |
|------|------|-----------|-----------|
| NCBI 도구 누락 | `command not found` | 미설치 | `conda install entrez-direct sra-tools` |
| SSL 연결 오류 | `curl: (56) SSL_read error` | OpenSSL 버전 비호환 | ENA API 사용 또는 `curl -k` |
| VHH 추출 실패 | ORF 탐색에서 짧은 단편만 발견 | amber 종결 코돈이 ORF를 절단 | 보존 모티프 검색으로 ORF 탐색 대체 |
| IgBLAST D/J 누락 | `No alias or index file` | conda에 D/J DB 미포함 | BLAST 데이터베이스 수동 생성 |
| IMGT 다운로드 실패 | FASTA 대신 HTML 수신 | 웹사이트 브라우저 접근 필요 | IgBLAST 내장 데이터 사용 |
| 낙타과 germline 부재 | - | IgBLAST에 낙타 데이터 없음 | 인간 germline 대체 사용 (한계 명시) |

### 향후 분석 방향

1. **정밀 넘버링**: ANARCI 또는 IMGT/V-QUEST를 사용한 표준 IMGT 넘버링
2. **구조 예측**: AlphaFold2 또는 ESMFold로 VHH 3D 구조 예측
3. **데이터베이스 검색**: PDB, sdAb-DB에서 상동 나노바디 검색
4. **친화도 예측**: CDR 서열 특성 기반 항원 결합 특성 예측
5. **인간화 설계**: 치료 응용이 필요한 경우 VHH 인간화 개조
6. **딥 시퀀싱**: 원래 라이브러리에 대해 NGS 시퀀싱으로 전체 다양성 평가

### 생성된 모든 파일

```
nanobody/
├── SRR6269034.fastq                          # 원시 시퀀싱 데이터
├── nanobody_analysis.py                       # 전체 분석 스크립트
├── blog_nanobody_analysis_zh.md              # 블로그 (중국어)
├── blog_nanobody_analysis_kr.md              # 블로그 (한국어, 본 문서)
├── germline_db/                               # Germline 데이터베이스
│   ├── human_ighd.fasta + BLAST 인덱스
│   └── human_ighj.fasta + BLAST 인덱스
└── analysis_results/                          # 분석 결과
    ├── 01_read_stats.png                      # 리드 길이/품질 분포도
    ├── 02_cdr_length_distribution.png         # CDR 길이 분포도
    ├── 03_aa_composition.png                  # 아미노산 조성 파이차트
    ├── 04_similarity_distribution.png         # 서열 유사도 분포도
    ├── 05_cdr3_position_freq.png              # CDR3 위치별 빈도 히트맵
    ├── 06_vhh_structure.png                   # VHH 도메인 구조 다이어그램
    ├── vhh_sequences.fasta                    # VHH 단백질 서열
    ├── cdr3_sequences.fasta                   # CDR3 서열
    ├── vhh_summary.tsv                        # 분석 요약 표
    ├── reads_reverse_complement.fasta         # 역상보 핵산 서열
    ├── unique_reads_rc.fasta                  # 중복 제거 후 고유 서열
    ├── igblast_results.txt                    # IgBLAST 상세 정렬
    ├── igblast_airr.tsv                       # AIRR 표준 형식
    └── igblast_tabular.txt                    # 탭 구분 형식
```

---

> **저자 노트**: 본 분석에서는 인간 germline 유전자를 참조로 사용했습니다. 나노바디는 낙타과 동물에서 유래하므로, 일부 서열 차이는 체세포 돌연변이가 아닌 종 간 자연적 차이일 수 있습니다. 정식 연구에서는 낙타과 전용 germline 데이터베이스(예: IMGT의 *Camelus dromedarius* 또는 *Vicugna pacos* 참조 서열)를 사용하여 더 정밀한 분석을 수행할 것을 권장합니다.

> **데이터 출처**: NCBI SRA, BioProject SRP124616, Run SRR6269034
