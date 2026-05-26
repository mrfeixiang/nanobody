#!/usr/bin/env python3
"""
=============================================================================
나노바디(Nanobody/VHH) 시퀀싱 데이터 분석 - 초보자용 케이스 스터디
=============================================================================
프로젝트: SRP124616 (SRR6269034)
데이터: PacBio CCS (Circular Consensus Sequencing) 리드
목적: 나노바디 라이브러리의 서열 분석 및 CDR(상보성 결정 영역) 다양성 파악

배경 지식:
-----------
나노바디(Nanobody)란?
  - 낙타과 동물(낙타, 라마, 알파카)에서 발견되는 중쇄 항체(HCAb)의
    가변 영역(VHH)만으로 이루어진 소형 항체 단편입니다.
  - 일반 항체(~150kDa)에 비해 매우 작고(~15kDa),
    안정적이며, 대장균에서 생산 가능합니다.

VHH 구조:
  FR1 - CDR1 - FR2 - CDR2 - FR3 - CDR3 - FR4
  (FR = Framework Region, CDR = Complementarity Determining Region)

  CDR3가 항원 결합 특이성을 결정하는 가장 중요한 영역입니다.

이 데이터의 특징:
  - 파지 디스플레이(Phage Display) 벡터에 클로닝된 VHH 서열
  - 벡터 내에 amber stop codon(TAG → *)이 포함됨
    (파지 디스플레이에서 흔히 사용되는 구조)
  - 서열이 역상보(reverse complement) 방향으로 읽힘
=============================================================================
"""

import os
import sys
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import re

# ===== 설정 =====
WORK_DIR = "/home/syslab/Desktop/nanobody"
FASTQ_FILE = os.path.join(WORK_DIR, "SRR6269034.fastq")
OUTPUT_DIR = os.path.join(WORK_DIR, "analysis_results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.rcParams['font.size'] = 10

print("=" * 70)
print("  나노바디(Nanobody) 시퀀싱 데이터 분석 케이스 스터디")
print("  프로젝트: SRP124616 | 데이터: SRR6269034")
print("=" * 70)

# =============================================================================
# 1단계: FASTQ 파일 읽기 및 기본 통계
# =============================================================================
print("\n" + "=" * 70)
print("  1단계: FASTQ 데이터 읽기 및 품질 확인")
print("=" * 70)
print("""
[설명] FASTQ 파일이란?
  시퀀싱 결과를 저장하는 표준 형식으로, 각 리드(read)는 4줄로 구성됩니다:
    1줄: @리드ID (헤더)
    2줄: 염기서열 (ATCG)
    3줄: + (구분자)
    4줄: 품질점수 (각 염기의 정확도를 나타내는 ASCII 문자)

  이 데이터는 PacBio CCS(Circular Consensus Sequencing)로 생산되었습니다.
  CCS는 같은 분자를 여러 번 읽어 정확도를 높인 리드입니다.
""")

records = list(SeqIO.parse(FASTQ_FILE, "fastq"))
print(f"  총 리드 수: {len(records)}개")

lengths = [len(r.seq) for r in records]
print(f"  리드 길이 범위: {min(lengths)} ~ {max(lengths)} bp")
print(f"  평균 리드 길이: {np.mean(lengths):.1f} bp")
print(f"  중앙값 리드 길이: {np.median(lengths):.1f} bp")

all_quals = []
mean_quals = []
for r in records:
    quals = r.letter_annotations["phred_quality"]
    all_quals.extend(quals)
    mean_quals.append(np.mean(quals))

print(f"\n  [품질 점수 (Phred Score)]")
print(f"  평균 품질 점수: {np.mean(all_quals):.1f}")
print(f"  리드별 평균 품질: {np.mean(mean_quals):.1f} (범위: {min(mean_quals):.1f} ~ {max(mean_quals):.1f})")
print("""
  Phred 품질 점수 해석:
    Q20 = 99% 정확도    (100개 중 1개 오류)
    Q30 = 99.9% 정확도  (1000개 중 1개 오류)
    Q40 = 99.99% 정확도 (10000개 중 1개 오류)
    Q90+ = PacBio CCS 최고 품질 ('~' = Q93)
""")

# 그림 1: 리드 길이 분포 및 품질 분포
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].hist(lengths, bins=20, color='steelblue', edgecolor='black', alpha=0.7)
axes[0].set_xlabel('Read Length (bp)')
axes[0].set_ylabel('Count')
axes[0].set_title('Read Length Distribution')
axes[0].axvline(np.mean(lengths), color='red', linestyle='--', label=f'Mean: {np.mean(lengths):.0f} bp')
axes[0].legend()

axes[1].hist(mean_quals, bins=20, color='forestgreen', edgecolor='black', alpha=0.7)
axes[1].set_xlabel('Mean Phred Quality Score')
axes[1].set_ylabel('Count')
axes[1].set_title('Per-Read Mean Quality Score')
axes[1].axvline(np.mean(mean_quals), color='red', linestyle='--', label=f'Mean: {np.mean(mean_quals):.1f}')
axes[1].legend()

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "01_read_stats.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  >> 그림 저장: analysis_results/01_read_stats.png")


# =============================================================================
# 2단계: VHH 서열 추출 (역상보 + 내부 stop codon 처리)
# =============================================================================
print("\n" + "=" * 70)
print("  2단계: VHH 서열 추출")
print("=" * 70)
print("""
[설명] 이 데이터의 서열 구조:
  이 시퀀싱 데이터는 파지 디스플레이(Phage Display) 벡터에 클로닝된
  나노바디 서열입니다.

  파지 디스플레이란?
    박테리오파지(세균을 감염시키는 바이러스) 표면에 나노바디를 발현시켜
    항원에 결합하는 나노바디를 선별하는 기술입니다.

  데이터 특징:
  1. 서열이 역상보(reverse complement) 방향으로 읽힘
     → 역상보 변환 후 번역해야 합니다
  2. 벡터 구조 내에 amber stop codon (TAG → *)이 포함
     → 대장균 suppressor strain에서는 읽히지만,
       단순 번역시 stop으로 인식됩니다
  3. VHH는 'QVQL' 또는 'EVQL'로 시작하고 'WGQG'로 끝남
""")

def extract_vhh_from_read(seq_record):
    """
    PacBio CCS 리드에서 VHH 서열을 추출합니다.
    역상보 변환 → 6-frame 번역 → VHH 모티프 검색
    """
    seq = seq_record.seq
    # 정방향과 역상보 모두 시도
    for strand, nuc in [("+", seq), ("-", seq.reverse_complement())]:
        for frame in range(3):
            prot = str(nuc[frame:].translate())
            # internal stop codon을 무시하고 VHH 시작 모티프 검색
            prot_no_stop = prot.replace("*", "X")  # stop을 X로 치환

            for start_motif in ['QVQL', 'EVQL', 'DVQL']:
                pos = prot_no_stop.find(start_motif)
                if pos == -1:
                    continue

                # WGQG (FR4 시작 모티프) 검색
                wg_pos = prot_no_stop.find('WGQG', pos)
                if wg_pos == -1:
                    wg_pos = prot_no_stop.find('WGRG', pos)
                if wg_pos == -1:
                    continue

                # VHH + 일부 FR4 추출 (WGQGTLVTVSS 등)
                end_pos = min(wg_pos + 11, len(prot_no_stop))
                vhh_seq = prot_no_stop[pos:end_pos]

                # 합리적 길이 확인 (VHH는 보통 110-140aa)
                if 100 <= len(vhh_seq) <= 160:
                    return {
                        'vhh_seq': vhh_seq,
                        'strand': strand,
                        'frame': frame + 1,
                        'start_pos': pos,
                        'has_internal_stop': '*' in prot[pos:end_pos],
                        'read_id': seq_record.id
                    }
    return None

# 모든 리드에서 VHH 추출
vhh_results = []
failed_reads = []
for r in records:
    result = extract_vhh_from_read(r)
    if result:
        vhh_results.append(result)
    else:
        failed_reads.append(r.id)

print(f"  VHH 추출 결과:")
print(f"  - 성공: {len(vhh_results)}개 / {len(records)}개 리드")
print(f"  - 실패: {len(failed_reads)}개 리드")

if vhh_results:
    # 방향 통계
    strand_counts = Counter(v['strand'] for v in vhh_results)
    frame_counts = Counter(f"{v['strand']}{v['frame']}" for v in vhh_results)
    stop_count = sum(1 for v in vhh_results if v['has_internal_stop'])

    print(f"\n  서열 방향 분포:")
    for s, c in strand_counts.items():
        direction = "정방향" if s == "+" else "역상보"
        print(f"    {s} ({direction}): {c}개")
    print(f"\n  Reading frame 분포:")
    for f, c in sorted(frame_counts.items()):
        print(f"    {f}: {c}개")
    print(f"\n  내부 stop codon 포함 서열: {stop_count}개")
    print(f"  (파지 디스플레이 벡터의 amber stop codon)")

    vhh_lengths = [len(v['vhh_seq']) for v in vhh_results]
    print(f"\n  VHH 길이: {min(vhh_lengths)} ~ {max(vhh_lengths)} aa (평균: {np.mean(vhh_lengths):.1f})")

    # 예시 VHH 출력
    print(f"\n  [예시 VHH 서열 (첫 번째 리드)]")
    ex = vhh_results[0]['vhh_seq']
    # 10aa씩 끊어서 보기 좋게 출력
    for i in range(0, len(ex), 60):
        print(f"    {i+1:3d}: {ex[i:i+60]}")


# =============================================================================
# 3단계: CDR 영역 추출
# =============================================================================
print("\n" + "=" * 70)
print("  3단계: CDR(상보성 결정 영역) 추출")
print("=" * 70)
print("""
[설명] VHH의 CDR 영역:
  VHH는 4개의 프레임워크(FR)와 3개의 CDR로 구성됩니다:

    FR1 ─ CDR1 ─ FR2 ─ CDR2 ─ FR3 ─ CDR3 ─ FR4

  CDR (Complementarity Determining Region):
    항원과 직접 접촉하는 루프 영역입니다.
    CDR1, CDR2: 항원 결합에 보조적 역할
    CDR3: 항원 결합 특이성의 핵심 영역 (가장 가변적)

  FR (Framework Region):
    VHH의 구조적 골격을 형성하는 보존된 영역입니다.

  CDR 경계 결정 방법:
    - Kabat 넘버링: 서열 가변성 기반
    - IMGT 넘버링: 구조 기반 (현재 표준)
    - Chothia 넘버링: 구조적 루프 기반

  여기서는 보존된 모티프 패턴을 이용한 간단한 방법을 사용합니다.
  (정확한 분석에는 ANARCI, IMGT/V-QUEST 등 전문 도구를 사용합니다)
""")

def extract_cdr_regions(vhh_seq):
    """
    VHH 서열에서 CDR 영역을 보존된 모티프 기반으로 추출합니다.

    주요 보존 모티프:
    - CDR1 시작: C 다음 (위치 ~23 근처) + 5잔기 뒤
    - CDR1 끝: W 앞 (위치 ~34 근처)
    - CDR2 시작: FR2의 보존된 패턴 후
    - CDR3 시작: 보존된 C (위치 ~92) 뒤, YYC 패턴
    - CDR3 끝: WGQG 앞
    """
    result = {
        'full_vhh': vhh_seq,
        'FR1': '', 'CDR1': '', 'FR2': '', 'CDR2': '',
        'FR3': '', 'CDR3': '', 'FR4': ''
    }

    # 1) FR4 끝 찾기: WGQG 또는 WGRG
    fr4_start = -1
    for motif in ['WGQG', 'WGRG', 'WGKG']:
        pos = vhh_seq.find(motif)
        if pos != -1:
            fr4_start = pos
            break

    if fr4_start == -1:
        return None

    result['FR4'] = vhh_seq[fr4_start:fr4_start + 11]

    # 2) CDR3 찾기: YYC 또는 YFC 패턴 (FR3 끝의 보존 모티프)
    cdr3_anchor = -1
    for pattern in ['YYC', 'YFC', 'YHC', 'YLC']:
        # FR3 끝 부근에서 검색 (위치 85-105)
        search_region = vhh_seq[80:fr4_start]
        pos = search_region.rfind(pattern)
        if pos != -1:
            cdr3_anchor = 80 + pos + len(pattern)
            break

    if cdr3_anchor == -1:
        # 대략적 위치 사용
        cdr3_anchor = fr4_start - 15 if fr4_start > 100 else 96

    result['CDR3'] = vhh_seq[cdr3_anchor:fr4_start]

    # 3) CDR1 찾기: 첫 번째 C 뒤 + 보존 패턴
    # 일반적으로 VHH 위치 26-35 근처
    first_c = vhh_seq.find('C', 20, 30)
    if first_c == -1:
        first_c = 22

    # CDR1은 C 이후 약 4잔기 뒤에서 시작
    cdr1_start = first_c + 4
    # W가 CDR1 끝을 표시 (위치 ~36)
    w_pos = vhh_seq.find('W', cdr1_start, cdr1_start + 20)
    if w_pos == -1:
        w_pos = cdr1_start + 8
    cdr1_end = w_pos

    result['FR1'] = vhh_seq[0:cdr1_start]
    result['CDR1'] = vhh_seq[cdr1_start:cdr1_end]

    # 4) CDR2 찾기
    # FR2는 CDR1 끝부터, CDR2는 약 위치 50-58
    fr2_end = cdr1_end + 15  # FR2는 약 15잔기
    # CDR2 끝은 약 위치 58-65
    # 보존된 R/K-F/L 패턴으로 FR3 시작 찾기
    cdr2_end = fr2_end + 8  # CDR2는 약 8잔기

    result['FR2'] = vhh_seq[cdr1_end:fr2_end]
    result['CDR2'] = vhh_seq[fr2_end:cdr2_end]
    result['FR3'] = vhh_seq[cdr2_end:cdr3_anchor]

    return result

# 모든 VHH에서 CDR 추출
vhh_data = []
for v in vhh_results:
    regions = extract_cdr_regions(v['vhh_seq'])
    if regions and regions['CDR3']:
        regions['read_id'] = v['read_id']
        vhh_data.append(regions)

print(f"  CDR 추출 결과:")
print(f"  - CDR 추출 성공: {len(vhh_data)}개 / {len(vhh_results)}개 VHH")

if vhh_data:
    # 예시 출력
    ex = vhh_data[0]
    print(f"\n  [예시 VHH 영역 구조 (첫 번째 서열)]")
    print(f"  ┌─────────┬───────┬─────────────────────────────────────────┐")
    print(f"  │ 영역    │ 길이  │ 서열                                    │")
    print(f"  ├─────────┼───────┼─────────────────────────────────────────┤")
    for region in ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']:
        seq = ex[region]
        marker = " ◀◀" if 'CDR' in region else ""
        print(f"  │ {region:7s} │ {len(seq):3d}aa │ {seq:40s}│{marker}")
    print(f"  └─────────┴───────┴─────────────────────────────────────────┘")
    print(f"  (◀◀ = 항원 결합 영역)")

    # CDR 길이 통계
    print(f"\n  [CDR 길이 통계]")
    for cdr in ['CDR1', 'CDR2', 'CDR3']:
        cdr_lens = [len(v[cdr]) for v in vhh_data if v[cdr]]
        if cdr_lens:
            print(f"  {cdr}: 평균 {np.mean(cdr_lens):.1f}aa (범위: {min(cdr_lens)} ~ {max(cdr_lens)})")


# =============================================================================
# 4단계: 서열 다양성 분석
# =============================================================================
print("\n" + "=" * 70)
print("  4단계: 서열 다양성 분석")
print("=" * 70)
print("""
[설명] 나노바디 라이브러리의 다양성:
  좋은 나노바디 라이브러리는 높은 서열 다양성을 가져야 합니다.
  특히 CDR3 영역의 다양성이 항원 결합 특이성을 결정합니다.

  다양성 지표:
  - 고유 서열(unique sequence) 비율: 높을수록 다양한 라이브러리
  - 싱글턴(singleton) 비율: 1번만 나타나는 서열의 비율
  - 클론 확장(clonal expansion): 같은 서열이 여러 번 나타남
    → 선별(selection) 과정에서 농축된 클론을 의미
""")

if vhh_data:
    # 전체 VHH 다양성
    full_seqs = [v['full_vhh'] for v in vhh_data]
    unique_full = set(full_seqs)
    full_counter = Counter(full_seqs)

    print(f"  [전체 VHH 서열 다양성]")
    print(f"  총 VHH 서열: {len(full_seqs)}개")
    print(f"  고유 서열: {len(unique_full)}개")
    print(f"  다양성 비율: {len(unique_full)/len(full_seqs)*100:.1f}%")

    singleton_full = sum(1 for c in full_counter.values() if c == 1)
    expanded_full = sum(1 for c in full_counter.values() if c > 1)
    print(f"  싱글턴 (1회): {singleton_full}개")
    print(f"  확장 클론 (2회+): {expanded_full}개")

    if expanded_full > 0:
        print(f"\n  [가장 많이 나타난 클론 (Top 5)]")
        for seq, count in full_counter.most_common(5):
            if count > 1:
                print(f"    {count}회: ...{seq[90:]}...")

    # CDR별 다양성
    print(f"\n  [CDR별 다양성]")
    print(f"  {'CDR':6s} {'총서열':>6s} {'고유':>6s} {'다양성':>8s}")
    print(f"  {'─'*30}")
    for cdr_name in ['CDR1', 'CDR2', 'CDR3']:
        cdr_seqs = [v[cdr_name] for v in vhh_data if v[cdr_name]]
        unique_cdr = set(cdr_seqs)
        diversity = len(unique_cdr) / len(cdr_seqs) * 100 if cdr_seqs else 0
        print(f"  {cdr_name:6s} {len(cdr_seqs):6d} {len(unique_cdr):6d} {diversity:7.1f}%")

    # CDR3 상세 분석
    cdr3_seqs = [v['CDR3'] for v in vhh_data if v['CDR3']]
    cdr3_counter = Counter(cdr3_seqs)

    print(f"\n  [CDR3 서열 상세 (모든 고유 서열)]")
    for seq, count in cdr3_counter.most_common():
        bar = '█' * count
        print(f"    {seq:25s} ({count}회) {bar}")

    # 그림 2: CDR 길이 분포
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    colors = ['#e74c3c', '#2ecc71', '#3498db']

    for idx, cdr_name in enumerate(['CDR1', 'CDR2', 'CDR3']):
        cdr_lens = [len(v[cdr_name]) for v in vhh_data if v[cdr_name]]
        if cdr_lens:
            axes[idx].hist(cdr_lens,
                          bins=range(min(cdr_lens), max(cdr_lens) + 2),
                          color=colors[idx], edgecolor='black', alpha=0.7,
                          align='left')
            axes[idx].set_xlabel('Length (aa)')
            axes[idx].set_ylabel('Count')
            axes[idx].set_title(f'{cdr_name} Length Distribution')
            axes[idx].axvline(np.mean(cdr_lens), color='black', linestyle='--',
                            label=f'Mean: {np.mean(cdr_lens):.1f}')
            axes[idx].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "02_cdr_length_distribution.png"), dpi=150, bbox_inches='tight')
    plt.close()
    print("\n  >> 그림 저장: analysis_results/02_cdr_length_distribution.png")


# =============================================================================
# 5단계: 아미노산 조성 분석
# =============================================================================
print("\n" + "=" * 70)
print("  5단계: 아미노산 조성 분석")
print("=" * 70)
print("""
[설명] 아미노산 특성:
  아미노산은 물리화학적 특성에 따라 그룹으로 나뉩니다:

  소수성(Hydrophobic): A, V, I, L, M, F, W, P
    → 단백질 내부에 묻히는 경향, 구조적 안정성에 기여
  극성(Polar): S, T, N, Q, Y, C
    → 수소결합 형성, 항원과의 상호작용에 중요
  양전하(Positive): R, H, K
    → 음전하 항원과의 정전기적 상호작용
  음전하(Negative): D, E
    → 양전하 항원과의 정전기적 상호작용
  글리신(Glycine): G
    → 가장 작은 아미노산, 구조적 유연성 제공
""")

AA_GROUPS = {
    'Hydrophobic': set('AVILMFWP'),
    'Polar': set('STNQYC'),
    'Positive': set('RHK'),
    'Negative': set('DE'),
    'Glycine': set('G'),
}

AA_GROUP_COLORS = {
    'Hydrophobic': '#f39c12',
    'Polar': '#2ecc71',
    'Positive': '#3498db',
    'Negative': '#e74c3c',
    'Glycine': '#95a5a6',
}

if vhh_data:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    regions_to_plot = ['CDR1', 'CDR2', 'CDR3', 'full_vhh']
    region_labels = ['CDR1', 'CDR2', 'CDR3', 'Full VHH']

    for idx, (region, label) in enumerate(zip(regions_to_plot, region_labels)):
        ax = axes[idx // 2][idx % 2]
        all_aa = ''.join([v[region] for v in vhh_data if v[region]])
        # X(internal stop)를 제거
        all_aa = all_aa.replace('X', '')
        aa_counts = Counter(all_aa)
        total = sum(aa_counts.values())

        if total > 0:
            group_fracs = {}
            for group, aas in AA_GROUPS.items():
                group_count = sum(aa_counts.get(aa, 0) for aa in aas)
                group_fracs[group] = group_count / total * 100

            sizes = list(group_fracs.values())
            labels_pie = list(group_fracs.keys())
            colors = [AA_GROUP_COLORS[g] for g in labels_pie]

            wedges, texts, autotexts = ax.pie(sizes, labels=labels_pie, colors=colors,
                                               autopct='%1.1f%%', startangle=90)
            ax.set_title(f'{label} (n={total} residues)')

    plt.suptitle('Amino Acid Composition by Region', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "03_aa_composition.png"), dpi=150, bbox_inches='tight')
    plt.close()
    print("  >> 그림 저장: analysis_results/03_aa_composition.png")

    # CDR3 아미노산 빈도 상세
    cdr3_all = ''.join([v['CDR3'] for v in vhh_data if v['CDR3']]).replace('X', '')
    cdr3_aa_counts = Counter(cdr3_all)
    total_cdr3_aa = sum(cdr3_aa_counts.values())

    print(f"\n  [CDR3 아미노산 빈도]")
    for aa, count in cdr3_aa_counts.most_common():
        bar_len = int(count / total_cdr3_aa * 50)
        bar = '█' * bar_len
        # 물리화학적 특성 표시
        group = "?"
        for g, aas in AA_GROUPS.items():
            if aa in aas:
                group = g[:3]
                break
        print(f"    {aa} ({group:3s}): {count:3d} ({count/total_cdr3_aa*100:5.1f}%) {bar}")


# =============================================================================
# 6단계: 서열 유사도 및 클러스터링
# =============================================================================
print("\n" + "=" * 70)
print("  6단계: 서열 유사도 분석")
print("=" * 70)
print("""
[설명] 서열 유사도(Sequence Identity):
  두 서열 간의 동일한 잔기 비율을 의미합니다.

  예시:  QVQLVESGG  vs  QVQLVESAG
         ████████       ████████
         → 8/9 = 88.9% 동일성

  높은 유사도 (>80%): 같은 생식세포 유전자(germline)에서 유래
  낮은 유사도 (<50%): 다른 계통의 항체

  나노바디 라이브러리에서:
  - 높은 전체 유사도 → 제한된 germline 다양성
  - CDR3 다양성이 높음 → 좋은 라이브러리
""")

if vhh_data:
    full_seqs = [v['full_vhh'] for v in vhh_data]
    unique_seqs = list(set(full_seqs))

    # 전체 VHH 쌍별 유사도
    if len(unique_seqs) > 1:
        similarities = []
        for i in range(len(unique_seqs)):
            for j in range(i + 1, len(unique_seqs)):
                s1, s2 = unique_seqs[i], unique_seqs[j]
                # 길이가 다르면 짧은 쪽에 맞춤
                min_len = min(len(s1), len(s2))
                matches = sum(a == b for a, b in zip(s1[:min_len], s2[:min_len]))
                sim = matches / min_len * 100
                similarities.append(sim)

        print(f"  [전체 VHH 쌍별 유사도]")
        print(f"  비교 쌍 수: {len(similarities)}")
        print(f"  평균 유사도: {np.mean(similarities):.1f}%")
        print(f"  범위: {min(similarities):.1f}% ~ {max(similarities):.1f}%")

        # CDR3 쌍별 유사도
        cdr3_unique = list(set(v['CDR3'] for v in vhh_data if v['CDR3']))
        cdr3_sims = []
        for i in range(len(cdr3_unique)):
            for j in range(i + 1, len(cdr3_unique)):
                s1, s2 = cdr3_unique[i], cdr3_unique[j]
                min_len = min(len(s1), len(s2))
                if min_len > 0:
                    matches = sum(a == b for a, b in zip(s1[:min_len], s2[:min_len]))
                    sim = matches / min_len * 100
                    cdr3_sims.append(sim)

        if cdr3_sims:
            print(f"\n  [CDR3 쌍별 유사도]")
            print(f"  비교 쌍 수: {len(cdr3_sims)}")
            print(f"  평균 유사도: {np.mean(cdr3_sims):.1f}%")
            print(f"  범위: {min(cdr3_sims):.1f}% ~ {max(cdr3_sims):.1f}%")

        # 그림: 유사도 분포
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        axes[0].hist(similarities, bins=20, color='purple', edgecolor='black', alpha=0.7)
        axes[0].set_xlabel('Pairwise Sequence Identity (%)')
        axes[0].set_ylabel('Count')
        axes[0].set_title('Full VHH Pairwise Identity')
        axes[0].axvline(np.mean(similarities), color='red', linestyle='--',
                       label=f'Mean: {np.mean(similarities):.1f}%')
        axes[0].legend()

        if cdr3_sims:
            axes[1].hist(cdr3_sims, bins=20, color='teal', edgecolor='black', alpha=0.7)
            axes[1].set_xlabel('Pairwise Sequence Identity (%)')
            axes[1].set_ylabel('Count')
            axes[1].set_title('CDR3 Pairwise Identity')
            axes[1].axvline(np.mean(cdr3_sims), color='red', linestyle='--',
                           label=f'Mean: {np.mean(cdr3_sims):.1f}%')
            axes[1].legend()

        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, "04_similarity_distribution.png"), dpi=150, bbox_inches='tight')
        plt.close()
        print("\n  >> 그림 저장: analysis_results/04_similarity_distribution.png")


# =============================================================================
# 7단계: CDR3 위치별 아미노산 빈도 히트맵
# =============================================================================
print("\n" + "=" * 70)
print("  7단계: CDR3 위치별 아미노산 빈도 (히트맵)")
print("=" * 70)
print("""
[설명] 위치별 아미노산 빈도 분석:
  CDR3의 각 위치에서 어떤 아미노산이 나타나는지 분석합니다.

  보존된 위치: 한 가지 아미노산이 지배적 → 구조/기능에 필수
  가변적 위치: 다양한 아미노산 → 항원 결합 다양성에 기여

  Shannon Entropy (H):
    H = -Σ p_i × log2(p_i)
    H = 0: 완전 보존 (한 가지 아미노산만 존재)
    H > 3: 매우 가변적 (다양한 아미노산)
""")

if vhh_data:
    # 가장 흔한 CDR3 길이의 서열만 사용
    cdr3_by_len = defaultdict(list)
    for v in vhh_data:
        if v['CDR3']:
            cdr3_by_len[len(v['CDR3'])].append(v['CDR3'])

    if cdr3_by_len:
        most_common_len = max(cdr3_by_len.keys(), key=lambda k: len(cdr3_by_len[k]))
        cdr3_subset = cdr3_by_len[most_common_len]

        print(f"  가장 빈번한 CDR3 길이: {most_common_len}aa ({len(cdr3_subset)}개 서열)")

        if len(cdr3_subset) >= 2:
            aa_list = sorted(set(''.join(cdr3_subset)))
            freq_matrix = np.zeros((len(aa_list), most_common_len))

            for seq in cdr3_subset:
                for pos, aa in enumerate(seq):
                    if aa in aa_list:
                        freq_matrix[aa_list.index(aa)][pos] += 1

            freq_matrix = freq_matrix / len(cdr3_subset) * 100

            # 히트맵
            fig, ax = plt.subplots(figsize=(max(10, most_common_len * 0.8),
                                             max(6, len(aa_list) * 0.4)))
            im = ax.imshow(freq_matrix, cmap='YlOrRd', aspect='auto')
            ax.set_yticks(range(len(aa_list)))
            ax.set_yticklabels(aa_list)
            ax.set_xticks(range(most_common_len))
            ax.set_xticklabels([str(i + 1) for i in range(most_common_len)])
            ax.set_xlabel('CDR3 Position')
            ax.set_ylabel('Amino Acid')
            ax.set_title(f'CDR3 Position-Specific AA Frequency (len={most_common_len}, n={len(cdr3_subset)})')
            plt.colorbar(im, ax=ax, label='Frequency (%)')

            for i in range(len(aa_list)):
                for j in range(most_common_len):
                    val = freq_matrix[i][j]
                    if val > 0:
                        color = 'white' if val > 50 else 'black'
                        ax.text(j, i, f'{val:.0f}', ha='center', va='center',
                               fontsize=7, color=color)

            plt.tight_layout()
            plt.savefig(os.path.join(OUTPUT_DIR, "05_cdr3_position_freq.png"), dpi=150, bbox_inches='tight')
            plt.close()
            print("  >> 그림 저장: analysis_results/05_cdr3_position_freq.png")

            # Shannon entropy
            print(f"\n  [CDR3 위치별 Shannon Entropy]")
            entropies = []
            for pos in range(most_common_len):
                col = freq_matrix[:, pos] / 100
                col = col[col > 0]
                entropy = -np.sum(col * np.log2(col)) if len(col) > 0 else 0
                entropies.append(entropy)
                if entropy < 0.5:
                    status = "매우 보존"
                elif entropy < 1.5:
                    status = "보존"
                elif entropy < 2.5:
                    status = "중간"
                else:
                    status = "가변"
                bar = '▓' * int(entropy * 5)
                print(f"    위치 {pos+1:2d}: H={entropy:.2f} [{status:5s}] {bar}")

    # 그림 6: VHH 전체 구조 다이어그램
    fig, ax = plt.subplots(figsize=(14, 4))

    # 모든 VHH의 평균 영역 길이 계산
    region_names = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']
    avg_lens = {}
    for region in region_names:
        lens = [len(v[region]) for v in vhh_data if v[region]]
        avg_lens[region] = np.mean(lens) if lens else 0

    # 막대 그래프로 영역 표시
    colors_region = {
        'FR1': '#a8d8ea', 'CDR1': '#e74c3c',
        'FR2': '#a8d8ea', 'CDR2': '#2ecc71',
        'FR3': '#a8d8ea', 'CDR3': '#3498db',
        'FR4': '#a8d8ea'
    }

    x_pos = 0
    for region in region_names:
        width = avg_lens[region]
        color = colors_region[region]
        rect = plt.Rectangle((x_pos, 0), width, 1, facecolor=color,
                             edgecolor='black', linewidth=1)
        ax.add_patch(rect)
        # 라벨
        ax.text(x_pos + width / 2, 0.5, f'{region}\n({width:.0f}aa)',
               ha='center', va='center', fontsize=9, fontweight='bold')
        x_pos += width

    ax.set_xlim(-2, x_pos + 2)
    ax.set_ylim(-0.5, 1.8)
    ax.set_xlabel('Position (amino acids)')
    ax.set_title('VHH Domain Structure (Average Region Lengths)')
    ax.set_yticks([])

    # 범례
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#a8d8ea', edgecolor='black', label='Framework (FR)'),
        Patch(facecolor='#e74c3c', edgecolor='black', label='CDR1'),
        Patch(facecolor='#2ecc71', edgecolor='black', label='CDR2'),
        Patch(facecolor='#3498db', edgecolor='black', label='CDR3'),
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "06_vhh_structure.png"), dpi=150, bbox_inches='tight')
    plt.close()
    print("  >> 그림 저장: analysis_results/06_vhh_structure.png")


# =============================================================================
# 8단계: FASTA 파일 및 요약 저장
# =============================================================================
print("\n" + "=" * 70)
print("  8단계: 결과 파일 저장")
print("=" * 70)

if vhh_data:
    # VHH 서열 FASTA
    fasta_records = []
    for i, v in enumerate(vhh_data):
        # X(internal stop)를 제거한 서열 저장
        clean_seq = v['full_vhh'].replace('X', '')
        record = SeqRecord(
            Seq(clean_seq),
            id=f"VHH_{i+1:03d}",
            description=f"CDR3={v['CDR3']} CDR3_len={len(v['CDR3'])}aa total_len={len(clean_seq)}aa"
        )
        fasta_records.append(record)

    fasta_path = os.path.join(OUTPUT_DIR, "vhh_sequences.fasta")
    SeqIO.write(fasta_records, fasta_path, "fasta")
    print(f"  VHH 단백질 서열: {fasta_path} ({len(fasta_records)}개)")

    # CDR3 서열 FASTA
    cdr3_records = []
    for i, v in enumerate(vhh_data):
        if v['CDR3']:
            record = SeqRecord(
                Seq(v['CDR3'].replace('X', '')),
                id=f"CDR3_{i+1:03d}",
                description=f"len={len(v['CDR3'])}aa"
            )
            cdr3_records.append(record)

    cdr3_path = os.path.join(OUTPUT_DIR, "cdr3_sequences.fasta")
    SeqIO.write(cdr3_records, cdr3_path, "fasta")
    print(f"  CDR3 서열: {cdr3_path} ({len(cdr3_records)}개)")

    # 요약 TSV
    summary_path = os.path.join(OUTPUT_DIR, "vhh_summary.tsv")
    with open(summary_path, 'w') as f:
        f.write("ID\tVHH_Length\tFR1\tCDR1\tCDR1_Len\tFR2\tCDR2\tCDR2_Len\tFR3\tCDR3\tCDR3_Len\tFR4\n")
        for i, v in enumerate(vhh_data):
            f.write(f"VHH_{i+1:03d}\t{len(v['full_vhh'])}\t"
                   f"{v['FR1']}\t{v['CDR1']}\t{len(v['CDR1'])}\t"
                   f"{v['FR2']}\t{v['CDR2']}\t{len(v['CDR2'])}\t"
                   f"{v['FR3']}\t{v['CDR3']}\t{len(v['CDR3'])}\t{v['FR4']}\n")
    print(f"  분석 요약: {summary_path}")

    # 핵산 서열 (역상보 후) FASTA 저장
    nt_records = []
    for i, r in enumerate(records):
        rc = r.seq.reverse_complement()
        nt_records.append(SeqRecord(rc, id=f"Read_{i+1:03d}_RC",
                                     description=f"reverse_complement of {r.id}"))
    nt_path = os.path.join(OUTPUT_DIR, "reads_reverse_complement.fasta")
    SeqIO.write(nt_records, nt_path, "fasta")
    print(f"  역상보 핵산 서열: {nt_path} ({len(nt_records)}개)")


# =============================================================================
# 최종 요약
# =============================================================================
print("\n" + "=" * 70)
print("  분석 완료 - 최종 요약")
print("=" * 70)
print(f"""
  ┌─────────────────────────────────────────────────────────┐
  │              데이터 요약                                  │
  ├─────────────────────────────────────────────────────────┤
  │  프로젝트: SRP124616                                     │
  │  데이터: SRR6269034 (PacBio CCS)                         │
  │  총 리드 수: {len(records):3d}개                                     │
  │  평균 리드 길이: {np.mean(lengths):.0f} bp                           │
  │  평균 품질: Q{np.mean(all_quals):.0f} (매우 높음)                     │
  └─────────────────────────────────────────────────────────┘
""")

if vhh_data:
    cdr3_seqs_all = [v['CDR3'] for v in vhh_data if v['CDR3']]
    print(f"""  ┌─────────────────────────────────────────────────────────┐
  │              VHH 분석 요약                                │
  ├─────────────────────────────────────────────────────────┤
  │  식별된 VHH: {len(vhh_data):3d}개 / {len(records)}개 리드                      │
  │  고유 VHH 서열: {len(set(v['full_vhh'] for v in vhh_data)):3d}개                                  │
  │  고유 CDR3: {len(set(cdr3_seqs_all)):3d}개                                      │
  │  CDR3 길이: {min(len(c) for c in cdr3_seqs_all)} ~ {max(len(c) for c in cdr3_seqs_all)} aa                              │
  │  CDR3 다양성: {len(set(cdr3_seqs_all))/len(cdr3_seqs_all)*100:.1f}%                                │
  └─────────────────────────────────────────────────────────┘
""")

print("""  [생성된 파일]
  analysis_results/
  ├── 01_read_stats.png              - 리드 길이/품질 분포
  ├── 02_cdr_length_distribution.png - CDR 길이 분포
  ├── 03_aa_composition.png          - 아미노산 조성 파이차트
  ├── 04_similarity_distribution.png - 서열 유사도 분포
  ├── 05_cdr3_position_freq.png      - CDR3 위치별 아미노산 빈도
  ├── 06_vhh_structure.png           - VHH 도메인 구조도
  ├── vhh_sequences.fasta            - VHH 전체 서열
  ├── cdr3_sequences.fasta           - CDR3 서열
  ├── vhh_summary.tsv                - 분석 요약 표
  └── reads_reverse_complement.fasta - 역상보 핵산 서열

  [후속 분석 제안]
  1. ANARCI/IMGT V-QUEST → 정확한 CDR 넘버링
  2. IgBLAST → germline 유전자 할당
  3. BLAST → 기존 나노바디 DB 검색 (sdAb-DB 등)
  4. AlphaFold2/ESMFold → 3D 구조 예측
  5. 계통수(phylogenetic tree) → 클론 진화 분석
  6. 더 많은 리드로 deep sequencing → 라이브러리 전체 다양성 파악
""")

print("=" * 70)
print("  케이스 스터디 분석 완료!")
print("=" * 70)
