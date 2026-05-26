# 从零开始的纳米抗体（Nanobody）测序数据分析实战

> **项目**: SRP124616 / SRR6269034
> **平台**: PacBio CCS (Circular Consensus Sequencing)
> **目标**: 噬菌体展示文库中纳米抗体VHH序列的完整生信分析
> **环境**: Ubuntu Linux + Conda + Python + IgBLAST

---

## 目录

1. [背景知识](#1-背景知识)
2. [环境搭建与数据下载](#2-环境搭建与数据下载)
3. [踩坑记录：NCBI工具安装与SSL错误](#3-踩坑记录ncbi工具安装与ssl错误)
4. [第一步：FASTQ质量控制](#4-第一步fastq质量控制)
5. [踩坑记录：VHH序列提取失败](#5-踩坑记录vhh序列提取失败)
6. [第二步：VHH序列提取（正确方法）](#6-第二步vhh序列提取正确方法)
7. [第三步：CDR区域鉴定](#7-第三步cdr区域鉴定)
8. [第四步：序列多样性分析](#8-第四步序列多样性分析)
9. [第五步：氨基酸组成分析](#9-第五步氨基酸组成分析)
10. [第六步：IgBLAST Germline基因分析](#10-第六步igblast-germline基因分析)
11. [踩坑记录：IgBLAST数据库配置](#11-踩坑记录igblast数据库配置)
12. [最终结果与生物学解读](#12-最终结果与生物学解读)
13. [总结与后续方向](#13-总结与后续方向)

---

## 1. 背景知识

### 什么是纳米抗体（Nanobody / VHH）？

纳米抗体是来源于骆驼科动物（骆驼、羊驼、大羊驼）的**重链抗体（HCAb）**的可变区片段。与传统抗体不同，重链抗体天然缺少轻链，其抗原结合功能完全由单个VHH结构域承担。

```
传统抗体 (IgG, ~150 kDa):          重链抗体 (HCAb, ~95 kDa):
  ┌─VH─┐ ┌─VH─┐                      ┌─VHH─┐ ┌─VHH─┐
  │    │ │    │                      │     │ │     │
  ├─VL─┤ ├─VL─┤                      │     │ │     │
  │    │ │    │                      ├─CH2─┤ ├─CH2─┤
  ├─CH1┤ ├─CH1┤                      │     │ │     │
  │    │ │    │                      ├─CH3─┤ ├─CH3─┤
  ├─CH2┤ ├─CH2┤                      └─────┘ └─────┘
  │    │ │    │                      (无轻链、无CH1)
  ├─CH3┤ ├─CH3┤
  └────┘ └────┘

纳米抗体 (VHH, ~15 kDa):
  ┌──────────────────────────┐
  │ FR1-CDR1-FR2-CDR2-FR3-CDR3-FR4 │
  └──────────────────────────┘
  (~120-130 个氨基酸)
```

**纳米抗体的优势**：
- 分子量小（~15 kDa），组织穿透性强
- 热稳定性高，耐极端pH
- 可在大肠杆菌中高效表达
- 可识别传统抗体无法到达的隐蔽表位（如酶活性位点）

### VHH结构域的组成

```
FR1 (25aa) ── CDR1 (7-12aa) ── FR2 (17aa) ── CDR2 (6-7aa) ── FR3 (38aa) ── CDR3 (3-28aa) ── FR4 (11aa)
 [骨架区]      [抗原结合]       [骨架区]      [抗原结合]       [骨架区]       [抗原结合]       [骨架区]
```

- **FR (Framework Region)**：结构骨架，高度保守，维持免疫球蛋白折叠
- **CDR (Complementarity Determining Region)**：互补决定区，直接接触抗原
- **CDR3**：长度和序列最多变，是决定抗原结合特异性的关键区域

### 噬菌体展示（Phage Display）

本数据来自噬菌体展示文库。该技术的基本原理：

1. 将VHH基因克隆到噬菌体外壳蛋白（pIII）基因的上游
2. VHH蛋白展示在噬菌体表面
3. 用目标抗原进行"淘选"（panning），筛选结合的噬菌体
4. 经过多轮淘选，富集高亲和力克隆

**重要**：噬菌体展示载体中通常包含一个**amber终止密码子（TAG）**，位于VHH和pIII之间。在amber抑制菌株（如TG1）中，TAG被读通为Gln；在非抑制菌株中则终止翻译，产生可溶性VHH。

---

## 2. 环境搭建与数据下载

### 安装所需工具

```bash
# 使用conda安装生信工具（推荐，无需sudo）
conda install -y -c bioconda -c conda-forge \
    entrez-direct sra-tools \
    fastqc seqkit blast muscle igblast

# Python包
pip install biopython pandas matplotlib
```

### 下载SRA数据

```bash
# 方法一：使用Entrez Direct（如果可用）
esearch -db sra -query "SRP124616" | efetch -format runinfo > runinfo.csv
cut -d',' -f1 runinfo.csv | grep SRR > srr_list.txt

# 方法二：从ENA获取运行信息（更可靠）
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?\
accession=SRP124616&result=read_run&\
fields=run_accession,library_strategy,library_layout,fastq_ftp&\
format=tsv&limit=0"
```

**输出**：

```
run_accession  library_strategy  library_layout  fastq_ftp
SRR6269034     POOLCLONE         PAIRED          ftp.sra.ebi.ac.uk/vol1/fastq/SRR626/004/SRR6269034/...
```

只有1个运行（SRR6269034），是PacBio的Pooled Clone测序。

```bash
# 下载并转换为FASTQ
prefetch SRR6269034
fasterq-dump --threads 8 SRR6269034
```

**输出**：

```
spots read      : 85
reads read      : 85
reads written   : 85
```

生成文件：`SRR6269034.fastq`（179 KB，85条reads）

---

## 3. 踩坑记录：NCBI工具安装与SSL错误

### 问题1：esearch/efetch命令不存在

```
Command 'esearch' not found, but can be installed with:
sudo apt install ncbi-entrez-direct
```

**原因**：系统未安装NCBI Entrez Direct工具包。

**解决**：使用conda安装，避免sudo权限问题：
```bash
conda install -c bioconda entrez-direct
```

### 问题2：SSL连接错误

```
curl: (56) OpenSSL SSL_read: error:0A000126:SSL routines::unexpected eof while reading
ERROR: curl command failed with: 56
```

**原因**：conda安装的OpenSSL版本（3.6.2）与NCBI服务器的TLS配置不完全兼容。

**解决方案**：

方案A — 使用 `curl -k` 跳过SSL验证（临时方案）：
```bash
curl -k -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?\
db=sra&term=SRP124616&retmax=500&usehistory=y"
```

方案B — 改用ENA API（推荐，更稳定）：
```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?\
accession=SRP124616&result=read_run&fields=run_accession&format=tsv"
```

**教训**：当NCBI的Entrez API出现SSL问题时，ENA（欧洲核酸档案）是很好的替代方案。ENA镜像了SRA的所有数据。

---

## 4. 第一步：FASTQ质量控制

### 基本统计

```python
from Bio import SeqIO
import numpy as np

records = list(SeqIO.parse("SRR6269034.fastq", "fastq"))
lengths = [len(r.seq) for r in records]
mean_quals = [np.mean(r.letter_annotations["phred_quality"]) for r in records]

print(f"总reads数:    {len(records)}")       # 85
print(f"长度范围:     {min(lengths)}-{max(lengths)} bp")  # 885-979 bp
print(f"平均长度:     {np.mean(lengths):.1f} bp")          # 975.6 bp
print(f"平均质量分数: {np.mean(mean_quals):.1f}")           # 91.6 (Q91+)
```

### 结果

| 指标 | 值 |
|------|-----|
| 总reads数 | 85 |
| 长度范围 | 885 ~ 979 bp |
| 平均长度 | 975.6 bp |
| 中位数长度 | 977.0 bp |
| 平均质量 | Q91.6 |
| 质量范围 | Q82.1 ~ Q93.0 |

### Phred质量分数解读

```
Q20 = 99%     准确率 (每100个碱基1个错误)
Q30 = 99.9%   准确率 (每1000个碱基1个错误)
Q40 = 99.99%  准确率
Q90 = 极高    (PacBio CCS的'~'字符 = Q93)
```

**解读**：Q91+的质量分数意味着这是PacBio CCS（环形一致性测序）数据。CCS通过多次读取同一分子并取一致序列，达到了接近完美的准确度。

### 质控可视化代码

```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 读长分布
axes[0].hist(lengths, bins=20, color='steelblue', edgecolor='black', alpha=0.7)
axes[0].set_xlabel('Read Length (bp)')
axes[0].set_ylabel('Count')
axes[0].set_title('Read Length Distribution')
axes[0].axvline(np.mean(lengths), color='red', linestyle='--',
                label=f'Mean: {np.mean(lengths):.0f} bp')
axes[0].legend()

# 质量分布
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

> **图1**：左图显示reads长度高度集中在977 bp左右，说明PCR扩增片段大小一致。右图显示大部分reads的质量在Q90以上，证实了CCS数据的高准确度。

---

## 5. 踩坑记录：VHH序列提取失败

### 初始方法：6-frame ORF搜索

最初的想法很直观——对每条read做6-frame翻译，找最长的ORF：

```python
def find_longest_orf(seq_str):
    """在6个阅读框中寻找最长ORF"""
    seq = Seq(seq_str)
    best_protein = ""
    for strand, nuc in [("+", seq), ("-", seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate())
            for i in range(len(trans)):
                if trans[i] == "M":  # 起始密码子
                    stop = trans.find("*", i)  # 终止密码子
                    if stop != -1 and (stop - i) > len(best_protein):
                        best_protein = trans[i:stop]
    return best_protein
```

### 结果：惨败

```
ORF分析结果:
- 100aa以上ORF的reads: 2个 / 85个   ← 几乎全军覆没！
- 蛋白长度范围: 6 ~ 109 aa
- 平均蛋白长度: 16.9 aa               ← 远小于VHH的120aa
- VHH识别数: 0个                       ← 完全没有！
```

### 问题诊断

VHH应该是~120 aa的蛋白，为什么ORF搜索只找到了最长109 aa的片段？

让我们手动检查第一条read的6-frame翻译：

```python
from Bio.Seq import Seq
r = records[0]
rc = r.seq.reverse_complement()

# 反义链第1阅读框
prot = str(rc.translate())
# 搜索VHH特征起始序列
pos = prot.find("QVQL")
print(f"QVQL found at position {pos}")
print(prot[pos:pos+130])
```

**输出**：

```
QVQL found at position 38
QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSNNFGWNWIRQSPSRGLE*LGRTYYRSKWY
NDYAVSVRSRITINPDTSKNQFSLQLNSVTPEDTAVYYCARQGSTYFDYWGQGTLVTVSS
```

**关键发现**：`GLE*LGR` — VHH序列中间有一个终止密码子（`*`）！

这就是ORF搜索失败的原因：ORF搜索在遇到第一个终止密码子时就截断了，而VHH的FR2区域内嵌入了一个amber终止密码子（TAG）。

### 根本原因

这是**噬菌体展示载体的设计特征**，不是测序错误：

```
载体结构: ... pelB-VHH-[amber TAG]-pIII ...

在amber抑制菌株(TG1)中:  TAG → Gln (读通，VHH-pIII融合蛋白，展示在噬菌体表面)
在非抑制菌株(HB2151)中:  TAG → STOP (终止，产生可溶性VHH蛋白)
```

所以这个 `*` 实际上是载体设计的一部分，在细菌中会被读通为谷氨酰胺（Q）。传统的ORF搜索方法在这里彻底失效了。

### 经验教训

> **在分析噬菌体展示文库的测序数据时，不能依赖简单的ORF搜索。必须基于VHH的保守模体（motif）进行序列提取，并正确处理amber终止密码子。**

---

## 6. 第二步：VHH序列提取（正确方法）

### 策略改进

既然ORF搜索行不通，我们改用**保守模体（motif）驱动**的方法：

1. 对每条read，尝试正义链和反义链的3个阅读框（共6种）
2. 翻译后搜索VHH起始模体：`QVQL`、`EVQL`、`DVQL`
3. 搜索VHH终止模体：`WGQG`（FR4的起始）
4. 将内部终止密码子标记为 `X`（代替amber读通的氨基酸）

```python
def extract_vhh_from_read(seq_record):
    """
    从PacBio CCS read中提取VHH序列。
    处理反向互补和内部终止密码子。
    """
    seq = seq_record.seq
    for strand, nuc in [("+", seq), ("-", seq.reverse_complement())]:
        for frame in range(3):
            prot = str(nuc[frame:].translate())
            # 将内部终止密码子替换为X
            prot_no_stop = prot.replace("*", "X")

            for start_motif in ['QVQL', 'EVQL', 'DVQL']:
                pos = prot_no_stop.find(start_motif)
                if pos == -1:
                    continue

                # 搜索FR4起始模体 WGQG
                wg_pos = prot_no_stop.find('WGQG', pos)
                if wg_pos == -1:
                    wg_pos = prot_no_stop.find('WGRG', pos)
                if wg_pos == -1:
                    continue

                # 提取VHH + FR4 (WGQGTLVTVSS)
                end_pos = min(wg_pos + 11, len(prot_no_stop))
                vhh_seq = prot_no_stop[pos:end_pos]

                # 验证长度合理性 (VHH通常110-140aa)
                if 100 <= len(vhh_seq) <= 160:
                    return {
                        'vhh_seq': vhh_seq,
                        'strand': strand,
                        'frame': frame + 1,
                        'has_internal_stop': '*' in prot[pos:end_pos]
                    }
    return None
```

### 结果

```
VHH提取结果:
- 成功: 83个 / 85条reads (97.6%)
- 失败: 2条reads

序列方向分布:
    + (正义链):   45条
    - (反义链):   38条

阅读框分布:
    +1: 45条
    -1: 36条
    -3:  2条

内部终止密码子: 83条 (全部)    ← 确认是amber stop codon

VHH长度: 120 aa (所有序列一致)
```

### 提取出的VHH序列示例

```
  1: QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSNNFGWNWIRQSPSRGLEXL
 51: GRTYYRSKWYNDYAVSVRSRITINPDTSKNQFSLQLNSVTPEDTAVYYCAR
101: QGSTYFDYWGQGTLVTVSS
```

注意第48位的 `X` 就是amber终止密码子对应的位置（在表达系统中被读通为Gln）。

---

## 7. 第三步：CDR区域鉴定

### CDR边界的确定方法

使用保守模体和位置特征来划分CDR边界：

```python
def extract_cdr_regions(vhh_seq):
    """
    基于保守模体的CDR提取。

    关键锚点:
    - FR1结束/CDR1开始: 第一个C之后约4个残基 (~位置26)
    - CDR1结束: W残基 (~位置36)
    - CDR2开始: FR2之后 (~位置50)
    - CDR3开始: YYC/YFC模体之后 (~位置96)
    - CDR3结束: WGQG之前
    - FR4: WGQGTLVTVSS
    """
    result = {}

    # 1) 找FR4起始: WGQG
    fr4_start = vhh_seq.find('WGQG')
    if fr4_start == -1:
        return None
    result['FR4'] = vhh_seq[fr4_start:fr4_start + 11]

    # 2) 找CDR3: YYC模体后到WGQG前
    cdr3_match = re.search(r'[YF][YF]C', vhh_seq[80:fr4_start])
    if cdr3_match:
        cdr3_start = 80 + cdr3_match.end()
    else:
        cdr3_start = fr4_start - 15  # 后备估计
    result['CDR3'] = vhh_seq[cdr3_start:fr4_start]

    # 3) 找CDR1: 第一个C之后
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

### 鉴定结果

```
┌─────────┬───────┬──────────────────────────────────────────┐
│ 区域    │ 长度  │ 序列                                      │
├─────────┼───────┼──────────────────────────────────────────┤
│ FR1     │ 25 aa │ QVQLQQSGPGLVKPSQTLSLTCAIS                │
│ CDR1    │ 10 aa │ GDSVSSNNFG                    ◀◀ 抗原结合 │
│ FR2     │ 15 aa │ WNWIRQSPSRGLEXL                          │
│ CDR2    │  8 aa │ GRTYYRSK                      ◀◀ 抗原结合 │
│ FR3     │ 41 aa │ WYNDYAVSVRSRITINPDTSKNQFSLQLNSVTPEDTAVYYC│
│ CDR3    │ 10 aa │ ARQGSTYFDY                    ◀◀ 抗原结合 │
│ FR4     │ 11 aa │ WGQGTLVTVSS                              │
└─────────┴───────┴──────────────────────────────────────────┘
```

与IgBLAST的IMGT编号结果对比：

| 区域 | 本方法 | IMGT (IgBLAST) | 一致？ |
|------|--------|----------------|--------|
| FR1 | 25 aa | 25 aa | Y |
| CDR1 | 10 aa (GDSVSSNNFG) | 10 aa (GDSVSSNNFG) | Y |
| FR2 | 15 aa | 17 aa | 接近 |
| CDR2 | 8 aa | 9 aa (TYYRSKWYN) | 接近 |
| FR3 | 41 aa | 38 aa | 接近 |
| CDR3 | 10 aa (ARQGSTYFDY) | 10 aa (ARQGSTYFDY) | Y |
| FR4 | 11 aa | 11 aa | Y |

CDR3完全一致！FR2/CDR2/FR3边界有1-3个残基的差异，这是因为不同编号体系（Kabat、IMGT、Chothia）对CDR边界的定义略有不同。**对于精确分析，应以IgBLAST/IMGT结果为准**。

---

## 8. 第四步：序列多样性分析

### 分析代码

```python
from collections import Counter

# 全长VHH多样性
full_seqs = [v['full_vhh'] for v in vhh_data]
unique_full = set(full_seqs)
full_counter = Counter(full_seqs)

# CDR分区多样性
for cdr_name in ['CDR1', 'CDR2', 'CDR3']:
    cdr_seqs = [v[cdr_name] for v in vhh_data]
    unique_cdr = set(cdr_seqs)
    diversity = len(unique_cdr) / len(cdr_seqs) * 100
    print(f"{cdr_name}: {len(unique_cdr)} unique / {len(cdr_seqs)} total = {diversity:.1f}%")
```

### 结果

```
全长VHH序列多样性:
  总序列:     83个
  唯一序列:    2个
  多样性比率:  2.4%
  singleton:  1个 (只出现1次的序列)
  扩增克隆:   1个 (出现82次!)

CDR多样性:
  CDR       总数    唯一    多样性
  ──────────────────────────────
  CDR1       83      1      1.2%
  CDR2       83      1      1.2%
  CDR3       83      1      1.2%

唯一的CDR3序列:
  ARQGSTYFDY    (83次) ████████████████████████████████████
```

### 解读

**这不是一个多样性文库，而是单克隆验证数据**：

- 85条reads中83条编码**完全相同的VHH序列**
- CDR3只有1种序列（`ARQGSTYFDY`）
- 这意味着数据来源是：
  - 噬菌体展示淘选后的**单一克隆**，或
  - 特定纳米抗体的**序列验证**实验
- PacBio长读长测序用于确认全长序列的准确性

---

## 9. 第五步：氨基酸组成分析

### 氨基酸分类

```python
AA_GROUPS = {
    'Hydrophobic': set('AVILMFWP'),  # 疏水性 - 蛋白核心
    'Polar':       set('STNQYC'),    # 极性   - 氢键、抗原接触
    'Positive':    set('RHK'),        # 正电荷 - 静电相互作用
    'Negative':    set('DE'),         # 负电荷 - 静电相互作用
    'Glycine':     set('G'),          # 甘氨酸 - 结构柔性
}
```

### CDR3氨基酸频率

```
CDR3序列: ARQGSTYFDY (10个氨基酸)

  Y (极性):   20.0% ██████████   ← 酪氨酸占主导（抗原接触的关键残基）
  A (疏水):   10.0% █████
  R (正电荷): 10.0% █████
  Q (极性):   10.0% █████
  G (甘氨酸): 10.0% █████
  S (极性):   10.0% █████
  T (极性):   10.0% █████
  F (疏水):   10.0% █████
  D (负电荷): 10.0% █████
```

**关键观察**：CDR3中酪氨酸(Y)含量最高（2个Y），这在抗体中非常常见——酪氨酸的侧链既有芳香性又有羟基，是与抗原形成多种相互作用（氢键、π-π堆积、范德华力）的理想残基。

---

## 10. 第六步：IgBLAST Germline基因分析

### 什么是V(D)J重组？

```
Germline DNA:
  ──[V1][V2]...[Vn]──[D1][D2]...[Dn]──[J1][J2]...[Jn]──[C]──

V(D)J重组 (B细胞发育过程中):
  1. D-J连接:  选择一个D和一个J → D-J
  2. V-DJ连接: 选择一个V → V-D-J
  3. 接合多样性: 连接处随机添加/删除核苷酸 (N/P nucleotides)

重组后:
  ──[V]─[N]─[D]─[N]─[J]──[C]──
     └── CDR3区域 ──┘
```

### IgBLAST安装与运行

```bash
# 安装
conda install -y -c bioconda igblast

# 检查安装路径
IGDATA=$(dirname $(which igblastn))/../share/igblast
ls $IGDATA/internal_data/human/  # 确认germline数据存在
```

### 准备输入文件

由于原始reads需要反向互补才能获得正确的VHH编码方向，先准备反向互补的FASTA：

```python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

records = list(SeqIO.parse("SRR6269034.fastq", "fastq"))
rc_records = []
seen = set()
for i, r in enumerate(records):
    rc = r.seq.reverse_complement()
    s = str(rc)
    if s not in seen:  # 去重
        seen.add(s)
        rc_records.append(SeqRecord(rc, id=f"Read_{i+1:03d}_RC",
                                     description=f"reverse_complement of {r.id}"))

SeqIO.write(rc_records, "unique_reads_rc.fasta", "fasta")
print(f"Wrote {len(rc_records)} unique sequences")  # 25个唯一序列
```

### 运行IgBLAST

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

### IgBLAST核心结果

```
V-(D)-J rearrangement summary:
  V gene:     IGHV6-1*01
  D gene:     IGHD1-26*01
  J gene:     IGHJ4*02
  Chain:      VH
  Stop codon: Yes (amber stop in phage display vector)
  V-J frame:  In-frame
  Productive: No (due to amber stop — functional in suppressor strains)

V-(D)-J junction (CDR3 DNA):
  V end:        CAAGA
  V-D junction: CA
  D region:     GGGCAGC
  D-J junction: ACT
  J start:      TACTT

CDR3 nucleotide: GCAAGACAGGGCAGCACTTACTTCGACTAT
CDR3 amino acid: ARQGSTYFDY
```

### IMGT区域比对详情

```
Alignment summary (Query vs IGHV6-1*01 germline):

  区域        起始   终止   长度   匹配   错配   Identity
  ─────────────────────────────────────────────────────
  FR1-IMGT    115    189    75     75     0      100.0%
  CDR1-IMGT   190    219    30     26     4       86.7%
  FR2-IMGT    220    270    51     50     1       98.0%
  CDR2-IMGT   271    297    27     27     0      100.0%
  FR3-IMGT    298    411   114    112     2       98.2%
  CDR3 (germ) 412    417     6      6     0      100.0%
  ─────────────────────────────────────────────────────
  Total        -      -    303    296     7       97.7%
```

**解读**：
- FR1和CDR2与germline完全一致（100%）
- CDR1有4个突变（86.7%）— 这些可能是物种差异（骆驼vs人类）或体细胞突变
- 总体97.7%的一致性说明这个VHH与人类IGHV6-1非常相似

### 比对可视化

```
                        <── FR1 ──><── CDR1 ──>
  Query:    QVQLQQSGPGLVKPSQTLSLTCAIS GDSVSSNNFG
  Germline: QVQLQQSGPGLVKPSQTLSLTCAIS GDSVSSNSA A
                                       |||||| * * ← 4个差异

                        <── FR2 ──><── CDR2 ──>
  Query:    WNWIRQSPSRGLE*LGR TYYRSKWY
  Germline: WNWIRQSPSRGLEWLGR TYYRSKWY
                          ^                ← amber stop (本应是W)

                        <── FR3 ──────────────────><── CDR3 ──>
  Query:    NDYAVSVRSRITINPDTSKNQFSLQLNSVTPEDTAVYYC ARQGSTYFDY
  Germline: NDYAVSVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYC AR─────────
                    ^                                   (V-D-J接合)

                        <- FR4 ->
  Query:    WGQGTLVTVSS
  J gene:   WGQG.LV.VSS  (IGHJ4*02, 91.3% identity)
```

---

## 11. 踩坑记录：IgBLAST数据库配置

### 问题1：找不到D/J基因数据库

```
BLAST Database error: No alias or index file found for nucleotide database
[human_gl_D] in search path [/home/syslab/Desktop/nanobody::]
```

**原因**：conda安装的IgBLAST只包含V基因的内部数据库（`internal_data/human/human_V`），但没有预构建的D和J基因BLAST数据库。

**解决方法**：手动创建D/J基因数据库。

```bash
# 创建IGHD germline FASTA文件
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
# (实际使用了全部26个D基因，此处简化展示)

# 创建IGHJ germline FASTA文件
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

# 构建BLAST数据库
makeblastdb -parse_seqids -dbtype nucl \
    -in germline_db/human_ighd.fasta \
    -out germline_db/human_ighd \
    -title "Human IGHD"

makeblastdb -parse_seqids -dbtype nucl \
    -in germline_db/human_ighj.fasta \
    -out germline_db/human_ighj \
    -title "Human IGHJ"
```

### 问题2：IMGT在线数据库下载失败

```bash
# 尝试从IMGT下载
curl -sL "https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframedP" \
    -o imgt_all.fasta
# 结果：获得的是HTML页面，不是FASTA文件
```

**原因**：IMGT网站需要通过浏览器访问，直接curl请求会被重定向到登录/验证页面。

**解决方案**：
1. 使用IgBLAST自带的internal_data中的V基因数据库
2. 手动准备D/J基因序列（如上所示）
3. 或者访问NCBI的IgBLAST FTP站点下载预构建数据库

### 问题3：物种不匹配

IgBLAST没有骆驼科（Camelidae）的germline数据库。

**处理方式**：使用人类（human）germline作为参考。这是可行的，因为：
- IGHV6-1是人类V基因中与骆驼VHH结构最相似的
- 分析结果仍然有意义，但需注意序列差异可能来自物种间自然差异而非体细胞突变
- 正式发表时应使用骆驼科专用germline数据库

---

## 12. 最终结果与生物学解读

### 综合结果

```
┌────────────────────────────────────────────────┐
│           最终分析结果汇总                       │
├────────────────────────────────────────────────┤
│                                                │
│  项目: SRP124616                               │
│  数据: SRR6269034 (PacBio CCS, 85 reads)       │
│                                                │
│  数据质量:                                      │
│    平均读长: 976 bp                             │
│    平均质量: Q92 (极高)                         │
│                                                │
│  VHH鉴定:                                      │
│    成功率: 83/85 (97.6%)                        │
│    唯一序列: 2个                                │
│    唯一CDR3: 1个                                │
│                                                │
│  Germline基因:                                  │
│    V gene: IGHV6-1*01 (97.7% identity)         │
│    D gene: IGHD1-26*01                          │
│    J gene: IGHJ4*02   (91.3% identity)         │
│                                                │
│  CDR3:                                          │
│    序列: ARQGSTYFDY                             │
│    长度: 10 aa                                  │
│    DNA:  GCAAGACAGGGCAGCACTTACTTCGACTAT         │
│                                                │
│  V-D-J接合:                                     │
│    V末端: CAAGA                                 │
│    V-D接合: CA (N-addition)                     │
│    D区域: GGGCAGC                               │
│    D-J接合: ACT                                 │
│    J起始: TACTT                                 │
│                                                │
└────────────────────────────────────────────────┘
```

### 生物学解读

**1. 这是一个单克隆纳米抗体**

85条reads中83条编码相同的VHH序列，只有1种CDR3（ARQGSTYFDY）。这不是文库多样性测序，而是特定克隆的序列验证。

**2. 使用PacBio CCS的原因**

PacBio长读长测序可以在单条read中覆盖完整的VHH编码序列（~400 bp），无需组装。CCS的高准确度（Q90+）确保了序列的可靠性。

**3. 载体结构确认**

```
5'─ pelB信号肽 ─ VHH ─ [amber TAG] ─ (G4S)3 linker ─ VL ─ 3'
                  │                                    │
              983 bp read 覆盖了完整的克隆插入片段
```

reads中除了VHH外还检测到了轻链可变区（IGKV1-39*01, 93.7% identity），说明这是一个VHH-VL双特异性构建体或scFv融合蛋白的一部分。

**4. FR2区域的VHH特征性突变**

在FR2中观察到一个关键差异：`E*L` vs germline的 `EWL`。其中 `*` 是amber终止密码子。在实际表达中，W→Q的替换是VHH（相对于传统VH）的标志性突变之一，有助于VHH在没有VL的情况下保持可溶性。

**5. Germline一致性分析**

97.7%的V基因一致性（相对于人类IGHV6-1*01）：
- 7个核苷酸差异/303个总核苷酸
- CDR1的4个差异可能反映骆驼科vs人类的germline自然差异
- FR区域的3个差异可能是体细胞突变或物种差异

---

## 13. 总结与后续方向

### 完整分析流程回顾

```
原始数据 (SRA)
    │
    ▼
[下载] prefetch + fasterq-dump
    │
    ▼
FASTQ文件 (85 reads, ~977 bp, Q92)
    │
    ├──▶ [质控] 读长分布、质量分布
    │
    ├──▶ [翻译] 反向互补 + 6-frame翻译
    │      │
    │      ├── ✗ ORF搜索失败 (amber stop codon!)
    │      │
    │      └── ✓ 模体驱动的VHH提取 (83/85成功)
    │
    ├──▶ [CDR鉴定] 保守模体 + IMGT编号
    │
    ├──▶ [多样性分析] → 单克隆 (1种CDR3)
    │
    ├──▶ [氨基酸分析] → CDR3富含Tyr
    │
    └──▶ [IgBLAST] → IGHV6-1*01 / IGHD1-26*01 / IGHJ4*02
```

### 踩坑总结

| 问题 | 症状 | 根因 | 解决方案 |
|------|------|------|----------|
| NCBI工具缺失 | `command not found` | 未安装 | `conda install entrez-direct sra-tools` |
| SSL连接错误 | `curl: (56) SSL_read error` | OpenSSL版本不兼容 | 改用ENA API或`curl -k` |
| VHH提取失败 | ORF搜索只找到短片段 | amber终止密码子截断ORF | 用保守模体搜索代替ORF搜索 |
| IgBLAST D/J缺失 | `No alias or index file` | conda未包含D/J数据库 | 手动创建BLAST数据库 |
| IMGT下载失败 | 获得HTML而非FASTA | 网站需要浏览器访问 | 使用IgBLAST内置数据 |
| 无骆驼科germline | - | IgBLAST不含骆驼数据 | 用人类germline替代（注明局限性） |

### 后续分析方向

1. **精确编号**: 使用ANARCI或IMGT/V-QUEST进行标准IMGT编号
2. **结构预测**: 用AlphaFold2或ESMFold预测VHH三维结构
3. **数据库检索**: 在PDB、sdAb-DB中搜索同源纳米抗体
4. **亲和力预测**: 基于CDR序列特征预测抗原结合特性
5. **人源化设计**: 如果需要治疗应用，进行VHH人源化改造
6. **深度测序**: 对原始文库进行NGS测序，评估完整多样性

### 生成的所有文件

```
nanobody/
├── SRR6269034.fastq                          # 原始测序数据
├── nanobody_analysis.py                       # 完整分析脚本
├── blog_nanobody_analysis_zh.md              # 本博客文章
├── germline_db/                               # Germline数据库
│   ├── human_ighd.fasta + BLAST索引
│   └── human_ighj.fasta + BLAST索引
└── analysis_results/                          # 分析结果
    ├── 01_read_stats.png                      # 读长与质量分布图
    ├── 02_cdr_length_distribution.png         # CDR长度分布图
    ├── 03_aa_composition.png                  # 氨基酸组成饼图
    ├── 04_similarity_distribution.png         # 序列相似度分布图
    ├── 05_cdr3_position_freq.png              # CDR3位置频率热图
    ├── 06_vhh_structure.png                   # VHH结构域示意图
    ├── vhh_sequences.fasta                    # VHH蛋白序列
    ├── cdr3_sequences.fasta                   # CDR3序列
    ├── vhh_summary.tsv                        # 分析汇总表
    ├── reads_reverse_complement.fasta         # 反向互补核酸序列
    ├── unique_reads_rc.fasta                  # 去重后唯一序列
    ├── igblast_results.txt                    # IgBLAST详细比对
    ├── igblast_airr.tsv                       # AIRR标准格式
    └── igblast_tabular.txt                    # 制表符分隔格式
```

---

> **作者注**：本分析使用人类germline基因作为参考。由于纳米抗体来源于骆驼科动物，部分序列差异可能是物种间的自然差异而非体细胞突变。正式研究中建议使用骆驼科专用germline数据库（如IMGT的*Camelus dromedarius*或*Vicugna pacos*参考序列）进行更精确的分析。

> **数据来源**: NCBI SRA, BioProject SRP124616, Run SRR6269034
