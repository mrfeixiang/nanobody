# 使用 Protenix v2 设计 EGFRvIII 纳米抗体：完整操作指南

## 目录

1. [背景知识](#1-背景知识)
2. [环境准备与安装](#2-环境准备与安装)
3. [准备 EGFRvIII 抗原信息](#3-准备-egfrviii-抗原信息)
4. [方案一：使用 Protenix v2 CLI 进行结构预测与纳米抗体-抗原对接](#4-方案一使用-protenix-v2-cli-进行结构预测与纳米抗体-抗原对接)
5. [方案二：使用 PXDesign 进行从头设计（De Novo Design）](#5-方案二使用-pxdesign-进行从头设计de-novo-design)
6. [方案三：使用 Protenix Web Server 在线操作](#6-方案三使用-protenix-web-server-在线操作)
7. [结果分析与验证](#7-结果分析与验证)
8. [完整工作流程示例](#8-完整工作流程示例)
9. [常见问题与注意事项](#9-常见问题与注意事项)
10. [参考资源](#10-参考资源)

---

## 1. 背景知识

### 1.1 什么是 EGFRvIII？

EGFRvIII（表皮生长因子受体变体III）是 EGFR 最常见的胞外结构域截短突变体，主要特征：

- **基因变异**：外显子2-7的框内缺失（in-frame deletion）
- **蛋白变化**：胞外结构域缺失268个氨基酸（野生型第6-273位），并在融合点处插入一个新的甘氨酸残基
- **功能特点**：组成型活化（constitutively active），不结合任何已知配体
- **临床意义**：在胶质母细胞瘤（GBM）等肿瘤中高频表达，是重要的治疗靶点
- **PDB结构**：8UKX（pH 7.0下的EGFRvIII胞外区晶体结构）；8UKV（纳米抗体34E5与EGFRvIII复合物结构）

### 1.2 什么是纳米抗体（Nanobody / VHH）？

纳米抗体是来源于骆驼科动物的重链抗体（HCAb）可变区（VHH），具有以下优势：

- 分子量小（~15 kDa），组织穿透性强
- 稳定性高，易于基因工程改造
- 可识别传统抗体难以到达的隐蔽表位
- 包含3个互补决定区（CDR1, CDR2, CDR3），CDR3通常最长且对抗原结合最关键

### 1.3 什么是 Protenix v2？

Protenix v2 是字节跳动（ByteDance）开发的开源生物分子结构预测系统，主要特性：

- **模型参数**：~464M，相比v1有显著提升
- **抗体-抗原预测**：DockQ > 0.23阈值下，相比v1提升9-13个百分点
- **效率提升**：仅5个seeds即可超越v1在1000个seeds下的性能
- **TFG模块**：Training-Free Guidance，在扩散采样中施加几何/物理约束
- **设计能力**：VHH-Fc设计中达到100%靶标级成功率，命中率高达48%

**PXDesign** 是基于 Protenix 基础模型构建的从头蛋白结合物设计套件，实验成功率20-73%。

---

## 2. 环境准备与安装

### 2.1 系统要求

| 项目 | 最低要求 | 推荐配置 |
|------|---------|---------|
| Python | >= 3.11 | 3.11 |
| GPU | NVIDIA GPU (6GB+) | A100 40GB / 80GB |
| CUDA | 11.8+ | 12.1+ |
| 内存 | 16 GB | 64 GB+ |
| 磁盘 | 50 GB | 200 GB+（含数据库） |

### 2.2 安装 Protenix

#### 方法一：pip 安装（推荐）

```bash
# 创建虚拟环境
conda create -n protenix python=3.11 -y
conda activate protenix

# 安装 protenix（使用官方 PyPI 源确保最新版本）
pip install --upgrade protenix --index-url https://pypi.org/simple

# 验证安装
protenix --help
```

#### 方法二：从 GitHub 源码安装

```bash
git clone https://github.com/bytedance/Protenix.git
cd Protenix
pip install -e .
```

#### 方法三：Docker 安装

```bash
# 从 Protenix 仓库构建 Docker 镜像
git clone https://github.com/bytedance/Protenix.git
cd Protenix
docker build -t protenix .

# 运行容器
docker run --gpus all -it -v $(pwd)/data:/data protenix bash
```

### 2.3 安装外部依赖

Docker 用户已预装以下工具，其他用户需手动安装：

```bash
# Ubuntu/Debian
sudo apt-get install kalign hmmer

# 或通过 conda
conda install -c bioconda kalign2 hmmer
```

### 2.4 验证安装

```bash
# 检查版本
pip show protenix
# 应显示 Version: 2.0.0 或更高

# 检查 CLI 可用命令
protenix --help
```

可用的主要 CLI 命令：

| 命令 | 功能 |
|------|------|
| `protenix pred` / `predict` | 模型推理/结构预测 |
| `protenix json` / `tojson` | PDB/CIF 文件转 JSON |
| `protenix msa` | 生成多序列比对 (MSA) |
| `protenix mt` / `msatemplate` | MSA + 模板搜索 |
| `protenix prep` / `inputprep` | 完整预处理流程 |

---

## 3. 准备 EGFRvIII 抗原信息

### 3.1 获取 EGFRvIII 序列

EGFRvIII 的序列通过以下方式构建：

1. 从 UniProt 下载 EGFR 野生型序列（UniProt ID: P00533）
2. 删除信号肽（第1-25位）
3. 删除第6-273位氨基酸（268个残基）
4. 在第5位后插入新的甘氨酸（G）残基

#### EGFRvIII 胞外域参考序列（基于 PDB 8UKX）

```
LEEKKGNYVVTDHGSCVRACGADSYEMEEDGVRKCKQLREEEECQKPGASFK
KSQLKKLEPLKNLLDKSLCDIGGNALEVVDKSGNATHDNIDRIIKYGLRCEPS
TLHYLPRSCPESCEGPCELVKGKGKSGKNQCEPCLPGKDSGEVCQPCGPKGE
LDGHQMFAYQADEEATCKMSSGGRIIDKNLCETDSTSGNSIVPEALQFYNDDL
ACPCEIGGHSECDHGKNLIEETFLIIEEKIPLCRKELKEKFQEALRMPEFEMDI
HFTNGRIRDLNQHCSYQRKL
```

> **注意**：以上为简化参考序列。实际使用时请以PDB 8UKX或UniProt数据库中的精确序列为准。

### 3.2 获取 EGFRvIII 结构

```bash
# 下载 EGFRvIII 胞外域晶体结构 (pH 7.0)
wget https://files.rcsb.org/download/8UKX.pdb -O EGFRvIII_8UKX.pdb

# 下载纳米抗体-EGFRvIII 复合物结构（作为参考）
wget https://files.rcsb.org/download/8UKV.pdb -O nanobody_EGFRvIII_8UKV.pdb
```

### 3.3 将 PDB 转换为 Protenix JSON 格式

```bash
protenix json --input ./EGFRvIII_8UKX.pdb --out_dir ./egfrviii_json --altloc first
```

### 3.4 确定关键表位残基

已知的 EGFRvIII 关键表位区域：

| 区域 | 描述 | 重要性 |
|------|------|--------|
| 融合接点（Junction） | 外显子1和外显子8融合处，含新甘氨酸 | EGFRvIII特异性表位 |
| Domain I | 截短后保留的N端区域 | 可及性好 |
| Domain IV | C端区域 | 已有纳米抗体（如34E5）靶向此域 |

> **提示**：PDB 8UKV 中的纳米抗体 34E5 靶向 EGFRvIII 的 Domain IV，该表位在未激活的完整 EGFR 中被遮蔽，因此具有 EGFRvIII 选择性。

---

## 4. 方案一：使用 Protenix v2 CLI 进行结构预测与纳米抗体-抗原对接

此方案适用于：**已有候选纳米抗体序列**，需要预测其与 EGFRvIII 的复合物结构。

### 4.1 准备输入 JSON 文件

创建文件 `nanobody_egfrviii_pred.json`：

```json
[
  {
    "name": "VHH_EGFRvIII_complex",
    "sequences": [
      {
        "proteinChain": {
          "sequence": "此处填入你的纳米抗体VHH序列",
          "count": 1,
          "id": ["H"]
        }
      },
      {
        "proteinChain": {
          "sequence": "此处填入EGFRvIII胞外域序列",
          "count": 1,
          "id": ["A"]
        }
      }
    ],
    "constraint": {
      "pocket": {
        "binder_chain": {
          "entity": 1,
          "copy": 1
        },
        "contact_residues": [
          {"entity": 2, "copy": 1, "position": 1},
          {"entity": 2, "copy": 1, "position": 2},
          {"entity": 2, "copy": 1, "position": 3},
          {"entity": 2, "copy": 1, "position": 4},
          {"entity": 2, "copy": 1, "position": 5}
        ],
        "max_distance": 6
      }
    }
  }
]
```

#### JSON 格式说明

| 字段 | 说明 |
|------|------|
| `name` | 任务名称，用于输出目录命名 |
| `sequences` | 序列列表，每个元素是一条链 |
| `proteinChain.sequence` | 氨基酸序列（标准20种 + X表示未知） |
| `proteinChain.count` | 该链的拷贝数 |
| `proteinChain.id` | 链ID标识符 |
| `constraint.pocket` | 口袋约束，指定结合界面 |
| `binder_chain` | 结合物链（此处为纳米抗体） |
| `contact_residues` | 抗原上的接触残基（表位残基） |
| `max_distance` | 最大接触距离（埃，通常6 A） |

#### 使用已知纳米抗体作为示例

以 PDB 8UKV 中的纳米抗体 34E5 为例：

```json
[
  {
    "name": "34E5_EGFRvIII_prediction",
    "sequences": [
      {
        "proteinChain": {
          "sequence": "此处从PDB 8UKV中提取34E5的VHH序列",
          "count": 1,
          "id": ["H"]
        }
      },
      {
        "proteinChain": {
          "sequence": "此处从PDB 8UKX中提取EGFRvIII序列",
          "count": 1,
          "id": ["A"]
        }
      }
    ]
  }
]
```

### 4.2 生成 MSA（多序列比对）

MSA 对提高预测精度至关重要：

```bash
# 为纳米抗体序列生成 MSA
protenix msa --input nanobody.fasta --out_dir ./msa_output --msa_server_mode protenix

# 或使用完整预处理流程（包含MSA + 模板搜索）
protenix prep --input nanobody_egfrviii_pred.json --out_dir ./prep_output
```

这会自动搜索 MSA 和模板，并更新输入 JSON 文件中的相关路径。

### 4.3 运行结构预测

```bash
# 基本预测（使用 Protenix v2 模型）
protenix pred \
  -i nanobody_egfrviii_pred.json \
  -o ./output_prediction \
  -n protenix-v2 \
  -s 101,102,103,104,105

# 启用模板
protenix pred \
  -i nanobody_egfrviii_pred.json \
  -o ./output_prediction \
  -n protenix-v2 \
  -s 101,102,103,104,105 \
  --use_template true

# 启用 Training-Free Guidance (TFG) 以增强物理合理性
protenix pred \
  -i nanobody_egfrviii_pred.json \
  -o ./output_prediction \
  -n protenix-v2 \
  -s 101,102,103,104,105 \
  --use_tfg_guidance true
```

#### 关键参数说明

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-i / --input` | 输入 JSON 文件路径 | 必填 |
| `-o / --output` | 输出目录 | 必填 |
| `-n / --model_name` | 模型名称 | `protenix-v2` |
| `-s / --seeds` | 随机种子（逗号分隔） | 101 |
| `--use_template` | 是否使用模板 | false |
| `--use_msa` | 是否使用 MSA | true |
| `--use_tfg_guidance` | 启用 TFG 物理约束引导 | false |
| `--dtype` | 精度 (bf16/fp32) | bf16 |
| `--enable_cache` | 启用缓存加速 | false |

### 4.4 输出文件

预测完成后，输出目录结构：

```
output_prediction/
└── VHH_EGFRvIII_complex/
    ├── 101/
    │   ├── VHH_EGFRvIII_complex_101_sample_0.cif     # 预测结构
    │   └── VHH_EGFRvIII_complex_101_summary_confidence_sample_0.json  # 置信度
    ├── 102/
    │   ├── ...
    └── ...
```

#### 置信度指标

| 指标 | 说明 | 理想值 |
|------|------|--------|
| `plddt` | 每残基置信度 | > 70（可靠），> 90（高置信） |
| `ptm` | 预测TM-score | > 0.5 |
| `iptm` | 界面预测TM-score | > 0.6（良好对接） |
| `ranking_score` | 综合排名分数 | 越高越好 |
| `has_clash` | 是否存在空间冲突 | false |

---

## 5. 方案二：使用 PXDesign 进行从头设计（De Novo Design）

此方案适用于：**尚无候选纳米抗体**，需要从头设计靶向 EGFRvIII 的全新结合物。

### 5.1 PXDesign 简介

PXDesign 是基于 Protenix 基础模型的从头蛋白结合物设计套件：

- **PXDesign-d**：基于 DiT 架构的扩散模型，用于骨架生成
- **ProteinMPNN**：为生成的骨架分配序列
- **Protenix + AF2-IG**：预测复合物结构并通过置信度分数筛选
- **Foldseek**：对通过筛选的设计进行聚类以保持多样性

### 5.2 设计流程概览

```
靶标抗原（EGFRvIII）+ 表位定义
         ↓
    [Step 1] 骨架生成（PXDesign-d 扩散模型）
    - 给定靶标条件，生成大量候选结合物骨架
         ↓
    [Step 2] 序列设计（ProteinMPNN）
    - 为每个骨架分配最优氨基酸序列
         ↓
    [Step 3] 结构预测与筛选（Protenix v2）
    - 预测复合物结构
    - 用 ipTM 和结构一致性筛选
         ↓
    [Step 4] 聚类与多样性选择（Foldseek）
    - 按结构相似性聚类
    - 从每个簇中选择代表候选物
         ↓
    最终候选纳米抗体序列与结构
```

### 5.3 使用 PXDesign Web Server

PXDesign 目前主要通过 Protenix 服务器提供：

**访问地址**：https://protenix-server.com/

#### 操作步骤：

1. **注册/登录** Protenix Server 账号

2. **选择 PXDesign 模式**
   - 在服务器主页找到 "Design" 或 "PXDesign" 入口

3. **上传靶标结构**
   - 上传 EGFRvIII 的 PDB 文件（8UKX.pdb）
   - 或直接输入 PDB ID: `8UKX`

4. **定义表位/结合区域**
   - 指定 EGFRvIII 上的目标结合残基
   - 可选择融合接点区域（EGFRvIII特异性）或 Domain IV 区域

5. **设置设计参数**
   - 结合物类型：选择 VHH / Nanobody
   - 生成数量：建议 100-1000 个候选物
   - 筛选阈值：ipTM, plddt 等

6. **提交任务并等待结果**

7. **下载并分析结果**

### 5.4 使用 PXDesign 开源代码（本地运行）

```bash
# 克隆 PXDesign 仓库
git clone https://github.com/bytedance/PXDesign.git
cd PXDesign

# 安装依赖
pip install -e .

# 查看帮助
python -m pxdesign --help
```

> **注意**：PXDesign 的开源代码和详细CLI文档请参考其 GitHub 仓库（https://github.com/bytedance/PXDesign）和官方文档（https://protenix.github.io/pxdesign/）。功能和接口可能随版本更新。

---

## 6. 方案三：使用 Protenix Web Server 在线操作

适合没有本地 GPU 资源的用户。

### 6.1 访问方式

| 平台 | 地址 | 说明 |
|------|------|------|
| Protenix 官方服务器 | https://protenix-server.com/ | 完整功能，含PXDesign |
| ProteinIQ | https://proteiniq.io/app/protenix | 第三方在线接口 |
| Tamarind Bio | https://app.tamarind.bio/protenix | 第三方在线接口 |

### 6.2 在线预测操作步骤

1. **打开** Protenix Server 网站
2. **创建新任务** → 选择 "Structure Prediction"
3. **输入序列**：
   - Chain 1：粘贴纳米抗体 VHH 序列
   - Chain 2：粘贴 EGFRvIII 胞外域序列
4. **选择模型**：选择 `protenix-v2`
5. **配置参数**：
   - Seeds 数量：建议 5 个以上
   - 启用 MSA 搜索
   - 可选启用模板搜索
6. **添加约束**（可选）：
   - 指定表位残基作为 pocket constraint
7. **提交任务**
8. **等待计算完成**（通常几分钟到几小时，取决于序列长度和队列）
9. **查看和下载结果**

---

## 7. 结果分析与验证

### 7.1 结构可视化

```bash
# 使用 PyMOL 查看预测结构
pymol output_prediction/VHH_EGFRvIII_complex/101/VHH_EGFRvIII_complex_101_sample_0.cif

# 使用 ChimeraX
chimerax output_prediction/VHH_EGFRvIII_complex/101/VHH_EGFRvIII_complex_101_sample_0.cif
```

#### PyMOL 分析脚本示例

```python
# pymol_analysis.py
from pymol import cmd

# 加载预测结构
cmd.load("VHH_EGFRvIII_complex_101_sample_0.cif", "complex")

# 分离链
cmd.select("nanobody", "chain H")
cmd.select("antigen", "chain A")

# 着色
cmd.color("cyan", "nanobody")
cmd.color("salmon", "antigen")

# 显示CDR区域（需根据实际编号调整）
cmd.select("cdr1", "chain H and resi 26-35")
cmd.select("cdr2", "chain H and resi 50-65")
cmd.select("cdr3", "chain H and resi 95-102")
cmd.color("yellow", "cdr1")
cmd.color("orange", "cdr2")
cmd.color("red", "cdr3")

# 显示界面
cmd.select("interface", "nanobody within 4.0 of antigen")
cmd.show("sticks", "interface")

# 保存图片
cmd.ray(2400, 2400)
cmd.png("nanobody_egfrviii_complex.png")
```

### 7.2 分析置信度分数

```python
import json

# 读取置信度文件
with open("VHH_EGFRvIII_complex_101_summary_confidence_sample_0.json") as f:
    confidence = json.load(f)

print(f"pLDDT (平均): {confidence.get('plddt', 'N/A')}")
print(f"pTM:          {confidence.get('ptm', 'N/A')}")
print(f"ipTM:         {confidence.get('iptm', 'N/A')}")
print(f"排名分数:     {confidence.get('ranking_score', 'N/A')}")
print(f"存在冲突:     {confidence.get('has_clash', 'N/A')}")

# 评估标准
iptm = confidence.get('iptm', 0)
if iptm > 0.8:
    print("界面预测质量：优秀")
elif iptm > 0.6:
    print("界面预测质量：良好")
elif iptm > 0.4:
    print("界面预测质量：一般")
else:
    print("界面预测质量：较差，建议重新设计")
```

### 7.3 对接质量评估（DockQ）

若有实验结构参考，可用 DockQ 评估：

```bash
pip install DockQ

# 计算 DockQ 分数
DockQ predicted_complex.pdb reference_complex.pdb
```

| DockQ 分数 | 质量等级 |
|-----------|---------|
| > 0.80 | 高质量 |
| 0.49 - 0.80 | 中等质量 |
| 0.23 - 0.49 | 可接受 |
| < 0.23 | 不合格 |

### 7.4 候选物排名策略

对于多个候选纳米抗体，建议按以下综合标准排名：

1. **ipTM > 0.6**：界面预测质量阈值
2. **pLDDT > 70**：结构整体置信度
3. **无空间冲突**：`has_clash = false`
4. **ranking_score**：综合排名分数排序
5. **结构多样性**：从不同聚类中选择代表

---

## 8. 完整工作流程示例

以下是一个从头到尾的完整操作示例：

### Step 1：设置工作目录

```bash
mkdir -p egfrviii_nanobody_design && cd egfrviii_nanobody_design
mkdir -p structures input output msa_data results
```

### Step 2：下载 EGFRvIII 结构

```bash
wget https://files.rcsb.org/download/8UKX.pdb -O structures/EGFRvIII_8UKX.pdb
wget https://files.rcsb.org/download/8UKV.pdb -O structures/nanobody_EGFRvIII_8UKV.pdb
```

### Step 3：转换为 JSON 格式

```bash
protenix json --input structures/EGFRvIII_8UKX.pdb --out_dir input/ --altloc first
```

### Step 4：准备输入文件

创建 `input/design_input.json`（以一个示例纳米抗体序列为例）：

```json
[
  {
    "name": "VHH_anti_EGFRvIII",
    "sequences": [
      {
        "proteinChain": {
          "sequence": "QVQLQESGGGLVQPGGSLRLSCAASGFTFSNYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDWRAYFDYWGQGTQVTVSS",
          "count": 1,
          "id": ["H"]
        }
      },
      {
        "proteinChain": {
          "sequence": "此处填入EGFRvIII胞外域完整序列",
          "count": 1,
          "id": ["A"]
        }
      }
    ],
    "constraint": {
      "pocket": {
        "binder_chain": {"entity": 1, "copy": 1},
        "contact_residues": [
          {"entity": 2, "copy": 1, "position": 1},
          {"entity": 2, "copy": 1, "position": 2},
          {"entity": 2, "copy": 1, "position": 3},
          {"entity": 2, "copy": 1, "position": 4},
          {"entity": 2, "copy": 1, "position": 5},
          {"entity": 2, "copy": 1, "position": 6}
        ],
        "max_distance": 6
      }
    }
  }
]
```

> **重要提示**：上述VHH序列为通用示例框架。`contact_residues` 中的 position 应根据你选择的 EGFRvIII 表位区域填写实际残基编号。

### Step 5：预处理（MSA + 模板搜索）

```bash
protenix prep --input input/design_input.json --out_dir msa_data/
```

### Step 6：运行结构预测

```bash
protenix pred \
  -i input/design_input.json \
  -o output/ \
  -n protenix-v2 \
  -s 101,102,103,104,105 \
  --use_template true \
  --use_tfg_guidance true
```

### Step 7：分析结果

```bash
# 查看输出文件
ls -la output/VHH_anti_EGFRvIII/

# 使用 Python 脚本批量分析置信度
python3 << 'EOF'
import json
import glob
import os

results = []
for conf_file in sorted(glob.glob("output/VHH_anti_EGFRvIII/*/VHH_anti_EGFRvIII_*_summary_confidence_*.json")):
    with open(conf_file) as f:
        data = json.load(f)
    seed = os.path.basename(os.path.dirname(conf_file))
    results.append({
        "seed": seed,
        "iptm": data.get("iptm", "N/A"),
        "ptm": data.get("ptm", "N/A"),
        "ranking_score": data.get("ranking_score", "N/A"),
    })

results.sort(key=lambda x: x.get("ranking_score", 0), reverse=True)

print(f"{'Seed':<8} {'ipTM':<8} {'pTM':<8} {'Ranking':<10}")
print("-" * 36)
for r in results:
    print(f"{r['seed']:<8} {r['iptm']:<8} {r['ptm']:<8} {r['ranking_score']:<10}")
EOF
```

### Step 8：可视化最佳结果

```bash
# 安装 PyMOL（如果尚未安装）
# conda install -c conda-forge pymol-open-source

# 查看排名最高的预测结构
pymol output/VHH_anti_EGFRvIII/101/VHH_anti_EGFRvIII_101_sample_0.cif
```

---

## 9. 常见问题与注意事项

### Q1: GPU 显存不足怎么办？

```bash
# 使用 Protenix-Mini 轻量模型
protenix pred --input input.json --model_name protenix_mini_default_v0.5.0

# 或使用 fp32 精度（更慢但更稳定）
protenix pred --input input.json -n protenix-v2 --dtype fp32
```

GPU 显存需求参考：
| 序列长度 (tokens) | 大致显存需求 |
|-------------------|-------------|
| 500 | ~6 GB |
| 1000 | ~15 GB |
| 2000 | ~40 GB |
| 4000 | ~78 GB |

### Q2: 如何提高预测准确性？

- **增加 seeds 数量**：从 5 个增加到 20-50 个
- **启用 MSA**：运行 `protenix prep` 预处理
- **启用模板**：添加 `--use_template true`
- **使用 TFG**：添加 `--use_tfg_guidance true`
- **使用口袋约束**：在 JSON 中指定已知表位残基

### Q3: 纳米抗体序列从哪里获取？

| 来源 | 说明 |
|------|------|
| 文献已报道 | 从发表论文中获取已知的抗EGFRvIII纳米抗体序列 |
| PDB 数据库 | 如 8UKV 中的 34E5 纳米抗体 |
| 免疫文库筛选 | 实验获得的候选序列 |
| 计算设计 | 通过 PXDesign 从头生成 |
| CDR 嫁接 | 从已知纳米抗体框架出发，替换 CDR 区域 |

### Q4: 如何选择最佳表位区域？

对于 EGFRvIII，建议优先考虑：

1. **融合接点区域**：EGFRvIII 独有的新表位，野生型 EGFR 上不存在，可实现高特异性
2. **Domain IV**：在 EGFRvIII 中暴露但在完整 EGFR 中被遮蔽的区域（如 34E5 纳米抗体的靶点）
3. **避免**靶向 EGFRvIII 和 EGFR 共有的保守区域（会导致脱靶毒性）

### Q5: constraint 中的 contact_residues 怎么确定？

```bash
# 方法1：从已知复合物结构中提取界面残基
python3 << 'SCRIPT'
from Bio.PDB import PDBParser, NeighborSearch

parser = PDBParser(QUIET=True)
structure = parser.get_structure("complex", "structures/nanobody_EGFRvIII_8UKV.pdb")

# 获取所有原子
nanobody_atoms = []
antigen_atoms = []
for chain in structure[0]:
    for residue in chain:
        for atom in residue:
            if chain.id == "H":  # 纳米抗体链
                nanobody_atoms.append(atom)
            else:  # 抗原链
                antigen_atoms.append(atom)

# 搜索界面残基
ns = NeighborSearch(antigen_atoms)
interface_residues = set()
for atom in nanobody_atoms:
    neighbors = ns.search(atom.coord, 5.0)
    for n in neighbors:
        res = n.get_parent()
        interface_residues.add(res.id[1])

print("抗原界面残基:", sorted(interface_residues))
SCRIPT
```

### Q6: 运行时间大概多久？

| 配置 | 5 seeds | 20 seeds |
|------|---------|----------|
| A100 80GB + 短序列 | ~5-15 分钟 | ~20-60 分钟 |
| A100 40GB + 中序列 | ~15-30 分钟 | ~1-2 小时 |
| V100 32GB + 短序列 | ~30-60 分钟 | ~2-4 小时 |

---

## 10. 参考资源

### 论文与文档

- **Protenix v2 技术报告**：[Protenix-v2: Broadening the Reach of Structure Prediction and Biomolecular Design (bioRxiv 2026)](https://www.biorxiv.org/content/10.64898/2026.04.10.717613v1)
- **PXDesign 论文**：[PXDesign: Fast, Modular, and Accurate De Novo Design of Protein Binders (bioRxiv 2025)](https://www.biorxiv.org/content/10.1101/2025.08.15.670450.full.pdf)
- **EGFRvIII 结构论文**：[Structural insights into the role and targeting of EGFRvIII (Structure 2024)](https://www.cell.com/structure/fulltext/S0969-2126(24)00195-3)
- **Protenix v1 论文**：[Protenix - Advancing Structure Prediction Through a Comprehensive AlphaFold3 Reproduction (bioRxiv 2025)](https://www.biorxiv.org/content/10.1101/2025.01.08.631967.full.pdf)

### 代码仓库

- **Protenix GitHub**：https://github.com/bytedance/Protenix
- **PXDesign GitHub**：https://github.com/bytedance/PXDesign
- **PXDesign 文档**：https://protenix.github.io/pxdesign/

### 在线服务

- **Protenix Server**：https://protenix-server.com/
- **ProteinIQ (第三方)**：https://proteiniq.io/app/protenix
- **PyPI 包**：https://pypi.org/project/protenix/

### PDB 结构

- **8UKX**：EGFRvIII 胞外区晶体结构 (pH 7.0)
- **8UKV**：纳米抗体 34E5 与 EGFRvIII 复合物晶体结构
- **8UKW**：EGFRvIII 胞外区晶体结构 (pH 5.0)

---

> **免责声明**：本指南基于截至2026年5月的公开信息编写。Protenix v2 和 PXDesign 的具体功能和接口可能随版本更新而变化，请始终参考官方最新文档。计算设计的纳米抗体候选物需经过实验验证（如 SPR 亲和力测定、细胞结合实验等）方可用于后续开发。
