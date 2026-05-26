# Nanobody (VHH) Sequencing Data Analysis

Complete bioinformatics pipeline for analyzing nanobody/VHH sequences from phage display libraries using PacBio CCS long-read sequencing data.

## Project

- **BioProject**: [SRP124616](https://www.ncbi.nlm.nih.gov/sra/?term=SRP124616)
- **Run**: SRR6269034 (PacBio CCS, 85 reads)
- **Organism**: Camelid VHH (nanobody)

## Key Results

| Item | Result |
|------|--------|
| V Gene | IGHV6-1*01 (97.7% identity) |
| D Gene | IGHD1-26*01 |
| J Gene | IGHJ4*02 (91.3% identity) |
| CDR3 | ARQGSTYFDY (10 aa) |
| Clone type | Monoclonal (single clone) |

## Analysis Pipeline

```
FASTQ → Quality Control → VHH Extraction → CDR Identification
     → Diversity Analysis → Amino Acid Composition → IgBLAST Germline Analysis
```

## Blog Posts (Case Study Tutorials)

- [Chinese (中文)](blog_nanobody_analysis_zh.md) - 从零开始的纳米抗体测序数据分析实战
- [Korean (한국어)](blog_nanobody_analysis_kr.md) - 처음부터 시작하는 나노바디 시퀀싱 데이터 분석 실전 가이드

Both blogs include:
- Background knowledge on nanobodies, VHH structure, and phage display
- Step-by-step analysis with real code
- 3 detailed troubleshooting records (SSL errors, amber stop codon, IgBLAST DB setup)
- Complete results with biological interpretation

## Generated Files

```
analysis_results/
├── 01_read_stats.png              - Read length & quality distribution
├── 02_cdr_length_distribution.png - CDR length distribution
├── 03_aa_composition.png          - Amino acid composition pie charts
├── 04_similarity_distribution.png - Sequence similarity distribution
├── 05_cdr3_position_freq.png      - CDR3 position-specific frequency heatmap
├── 06_vhh_structure.png           - VHH domain structure diagram
├── vhh_sequences.fasta            - VHH protein sequences
├── cdr3_sequences.fasta           - CDR3 sequences
├── vhh_summary.tsv                - Analysis summary table
├── igblast_results.txt            - IgBLAST detailed alignments
├── igblast_airr.tsv               - AIRR standard format
└── igblast_tabular.txt            - Tabular format
```

## Requirements

```bash
conda install -c bioconda -c conda-forge entrez-direct sra-tools igblast blast
pip install biopython pandas matplotlib
```

## Usage

```bash
# Download data
prefetch SRR6269034
fasterq-dump --threads 8 SRR6269034

# Run analysis
python nanobody_analysis.py
```

## License

MIT
