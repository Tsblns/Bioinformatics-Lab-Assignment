# Window 3 Method Notes: CNV + Mutation Preprocessing

Generated: 2026-05-02 20:01:52

## Actual input files used
- master_samples_6omics.csv: `C:\Users\zebli\Desktop\大三下\生物信息\数据清点 - 样本对齐 - 共同样本表 - 后续输入规范\window1_outputs\master_samples_6omics.csv`
- sample_presence_6omics.csv: `C:\Users\zebli\Desktop\大三下\生物信息\数据清点 - 样本对齐 - 共同样本表 - 后续输入规范\window1_outputs\sample_presence_6omics.csv`
- preprocess_rules.md: `C:\Users\zebli\Desktop\大三下\生物信息\数据清点 - 样本对齐 - 共同样本表 - 后续输入规范\window1_outputs\preprocess_rules.md`
- window_input_spec.md: `C:\Users\zebli\Desktop\大三下\生物信息\数据清点 - 样本对齐 - 共同样本表 - 后续输入规范\window1_outputs\window_input_spec.md`
- CNV GISTIC thresholded: `C:\Users\zebli\Desktop\大三下\生物信息\数据清点 - 样本对齐 - 共同样本表 - 后续输入规范\outputs_window3\extracted_inputs\原始实验数据\TCGA.STAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes`
- mutation gene-level partial data: `C:\Users\zebli\Desktop\大三下\生物信息\数据清点 - 样本对齐 - 共同样本表 - 后续输入规范\outputs_window3\extracted_inputs\原始实验数据\mc3_gene_level_STAD_mc3_gene_level.txt`
- 2txt innovation notes: `C:\Users\zebli\Desktop\大三下\生物信息\数据清点 - 样本对齐 - 共同样本表 - 后续输入规范\2txt创新点改进.md`

## Archive handling
- extracted gz: C:\Users\zebli\Desktop\大三下\生物信息\数据清点 - 样本对齐 - 共同样本表 - 后续输入规范\原始实验数据\mc3_gene_level_STAD_mc3_gene_level.txt.gz -> C:\Users\zebli\Desktop\大三下\生物信息\数据清点 - 样本对齐 - 共同样本表 - 后续输入规范\outputs_window3\extracted_inputs\原始实验数据\mc3_gene_level_STAD_mc3_gene_level.txt
- extracted gz: C:\Users\zebli\Desktop\大三下\生物信息\数据清点 - 样本对齐 - 共同样本表 - 后续输入规范\原始实验数据\TCGA.STAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz -> C:\Users\zebli\Desktop\大三下\生物信息\数据清点 - 样本对齐 - 共同样本表 - 后续输入规范\outputs_window3\extracted_inputs\原始实验数据\TCGA.STAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes

## Final common samples
- Final common primary tumor samples: 288
- Reason/check: 与预期 288 一致。

## CNV
- Raw matrix after sample alignment: 288 samples x 24776 genes
- Missing values before median fill: 0
- Filtered matrix: 288 samples x 3000 genes
- Feature filter: non-zero variance genes, top 3000 by variance or all available if fewer than 3000
- PCA final dimension: 30
- PCA comparison:
- 20 PCs: cumulative explained variance = 0.8275
- 30 PCs: cumulative explained variance = 0.8854
- 50 PCs: cumulative explained variance = 0.9382

- CNV was not log-transformed because GISTIC thresholded calls are discrete copy-number states, usually -2, -1, 0, 1, 2; a log transform would distort the ordinal event encoding.

## Mutation
- Raw matrix after sample alignment: 288 samples x 40542 genes
- Missing values before binary fill: 0
- All-zero genes: 23889
- Baseline after removing all-zero genes: 288 samples x 16653 genes
- Frequency >= 1% version: 288 samples x 11149 genes
- Final SVD embedding dimension: 8
- Mutation was not z-scored because it is a 0/1 event matrix, not a continuous expression-like abundance matrix; z-score would obscure the binary event meaning.

## Lightweight enhanced version
- Enhanced method: CNV-NMF on selected CNV features after shifting by +2
- Minimum value after shift: 0.0000
- Final NMF rank: 15
- CNV-NMF was prioritized because CNV has genome-wide discrete copy-number features suitable for a complementary nonnegative event-pattern embedding; mutation remains binary/sparse, and pathway aggregation was intentionally left out of this baseline-focused window.

## Recommended inputs for Window 5
- Baseline CNV: `outputs_window3/cnv/cnv_pca_6omics.csv`
- Baseline mutation: `outputs_window3/mutation/mutation_embedding_6omics.csv`
- Optional enhanced CNV: `outputs_window3/enhanced/cnv_nmf_embedding_6omics.csv`

## Direct output files
- `outputs_window3/cnv/cnv_processed_6omics.csv`
- `outputs_window3/cnv/cnv_pca_6omics.csv`
- `outputs_window3/cnv/cnv_feature_variance_top.xlsx`
- `outputs_window3/figures/cnv_pca_explained_variance.png`
- `outputs_window3/mutation/mutation_processed_6omics.csv`
- `outputs_window3/mutation/mutation_processed_freq1pct.csv`
- `outputs_window3/mutation/mutation_embedding_6omics.csv`
- `outputs_window3/mutation/mutation_frequency_summary.xlsx`
- `outputs_window3/enhanced/cnv_nmf_embedding_6omics.csv`
- `outputs_window3/enhanced/cnv_nmf_rank_comparison.xlsx`
