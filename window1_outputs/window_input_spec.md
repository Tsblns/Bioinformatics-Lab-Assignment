# 后续窗口统一输入说明

## 统一样本 ID 列
- `sample_barcode`: 标准样本 ID，格式示例 `TCGA-3M-AB46-01`
- `patient_barcode`: 病人级 ID，格式示例 `TCGA-3M-AB46`

## 主分析默认输入
- 主分析共同样本：`e:\swxxx\window1_outputs\master_samples_4omics.csv`
- 适用窗口：RNA/miRNA/RPPA、CNV、整合聚类、临床/生存验证、监督模型。

## 扩展分析输入
- 六组学共同样本：`e:\swxxx\window1_outputs\master_samples_6omics.csv`
- 适用窗口：methylation、mutation 扩展流程，或六组学敏感性分析。

## 输出矩阵规范
1. 所有清洗后矩阵统一为“行 = 样本，列 = 特征或主成分”。
2. 第一列固定为 `sample_barcode`；第二列推荐保留 `patient_barcode`。
3. 聚类窗口只接收各组学清洗后矩阵或降维表示，不允许各窗口自行重新取交集。
4. 中间结果统一保存到各自窗口输出目录，并同步记录随机种子和参数。

## 文件命名建议
- `rna_processed.csv`
- `mirna_processed.csv`
- `rppa_processed.csv`
- `cnv_processed.csv`
- `mutation_processed.csv`
- `methylation_processed.csv`
- `cluster_labels.csv`
- `clinical_cluster_association.xlsx`
- `classifier_metrics.xlsx`
