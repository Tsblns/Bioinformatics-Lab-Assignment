# 窗口3：CNV + mutation 预处理说明

生成时间：2026-05-02  
主脚本：`run_window3.py`  
输出目录：`outputs_window3/`

本窗口的任务是接着窗口1已经完成的“数据清点、样本对齐、共同样本表、输入规范”继续做 CNV 和 mutation 的预处理，形成后续窗口5可以直接使用的输入矩阵。本窗口没有重做窗口1，没有做聚类，也没有做 clinical / survival 分析。

## 一、实际使用的输入文件

本窗口只使用窗口1给出的共同样本表进行样本限制和排序。

| 类型 | 实际使用文件 |
|---|---|
| 共同样本表 | `window1_outputs/master_samples_6omics.csv` |
| 六组学样本存在表 | `window1_outputs/sample_presence_6omics.csv` |
| 预处理规则 | `window1_outputs/preprocess_rules.md` |
| 后续窗口输入规范 | `window1_outputs/window_input_spec.md` |
| CNV 原始文件 | `outputs_window3/extracted_inputs/原始实验数据/TCGA.STAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes` |
| mutation 原始文件 | `outputs_window3/extracted_inputs/原始实验数据/mc3_gene_level_STAD_mc3_gene_level.txt` |
| 创新点说明 | `2txt创新点改进.md` |

原始 CNV 和 mutation 文件来自 `原始实验数据/` 中的 `.gz` 压缩文件，脚本已自动解压到 `outputs_window3/extracted_inputs/` 后再读取。

## 二、共同样本处理

本窗口严格使用 `master_samples_6omics.csv` 作为共同样本依据。

- 最终共同样本数：`288`
- 样本类型：TCGA 原发肿瘤样本，样本条码末尾为 `-01`
- 输出矩阵统一格式：行 = 样本，列 = 特征或低维成分
- 输出矩阵前两列统一保留：`sample_barcode`、`patient_barcode`

本次最终样本数正好为 288，符合项目预期。

## 三、CNV 预处理方法和步骤

### 1. 读取和方向识别

CNV 原始文件是 GISTIC thresholded by genes 数据，原始结构为“行 = 基因，列 = 样本”，数值主要为 `-2, -1, 0, 1, 2`。脚本自动识别样本列和基因行，并转换为后续统一格式：“行 = 样本，列 = 基因”。

### 2. 样本对齐

只保留 `master_samples_6omics.csv` 中的 288 个共同样本，并按 master 表中的样本顺序重新排序。

### 3. 基础检查

本次 CNV 对齐后的基础情况：

- 样本数：`288`
- 原始基因数：`24776`
- 缺失值数：`0`
- 数据类型：GISTIC 离散阈值 CNV

CNV 没有做 log 转换。原因是 GISTIC thresholded 数据本身不是连续表达量，而是离散拷贝数状态；对 `-2, -1, 0, 1, 2` 做 log 会破坏其事件含义。

### 4. 特征筛选

使用默认 baseline 策略：计算每个基因在 288 个共同样本中的方差，去掉零方差基因，按方差从高到低排序，保留 Top 3000 个 CNV 基因。如果可用基因少于 3000，则保留全部非零方差基因。

本次结果：

- 筛选后 CNV 特征数：`3000`
- CNV 主特征矩阵：`288 x 3000`

### 5. 标准化和 PCA

筛选后的 CNV 特征先做标准化，再做 PCA。这里的标准化只用于 PCA 的数值尺度平衡，不改变原始 CNV baseline 矩阵的离散含义。

PCA 尝试了 20、30、50 维：

| PCA 维度 | 累计解释率 |
|---:|---:|
| 20 | 0.8275 |
| 30 | 0.8854 |
| 50 | 0.9382 |

本窗口默认采用 30 维作为 CNV 主输出，因为 30 维可行，并且能在较紧凑维度下保留约 88.54% 的累计解释率。

### 6. CNV 输出文件

| 文件 | 内容 | 用途 |
|---|---|---|
| `outputs_window3/cnv/cnv_processed_6omics.csv` | CNV 方差筛选后的 3000 基因矩阵，含样本 ID | CNV baseline 特征备查或备用输入 |
| `outputs_window3/cnv/cnv_pca_6omics.csv` | CNV PCA 30 维矩阵，含样本 ID | 推荐给窗口5作为 baseline CNV 输入 |
| `outputs_window3/cnv/cnv_feature_variance_top.xlsx` | CNV 基因方差排序和筛选记录 | 记录特征选择依据 |
| `outputs_window3/figures/cnv_pca_explained_variance.png` | PCA 累计解释率图 | 说明 PCA 维度选择 |
| `outputs_window3/logs/cnv_value_distribution.csv` | CNV 数值分布 | 检查是否符合 GISTIC thresholded 数据特征 |
| `outputs_window3/logs/cnv_pca_explained_variance.csv` | PCA 维度比较表 | 记录 20、30、50 维解释率 |

## 四、mutation 预处理方法和步骤

### 1. 读取和方向识别

mutation 原始文件为基因层面的突变矩阵，脚本识别后转换为“行 = 样本，列 = 基因，数值 = 0/1 突变事件”。当前 mutation 文件虽然是“部分数据”版本，但本次实际读入后仍包含较多基因。脚本按实际文件内容处理，不硬编码维度。

### 2. 样本对齐

只保留 `master_samples_6omics.csv` 中的 288 个共同样本，并按 master 表顺序重排。

### 3. 基础检查

本次 mutation 对齐后的基础情况：

- 样本数：`288`
- 原始基因数：`40542`
- 缺失值数：`0`
- 全 0 基因数：`23889`

mutation 没有做 z-score。原因是 mutation 是 0/1 二值事件矩阵，不是连续表达矩阵；z-score 会让二值突变事件失去清晰的“是否突变”含义。

### 4. baseline 策略

baseline 处理策略：将 mutation 矩阵统一为 0/1，计算每个基因的突变频率，去掉所有共同样本中全 0 的基因，保留剩余全部非零突变基因。

本次 baseline 结果：

- 去全 0 后 mutation 基因数：`16653`
- baseline mutation 矩阵：`288 x 16653`

### 5. 频率 1% 紧凑版本

另外生成一个较紧凑版本用于对比：筛选条件为突变频率 `>= 1%`，结果基因数为 `11149`。该文件只是对比输出，不替代 baseline 主输出。

### 6. TruncatedSVD / SVD embedding

在 baseline mutation 矩阵上尝试低维 embedding：尝试 5 维、8 维，本次选择 8 维。由于 mutation 是稀疏 0/1 矩阵，使用 SVD embedding 可以给后续窗口提供更紧凑的 mutation 表示，同时避免把二值矩阵错误当成连续表达量做 z-score。

### 7. mutation 输出文件

| 文件 | 内容 | 用途 |
|---|---|---|
| `outputs_window3/mutation/mutation_processed_6omics.csv` | 去全 0 后的 baseline mutation 0/1 矩阵，含样本 ID | mutation baseline 原始特征备查或备用输入 |
| `outputs_window3/mutation/mutation_processed_freq1pct.csv` | 突变频率 >= 1% 的紧凑 mutation 矩阵，含样本 ID | 对比分析备用 |
| `outputs_window3/mutation/mutation_embedding_6omics.csv` | mutation SVD 8 维 embedding，含样本 ID | 推荐给窗口5作为 baseline mutation 输入 |
| `outputs_window3/mutation/mutation_frequency_summary.xlsx` | 每个基因突变频率、突变样本数、是否保留 | mutation 特征筛选记录 |

## 五、轻量增强版：CNV-NMF

根据任务要求，本窗口优先做 CNV 的轻量增强版，而不是 mutation 通路聚合。

### 1. 为什么选择 CNV-NMF

选择 CNV-NMF 的原因：CNV 是全基因组范围的离散拷贝数事件，适合提取拷贝数改变模式；mutation 当前不适合在本窗口做复杂通路聚合，因为这会引入额外注释库和更强的生物学假设；本窗口目标是先完成稳健 baseline，再提供一个轻量增强输入，因此 CNV-NMF 更符合任务边界。

### 2. NMF 输入构建

NMF 要求输入非负矩阵，但原始 CNV 取值包含负值。因此使用以下方式平移：

```text
X_nonnegative = X_cnv_selected + 2
```

检查结果：平移后最小值为 `0.0000`，满足 NMF 非负输入要求。

### 3. rank 比较

| rank | 迭代次数 | reconstruction error | relative error |
|---:|---:|---:|---:|
| 10 | 300 | 400.2218 | 0.1902 |
| 15 | 300 | 342.0692 | 0.1626 |
| 20 | 300 | 317.0485 | 0.1507 |

虽然 rank 20 的重构误差更低，但默认选择 rank 15，原因是 rank 15 在压缩程度和模式表达能力之间更平衡，也符合任务中“默认优先保留 15 维版本”的要求。

### 4. 增强版输出文件

| 文件 | 内容 | 用途 |
|---|---|---|
| `outputs_window3/enhanced/cnv_nmf_embedding_6omics.csv` | CNV-NMF rank 15 样本低维表示，含样本 ID | 可选增强输入，推荐窗口5作为 optional enhanced |
| `outputs_window3/enhanced/cnv_nmf_rank_comparison.xlsx` | NMF rank 10、15、20 的重构误差比较 | 记录增强版参数选择依据 |

## 六、提供给后续窗口继续做的文件

### 1. 推荐 baseline 输入

| 组学 | 推荐文件 | 维度说明 |
|---|---|---|
| CNV | `outputs_window3/cnv/cnv_pca_6omics.csv` | 288 个样本，30 个 CNV PC，另含 2 个 ID 列 |
| mutation | `outputs_window3/mutation/mutation_embedding_6omics.csv` | 288 个样本，8 个 mutation SVD 成分，另含 2 个 ID 列 |

### 2. 可选增强输入

| 类型 | 推荐文件 | 维度说明 |
|---|---|---|
| CNV-NMF | `outputs_window3/enhanced/cnv_nmf_embedding_6omics.csv` | 288 个样本，15 个 CNV-NMF 成分，另含 2 个 ID 列 |

### 3. 备用或审计文件

| 文件 | 说明 |
|---|---|
| `outputs_window3/cnv/cnv_processed_6omics.csv` | CNV 方差筛选后的 3000 基因矩阵 |
| `outputs_window3/mutation/mutation_processed_6omics.csv` | mutation 去全 0 后的 0/1 baseline 矩阵 |
| `outputs_window3/mutation/mutation_processed_freq1pct.csv` | mutation 突变频率 >= 1% 的紧凑矩阵 |
| `outputs_window3/cnv/cnv_feature_variance_top.xlsx` | CNV 特征筛选记录 |
| `outputs_window3/mutation/mutation_frequency_summary.xlsx` | mutation 频率统计记录 |
| `outputs_window3/figures/cnv_pca_explained_variance.png` | CNV PCA 解释率图 |

## 七、本窗口的主要贡献

本窗口完成了以下贡献：

1. 自动定位并读取窗口1输出和原始 CNV / mutation 文件。
2. 自动解压 CNV 和 mutation 的 `.gz` 文件，并保留解压后的可追溯输入。
3. 严格使用 `master_samples_6omics.csv`，没有重新取交集，也没有重做窗口1。
4. 将 CNV 和 mutation 都整理为统一格式：行 = 样本，列 = 特征。
5. 保证样本顺序与 master 表一致，便于后续窗口直接拼接。
6. 对 CNV 完成方差筛选、标准化、PCA，并输出推荐的 30 维表示。
7. 对 mutation 保留 0/1 事件含义，去掉全 0 基因，输出 baseline 矩阵和 8 维 SVD embedding。
8. 额外完成 CNV-NMF 轻量增强版，给后续窗口提供可选输入。
9. 生成筛选记录、解释率图、rank 比较表和日志，保证处理过程可复核。

结论：窗口3“CNV + mutation 预处理”内容已经完成，产物可以交给后续窗口继续使用。
