from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List

import pandas as pd
## 记得修改路径
BASE_DIR = r"e:\swxxx"
DATA_DIR = os.path.join(BASE_DIR, "生物信息学实验数据")
OUT_DIR = os.path.join(BASE_DIR, "window1_outputs")
MISSING_DIR = os.path.join(OUT_DIR, "missingness")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(MISSING_DIR, exist_ok=True)


@dataclass
class DatasetConfig:
    key: str
    display_name: str
    omics_type: str
    path: str
    format: str  # matrix_sample_columns | table_sample_rows
    compression: str | None = None
    sample_col: str | None = None
    notes: str = ""


DATASETS: List[DatasetConfig] = [
    DatasetConfig(
        key="rna",
        display_name="TCGA.STAD.sampleMap_HiSeqV2.gz",
        omics_type="RNA",
        path=os.path.join(DATA_DIR, "TCGA.STAD.sampleMap_HiSeqV2.gz"),
        format="matrix_sample_columns",
        compression="gzip",
        notes="RSEM normalized expression; assignment notes indicate log2(x+1) already applied.",
    ),
    DatasetConfig(
        key="mirna",
        display_name="TCGA.STAD.sampleMap_miRNA_HiSeq_gene.gz",
        omics_type="miRNA",
        path=os.path.join(DATA_DIR, "TCGA.STAD.sampleMap_miRNA_HiSeq_gene.gz"),
        format="matrix_sample_columns",
        compression="gzip",
        notes="miRNA expression; assignment notes indicate log2(total_RPM+1) already applied.",
    ),
    DatasetConfig(
        key="cnv",
        display_name="TCGA.STAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz",
        omics_type="CNV",
        path=os.path.join(DATA_DIR, "TCGA.STAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"),
        format="matrix_sample_columns",
        compression="gzip",
        notes="GISTIC thresholded CNV with discrete values -2,-1,0,1,2.",
    ),
    DatasetConfig(
        key="methylation",
        display_name="TCGA.STAD.sampleMap_HumanMethylation450.gz",
        omics_type="methylation",
        path=os.path.join(DATA_DIR, "TCGA.STAD.sampleMap_HumanMethylation450.gz"),
        format="matrix_sample_columns",
        compression="gzip",
        notes="HumanMethylation450 beta values; high-dimensional methylation matrix.",
    ),
    DatasetConfig(
        key="rppa",
        display_name="TCGA.STAD.sampleMap_RPPA.gz",
        omics_type="RPPA",
        path=os.path.join(DATA_DIR, "TCGA.STAD.sampleMap_RPPA.gz"),
        format="matrix_sample_columns",
        compression="gzip",
        notes="Protein expression data; already normalized by source pipeline.",
    ),
    DatasetConfig(
        key="mutation",
        display_name="mc3_gene_level_STAD_mc3_gene_level.txt.gz",
        omics_type="mutation",
        path=os.path.join(DATA_DIR, "mc3_gene_level_STAD_mc3_gene_level.txt.gz"),
        format="matrix_sample_columns",
        compression="gzip",
        notes="Gene-level non-silent mutation matrix with binary 0/1 style encoding.",
    ),
    DatasetConfig(
        key="clinical",
        display_name="TCGA.STAD.sampleMap_STAD_clinicalMatrix",
        omics_type="clinical",
        path=os.path.join(DATA_DIR, "TCGA.STAD.sampleMap_STAD_clinicalMatrix"),
        format="table_sample_rows",
        sample_col="sampleID",
        notes="Clinical matrix; use only for post-clustering validation, not for clustering input.",
    ),
    DatasetConfig(
        key="survival",
        display_name="survival_STAD_survival.txt",
        omics_type="survival",
        path=os.path.join(DATA_DIR, "survival_STAD_survival.txt"),
        format="table_sample_rows",
        sample_col="sample",
        notes="Survival endpoints include OS, DSS, DFI, PFI; use only for post-clustering validation.",
    ),
]


def normalize_sample_barcode(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    parts = text.split("-")
    if len(parts) < 4 or parts[0] != "TCGA":
        return None
    return "-".join(parts[:3] + [parts[3][:2]])


def patient_barcode_from_sample(sample_barcode: str | None) -> str | None:
    if not sample_barcode:
        return None
    parts = sample_barcode.split("-")
    if len(parts) < 3:
        return None
    return "-".join(parts[:3])


def sample_type_from_sample(sample_barcode: str | None) -> str | None:
    if not sample_barcode:
        return None
    parts = sample_barcode.split("-")
    if len(parts) < 4:
        return None
    return parts[3][:2]


summary_rows: List[Dict[str, object]] = []
sample_filter_rows: List[Dict[str, object]] = []
overlap_rows: List[Dict[str, object]] = []
dataset_samples: Dict[str, set[str]] = {}
dataset_patients: Dict[str, set[str]] = {}


def process_matrix_sample_columns(cfg: DatasetConfig) -> None:
    print(f"Processing matrix dataset: {cfg.key}")
    chunk_iter = pd.read_csv(
        cfg.path,
        sep="\t",
        compression=cfg.compression,
        chunksize=2000,
        low_memory=False,
    )

    n_features = 0
    total_missing = 0
    total_cells = 0
    tumor_sample_ids: List[str] | None = None
    tumor_patient_ids: List[str] | None = None
    tumor_mask = None
    raw_sample_ids: List[str] | None = None
    feature_axis = None
    sample_missing_counts = None
    feature_missing_path = os.path.join(MISSING_DIR, f"{cfg.key}_feature_missing.csv")
    if os.path.exists(feature_missing_path):
        os.remove(feature_missing_path)

    for chunk_idx, chunk in enumerate(chunk_iter):
        if chunk_idx == 0:
            feature_axis = chunk.columns[0]
            raw_sample_ids = [str(col) for col in chunk.columns[1:]]
            normalized = [normalize_sample_barcode(col) for col in raw_sample_ids]
            tumor_mask = [sample_type_from_sample(sample) == "01" for sample in normalized]
            tumor_sample_ids = [sample for sample, keep in zip(normalized, tumor_mask) if keep and sample]
            tumor_patient_ids = [patient_barcode_from_sample(sample) for sample in tumor_sample_ids]
            sample_missing_counts = pd.Series(0, index=tumor_sample_ids, dtype="int64")

        feature_ids = chunk.iloc[:, 0].astype(str)
        data = chunk.iloc[:, 1:]
        tumor_data = data.loc[:, tumor_mask]
        missing_mask = tumor_data.isna()

        n_features += len(chunk)
        total_missing += int(missing_mask.to_numpy().sum())
        total_cells += tumor_data.shape[0] * tumor_data.shape[1]
        sample_missing_counts = sample_missing_counts.add(missing_mask.sum(axis=0), fill_value=0).astype("int64")

        feature_missing_df = pd.DataFrame(
            {
                "feature_id": feature_ids,
                "missing_count": missing_mask.sum(axis=1).astype("int64").to_numpy(),
                "missing_rate": missing_mask.mean(axis=1).to_numpy(),
            }
        )
        feature_missing_df.to_csv(
            feature_missing_path,
            mode="a",
            index=False,
            header=(chunk_idx == 0),
        )

    assert raw_sample_ids is not None and tumor_sample_ids is not None and tumor_patient_ids is not None and sample_missing_counts is not None

    sample_missing_df = pd.DataFrame(
        {
            "sample_barcode": tumor_sample_ids,
            "patient_barcode": tumor_patient_ids,
            "missing_count": sample_missing_counts.values,
            "missing_rate": sample_missing_counts.values / max(n_features, 1),
        }
    )
    sample_missing_df.to_csv(os.path.join(MISSING_DIR, f"{cfg.key}_sample_missing.csv"), index=False)

    unique_tumor_samples = list(dict.fromkeys(tumor_sample_ids))
    unique_tumor_patients = list(dict.fromkeys([pid for pid in tumor_patient_ids if pid]))
    dataset_samples[cfg.key] = set(unique_tumor_samples)
    dataset_patients[cfg.key] = set(unique_tumor_patients)

    summary_rows.append(
        {
            "dataset_name": cfg.display_name,
            "omics_type": cfg.omics_type,
            "n_samples_raw": len(raw_sample_ids),
            "n_features_raw": n_features,
            "sample_axis": "columns",
            "feature_axis": feature_axis,
            "sample_id_example": raw_sample_ids[0] if raw_sample_ids else None,
            "missing_rate_overall": total_missing / total_cells if total_cells else 0,
            "notes": cfg.notes,
        }
    )
    sample_filter_rows.append(
        {
            "dataset_name": cfg.display_name,
            "omics_type": cfg.omics_type,
            "raw_sample_entries": len(raw_sample_ids),
            "tumor_sample_entries": len(tumor_sample_ids),
            "normal_or_non_primary_removed": len(raw_sample_ids) - len(tumor_sample_ids),
            "unique_tumor_samples": len(unique_tumor_samples),
            "unique_tumor_patients": len(unique_tumor_patients),
            "feature_count": n_features,
        }
    )


def process_table_sample_rows(cfg: DatasetConfig) -> None:
    print(f"Processing row-oriented table: {cfg.key}")
    df = pd.read_csv(cfg.path, sep="\t", low_memory=False)
    sample_col = cfg.sample_col
    assert sample_col is not None

    df["sample_barcode"] = df[sample_col].map(normalize_sample_barcode)
    df["patient_barcode"] = df["sample_barcode"].map(patient_barcode_from_sample)
    df["sample_type_code"] = df["sample_barcode"].map(sample_type_from_sample)

    tumor_df = df[df["sample_type_code"] == "01"].copy()
    dataset_samples[cfg.key] = set(tumor_df["sample_barcode"].dropna().astype(str).unique())
    dataset_patients[cfg.key] = set(tumor_df["patient_barcode"].dropna().astype(str).unique())

    exclude_cols = {sample_col, "sample_barcode", "patient_barcode", "sample_type_code"}
    data_cols = [col for col in tumor_df.columns if col not in exclude_cols]
    feature_missing = tumor_df[data_cols].isna().mean(axis=0)
    sample_missing = tumor_df[data_cols].isna().mean(axis=1)

    pd.DataFrame(
        {
            "feature_id": feature_missing.index,
            "missing_rate": feature_missing.values,
            "missing_count": tumor_df[data_cols].isna().sum(axis=0).values,
        }
    ).to_csv(os.path.join(MISSING_DIR, f"{cfg.key}_feature_missing.csv"), index=False)

    pd.DataFrame(
        {
            "sample_barcode": tumor_df["sample_barcode"],
            "patient_barcode": tumor_df["patient_barcode"],
            "missing_rate": sample_missing.values,
            "missing_count": tumor_df[data_cols].isna().sum(axis=1).values,
        }
    ).to_csv(os.path.join(MISSING_DIR, f"{cfg.key}_sample_missing.csv"), index=False)

    summary_rows.append(
        {
            "dataset_name": cfg.display_name,
            "omics_type": cfg.omics_type,
            "n_samples_raw": int(df.shape[0]),
            "n_features_raw": int(len(data_cols)),
            "sample_axis": "rows",
            "feature_axis": "columns",
            "sample_id_example": df[sample_col].iloc[0] if not df.empty else None,
            "missing_rate_overall": float(tumor_df[data_cols].isna().mean().mean()) if data_cols else 0,
            "notes": cfg.notes,
        }
    )
    sample_filter_rows.append(
        {
            "dataset_name": cfg.display_name,
            "omics_type": cfg.omics_type,
            "raw_sample_entries": int(df.shape[0]),
            "tumor_sample_entries": int(tumor_df.shape[0]),
            "normal_or_non_primary_removed": int(df.shape[0] - tumor_df.shape[0]),
            "unique_tumor_samples": int(tumor_df["sample_barcode"].nunique()),
            "unique_tumor_patients": int(tumor_df["patient_barcode"].nunique()),
            "feature_count": int(len(data_cols)),
        }
    )


for dataset in DATASETS:
    if dataset.format == "matrix_sample_columns":
        process_matrix_sample_columns(dataset)
    else:
        process_table_sample_rows(dataset)


main_omics = ["rna", "mirna", "cnv", "rppa"]
all_omics = ["rna", "mirna", "cnv", "methylation", "rppa", "mutation"]

common_6 = sorted(set.intersection(*(dataset_samples[key] for key in all_omics)))
common_4 = sorted(set.intersection(*(dataset_samples[key] for key in main_omics)))
clinical_overlap_6 = sorted(set(common_6) & dataset_samples["clinical"])
survival_overlap_6 = sorted(set(common_6) & dataset_samples["survival"])
clinical_overlap_4 = sorted(set(common_4) & dataset_samples["clinical"])
survival_overlap_4 = sorted(set(common_4) & dataset_samples["survival"])

master_6_df = pd.DataFrame(
    {
        "sample_barcode": common_6,
        "patient_barcode": [patient_barcode_from_sample(x) for x in common_6],
    }
)
master_4_df = pd.DataFrame(
    {
        "sample_barcode": common_4,
        "patient_barcode": [patient_barcode_from_sample(x) for x in common_4],
    }
)

master_6_df.to_csv(os.path.join(OUT_DIR, "master_samples_6omics.csv"), index=False)
master_4_df.to_csv(os.path.join(OUT_DIR, "master_samples_4omics.csv"), index=False)

sample_presence_df = master_6_df.copy()
for key in all_omics + ["clinical", "survival"]:
    sample_presence_df[f"present_in_{key}"] = sample_presence_df["sample_barcode"].isin(dataset_samples[key])

clinical_survival_overlap_df = pd.DataFrame(
    {
        "analysis_scope": ["6omics", "6omics", "4omics", "4omics"],
        "dataset": ["clinical", "survival", "clinical", "survival"],
        "overlap_sample_count": [
            len(clinical_overlap_6),
            len(survival_overlap_6),
            len(clinical_overlap_4),
            len(survival_overlap_4),
        ],
    }
)
clinical_survival_overlap_df.to_csv(os.path.join(OUT_DIR, "clinical_survival_overlap.csv"), index=False)
sample_presence_df.to_csv(os.path.join(OUT_DIR, "sample_presence_6omics.csv"), index=False)

for key in DATASETS:
    overlap_rows.append(
        {
            "dataset_name": key.key,
            "tumor_sample_count": len(dataset_samples[key.key]),
            "tumor_patient_count": len(dataset_patients[key.key]),
        }
    )

overlap_rows.extend(
    [
        {"dataset_name": "common_6omics", "tumor_sample_count": len(common_6), "tumor_patient_count": len({patient_barcode_from_sample(x) for x in common_6})},
        {"dataset_name": "common_4omics", "tumor_sample_count": len(common_4), "tumor_patient_count": len({patient_barcode_from_sample(x) for x in common_4})},
        {"dataset_name": "clinical_overlap_with_6omics", "tumor_sample_count": len(clinical_overlap_6), "tumor_patient_count": len({patient_barcode_from_sample(x) for x in clinical_overlap_6})},
        {"dataset_name": "survival_overlap_with_6omics", "tumor_sample_count": len(survival_overlap_6), "tumor_patient_count": len({patient_barcode_from_sample(x) for x in survival_overlap_6})},
        {"dataset_name": "clinical_overlap_with_4omics", "tumor_sample_count": len(clinical_overlap_4), "tumor_patient_count": len({patient_barcode_from_sample(x) for x in clinical_overlap_4})},
        {"dataset_name": "survival_overlap_with_4omics", "tumor_sample_count": len(survival_overlap_4), "tumor_patient_count": len({patient_barcode_from_sample(x) for x in survival_overlap_4})},
    ]
)

summary_df = pd.DataFrame(summary_rows)
filter_df = pd.DataFrame(sample_filter_rows)
overlap_df = pd.DataFrame(overlap_rows)

with pd.ExcelWriter(os.path.join(OUT_DIR, "dataset_qc_summary.xlsx"), engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="dataset_summary", index=False)
    filter_df.to_excel(writer, sheet_name="sample_filtering", index=False)
    overlap_df.to_excel(writer, sheet_name="overlap_summary", index=False)

preprocess_rules = f"""# 预处理统一规则

## 为什么先统一共同样本，再处理缺失
多组学整合的核心是同一批病人的不同组学观测，而不是不同组学的共同特征名。不同组学的特征空间天然不同：RNA 是基因表达，miRNA 是 miRNA 表达，CNV 是拷贝数，甲基化是探针 beta 值，RPPA 是蛋白，突变是 0/1 事件。先对齐共同样本，才能保证后续缺失处理、标准化、降维和聚类都在同一批病人上进行，避免把不属于共同分析集的样本噪声提前写进结果。

## 全组统一规则
1. 只保留原发肿瘤样本 `-01`，去掉正常样本 `-11`。
2. 先统一共同样本，再删高缺失特征，再填剩余少量缺失，最后标准化。
3. `clinical` 和 `survival` 只能做后验验证，不能参与聚类。
4. 所有后续矩阵统一格式为：行 = 样本，列 = 特征或主成分。
5. 所有人固定随机种子，并保留中间结果文件。

## 各组学预处理建议
### RNA
- 已做 `log2(x+1)`，不要重复 log。
- 先做低表达 / 低方差过滤。
- 再做 z-score 标准化。
- 后续建议 PCA。

### miRNA
- 已做 `log2(total_RPM+1)`，不要重复 log。
- 先删高缺失 miRNA，再对剩余少量缺失填补。
- 再做 z-score。
- 后续建议 PCA。

### RPPA
- 缺失通常较少，可中位数填充少量缺失。
- 再做 z-score。
- 后续建议 PCA。

### CNV
- GISTIC 阈值数据，取值一般为 `-2,-1,0,1,2`。
- 不做 log。
- 先做低方差过滤。
- 后续建议 PCA。

### mutation
- 0/1 二值矩阵，不按连续表达矩阵处理。
- 先去全 0 和极低频突变基因。
- 一般不做 z-score。
- 若维度仍高，可考虑 TruncatedSVD。

### methylation
- beta 值范围通常在 `0~1`。
- 先删高缺失探针。
- 再按高方差 / MAD / Top N 进行筛选。
- 后续再降维，不建议全量直接拼接进聚类。

## 主分析与扩展分析
- 主分析共同样本文件：`master_samples_4omics.csv`，用于 RNA + miRNA + CNV + RPPA。
- 扩展分析共同样本文件：`master_samples_6omics.csv`，用于六组学可替换流程或敏感性分析。
"""

with open(os.path.join(OUT_DIR, "preprocess_rules.md"), "w", encoding="utf-8") as f:
    f.write(preprocess_rules)

window_input_spec = f"""# 后续窗口统一输入说明

## 统一样本 ID 列
- `sample_barcode`: 标准样本 ID，格式示例 `TCGA-3M-AB46-01`
- `patient_barcode`: 病人级 ID，格式示例 `TCGA-3M-AB46`

## 主分析默认输入
- 主分析共同样本：`{os.path.join(OUT_DIR, 'master_samples_4omics.csv')}`
- 适用窗口：RNA/miRNA/RPPA、CNV、整合聚类、临床/生存验证、监督模型。

## 扩展分析输入
- 六组学共同样本：`{os.path.join(OUT_DIR, 'master_samples_6omics.csv')}`
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
"""

with open(os.path.join(OUT_DIR, "window_input_spec.md"), "w", encoding="utf-8") as f:
    f.write(window_input_spec)

print("Done.")
print(f"6-omics common tumor samples: {len(common_6)}")
print(f"4-omics common tumor samples: {len(common_4)}")
print(f"Clinical overlap with 6-omics: {len(clinical_overlap_6)}")
print(f"Survival overlap with 6-omics: {len(survival_overlap_6)}")
print(f"Clinical overlap with 4-omics: {len(clinical_overlap_4)}")
print(f"Survival overlap with 4-omics: {len(survival_overlap_4)}")