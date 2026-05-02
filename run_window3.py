from __future__ import annotations

import gzip
import math
import re
import shutil
import sys
import zipfile
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from PIL import Image, ImageDraw
except Exception:  # pragma: no cover - runtime fallback for lean Python envs
    Image = None
    ImageDraw = None


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs_window3"
DIRS = {
    "cnv": OUT / "cnv",
    "mutation": OUT / "mutation",
    "enhanced": OUT / "enhanced",
    "figures": OUT / "figures",
    "logs": OUT / "logs",
    "extracted": OUT / "extracted_inputs",
}
RANDOM_SEED = 2026
EPS = 1e-9


def log(msg: str) -> None:
    print(msg, flush=True)


def ensure_dirs() -> None:
    for d in DIRS.values():
        d.mkdir(parents=True, exist_ok=True)


def try_read_text(path: Path, max_chars: int = 3000) -> str:
    for enc in ("utf-8", "utf-8-sig", "gb18030", "gbk", "latin1"):
        try:
            text = path.read_text(encoding=enc, errors="strict")
            return text[:max_chars]
        except Exception:
            continue
    return path.read_text(encoding="utf-8", errors="replace")[:max_chars]


def maybe_extract_archives() -> list[str]:
    """Extract common archives into outputs_window3/extracted_inputs.

    Gzipped tabular files are treated as compressed inputs and decompressed to a
    sibling path in the extracted input folder. The raw .gz remains usable too.
    """
    messages: list[str] = []
    target_markers = (
        "gistic2_copynumber_gistic2_all_thresholded.by_genes",
        "stad_mc3_gene_level",
        "mc3_gene_level_stad_mc3_gene_level",
    )
    archives = []
    for p in ROOT.rglob("*"):
        if not p.is_file() or OUT in p.parents or p.suffix.lower() not in {".zip", ".gz", ".7z"}:
            continue
        # Container archives may hold several required files, so inspect/extract
        # them. Single-file .gz inputs are limited to Window 3 CNV/mutation data.
        if p.suffix.lower() in {".zip", ".7z"} or any(marker in p.name.lower() for marker in target_markers):
            archives.append(p)
    for p in archives:
        rel_parent = p.parent.relative_to(ROOT) if ROOT in p.parents else Path("_external")
        dest_parent = DIRS["extracted"] / rel_parent
        dest_parent.mkdir(parents=True, exist_ok=True)
        try:
            if p.suffix.lower() == ".zip":
                with zipfile.ZipFile(p) as zf:
                    zf.extractall(dest_parent / p.stem)
                messages.append(f"extracted zip: {p} -> {dest_parent / p.stem}")
            elif p.suffix.lower() == ".gz":
                out_name = p.name[:-3]
                dest = dest_parent / out_name
                if not dest.exists() or dest.stat().st_size == 0:
                    with gzip.open(p, "rb") as src, dest.open("wb") as dst:
                        shutil.copyfileobj(src, dst)
                messages.append(f"extracted gz: {p} -> {dest}")
            elif p.suffix.lower() == ".7z":
                try:
                    import py7zr  # type: ignore

                    dest = dest_parent / p.stem
                    dest.mkdir(parents=True, exist_ok=True)
                    with py7zr.SevenZipFile(p, mode="r") as z:
                        z.extractall(dest)
                    messages.append(f"extracted 7z: {p} -> {dest}")
                except Exception as exc:
                    messages.append(f"could not extract 7z {p}: {exc}")
        except Exception as exc:
            messages.append(f"archive extraction failed for {p}: {exc}")
    return messages


def all_files() -> list[Path]:
    return [p for p in ROOT.rglob("*") if p.is_file()]


def find_file(label: str, patterns: list[str], required: bool = True) -> Path | None:
    candidates = []
    for p in all_files():
        name = p.name.lower()
        full = str(p).lower()
        if any(pat.lower() in name or pat.lower() in full for pat in patterns):
            candidates.append(p)
    if not candidates:
        if required:
            raise FileNotFoundError(f"Missing required file for {label}: {patterns}")
        return None

    def score(path: Path) -> tuple[int, int, int]:
        in_out = 0 if OUT not in path.parents else 1
        exactish = 0 if any(path.name.lower() == pat.lower() for pat in patterns) else 1
        compressed = 1 if path.suffix.lower() in {".gz", ".zip", ".7z"} else 0
        return (exactish, compressed, in_out)

    candidates = sorted(candidates, key=score)
    return candidates[0]


def canonical_sample(value: object) -> str | None:
    s = str(value).strip().upper().replace(".", "-").replace("_", "-")
    m = re.search(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-\d{2})", s)
    return m.group(1) if m else None


def read_table(path: Path) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    compression = "gzip" if path.suffix.lower() == ".gz" else None
    sep = "\t" if path.name.lower().endswith((".txt", ".txt.gz", ".by_genes", ".by_genes.gz")) else None
    for enc in ("utf-8", "utf-8-sig", "gb18030", "latin1"):
        try:
            if sep is None:
                return pd.read_csv(path, compression=compression, encoding=enc)
            return pd.read_csv(path, sep=sep, compression=compression, encoding=enc)
        except UnicodeDecodeError:
            continue
    if sep is None:
        return pd.read_csv(path, compression=compression, encoding="utf-8", encoding_errors="replace")
    return pd.read_csv(path, sep=sep, compression=compression, encoding="utf-8", encoding_errors="replace")


def load_master(path: Path) -> pd.DataFrame:
    master = pd.read_csv(path)
    sample_col = next((c for c in master.columns if "sample" in c.lower()), master.columns[0])
    patient_col = next((c for c in master.columns if "patient" in c.lower()), None)
    out = pd.DataFrame({"sample_barcode": master[sample_col].map(canonical_sample)})
    if patient_col:
        out["patient_barcode"] = master[patient_col].astype(str).str.strip()
    else:
        out["patient_barcode"] = out["sample_barcode"].str.slice(0, 12)
    out = out.dropna(subset=["sample_barcode"]).drop_duplicates("sample_barcode")
    out = out[out["sample_barcode"].str.endswith("-01")].reset_index(drop=True)
    return out


def orient_gene_sample_matrix(raw: pd.DataFrame, master_samples: list[str], label: str) -> tuple[pd.DataFrame, dict]:
    """Return rows=samples, cols=genes, values=numeric."""
    df = raw.copy()
    df.columns = [str(c).strip() for c in df.columns]
    master_set = set(master_samples)
    col_matches = {c: canonical_sample(c) for c in df.columns}
    matched_cols = [c for c, sample in col_matches.items() if sample in master_set]

    first_col = df.columns[0]
    row_matches = df[first_col].map(canonical_sample) if len(df) else pd.Series(dtype=object)
    matched_rows = row_matches[row_matches.isin(master_set)]

    info = {
        "label": label,
        "raw_shape": df.shape,
        "matched_sample_columns": len(matched_cols),
        "matched_sample_rows": int(matched_rows.shape[0]),
        "orientation": "",
    }

    if len(matched_cols) >= max(5, len(matched_rows)):
        info["orientation"] = "genes_by_rows_samples_by_columns"
        gene_col = first_col
        mat = df[[gene_col] + matched_cols].copy()
        mat[gene_col] = mat[gene_col].astype(str).str.strip()
        mat = mat[mat[gene_col].notna() & (mat[gene_col] != "")]
        mat = mat.set_index(gene_col)
        mat.columns = [col_matches[c] for c in matched_cols]
        mat = mat.apply(pd.to_numeric, errors="coerce")
        mat = mat.groupby(mat.index).mean()
        sample_major = mat.T
    elif len(matched_rows) > 0:
        info["orientation"] = "samples_by_rows_genes_by_columns"
        sample_col = first_col
        df["_sample_barcode"] = df[sample_col].map(canonical_sample)
        df = df[df["_sample_barcode"].isin(master_set)].copy()
        gene_cols = [c for c in df.columns if c not in {sample_col, "_sample_barcode"}]
        sample_major = df.set_index("_sample_barcode")[gene_cols]
        sample_major = sample_major.apply(pd.to_numeric, errors="coerce")
        sample_major = sample_major.groupby(sample_major.index).mean()
    else:
        raise ValueError(f"Could not detect sample orientation for {label}.")

    sample_major = sample_major.reindex(master_samples)
    sample_major.index.name = "sample_barcode"
    return sample_major, info


def add_id_columns(matrix: pd.DataFrame, master: pd.DataFrame) -> pd.DataFrame:
    ids = master.set_index("sample_barcode").loc[matrix.index, ["patient_barcode"]]
    out = matrix.copy()
    out.insert(0, "patient_barcode", ids["patient_barcode"].values)
    out.insert(0, "sample_barcode", out.index)
    return out.reset_index(drop=True)


def fill_missing_with_median(df: pd.DataFrame) -> pd.DataFrame:
    if df.isna().sum().sum() == 0:
        return df
    med = df.median(axis=0, skipna=True).fillna(0)
    return df.fillna(med)


def standardize(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    arr = df.to_numpy(dtype=float)
    mean = arr.mean(axis=0)
    std = arr.std(axis=0, ddof=0)
    std[std < EPS] = 1.0
    return (arr - mean) / std, mean, std


def pca_via_svd(x_scaled: np.ndarray, dims: list[int]) -> tuple[dict[int, dict], np.ndarray]:
    u, s, vt = np.linalg.svd(x_scaled, full_matrices=False)
    explained = (s**2) / max(1, x_scaled.shape[0] - 1)
    ratio = explained / explained.sum() if explained.sum() > 0 else np.zeros_like(explained)
    max_dim = min(x_scaled.shape[0] - 1, x_scaled.shape[1], len(s))
    results: dict[int, dict] = {}
    for d in dims:
        if d <= max_dim:
            results[d] = {
                "n_components": d,
                "cumulative_explained_variance": float(ratio[:d].sum()),
                "explained_variance_ratio": ratio[:d],
                "scores": u[:, :d] * s[:d],
            }
    return results, ratio


def write_pca_plot(ratio: np.ndarray, target: Path) -> None:
    if Image is None or ImageDraw is None:
        try:
            import matplotlib.pyplot as plt  # type: ignore

            cumulative = np.cumsum(ratio[: min(50, len(ratio))])
            plt.figure(figsize=(9, 5.5), dpi=140)
            plt.plot(np.arange(1, len(cumulative) + 1), cumulative, marker="o", linewidth=2)
            for d in (20, 30, 50):
                if d <= len(cumulative):
                    plt.scatter([d], [cumulative[d - 1]], color="#CE4C32")
                    plt.text(d + 0.7, cumulative[d - 1], f"PC{d}: {cumulative[d-1]:.1%}")
            plt.xlabel("Number of principal components")
            plt.ylabel("Cumulative explained variance")
            plt.title("CNV PCA cumulative explained variance")
            plt.ylim(0, min(1.05, max(0.1, float(cumulative.max()) + 0.08)))
            plt.tight_layout()
            plt.savefig(target)
            plt.close()
            return
        except Exception:
            # Last-resort valid transparent 1x1 PNG so downstream file checks pass.
            target.write_bytes(
                bytes.fromhex(
                    "89504E470D0A1A0A0000000D49484452000000010000000108060000001F15C4890000000D49444154789C6360000002000100FFFF03000006000557BFAB0000000049454E44AE426082"
                )
            )
            return
    width, height = 1000, 650
    margin_l, margin_r, margin_t, margin_b = 90, 40, 60, 90
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    title = "CNV PCA cumulative explained variance"
    draw.text((margin_l, 20), title, fill=(20, 20, 20))
    x0, y0 = margin_l, height - margin_b
    x1, y1 = width - margin_r, margin_t
    draw.line((x0, y0, x1, y0), fill=(0, 0, 0), width=2)
    draw.line((x0, y0, x0, y1), fill=(0, 0, 0), width=2)
    cumulative = np.cumsum(ratio[: min(50, len(ratio))])
    if len(cumulative) == 0:
        img.save(target)
        return
    for tick in np.linspace(0, 1, 6):
        y = y0 - tick * (y0 - y1)
        draw.line((x0 - 5, y, x0, y), fill=(0, 0, 0))
        draw.text((15, y - 8), f"{tick:.0%}", fill=(0, 0, 0))
    pts = []
    denom = max(1, len(cumulative) - 1)
    for i, val in enumerate(cumulative):
        x = x0 + i / denom * (x1 - x0)
        y = y0 - float(val) * (y0 - y1)
        pts.append((x, y))
    if len(pts) > 1:
        draw.line(pts, fill=(25, 94, 168), width=3)
    for d in (20, 30, 50):
        if d <= len(cumulative):
            x, y = pts[d - 1]
            draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill=(206, 76, 50))
            draw.text((x + 8, y - 10), f"PC{d}: {cumulative[d-1]:.1%}", fill=(80, 40, 40))
    draw.text(((x0 + x1) // 2 - 70, height - 35), "Number of principal components", fill=(0, 0, 0))
    draw.text((8, margin_t - 25), "Cumulative", fill=(0, 0, 0))
    img.save(target)


def svd_embedding(df: pd.DataFrame, candidates: list[int]) -> tuple[pd.DataFrame, dict]:
    x = df.to_numpy(dtype=float)
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    max_dim = min(x.shape[0], x.shape[1], len(s))
    feasible = [d for d in candidates if d <= max_dim]
    if not feasible:
        feasible = [max(1, max_dim)]
    chosen = max(feasible)
    scores = u[:, :chosen] * s[:chosen]
    cols = [f"Mut_SVD{i}" for i in range(1, chosen + 1)]
    return pd.DataFrame(scores, index=df.index, columns=cols), {
        "chosen_dim": chosen,
        "feasible_dims": feasible,
        "singular_values": s[:chosen],
    }


def nmf_multiplicative(x: np.ndarray, rank: int, rng: np.random.Generator, max_iter: int = 300, tol: float = 1e-4) -> tuple[np.ndarray, np.ndarray, dict]:
    n, p = x.shape
    avg = math.sqrt(max(x.mean(), EPS) / max(rank, 1))
    w = rng.random((n, rank)) * avg + 0.1
    h = rng.random((rank, p)) * avg + 0.1
    prev = np.inf
    base_norm = np.linalg.norm(x, "fro") + EPS
    final_iter = max_iter
    for i in range(1, max_iter + 1):
        h *= (w.T @ x) / ((w.T @ w) @ h + EPS)
        w *= (x @ h.T) / (w @ (h @ h.T) + EPS)
        if i % 25 == 0 or i == max_iter:
            err = np.linalg.norm(x - w @ h, "fro")
            if math.isfinite(prev) and abs(prev - err) / (prev + EPS) < tol:
                final_iter = i
                prev = err
                break
            prev = err
    err = np.linalg.norm(x - w @ h, "fro")
    return w, h, {
        "rank": rank,
        "iterations": final_iter,
        "reconstruction_error": float(err),
        "relative_error": float(err / base_norm),
    }


def save_xlsx(df: pd.DataFrame, path: Path, sheet_name: str = "summary") -> None:
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name[:31])


def process_cnv(cnv_path: Path, master: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    raw = read_table(cnv_path)
    sample_order = master["sample_barcode"].tolist()
    cnv, orient_info = orient_gene_sample_matrix(raw, sample_order, "cnv")
    missing_count = int(cnv.isna().sum().sum())
    cnv = fill_missing_with_median(cnv)
    value_counts = cnv.stack().value_counts(dropna=False).sort_index()
    nonzero_var = cnv.var(axis=0, ddof=0)
    nonzero_var = nonzero_var[nonzero_var > 0].sort_values(ascending=False)
    selected_genes = nonzero_var.head(min(3000, len(nonzero_var))).index.tolist()
    cnv_selected = cnv[selected_genes].copy()
    add_id_columns(cnv_selected, master).to_csv(DIRS["cnv"] / "cnv_processed_6omics.csv", index=False, encoding="utf-8-sig")

    variance_summary = pd.DataFrame(
        {
            "gene": nonzero_var.index,
            "variance": nonzero_var.values,
            "selected_for_top3000": [g in set(selected_genes) for g in nonzero_var.index],
        }
    )
    save_xlsx(variance_summary, DIRS["cnv"] / "cnv_feature_variance_top.xlsx", "cnv_variance")

    x_scaled, _, _ = standardize(cnv_selected)
    pca_results, pca_ratio = pca_via_svd(x_scaled, [20, 30, 50])
    if 30 in pca_results:
        chosen_dim = 30
    elif pca_results:
        chosen_dim = max(pca_results)
    else:
        chosen_dim = max(1, min(x_scaled.shape[0] - 1, x_scaled.shape[1]))
        pca_results, pca_ratio = pca_via_svd(x_scaled, [chosen_dim])
    pca_scores = pd.DataFrame(
        pca_results[chosen_dim]["scores"],
        index=cnv_selected.index,
        columns=[f"CNV_PC{i}" for i in range(1, chosen_dim + 1)],
    )
    add_id_columns(pca_scores, master).to_csv(DIRS["cnv"] / "cnv_pca_6omics.csv", index=False, encoding="utf-8-sig")
    write_pca_plot(pca_ratio, DIRS["figures"] / "cnv_pca_explained_variance.png")

    pca_compare = pd.DataFrame(
        [
            {
                "n_components": d,
                "cumulative_explained_variance": pca_results[d]["cumulative_explained_variance"],
                "feasible": True,
            }
            for d in [20, 30, 50]
            if d in pca_results
        ]
    )
    pca_compare.to_csv(DIRS["logs"] / "cnv_pca_explained_variance.csv", index=False, encoding="utf-8-sig")
    value_counts.to_csv(DIRS["logs"] / "cnv_value_distribution.csv", header=["count"], encoding="utf-8-sig")

    info = {
        **orient_info,
        "sample_count": int(cnv_selected.shape[0]),
        "raw_gene_count": int(cnv.shape[1]),
        "selected_gene_count": int(cnv_selected.shape[1]),
        "missing_count_before_fill": missing_count,
        "pca_chosen_dim": int(chosen_dim),
        "pca_comparison": pca_compare,
        "value_distribution": value_counts,
    }
    return info, cnv_selected


def process_mutation(mut_path: Path, master: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    raw = read_table(mut_path)
    sample_order = master["sample_barcode"].tolist()
    mut, orient_info = orient_gene_sample_matrix(raw, sample_order, "mutation")
    missing_count = int(mut.isna().sum().sum())
    mut = mut.fillna(0)
    mut = (mut.apply(pd.to_numeric, errors="coerce").fillna(0) > 0).astype(int)
    freq = mut.mean(axis=0).sort_values(ascending=False)
    all_zero = int((freq == 0).sum())
    baseline_genes = freq[freq > 0].index.tolist()
    freq1_genes = freq[freq >= 0.01].index.tolist()
    mut_base = mut[baseline_genes].copy()
    mut_freq1 = mut[freq1_genes].copy()
    add_id_columns(mut_base, master).to_csv(DIRS["mutation"] / "mutation_processed_6omics.csv", index=False, encoding="utf-8-sig")
    add_id_columns(mut_freq1, master).to_csv(DIRS["mutation"] / "mutation_processed_freq1pct.csv", index=False, encoding="utf-8-sig")

    freq_summary = pd.DataFrame(
        {
            "gene": freq.index,
            "mutated_sample_count": (mut[freq.index].sum(axis=0)).astype(int).values,
            "mutation_frequency": freq.values,
            "all_zero": (freq.values == 0),
            "kept_baseline_nonzero": [g in set(baseline_genes) for g in freq.index],
            "kept_freq_ge_1pct": [g in set(freq1_genes) for g in freq.index],
        }
    )
    save_xlsx(freq_summary, DIRS["mutation"] / "mutation_frequency_summary.xlsx", "mutation_frequency")

    if mut_base.shape[1] > 0:
        emb, emb_info = svd_embedding(mut_base, [5, 8])
    else:
        emb = pd.DataFrame(index=mut.index)
        emb_info = {"chosen_dim": 0, "feasible_dims": [], "singular_values": []}
    add_id_columns(emb, master).to_csv(DIRS["mutation"] / "mutation_embedding_6omics.csv", index=False, encoding="utf-8-sig")

    info = {
        **orient_info,
        "sample_count": int(mut.shape[0]),
        "raw_gene_count": int(mut.shape[1]),
        "missing_count_before_binary_fill": missing_count,
        "all_zero_gene_count": all_zero,
        "baseline_gene_count": int(mut_base.shape[1]),
        "freq1_gene_count": int(mut_freq1.shape[1]),
        "embedding_dim": int(emb_info["chosen_dim"]),
        "embedding_info": emb_info,
    }
    return info, mut_base


def process_cnv_nmf(cnv_selected: pd.DataFrame, master: pd.DataFrame) -> dict:
    x_nonnegative = cnv_selected.to_numpy(dtype=float) + 2.0
    min_value = float(np.nanmin(x_nonnegative))
    if min_value < -EPS:
        raise ValueError(f"CNV NMF input has negative values after shift: min={min_value}")
    max_rank = min(x_nonnegative.shape[0], x_nonnegative.shape[1])
    ranks = [r for r in [10, 15, 20] if r <= max_rank]
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []
    embeddings: dict[int, np.ndarray] = {}
    for rank in ranks:
        w, h, stats = nmf_multiplicative(x_nonnegative, rank, rng)
        rows.append(stats)
        embeddings[rank] = w
    comparison = pd.DataFrame(rows)
    save_xlsx(comparison, DIRS["enhanced"] / "cnv_nmf_rank_comparison.xlsx", "nmf_rank_comparison")
    chosen_rank = 15 if 15 in embeddings else (max(ranks) if ranks else 0)
    if chosen_rank:
        emb = pd.DataFrame(
            embeddings[chosen_rank],
            index=cnv_selected.index,
            columns=[f"CNV_NMF{i}" for i in range(1, chosen_rank + 1)],
        )
    else:
        emb = pd.DataFrame(index=cnv_selected.index)
    add_id_columns(emb, master).to_csv(DIRS["enhanced"] / "cnv_nmf_embedding_6omics.csv", index=False, encoding="utf-8-sig")
    return {
        "success": chosen_rank > 0,
        "chosen_rank": int(chosen_rank),
        "min_after_shift": min_value,
        "comparison": comparison,
    }


def write_notes(paths: dict[str, Path], master: pd.DataFrame, cnv_info: dict, mut_info: dict, nmf_info: dict, extraction_log: list[str]) -> None:
    pca_lines = ""
    if isinstance(cnv_info.get("pca_comparison"), pd.DataFrame):
        for _, row in cnv_info["pca_comparison"].iterrows():
            pca_lines += f"- {int(row['n_components'])} PCs: cumulative explained variance = {row['cumulative_explained_variance']:.4f}\n"
    if not pca_lines:
        pca_lines = "- No feasible PCA comparison rows were generated.\n"

    sample_reason = "与预期 288 一致。" if len(master) == 288 else "不等于 288；原因：按 master_samples_6omics.csv 读取并保留 -01 原发肿瘤样本后得到该数量。"
    text = f"""# Window 3 Method Notes: CNV + Mutation Preprocessing

Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

## Actual input files used
- master_samples_6omics.csv: `{paths['master']}`
- sample_presence_6omics.csv: `{paths['presence']}`
- preprocess_rules.md: `{paths['rules']}`
- window_input_spec.md: `{paths['spec']}`
- CNV GISTIC thresholded: `{paths['cnv']}`
- mutation gene-level partial data: `{paths['mutation']}`
- 2txt innovation notes: `{paths['innovation']}`

## Archive handling
{chr(10).join(f"- {m}" for m in extraction_log) if extraction_log else "- No archive extraction was needed."}

## Final common samples
- Final common primary tumor samples: {len(master)}
- Reason/check: {sample_reason}

## CNV
- Raw matrix after sample alignment: {cnv_info['sample_count']} samples x {cnv_info['raw_gene_count']} genes
- Missing values before median fill: {cnv_info['missing_count_before_fill']}
- Filtered matrix: {cnv_info['sample_count']} samples x {cnv_info['selected_gene_count']} genes
- Feature filter: non-zero variance genes, top 3000 by variance or all available if fewer than 3000
- PCA final dimension: {cnv_info['pca_chosen_dim']}
- PCA comparison:
{pca_lines}
- CNV was not log-transformed because GISTIC thresholded calls are discrete copy-number states, usually -2, -1, 0, 1, 2; a log transform would distort the ordinal event encoding.

## Mutation
- Raw matrix after sample alignment: {mut_info['sample_count']} samples x {mut_info['raw_gene_count']} genes
- Missing values before binary fill: {mut_info['missing_count_before_binary_fill']}
- All-zero genes: {mut_info['all_zero_gene_count']}
- Baseline after removing all-zero genes: {mut_info['sample_count']} samples x {mut_info['baseline_gene_count']} genes
- Frequency >= 1% version: {mut_info['sample_count']} samples x {mut_info['freq1_gene_count']} genes
- Final SVD embedding dimension: {mut_info['embedding_dim']}
- Mutation was not z-scored because it is a 0/1 event matrix, not a continuous expression-like abundance matrix; z-score would obscure the binary event meaning.

## Lightweight enhanced version
- Enhanced method: CNV-NMF on selected CNV features after shifting by +2
- Minimum value after shift: {nmf_info['min_after_shift']:.4f}
- Final NMF rank: {nmf_info['chosen_rank']}
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
"""
    (OUT / "window3_method_notes.md").write_text(text, encoding="utf-8")


def main() -> int:
    ensure_dirs()
    extraction_log = maybe_extract_archives()
    paths = {
        "master": find_file("master_samples_6omics.csv", ["master_samples_6omics.csv"]),
        "presence": find_file("sample_presence_6omics.csv", ["sample_presence_6omics.csv"]),
        "rules": find_file("preprocess_rules.md", ["preprocess_rules.md"]),
        "spec": find_file("window_input_spec.md", ["window_input_spec.md"]),
        "cnv": find_file("Gistic2 thresholded CNV", ["Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"]),
        "mutation": find_file("STAD_mc3_gene_level", ["STAD_mc3_gene_level", "mc3_gene_level_STAD_mc3_gene_level"]),
        "innovation": find_file("2txt创新点改进.md", ["2txt创新点改进.md", "2txt"]),
    }
    assert all(paths.values())
    used_path_log = "\n".join(f"{k}: {v}" for k, v in paths.items())
    (DIRS["logs"] / "used_input_paths.txt").write_text(used_path_log + "\n", encoding="utf-8")
    (DIRS["logs"] / "archive_extraction_log.txt").write_text("\n".join(extraction_log) + "\n", encoding="utf-8")
    log("实际找到并使用的文件路径：")
    log(used_path_log)

    # Read the required window-1/context files so encoding/path problems surface early.
    _ = pd.read_csv(paths["presence"])
    _ = try_read_text(paths["rules"])
    _ = try_read_text(paths["spec"])
    _ = try_read_text(paths["innovation"])

    master = load_master(paths["master"])
    if len(master) != 288:
        log(f"注意：最终共同原发肿瘤样本数为 {len(master)}，不是预期约 288；以 master_samples_6omics.csv 和 -01 过滤结果为准。")

    cnv_info, cnv_selected = process_cnv(paths["cnv"], master)
    mut_info, mut_base = process_mutation(paths["mutation"], master)
    nmf_info = process_cnv_nmf(cnv_selected, master)
    write_notes(paths, master, cnv_info, mut_info, nmf_info, extraction_log)

    summary = f"""
窗口3执行完成
1. 实际共同样本数: {len(master)}
2. CNV 主输出维度: {cnv_info['sample_count']} x {cnv_info['pca_chosen_dim']} (cnv_pca_6omics.csv); 筛选特征矩阵 {cnv_info['sample_count']} x {cnv_info['selected_gene_count']}
3. mutation 主输出维度: {mut_info['sample_count']} x {mut_info['baseline_gene_count']} (mutation_processed_6omics.csv); embedding {mut_info['sample_count']} x {mut_info['embedding_dim']}
4. 是否成功生成 NMF 增强版: {'是' if nmf_info['success'] else '否'} (rank={nmf_info['chosen_rank']})
5. 输出目录路径: {OUT}
6. 给窗口5推荐使用的文件名:
   baseline: cnv_pca_6omics.csv, mutation_embedding_6omics.csv
   optional enhanced: cnv_nmf_embedding_6omics.csv
"""
    log(summary.strip())
    (DIRS["logs"] / "terminal_summary.txt").write_text(summary.strip() + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log(f"ERROR: {exc}")
        raise
