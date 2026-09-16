# -*- coding: utf-8 -*-
"""
utils.py —— 共享工具函数
============================================================
被各脚本复用的通用函数：路径注入、日志、样本标签读取、绘图样式、签名保存/加载等。
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # 无显示环境下静默出图，避免 GUI 报错
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

# 让脚本无论从哪里启动都能 import 到同目录的 config
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402


# ---------------------------------------------------------------------------
# 日志
# ---------------------------------------------------------------------------
def log(msg: str):
    """带时间戳的简易日志。"""
    import datetime
    print(f"[{datetime.datetime.now():%H:%M:%S}] {msg}", flush=True)


# ---------------------------------------------------------------------------
# 样本标签
# ---------------------------------------------------------------------------
def get_sample_labels(samples) -> pd.Series:
    """
    由样本 ID 列表生成药敏标签 Series。
    samples: 可迭代的样本 ID（类器官名）。
    返回:  index=样本ID, value=0/1 标签的 Series（未在 config.SAMPLE_METADATA 中的样本会报错）。
    """
    samples = pd.Series(list(samples)).astype(str)
    missing = [s for s in samples.unique() if s not in config.SAMPLE_METADATA]
    if missing:
        raise KeyError(
            f"以下样本在 config.SAMPLE_METADATA 中缺失，请补全: {missing}"
        )
    return samples.map(config.SAMPLE_METADATA)


# ---------------------------------------------------------------------------
# 绘图样式
# ---------------------------------------------------------------------------
def set_style():
    """统一 matplotlib 风格。"""
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "figure.figsize": (7, 5),
    })


def savefig(fig, path, tight=True):
    """统一保存图：建目录 + bbox_inches='tight'。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight" if tight else None, dpi=200)
    plt.close(fig)
    log(f"已保存图: {path}")


# ---------------------------------------------------------------------------
# 基因集打分（decoupler 封装，样本级 vs 细胞级）
# ---------------------------------------------------------------------------
def score_genesets(adata, genesets, method="aucell", use_raw=False):
    """
    用 decoupler 对 adata 做基因集打分。
    返回 (scores, acts)：scores 为 adata.obsm[f'{method}_estimate'] 型矩阵，
    实际返回为行=细胞、列=基因集 的 DataFrame。
    """
    import decoupler as dc
    # 把 dict -> decoupler 需要的长表格式
    gs = pd.DataFrame(
        [(name, gene) for name, genes in genesets.items() for gene in genes],
        columns=["geneset", "gene"],
    )
    if method == "aucell":
        # 基因名必须在 adata.var_names 中，decoupler 会自动处理
        dc.run_aucell(adata, gs, use_raw=use_raw, min_n=config.FEATURE["aucell_min_genes_pct"])
        scores = adata.obsm["aucell_estimate"]
    elif method == "ssgsea":
        dc.run_gsva(adata, gs, use_raw=use_raw, verbose=False)
        scores = adata.obsm["gsva_estimate"]
    else:
        raise ValueError(f"未知打分法: {method}")
    scores = pd.DataFrame(scores, index=adata.obs_names, columns=list(genesets.keys()))
    return scores


# ---------------------------------------------------------------------------
# 签名保存 / 加载（05/06 复用）
# ---------------------------------------------------------------------------
def save_signature(signature: dict, path):
    """
    保存基因签名到 joblib。
    signature: {"name": {"genes": [...], "weights": [...], "direction": 1/-1}, ...}
    """
    import joblib
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(signature, path)
    log(f"已保存签名: {path}")


def load_signature(path) -> dict:
    import joblib
    return joblib.load(path)


# ---------------------------------------------------------------------------
# 伪 bulk 聚合（decoupler 封装）
# ---------------------------------------------------------------------------
def make_pseudobulk(adata, sample_key=None, min_cells=10):
    """
    将单细胞 adata 按样本聚合为伪 bulk 表达矩阵。
    返回 (pseudobulk_adata, sample_ids)。
    """
    import decoupler as dc
    sample_key = sample_key or config.SAMPLE_KEY
    adata = adata.copy()
    # 过滤细胞数过少的样本
    keep = adata.obs.groupby(sample_key)[sample_key].transform("size") >= min_cells
    adata = adata[keep].copy()
    pb = dc.get_pseudobulk(
        adata, sample_col=sample_key, groups_col=None,
        min_prop=0.0, min_smpls=0, min_cells=min_cells, min_counts=0,
    )
    return pb
