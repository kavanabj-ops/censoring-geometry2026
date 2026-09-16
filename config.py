# -*- coding: utf-8 -*-
"""
config.py —— 全局配置中心
============================================================
集中管理所有路径、样本元数据、基因集与超参数。
各分析脚本通过 `import config` 引用本模块，避免路径/参数在各脚本间散落、不一致。

【用户必读】
1. 把所有 "TODO" 处的路径、样本标签、基因集按实际情况填好。
2. 样本量很小（类器官十几个~几十个），建模阶段务必用简单模型 + 交叉验证，
   本框架已内置嵌套 CV / LOOCV，切勿自行改成深度学习。
3. 跨「单细胞 → bulk」验证走"基因签名迁移"，而非完整模型迁移，详见 06。
"""
from pathlib import Path

# ============================================================================
# 一、目录结构（自动创建，无需手建）
# ============================================================================
BASE_DIR = Path(r"D:\paper\02_analysis")
DATA_DIR       = BASE_DIR / "data"
RAW_DIR        = DATA_DIR / "raw"        # 原始数据：h5ad / 参考 / 外部 bulk
PROCESSED_DIR  = DATA_DIR / "processed"  # 各步骤产出的中间 h5ad
FEATURE_DIR    = DATA_DIR / "features"   # 样本级特征矩阵
MODEL_DIR      = BASE_DIR / "models"     # 训练好的模型 + 签名
RESULT_DIR     = BASE_DIR / "results"    # 指标、表格、SHAP 值
FIG_DIR        = BASE_DIR / "figures"    # 所有图

for _d in [RAW_DIR, PROCESSED_DIR, FEATURE_DIR, MODEL_DIR, RESULT_DIR, FIG_DIR]:
    _d.mkdir(parents=True, exist_ok=True)

# ============================================================================
# 二、输入数据路径（TODO：填真实路径）
# ============================================================================
# GSE312260 结直肠癌类器官 scRNA h5ad（已合并 CRC5+CRC11，基因名用 ENSG 以对齐参考）
GSE312260_H5AD = Path(r"D:\paper\data\GSE312260_merged_ensg.h5ad")

# Tabula Sapiens 2.0 肠道参考（精简版：子采样 2.5万细胞 + HVG 4000，ENSG，含 broad_cell_class/donor_id）
TABULA_SAPIENS_GUT_H5AD = Path(r"D:\paper\data\Tabula_gut_reference_slim.h5ad")

# TODO: Kong 2020 (Nat Commun) CRC bulk 转录组 + IC50 表（外部验证）
#       expr 行为基因、列为样本；ic50 为 "样本 -> 伊立替康 IC50" 两列
KONG2020_EXPR_CSV = RAW_DIR / "Kong2020_expr.csv"
KONG2020_IC50_CSV = RAW_DIR / "Kong2020_ic50.csv"

# ============================================================================
# 三、样本元数据（TODO：补齐全部类器官 → 药敏标签）
# ============================================================================
# 药敏标签约定：0 = 敏感，1 = 耐药。
# GSE312260 中 CRC5 为伊立替康耐药、CRC11 为敏感（示例），需按数据集实际情况补全。
SAMPLE_METADATA = {
    "CRC5":  1,   # 耐药
    "CRC11": 0,   # 敏感
    # TODO: 继续补充其它类器官样本，例如：
    # "CRC1": 0,
    # "CRC7": 1,
}
RESISTANT_LABEL = 1   # 标签中代表"耐药"的值
SENSITIVE_LABEL = 0   # 标签中代表"敏感"的值

# h5ad 的 obs 里标记样本来源的列名
SAMPLE_KEY = "sample"          # 每个细胞属于哪个类器官样本
# TODO: 若 h5ad 中样本列名不同，改这里；例如 "donor_id" / "orig.ident"

# ============================================================================
# 四、细胞注释相关（01_annotation.py）
# ============================================================================
# Tabula Sapiens 参考数据里"细胞类型"注释列名（用粗粒度大类，14 类，训练更稳）
REFERENCE_LABEL_KEY = "broad_cell_class"
REFERENCE_BATCH_KEY = "donor_id"                 # 参考数据的批次/供体列（scANVI 需要）

# 参考 → 目标粗粒度细胞大类映射（Tabula Sapiens broad_cell_class 的 14 类 → 5 大类）
MAJOR_CELLTYPE_MAP = {
    # —— 上皮 ——
    "intestinal epithelial cell": "Epithelial",
    "glandular epithelial cell": "Epithelial",
    "endo-epithelial cell":       "Epithelial",
    "stem cell":                  "Epithelial",
    # —— 间质 ——
    "fibroblast":                 "Stromal",
    "contractile cell":           "Stromal",
    # —— 免疫 ——
    "t cell":                     "Immune",
    "lymphocyte of b lineage":    "Immune",
    "myeloid leukocyte":          "Immune",
    "dendritic cell":             "Immune",
    "granulocyte":                "Immune",
    "innate lymphoid cell":       "Immune",
    # —— 内皮 ——
    "endothelial cell":           "Endothelial",
    # —— 胶质 ——
    "glial cell":                 "Glial",
}

# ============================================================================
# 五、耐药相关基因集（03_feature_engineering.py 打分用）
# ============================================================================
# TODO: 按伊立替康耐药机制补充/修正基因集（可参考 MSigDB、文献）
RESISTANCE_GENESETS = {
    # 药物外排泵（ABCB1/ABCG2 等介导 SN-38 外排）
    "DRUG_EFFLUX": ["ABCB1", "ABCG2", "ABCC1", "ABCC2", "ABCC5", "ABCC10"],
    # 拓扑异构酶 I 靶点（伊立替康活性代谢物 SN-38 的靶点）
    "TOP1_TARGET": ["TOP1", "TOP1MT"],
    # DNA 损伤修复（修复 SN-38 引起的 DNA 断裂）
    "DNA_REPAIR": ["MSH2", "MLH1", "RAD51", "XRCC1", "ERCC1", "MRE11", "NBN", "RAD50", "BRCA1"],
    # 上皮-间质转化 EMT（常与耐药相关）
    "EMT": ["VIM", "CDH2", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "FN1", "CDH1"],
    # 干性 / 类器官干样细胞
    "STEMNESS": ["LGR5", "PROM1", "ALDH1A1", "ASCL2", "OLFM4", "SMOC2"],
    # 抗凋亡 / 生存
    "APOPTOSIS": ["BCL2", "BCL2L1", "MCL1", "BIRC5", "BAX", "BAK1"],
    # SN-38 葡萄糖醛酸化（UGT1A 失活 SN-38）
    "UGT_GLUCURONIDATION": ["UGT1A1", "UGT1A9"],
}

# ============================================================================
# 六、预处理 / 聚类超参数（00_preprocess.py）
# ============================================================================
PREPROCESS = {
    "min_genes": 200,       # 每个细胞至少表达的基因数
    "max_genes": 8000,      # 每个细胞最多表达基因数（过高可能是双峰）
    "max_pct_mt": 20.0,     # 线粒体比例上限（%）
    "min_cells": 3,         # 每个基因至少在多少个细胞中表达
    "target_sum": 1e4,      # 标准化目标总 counts
    "n_hvg": 2000,          # 高变基因数
    "n_pcs": 30,            # PCA 主成分数
    "n_neighbors": 15,      # 图邻居数
    "leiden_resolution": 1.0,  # Leiden 聚类分辨率
    "mt_prefix": "MT-",     # 线粒体基因前缀（人）
    "ribo_prefix": ("RPS", "RPL"),  # 核糖体基因前缀
    "scrublet_expected_doublet_rate": 0.06,  # scrublet 预期双峰率
}

# ============================================================================
# 七、差异分析 / 特征工程 / 建模超参数
# ============================================================================
DIFF = {
    "min_cells_per_group": 10,  # 每个细胞类型在每个条件里最少细胞数，低于则跳过
    "logfc_threshold": 0.5,     # 差异基因 logFC 阈值
    "pval_threshold": 0.05,     # 差异基因显著性阈值
    "n_top_deg": 50,            # 每个细胞类型保留的 top 差异基因数（构建签名用）
    "pseudobulk_min_cells": 20, # 伪 bulk 聚合所需最少细胞数
}

FEATURE = {
    "gene_set_method": "aucell",  # 基因集打分法：aucell（推荐）/ ssgsea
    "aucell_min_genes_pct": 0.05, # AUCell 通路内最小基因覆盖比例
}

MODEL = {
    "outer_cv_splits": 5,       # 外层 CV 折数（样本太少时脚本会自动退化为 LOOCV）
    "inner_cv_splits": 3,       # 内层（超参搜索）折数
    "n_repeats": 20,            # 重复分层 CV 次数（小样本稳定性评估）
    "random_state": 42,
    "test_size": 0.2,           # 留出测试集比例（Holdout 评估）
    "n_permutations": 100,      # 置换检验次数（评估是否显著优于随机）
}
