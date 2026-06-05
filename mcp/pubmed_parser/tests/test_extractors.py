"""
test_extractors.py — 变异信息提取单元测试

测试对象：pubmed_client.extractors 模块
  - extract_pathogenicity()  → 从句子中提取致病性评级
  - extract_zygosity()        → 从句子中提取合子状态（纯合/杂合/复合杂合等）
  - extract_inheritance()     → 从句子中提取遗传模式（显性/隐性等）
  - extract_patient_phenotypes() → 从句子中提取患者表型

测试数据：
  - cases/*.json — 各函数的测试用例文档（供人类阅读，非参数化运行）
"""

import sys
from pathlib import Path

_PACKAGE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PACKAGE_ROOT))

import pytest
from pubmed_client import (
    extract_pathogenicity,
    extract_zygosity,
    extract_inheritance,
    extract_patient_phenotypes,
)

# cases/*.json 是参考测试数据文件（供人类阅读）


# =============================================================================
# TestExtractPathogenicity — 致病性提取
# =============================================================================

class TestExtractPathogenicity:
    """
    测试 extract_pathogenicity() 函数

    输入：
      - sentences : list[str]，变异相关句子列表
      - full_text_lower : str，全文小写版本（用于 fallback 搜索）
    输出：致病性字符串，如 "致病 (pathogenic)" 或 "意义不明 (VUS)"

    策略：
      优先从 sentences 中搜索致病性关键词；
      若未匹配，则 fallback 到 full_text_lower 中搜索（且需有 c./p. 上下文）
    """

    def test_pathogenic_keyword(self):
        """sentences 中包含 'pathogenic' → 直接匹配"""
        sentences = ["This variant is pathogenic."]
        result = extract_pathogenicity(sentences, "")
        assert "pathogenic" in result.lower()

    def test_vus(self):
        """sentences 中包含 'VUS' → 匹配为意义不明"""
        sentences = ["The variant is a VUS."]
        result = extract_pathogenicity(sentences, "")
        assert "VUS" in result or "不确定" in result

    def test_empty_sentences(self):
        """sentences 为空时 → fallback 到 full_text_lower 搜索，结果非空"""
        result = extract_pathogenicity([], "no variant mentioned here")
        assert result != ""


class TestExtractPathogenicityFromCases:
    """
    ⚠️ informational only — cases/extract_pathogenicity_cases.json 供人类阅读。

    此类的参数化测试已禁用，原因：JSON case 结构无法 1:1 映射到函数签名。
    如需添加新测试用例，请在上面的 TestExtractPathogenicity 中添加显式测试方法。
    """


# =============================================================================
# TestExtractZygosity — 合子状态提取
# =============================================================================

class TestExtractZygosity:
    """
    测试 extract_zygosity() 函数

    输入：
      - sentences : list[str]，变异相关句子列表
      - full_text_lower : str，全文小写版本
      - target_cdna / target_protein : 目标变异（用于排除其他变异的干扰）
      - tables : 表格数据（可选）
      - keywords : build_variant_keywords 输出（可选，用于判断"是否包含目标变异"）

    输出：合子状态字符串，如 "复合杂合 (compound heterozygous)"

    策略：
      1. 仅在同时包含目标变异的句子中匹配（排除其他变异的干扰）
      2. "in trans" 证据存在时 → 排除纯合，确认复合杂合
      3. 表格行中的 Het/Hom 标记优先
    """

    def test_compound_heterozygous(self):
        """明确提到 'compound heterozygous' → 复合杂合"""
        sentences = ["The patient had compound heterozygous mutations."]
        result = extract_zygosity(
            sentences, "",
            target_cdna="c.2279C>T", target_protein="p.Thr760Met"
        )
        assert "compound" in result.lower() or "复合杂合" in result

    def test_homozygous(self):
        """明确提到 'homozygous' → 纯合"""
        sentences = ["Homozygous variant was detected."]
        result = extract_zygosity(sentences, "")
        assert "hom" in result.lower() or "纯合" in result

    def test_heterozygous(self):
        """明确提到 'heterozygous' → 杂合"""
        sentences = ["Heterozygous for the variant."]
        result = extract_zygosity(sentences, "")
        assert "het" in result.lower() or "杂合" in result

    def test_in_trans_corects_to_compound(self):
        """句子同时提到 'homozygous' 和 'in trans' → in trans 否定纯合，应修正为复合杂合"""
        sentences = [
            "The variant was homozygous. We could not formally confirm the in trans position."
        ]
        result = extract_zygosity(sentences, "", target_cdna="c.2279C>T")
        assert "compound" in result.lower() or "复合杂合" in result


class TestExtractZygosityFromCases:
    """
    ⚠️ informational only — cases/extract_zygosity_cases.json 供人类阅读。

    原因同上。
    """


# =============================================================================
# TestExtractInheritance — 遗传模式提取
# =============================================================================

class TestExtractInheritance:
    """
    测试 extract_inheritance() 函数

    输入：
      - sentences : list[str]，变异相关句子列表
      - full_text_lower : str，全文小写版本
    输出：遗传模式字符串，如 "常染色体隐性遗传" 或 "常染色体显性遗传"
    """

    def test_autosomal_recessive(self):
        """提到 'autosomal recessive' → 常染色体隐性遗传"""
        sentences = ["Inheritance is autosomal recessive."]
        result = extract_inheritance(sentences, "")
        assert "隐性" in result or "recessive" in result.lower()

    def test_autosomal_dominant(self):
        """提到 'autosomal dominant' → 常染色体显性遗传"""
        sentences = ["Autosomal dominant inheritance."]
        result = extract_inheritance(sentences, "")
        assert "显性" in result or "dominant" in result.lower()

    def test_de_novo(self):
        """提到 'De novo' → 新发突变"""
        sentences = ["De novo mutation."]
        result = extract_inheritance(sentences, "")
        assert "novo" in result.lower()


# =============================================================================
# TestExtractPatientPhenotypes — 患者表型提取
# =============================================================================

class TestExtractPatientPhenotypes:
    """
    测试 extract_patient_phenotypes() 函数

    输入：
      - sentences : list[str]，变异相关句子列表
      - tables : 表格数据（可选）
      - keywords : 关键词（可选）
    输出：list[str]，表型关键词列表
    """

    def test_phenotype_from_keyword(self):
        """句子中包含表型关键词（如 retinitis pigmentosa）→ 提取成功"""
        sentences = ["Patient had retinitis pigmentosa (RP)."]
        result = extract_patient_phenotypes(sentences)
        assert len(result) >= 1