"""
test_integration.py — 端到端集成测试

测试对象：server.py 模块（整个解析流程）
  - extract_variant_info() → 从 PDF 文件提取目标变异信息
  - analyze_variant() → 从文本（无 PDF）提取目标变异信息

测试策略：
  使用真实 PDF 文件，验证从文本提取 → 关键词构建 → 句子查找 →
  信息提取的完整流程是否正确。
"""

import sys
from pathlib import Path

_PACKAGE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PACKAGE_ROOT))


# =============================================================================
# TestExtractVariantInfo — PDF 端到端测试
# =============================================================================

class TestExtractVariantInfo:
    """
    测试 extract_variant_info() 完整流程

    输入：pdf_path, cdna, protein, gene
    输出：dict，包含以下字段：
          变异提及 | 基因 | 变异类型 | 致病性 | 合子状态 | 临床表型 等
    """

    def test_cftr_t760m_variant_mention(self, pdf_33374015):
        """验证能从 CFTR T760M 论文中检测到 p.Thr760Met 变异提及"""
        from server import extract_variant_info

        result = extract_variant_info(
            pdf_path=str(pdf_33374015),
            cdna="c.2279C>T",
            protein="p.Thr760Met",
            gene="CFTR",
        )
        assert result["变异提及"] is True
        assert result["基因"] == "CFTR"

    def test_cftr_t760m_variant_type_missense(self, pdf_33374015):
        """验证变异类型推断为错义突变"""
        from server import extract_variant_info

        result = extract_variant_info(
            pdf_path=str(pdf_33374015),
            cdna="c.2279C>T",
            protein="p.Thr760Met",
            gene="CFTR",
        )
        assert result["变异类型"] == "错义突变 (missense)", \
            f"实际结果: {result['变异类型']}"

    def test_cftr_t760m_zygosity_compound_heterozygous(self, pdf_33374015):
        """验证合子状态推断为复合杂合"""
        from server import extract_variant_info

        result = extract_variant_info(
            pdf_path=str(pdf_33374015),
            cdna="c.2279C>T",
            protein="p.Thr760Met",
            gene="CFTR",
        )
        assert "复合杂合" in result["合子状态"] or "compound" in result["合子状态"].lower()

    def test_paper_without_variant(self, pdf_34426522):
        """验证目标变异不存在时，变异提及=False"""
        from server import extract_variant_info

        result = extract_variant_info(
            pdf_path=str(pdf_34426522),
            cdna="c.9999G>A",
            protein="p.X9999X",
            gene="CFTR",
        )
        assert result["变异提及"] is False


# =============================================================================
# TestAnalyzeVariant — 纯文本端到端测试（无需 PDF）
# =============================================================================

class TestAnalyzeVariant:
    """
    测试 analyze_variant() 纯文本流程

    输入：text（论文文本）, cdna, protein, gene
    输出：与 extract_variant_info 相同的 dict 结构
    """

    def test_analyze_variant_from_text(self):
        """验证从文本中提取变异信息的完整流程"""
        from server import analyze_variant

        text = """
        We identified a pathogenic c.1166G>A (p.Arg389His) variant in
        compound heterozygous state. The patient was homozygous for
        the second variant. Inheritance is autosomal recessive.
        """
        result = analyze_variant(text, cdna="c.1166G>A", protein="p.Arg389His", gene="GENE")
        assert result["变异提及"] is True
        assert "pathogenic" in result["致病性"].lower()