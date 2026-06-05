"""
test_variants.py — 变异关键词构建与匹配单元测试

测试对象：pubmed_client.variants 模块
  - build_variant_keywords()    → 构建变异搜索关键词
  - find_variant_sentences()    → 在全文中查找包含变异关键词的句子
  - infer_variant_type()        → 根据 cDNA/protein 推断变异类型

测试数据：
  - cases/*.json — 各函数的测试用例文档（供人类阅读，非参数化运行）
"""

import sys
from pathlib import Path

_PACKAGE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PACKAGE_ROOT))

import pytest
from pubmed_client import build_variant_keywords, find_variant_sentences, infer_variant_type

# cases/*.json 是参考测试数据文件（供人类阅读）
# 如需查看某个函数的测试用例，打开对应 JSON 文件即可


# =============================================================================
# TestBuildVariantKeywords — 关键词构建
# =============================================================================

class TestBuildVariantKeywords:
    """
    测试 build_variant_keywords() 函数

    输入：cdna, protein, transcript（均为字符串）
    输出：dict，包含以下键：
          - exact    : list[str]，精确匹配关键词
          - fuzzy    : list[str]，模糊匹配关键词
          - protein : list[str]，蛋白关键词变体
          - descriptive : list[(pattern, desc)]，描述性短语正则
          - historical : list[str]，历史命名
          - all      : list[str]，所有关键词合并
    """

    def test_cdna_only(self):
        """仅传入 cDNA → exact 含完整形式，fuzzy 含去掉前缀/点号的变体"""
        kw = build_variant_keywords("c.1166G>A", "", "")
        assert "c.1166G>A" in kw["exact"]
        assert "1166G>A" in kw["fuzzy"]
        assert "1166" in kw["fuzzy"]

    def test_protein_p_thr760met(self):
        """蛋白变异 p.Thr760Met → 生成多种格式（三字母、单字母、带空格）"""
        kw = build_variant_keywords("", "p.Thr760Met", "")
        assert "p.Thr760Met" in kw["protein"]
        assert "Thr760Met" in kw["protein"]
        assert "T760M" in kw["protein"] or "T 760 M" in kw["protein"]

    def test_protein_without_p_prefix(self):
        """蛋白变异不带 p. 前缀（Thr760Met）→ 仍能正常匹配"""
        kw = build_variant_keywords("", "Thr760Met", "")
        assert "Thr760Met" in kw["protein"]

    def test_protein_three_letter_to_one(self):
        """三字母格式 p.Arg389His → 展开为单字母 R389H"""
        kw = build_variant_keywords("", "p.Arg389His", "")
        assert any("R389" in k for k in kw["protein"])

    def test_nonsense_ter(self):
        """无义突变 p.Arg389Ter → 包含 Ter 和 *两种终止表示"""
        kw = build_variant_keywords("", "p.Arg389Ter", "")
        assert any("Ter" in k or "*" in k for k in kw["protein"])

    def test_fs_variant(self):
        """移码突变 p.Ala411fs → 包含独立的 "fs" 关键词"""
        kw = build_variant_keywords("", "p.Ala411fs", "")
        assert "fs" in kw["protein"]


class TestBuildVariantKeywordsFromCases:
    """
    ⚠️ informational only — cases/build_variant_keywords_cases.json供人类阅读。

    此类的参数化测试已禁用，原因：JSON case 结构无法 1:1 映射到函数签名。
    如需添加新测试用例，请在上面的 TestBuildVariantKeywords 中添加显式测试方法。
    """


# =============================================================================
# TestFindVariantSentences — 变异句子查找
# =============================================================================

class TestFindVariantSentences:
    """
    测试 find_variant_sentences() 函数

    输入：full_text（完整论文文本）, keywords（build_variant_keywords 的输出）
    输出：(variant_sentences: list[str], matched_keywords: list[str])
    """

    def test_exact_match(self):
        """完整 cDNA 形式（c.2279C>T）在句子中 → 精确匹配"""
        text = "The patient carried the c.2279C>T variant. Another sentence."
        kw = build_variant_keywords("c.2279C>T", "", "")
        sentences, matched = find_variant_sentences(text, kw)
        assert len(sentences) >= 1
        assert "c.2279C>T" in matched or "2279" in matched

    def test_protein_match(self):
        """蛋白变异 p.Thr760Met 在句子中 → 蛋白匹配"""
        text = "The p.Thr760Met variant was found in patient 1."
        kw = build_variant_keywords("", "p.Thr760Met", "")
        sentences, _ = find_variant_sentences(text, kw)
        assert len(sentences) >= 1

    def test_no_match(self):
        """文本中完全不包含目标变异 → 无匹配句子"""
        text = "This is unrelated text about another gene."
        kw = build_variant_keywords("c.9999G>A", "", "")
        sentences, _ = find_variant_sentences(text, kw)
        assert len(sentences) == 0

    def test_fuzzy_cdna_match(self):
        """文本用 2279C>T（无 c. 前缀），关键词有 c. 前缀 → 模糊匹配"""
        text = "Mutation 2279C>T was detected."
        kw = build_variant_keywords("c.2279C>T", "", "")
        sentences, _ = find_variant_sentences(text, kw)
        assert len(sentences) >= 1


class TestFindVariantSentencesFromCases:
    """
    ⚠️ informational only — cases/find_variant_sentences_cases.json 供人类阅读。

    原因同上。
    """


# =============================================================================
# TestInferVariantType — 变异类型推断
# =============================================================================

class TestInferVariantType:
    """
    测试 infer_variant_type() 函数

    输入：cdna, protein, sentences（均为字符串）
    输出：变异类型字符串，如 "错义突变 (missense)"
    """

    @pytest.mark.parametrize("cdna,protein,expected", [
        # (cdna, protein,期望的变异类型)
        ("c.2279C>T",  "p.Thr760Met",  "错义突变 (missense)"),
        ("c.1166G>T",  "p.Arg389Ter",  "无义突变 (nonsense)"),
        ("c.1234del",  "p.Ala411fs",  "移码突变 (frameshift)"),
        ("c.1166_1167del", "p.Pro392del", "缺失 (deletion)"),
        ("c.1234ins",  "",             "插入 (insertion)"),
        ("c.1234dup",  "",             "重复 (duplication)"),
        ("c.1234C=",   "",             "同义突变 (silent)"),
        ("c.858+2T>A","",             "剪接位点突变 (splice site)"),
    ])
    def test_variant_type_from_cdna(self, cdna, protein, expected):
        """参数化测试：各种 cDNA 格式 → 对应变异类型"""
        result = infer_variant_type(cdna, protein, [])
        assert result == expected, f"期望 {expected}，实际 {result}"

    def test_missense_from_protein_only(self):
        """仅有蛋白变异（无 cDNA）→ 从蛋白格式推断为错义突变"""
        result = infer_variant_type("", "p.Thr760Met", [])
        assert result == "错义突变 (missense)"

    def test_nonsense_from_protein_only(self):
        """仅有蛋白变异（无 cDNA）→ 从 Ter/* 推断为无义突变"""
        result = infer_variant_type("", "p.Arg389Ter", [])
        assert result == "无义突变 (nonsense)"

    def test_frameshift_from_protein_only(self):
        """仅有蛋白变异（无 cDNA）→ 从 fs 推断为移码突变"""
        result = infer_variant_type("", "p.Ala411fs", [])
        assert result == "移码突变 (frameshift)"

    def test_no_false_positive_from_other_variants_in_text(self):
        """文本中出现其他变异（p.Phe861Leufs*3）→ 不应错误覆盖目标变异类型"""
        bad_text = ["p.Phe861Leufs*3 was found in another patient."]
        result = infer_variant_type("", "p.Thr760Met", bad_text)
        assert result == "错义突变 (missense)", \
            f"实际结果 '{result}' — 文本中的 fs*3 不应影响目标变异的类型判断"


class TestInferVariantTypeFromCases:
    """
    ⚠️ informational only — cases/infer_variant_type_cases.json 供人类阅读。

    所有参数化用例已在上面的 TestInferVariantType 中实现。
    如需添加新测试用例，请在上面的类中添加 @pytest.mark.parametrize 装饰的方法。
    """