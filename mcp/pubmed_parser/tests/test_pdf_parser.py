"""
test_pdf_parser.py — 原始单文件测试套件（保留）

本文件保留原有的单文件格式，便于快速运行全部测试。
所有 fixtures来自 conftest.py，不写死具体 PDF 文件名。

新增泛化测试：
  - 使用 a_pdf 参数化 fixture，遍历 data/pdf/ 下所有 PDF
  - 使用 pdf_dir fixture 支持自定义过滤逻辑

运行：
  .venv/bin/python -m pytest tests/test_pdf_parser.py -v
  .venv/bin/python -m pytest tests/test_pdf_parser.py -k "keyword" -v
"""

import sys
from pathlib import Path

_PACKAGE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PACKAGE_ROOT))

from pubmed_client import extract_text_from_pdf, extract_tables_from_pdf


# ═══════════════════════════════════════════════════════════════════════════════
# 泛化测试：遍历目录下所有 PDF（a_pdf 参数化 fixture）
# ═══════════════════════════════════════════════════════════════════════════════

class TestAnyPDF:
    """
    使用 a_pdf 参数化 fixture，自动遍历 data/pdf/ 下所有 PDF。

    测试目标：验证 extract_text_from_pdf / extract_tables_from_pdf
    在任意 PDF 上都能正常工作（返回正确类型，不抛异常）。
    """

    def test_extract_text_returns_string(self, a_pdf):
        """任意 PDF → extract_text_from_pdf 返回字符串"""
        text = extract_text_from_pdf(str(a_pdf))
        assert isinstance(text, str)

    def test_extract_tables_returns_list(self, a_pdf):
        """任意 PDF → extract_tables_from_pdf 返回 list"""
        tables = extract_tables_from_pdf(str(a_pdf))
        assert isinstance(tables, list)

    def test_extract_tables_cells_are_strings(self, a_pdf):
        """任意 PDF → 表格单元格均为字符串"""
        tables = extract_tables_from_pdf(str(a_pdf))
        for t in tables:
            for row in t["rows"]:
                for cell in row:
                    assert isinstance(cell, str), \
                        f"单元格类型错误：期望 str，实际 {type(cell)}"


# ═══════════════════════════════════════════════════════════════════════════════
# 特定 PDF 测试（通过 conftest.py 的 pdf_33374015 / pdf_34426522 fixtures）
# 这些 fixtures 按 PMID 动态查找文件名，不写死路径
# ═══════════════════════════════════════════════════════════════════════════════

class TestExtractTextFromPDF:
    """test extract_text_from_pdf（使用特定 PDF fixtures）"""

    def test_returns_string(self, pdf_33374015):
        text = extract_text_from_pdf(str(pdf_33374015))
        assert isinstance(text, str)
        assert len(text) > 1000, f"Text too short ({len(text)} chars)"

    def test_contains_gene_name(self, pdf_33374015):
        text = extract_text_from_pdf(str(pdf_33374015))
        assert "CFTR" in text or "cftr" in text.lower()

    def test_contains_variant_keyword(self, pdf_33374015):
        text = extract_text_from_pdf(str(pdf_33374015))
        assert "760" in text, "Position 760 not found in text"

    def test_empty_for_nonexistent_file(self):
        result = extract_text_from_pdf("/nonexistent/file.pdf")
        assert result == ""


class TestExtractTablesFromPDF:
    """test extract_tables_from_pdf（使用特定 PDF fixtures）"""

    def test_returns_list(self, pdf_33374015):
        tables = extract_tables_from_pdf(str(pdf_33374015))
        assert isinstance(tables, list)

    def test_table_has_rows(self, pdf_33374015):
        tables = extract_tables_from_pdf(str(pdf_33374015))
        assert len(tables) >= 1, "Should find at least 1 table"
        assert "rows" in tables[0]
        assert len(tables[0]["rows"]) >= 1

    def test_table_row_cells_are_strings(self, pdf_33374015):
        tables = extract_tables_from_pdf(str(pdf_33374015))
        row = tables[0]["rows"][0]
        for cell in row:
            assert isinstance(cell, str), f"Cell should be str, got {type(cell)}"