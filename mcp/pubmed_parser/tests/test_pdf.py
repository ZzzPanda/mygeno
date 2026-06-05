"""
test_pdf.py — PDF 文本和表格提取单元测试

测试对象：pubmed_client.pdf 模块
  - extract_text_from_pdf(pdf_path)  → 从 PDF 提取全文文本
  - extract_tables_from_pdf(pdf_path) → 从 PDF 提取表格（返回 list[dict]）

测试策略：
  - 使用真实 PDF 文件（通过 conftest.py 的 fixtures）
  - 每个测试函数验证一个具体的行为/边界条件
"""

import sys
from pathlib import Path

# 将 pubmed_client 包加入 sys.path
_PACKAGE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PACKAGE_ROOT))

import pytest
from pubmed_client import extract_text_from_pdf, extract_tables_from_pdf


class TestExtractTextFromPDF:
    """
    测试 extract_text_from_pdf() 函数

    输入：PDF 文件路径（字符串）
    输出：PDF 全文内容（字符串）
    """

    def test_returns_string(self, pdf_33374015):
        """验证函数返回类型为字符串"""
        text = extract_text_from_pdf(str(pdf_33374015))
        assert isinstance(text, str)
        assert len(text) > 1000, f"文本过短（{len(text)} 字符），可能提取失败"

    def test_contains_gene_name(self, pdf_33374015):
        """验证 PDF 中包含 CFTR 基因名（大小写不敏感）"""
        text = extract_text_from_pdf(str(pdf_33374015))
        assert "CFTR" in text or "cftr" in text.lower()

    def test_contains_variant_keyword(self, pdf_33374015):
        """验证 PDF 中包含变异位置 760（T760M 论文）"""
        text = extract_text_from_pdf(str(pdf_33374015))
        assert "760" in text, "PDF 中未找到位置 760，可能不是正确的 T760M 论文"

    def test_empty_for_nonexistent_file(self):
        """边界条件：文件不存在时返回空字符串（不抛出异常）"""
        result = extract_text_from_pdf("/nonexistent/file.pdf")
        assert result == ""


class TestExtractTablesFromPDF:
    """
    测试 extract_tables_from_pdf() 函数

    输入：PDF 文件路径（字符串）
    输出：list[dict]，每个 dict 包含：
          - id : str，表格 ID（如 "page_1"）
          - caption : str，表格标题（可能为空）
          - rows : list[list[str]]，每行是单元格字符串列表
    """

    def test_returns_list(self, pdf_33374015):
        """验证函数返回类型为 list"""
        tables = extract_tables_from_pdf(str(pdf_33374015))
        assert isinstance(tables, list)

    def test_table_has_rows(self, pdf_33374015):
        """验证每个表格 dict 包含 rows 字段，且至少有一行"""
        tables = extract_tables_from_pdf(str(pdf_33374015))
        assert len(tables) >= 1, "应该找到至少 1 个表格"
        assert "rows" in tables[0]
        assert len(tables[0]["rows"]) >= 1

    def test_table_row_cells_are_strings(self, pdf_33374015):
        """验证表格每行的每个单元格都是字符串（而非 None 或其他类型）"""
        tables = extract_tables_from_pdf(str(pdf_33374015))
        row = tables[0]["rows"][0]
        for cell in row:
            assert isinstance(cell, str), f"单元格应为 str，实际为 {type(cell)}"