"""
test_pdf.py — PDF 文本和表格提取单元测试

测试对象：pubmed_client.pdf 模块
  - extract_text_from_pdf(pdf_path)  → 从 PDF 提取全文文本
  - extract_tables_from_pdf(pdf_path) → 从 PDF 提取表格（返回 list[dict]）

测试策略：
  - 使用 conftest.py 的 a_pdf 参数化 fixture，自动遍历 data/pdf/ 下所有 PDF
  - 每个测试函数验证一个具体的行为/边界条件
  - 不依赖任何特定 PDF 文件名或内容，只验证函数签名的正确性
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

    注意：使用 a_pdf 参数化 fixture，自动遍历目录下所有 PDF，
          测试在任何一个 PDF 上失败都会导致整体失败（快速反馈）。
    """

    def test_returns_string(self, a_pdf):
        """验证函数返回类型为字符串"""
        text = extract_text_from_pdf(str(a_pdf))
        assert isinstance(text, str)

    def test_non_empty(self, a_pdf):
        """验证从 PDF 提取的文本非空"""
        text = extract_text_from_pdf(str(a_pdf))
        assert len(text) > 0, f"PDF {a_pdf.name} 提取结果为空"

    def test_contains_text_blocks(self, a_pdf):
        """验证提取的文本包含换行符（说明是分页提取的，不是单行拼接）"""
        text = extract_text_from_pdf(str(a_pdf))
        assert "\n" in text, f"PDF {a_pdf.name} 文本无双换行符，可能未分页提取"

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

    注意：表格提取能力因 PDF 而异，使用 a_pdf 遍历所有 PDF，
          只要有一个 PDF 能提取到表格就算基本通过。
    """

    def test_returns_list(self, a_pdf):
        """验证函数返回类型为 list"""
        tables = extract_tables_from_pdf(str(a_pdf))
        assert isinstance(tables, list)

    def test_table_structure(self, a_pdf):
        """验证每个表格 dict 包含必要字段：id, rows"""
        tables = extract_tables_from_pdf(str(a_pdf))
        #允许 PDF 没有表格（返回空 list），但不接受结构缺失
        for t in tables:
            assert "id" in t, f"表格缺少 'id' 字段: {t}"
            assert "rows" in t, f"表格缺少 'rows' 字段: {t}"
            assert isinstance(t["rows"], list), f"rows 应为 list，实际 {type(t['rows'])}"

    def test_table_row_cells_are_strings(self, a_pdf):
        """验证表格每行的每个单元格都是字符串（而非 None 或其他类型）"""
        tables = extract_tables_from_pdf(str(a_pdf))
        for t in tables:
            for row in t["rows"]:
                for cell in row:
                    assert isinstance(cell, str), \
                        f"PDF {a_pdf.name} 表格 {t['id']} 单元格类型错误：期望 str，实际 {type(cell)}"