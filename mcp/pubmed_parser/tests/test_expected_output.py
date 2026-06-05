"""
test_expected_output.py — 预期输出对比测试（informational）

测试策略：
  将 extract_variant_info() 的实际输出与 data/pdf/variant_results.json
  中的 ground-truth 进行对比，打印逐字段 diff，供人工审核。

运行方式（查看打印输出）：
  .venv/bin/python -m pytest tests/test_expected_output.py -v -s

注意：
  这类测试不以 PASS/FAIL 为判定标准，而是打印 diff 供人类判断。
  所有字段均打印，无论是否匹配。
"""

import sys
from pathlib import Path

_PACKAGE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PACKAGE_ROOT))

import json
import pytest
from server import extract_variant_info

# Ground-truth 文件路径（项目根目录的 data/pdf/）
DATA_DIR = Path(__file__).resolve().parent.parent.parent.parent / "data" / "pdf"
EXPECTED_FILE = DATA_DIR / "variant_results.json"


@pytest.fixture
def expected_results():
    """
    从 variant_results.json 加载 ground-truth 数据。

    文件格式：list[dict]，每个 dict 的 "文件" 字段存储 PDF 文件名。
    加载后转为 dict，以文件名为 key，方便按文件名查找。
    """
    if not EXPECTED_FILE.exists():
        pytest.skip(f"Ground-truth 文件未找到: {EXPECTED_FILE}")
    with open(EXPECTED_FILE, encoding="utf-8") as f:
        return {r["文件"]: r for r in json.load(f)}


class TestExpectedOutput:
    """
    对比实际输出与 ground-truth，打印逐字段 diff。

    ground-truth 来源：data/pdf/variant_results.json
    实际输出来源：server.extract_variant_info()

    打印格式：
      ✓ 字段名: expected=... | actual=...   （匹配）
      ✗ 字段名: expected=... | actual=...   （不匹配）
    """

    def test_33374015_vs_expected(self, pdf_33374015, expected_results):
        """PMID:33374015.pdf（CFTR T760M 论文）实际输出 vs 预期输出"""
        expected = expected_results.get("PMID:33374015.pdf", {})

        # 从 ground-truth 的匹配关键词中提取 cdna/protein（用于调用 extract_variant_info）
        cdna, protein = "", ""
        for kw in expected.get("匹配关键词", []):
            if not cdna and kw.startswith("c."):
                cdna = kw
            if not protein and kw.startswith("p."):
                protein = kw

        result = extract_variant_info(
            str(pdf_33374015), cdna=cdna, protein=protein, gene=expected.get("基因", "")
        )

        print("\n" + "=" * 60)
        print("PMID:33374015.pdf — actual vs expected")
        print("=" * 60)
        for field, exp_val in expected.items():
            if exp_val in ("", None):
                continue
            act_val = result.get(field)
            match = "✓" if exp_val == act_val else "✗"
            print(f"  {match} {field}:")
            print(f"      expected: {str(exp_val)[:80]}")
            print(f"      actual:   {str(act_val)[:80]}")

    def test_34426522_vs_expected(self, pdf_34426522, expected_results):
        """PMID:34426522.pdf（土耳其人群群体遗传学论文）实际输出 vs 预期输出"""
        expected = expected_results.get("PMID:34426522.pdf", {})

        result = extract_variant_info(
            str(pdf_34426522), cdna="", protein="", gene=expected.get("基因", "")
        )

        print("\n" + "=" * 60)
        print("PMID:34426522.pdf — actual vs expected")
        print("=" * 60)
        for field, exp_val in expected.items():
            if exp_val in ("", None):
                continue
            act_val = result.get(field)
            match = "✓" if exp_val == act_val else "✗"
            print(f"  {match} {field}:")
            print(f"      expected: {str(exp_val)[:80]}")
            print(f"      actual:   {str(act_val)[:80]}")