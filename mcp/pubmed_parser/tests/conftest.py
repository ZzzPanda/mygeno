"""
conftest.py — pytest 共享 fixtures

本文件定义所有测试文件共用的 pytest fixtures。

设计原则：
  - 所有 fixtures 从同一个根目录（DATA_DIR）派生，不写死具体文件名
  - 集成测试 fixtures（需要真实 PDF 内容）从 DATA_DIR 中动态发现
  - 单元测试 fixtures（只需函数签名）使用内存合成数据，不依赖外部文件

PDF +预期输出配对规范：
  每个 PDF 文件旁边应有一个同名 .expected.json 文件，如：
    PMID:33374015.pdf  →  PMID:33374015.expected.json
  conftest.py 提供以下配对查找 fixtures：
    pdf_with_expected(pdf) → tuple(Path(pdf), Path(expected_json) | None)
    expected_for(pdf_path) → dict | None（从 conftest.py 调用）

路径约定：
  - PDF 文件所在目录：项目根目录 / data / pdf/
  - 每个 fixture 在目录不存在或无 PDF 时自动跳过测试（pytest.skip）

用法示例：
  # 遍历目录下所有 PDF 及其配对预期输出
  @pytest.fixture(params=all_pdfs())
  def a_pdf_with_expected(request):
      pdf = request.param
      expected = load_expected(pdf)
      if expected is None:
          pytest.skip(f"无配对预期文件: {pdf.with_suffix('.expected.json')}")
      return pdf, expected

  # 加载某个特定 PDF 的预期输出（用于 test_pdf_parser.py）
  def test_foo(pdf_33374015, expected_33374015):
      ...
"""

import sys
from pathlib import Path
import json

# 将 pubmed_client 包加入 sys.path
_PACKAGE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PACKAGE_ROOT))

import pytest

# ── 目录级 fixtures ─────────────────────────────────────────────────────────

# PDF 文件所在根目录
DATA_DIR = Path(__file__).resolve().parent.parent.parent.parent / "data" / "pdf"


@pytest.fixture
def pdf_dir():
    """返回 PDF 文件所在目录路径（Path 对象）。

    用途：测试需要自行过滤/扫描 PDF 文件时使用。
    """
    if not DATA_DIR.exists():
        pytest.skip(f"PDF 目录未找到: {DATA_DIR}")
    return DATA_DIR


def all_pdfs():
    """返回 DATA_DIR 下所有 .pdf 文件的 Path 列表，按文件名排序。"""
    if not DATA_DIR.exists():
        return []
    return sorted(DATA_DIR.glob("*.pdf"))


def load_expected(pdf_path: Path) -> dict | None:
    """加载与 pdf_path配对的 .expected.json 文件。

    配对规则：pdf 文件名.pdf → pdf 文件名.expected.json
    返回 dict 或 None（文件不存在时返回 None，不抛异常）。
    """
    expected_path = pdf_path.with_suffix(".expected.json")
    if not expected_path.exists():
        return None
    with open(expected_path, encoding="utf-8") as f:
        return json.load(f)


# ── 参数化 fixture：遍历目录下所有 PDF ─────────────────────────────────────────

@pytest.fixture(params=all_pdfs())
def a_pdf(request):
    """遍历 DATA_DIR 下所有 PDF 文件，逐个传给测试函数。

    用法示例（直接使用，无需额外代码）：
        class TestFoo:
            def test_bar(self, a_pdf):
                # a_pdf 是 Path 对象，代表当前循环到的 PDF 文件
                text = extract_text_from_pdf(str(a_pdf))
                assert len(text) > 0

    等价于在类或模块级别：
        @pytest.mark.parametrize("a_pdf", all_pdfs(), indirect=True)
    """
    return request.param


# ── 参数化 fixture：遍历目录下所有 PDF 及其配对预期输出 ────────────────────

@pytest.fixture(params=all_pdfs())
def a_pdf_with_expected(request):
    """遍历 DATA_DIR 下所有 PDF及其配对 .expected.json 文件。

    返回 tuple：(pdf_path: Path, expected: dict | None)
    - 如果配对 .expected.json存在 → expected 是 dict
    - 如果不存在 → expected 为 None，测试自动跳过

    用法示例：
        def test_variant_extraction(self, a_pdf_with_expected):
            pdf, expected = a_pdf_with_expected
            if expected is None:
                pytest.skip("无配对预期文件")
            result = extract_variant_info(str(pdf), ...)
            assert result["变异提及"] == expected["变异提及"]
    """
    pdf = request.param
    expected = load_expected(pdf)
    if expected is None:
        pytest.skip(f"无配对预期文件: {pdf.with_suffix('.expected.json')}")
    return pdf, expected


# ── 固定 ID 的 fixtures（按 PMID动态查找，不写死文件名）─────────────────────

@pytest.fixture
def pdf_33374015():
    """CFTR T760M 变异论文 PDF路径。

    实现：动态查找文件名中包含 "33374015" 的 PDF，不写死具体文件名。
    用途：测试能提取到变异信息的完整流程。
    """
    if not DATA_DIR.exists():
        pytest.skip(f"PDF 目录未找到: {DATA_DIR}")
    matches = list(DATA_DIR.glob("*33374015*.pdf"))
    if not matches:
        pytest.skip(f"未找到包含 '33374015' 的 PDF，目录内容: {list(DATA_DIR.iterdir())}")
    return matches[0]


@pytest.fixture
def pdf_34426522():
    """土耳其人群群体遗传学论文 PDF 路径（无目标变异）。

    实现：动态查找文件名中包含 "34426522" 的 PDF，不写死具体文件名。
    用途：测试目标变异不存在时的行为（变异提及=False）。
    """
    if not DATA_DIR.exists():
        pytest.skip(f"PDF 目录未找到: {DATA_DIR}")
    matches = list(DATA_DIR.glob("*34426522*.pdf"))
    if not matches:
        pytest.skip(f"未找到包含 '34426522' 的 PDF，目录内容: {list(DATA_DIR.iterdir())}")
    return matches[0]


# ── 预期输出 fixtures ──────────────────────────────────────────────────────────

@pytest.fixture
def expected_33374015(pdf_33374015):
    """加载 PMID:33374015.pdf 的配对预期输出文件。

    返回 dict 或 None（文件不存在时）。
    """
    return load_expected(pdf_33374015)


@pytest.fixture
def expected_34426522(pdf_34426522):
    """加载 PMID:34426522.pdf 的配对预期输出文件。

    返回 dict 或 None（文件不存在时）。
    """
    return load_expected(pdf_34426522)