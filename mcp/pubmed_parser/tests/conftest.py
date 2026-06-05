"""
conftest.py — pytest 共享 fixtures

本文件定义所有测试文件共用的 pytest fixtures，
所有测试模块都从这里导入 PDF 路径等共享资源。

路径约定：
  - PDF 文件位于项目根目录的 data/pdf/ 下
  - 每个 fixture 在文件不存在时自动跳过测试（pytest.skip）
"""

import sys
from pathlib import Path

# 将 pubmed_client 包加入 sys.path，以便各测试模块直接 import
_PACKAGE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PACKAGE_ROOT))

import pytest

# PDF 文件所在目录（项目根目录的 data/pdf/）
DATA_DIR = Path(__file__).resolve().parent.parent.parent.parent / "data" / "pdf"

# 具体 PDF 文件路径
#   - PMID:33374015.pdf : CFTR T760M 变异论文（用于有变异提及的集成测试）
#   - PMID:34426522.pdf : 土耳其人群群体遗传学论文（用于无变异提及的对照测试）
PDF_33374015 = DATA_DIR / "PMID:33374015.pdf"
PDF_34426522 = DATA_DIR / "PMID:34426522.pdf"


@pytest.fixture
def pdf_33374015():
    """CFTR T760M 变异论文 PDF路径。

    用途：测试能提取到变异信息的完整流程。
    预期：PDF 中包含 p.Thr760Met 等变异提及。
    """
    if not PDF_33374015.exists():
        pytest.skip(f"PDF 未找到: {PDF_33374015}")
    return PDF_33374015


@pytest.fixture
def pdf_34426522():
    """土耳其人群群体遗传学论文 PDF 路径（无目标变异）。

    用途：测试目标变异不存在时的行为（变异提及=False）。
    预期：PDF 中不包含 c.9999G>A / p.X9999X 等测试用变异。
    """
    if not PDF_34426522.exists():
        pytest.skip(f"PDF 未找到: {PDF_34426522}")
    return PDF_34426522