"""
PDF variant extraction test suite.

Structure:
  Unit tests — test each function in isolation with controlled inputs
  Integration tests — test the full pipeline with real PDFs

Run:
  .venv/bin/python -m pytest tests/test_pdf_parser.py -v
  .venv/bin/python -m pytest tests/test_pdf_parser.py -k "keyword" -v   # run subset
"""

import sys
from pathlib import Path

# Add package root so unit tests can `from pubmed_client import ...`
_PACKAGE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_PACKAGE_ROOT))

import json
import pytest

# ── Paths ───────────────────────────────────────────────────────────────────────

DATA_DIR = Path(__file__).resolve().parent.parent.parent.parent / "data" / "pdf"
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
EXPECTED_FILE = DATA_DIR / "variant_results.json"

PDF_33374015 = DATA_DIR / "PMID:33374015.pdf"
PDF_34426522 = DATA_DIR / "PMID:34426522.pdf"

# ── Fixtures ───────────────────────────────────────────────────────────────────

@pytest.fixture
def pdf_33374015():
    """Path to CFTR T760M paper PDF."""
    if not PDF_33374015.exists():
        pytest.skip(f"PDF not found: {PDF_33374015}")
    return PDF_33374015

@pytest.fixture
def pdf_34426522():
    """Path to Turkish population genetics paper PDF."""
    if not PDF_34426522.exists():
        pytest.skip(f"PDF not found: {PDF_34426522}")
    return PDF_34426522

@pytest.fixture
def expected_results():
    """Load expected results from JSON file."""
    if not EXPECTED_FILE.exists():
        pytest.skip(f"Expected results not found: {EXPECTED_FILE}")
    with open(EXPECTED_FILE, encoding="utf-8") as f:
        return {r["文件"]: r for r in json.load(f)}


# ═══════════════════════════════════════════════════════════════════════════════
# UNIT TESTS — pdf.py
# ═══════════════════════════════════════════════════════════════════════════════

class TestExtractTextFromPDF:
    """test extract_text_from_pdf"""

    def test_returns_string(self, pdf_33374015):
        from pubmed_client import extract_text_from_pdf
        text = extract_text_from_pdf(str(pdf_33374015))
        assert isinstance(text, str)
        assert len(text) > 1000, f"Text too short ({len(text)} chars)"

    def test_contains_gene_name(self, pdf_33374015):
        from pubmed_client import extract_text_from_pdf
        text = extract_text_from_pdf(str(pdf_33374015))
        assert "CFTR" in text or "cftr" in text.lower()

    def test_contains_variant_keyword(self, pdf_33374015):
        from pubmed_client import extract_text_from_pdf
        text = extract_text_from_pdf(str(pdf_33374015))
        # T760M paper should mention Thr760Met or p.Thr760Met
        assert "760" in text, "Position 760 not found in text"

    def test_empty_for_nonexistent_file(self):
        from pubmed_client import extract_text_from_pdf
        result = extract_text_from_pdf("/nonexistent/file.pdf")
        assert result == ""


class TestExtractTablesFromPDF:
    """test extract_tables_from_pdf"""

    def test_returns_list(self, pdf_33374015):
        from pubmed_client import extract_tables_from_pdf
        tables = extract_tables_from_pdf(str(pdf_33374015))
        assert isinstance(tables, list)

    def test_table_has_rows(self, pdf_33374015):
        from pubmed_client import extract_tables_from_pdf
        tables = extract_tables_from_pdf(str(pdf_33374015))
        assert len(tables) >= 1, "Should find at least 1 table"
        assert "rows" in tables[0]
        assert len(tables[0]["rows"]) >= 1

    def test_table_row_cells_are_strings(self, pdf_33374015):
        from pubmed_client import extract_tables_from_pdf
        tables = extract_tables_from_pdf(str(pdf_33374015))
        row = tables[0]["rows"][0]
        for cell in row:
            assert isinstance(cell, str), f"Cell should be str, got {type(cell)}"


# ═══════════════════════════════════════════════════════════════════════════════
# UNIT TESTS — variants.py
# ═══════════════════════════════════════════════════════════════════════════════

class TestBuildVariantKeywords:
    """test build_variant_keywords"""

    def test_cdna_only(self):
        from pubmed_client import build_variant_keywords
        kw = build_variant_keywords("c.1166G>A", "", "")
        assert "c.1166G>A" in kw["exact"]
        assert "1166G>A" in kw["fuzzy"]
        assert "1166" in kw["fuzzy"]

    def test_protein_p_thr760met(self):
        from pubmed_client import build_variant_keywords
        kw = build_variant_keywords("", "p.Thr760Met", "")
        assert "p.Thr760Met" in kw["protein"]
        assert "Thr760Met" in kw["protein"]
        assert "T760M" in kw["protein"] or "T 760 M" in kw["protein"]

    def test_protein_without_p_prefix(self):
        from pubmed_client import build_variant_keywords
        kw = build_variant_keywords("", "Thr760Met", "")
        assert "Thr760Met" in kw["protein"]

    def test_protein_three_letter_to_one(self):
        from pubmed_client import build_variant_keywords
        kw = build_variant_keywords("", "p.Arg389His", "")
        # Should expand to R389H
        assert any("R389" in k for k in kw["protein"])

    def test_nonsense_ter(self):
        from pubmed_client import build_variant_keywords
        kw = build_variant_keywords("", "p.Arg389Ter", "")
        assert any("Ter" in k or "*" in k for k in kw["protein"])

    def test_fs_variant(self):
        from pubmed_client import build_variant_keywords
        kw = build_variant_keywords("", "p.Ala411fs", "")
        assert "fs" in kw["protein"]


class TestFindVariantSentences:
    """test find_variant_sentences"""

    def test_exact_match(self):
        from pubmed_client import build_variant_keywords, find_variant_sentences
        text = "The patient carried the c.2279C>T variant. Another sentence."
        kw = build_variant_keywords("c.2279C>T", "", "")
        sentences, matched = find_variant_sentences(text, kw)
        assert len(sentences) >= 1
        assert "c.2279C>T" in matched or "2279" in matched

    def test_protein_match(self):
        from pubmed_client import build_variant_keywords, find_variant_sentences
        text = "The p.Thr760Met variant was found in patient 1."
        kw = build_variant_keywords("", "p.Thr760Met", "")
        sentences, _ = find_variant_sentences(text, kw)
        assert len(sentences) >= 1

    def test_no_match(self):
        from pubmed_client import build_variant_keywords, find_variant_sentences
        text = "This is unrelated text about another gene."
        kw = build_variant_keywords("c.9999G>A", "", "")
        sentences, _ = find_variant_sentences(text, kw)
        assert len(sentences) == 0

    def test_fuzzy_cdna_match(self):
        from pubmed_client import build_variant_keywords, find_variant_sentences
        # Text uses c.2279C>T but keyword has c. prefix
        text = "Mutation 2279C>T was detected."
        kw = build_variant_keywords("c.2279C>T", "", "")
        sentences, _ = find_variant_sentences(text, kw)
        assert len(sentences) >= 1


class TestInferVariantType:
    """test infer_variant_type"""

    @pytest.mark.parametrize("cdna,protein,expected", [
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
        from pubmed_client import infer_variant_type
        result = infer_variant_type(cdna, protein, [])
        assert result == expected, f"Expected {expected}, got {result}"

    def test_missense_from_protein_only(self):
        from pubmed_client import infer_variant_type
        # Without cDNA, missense is inferred from protein
        result = infer_variant_type("", "p.Thr760Met", [])
        assert result == "错义突变 (missense)"

    def test_nonsense_from_protein_only(self):
        from pubmed_client import infer_variant_type
        result = infer_variant_type("", "p.Arg389Ter", [])
        assert result == "无义突变 (nonsense)"

    def test_frameshift_from_protein_only(self):
        from pubmed_client import infer_variant_type
        result = infer_variant_type("", "p.Ala411fs", [])
        assert result == "移码突变 (frameshift)"

    def test_no_false_positive_from_other_variants_in_text(self):
        from pubmed_client import infer_variant_type
        # Text contains p.Phe861Leufs*3 but we search for p.Thr760Met
        bad_text = ["p.Phe861Leufs*3 was found in another patient."]
        result = infer_variant_type("", "p.Thr760Met", bad_text)
        assert result == "错义突变 (missense)", f"Got '{result}' — fs*3 in text should not affect missense"


# ═══════════════════════════════════════════════════════════════════════════════
# UNIT TESTS — extractors.py
# ═══════════════════════════════════════════════════════════════════════════════

class TestExtractPathogenicity:
    """test extract_pathogenicity"""

    def test_pathogenic_keyword(self):
        from pubmed_client import extract_pathogenicity
        sentences = ["This variant is pathogenic."]
        result = extract_pathogenicity(sentences, "")
        assert "pathogenic" in result.lower()

    def test_vus(self):
        from pubmed_client import extract_pathogenicity
        sentences = ["The variant is a VUS."]
        result = extract_pathogenicity(sentences, "")
        assert "VUS" in result or "不确定" in result

    def test_empty_sentences(self):
        from pubmed_client import extract_pathogenicity
        result = extract_pathogenicity([], "no variant mentioned here")
        assert result != ""


class TestExtractZygosity:
    """test extract_zygosity"""

    def test_compound_heterozygous(self):
        from pubmed_client import extract_zygosity
        sentences = ["The patient had compound heterozygous mutations."]
        result = extract_zygosity(sentences, "", target_cdna="c.2279C>T", target_protein="p.Thr760Met")
        assert "compound" in result.lower() or "复合杂合" in result

    def test_homozygous(self):
        from pubmed_client import extract_zygosity
        sentences = ["Homozygous variant was detected."]
        result = extract_zygosity(sentences, "")
        assert "hom" in result.lower() or "纯合" in result

    def test_heterozygous(self):
        from pubmed_client import extract_zygosity
        sentences = ["Heterozygous for the variant."]
        result = extract_zygosity(sentences, "")
        assert "het" in result.lower() or "杂合" in result

    def test_in_trans_corects_to_compound(self):
        from pubmed_client import extract_zygosity
        # "in trans" with homozygous should become compound heterozygous
        sentences = [
            "The variant was homozygous. We could not formally confirm the in trans position."
        ]
        result = extract_zygosity(sentences, "", target_cdna="c.2279C>T")
        assert "compound" in result.lower() or "复合杂合" in result


class TestExtractInheritance:
    """test extract_inheritance"""

    def test_autosomal_recessive(self):
        from pubmed_client import extract_inheritance
        sentences = ["Inheritance is autosomal recessive."]
        result = extract_inheritance(sentences, "")
        assert "隐性" in result or "recessive" in result.lower()

    def test_autosomal_dominant(self):
        from pubmed_client import extract_inheritance
        sentences = ["Autosomal dominant inheritance."]
        result = extract_inheritance(sentences, "")
        assert "显性" in result or "dominant" in result.lower()

    def test_de_novo(self):
        from pubmed_client import extract_inheritance
        sentences = ["De novo mutation."]
        result = extract_inheritance(sentences, "")
        assert "novo" in result.lower()


class TestExtractPatientPhenotypes:
    """test extract_patient_phenotypes"""

    def test_phenotype_from_keyword(self):
        from pubmed_client import extract_patient_phenotypes
        sentences = ["Patient had retinitis pigmentosa (RP)."]
        result = extract_patient_phenotypes(sentences)
        assert len(result) >= 1


# ═══════════════════════════════════════════════════════════════════════════════
# INTEGRATION TESTS — server.py
# ═══════════════════════════════════════════════════════════════════════════════

class TestExtractVariantInfo:
    """Integration test: full pipeline from PDF"""

    def test_cftr_t760m_variant_mention(self, pdf_33374015):
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
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
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from server import extract_variant_info

        result = extract_variant_info(
            pdf_path=str(pdf_33374015),
            cdna="c.2279C>T",
            protein="p.Thr760Met",
            gene="CFTR",
        )
        assert result["变异类型"] == "错义突变 (missense)", f"Got {result['变异类型']}"

    def test_cftr_t760m_zygosity_compound_heterozygous(self, pdf_33374015):
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from server import extract_variant_info

        result = extract_variant_info(
            pdf_path=str(pdf_33374015),
            cdna="c.2279C>T",
            protein="p.Thr760Met",
            gene="CFTR",
        )
        assert "复合杂合" in result["合子状态"] or "compound" in result["合子状态"].lower()

    def test_paper_without_variant(self, pdf_34426522):
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from server import extract_variant_info

        result = extract_variant_info(
            pdf_path=str(pdf_34426522),
            cdna="c.9999G>A",
            protein="p.X9999X",
            gene="CFTR",
        )
        assert result["变异提及"] is False


class TestAnalyzeVariant:
    """Integration test: text-based variant analysis (no PDF)"""

    def test_analyze_variant_from_text(self):
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from server import analyze_variant

        text = """
        We identified a pathogenic c.1166G>A (p.Arg389His) variant in
        compound heterozygous state. The patient was homozygous for
        the second variant. Inheritance is autosomal recessive.
        """
        result = analyze_variant(text, cdna="c.1166G>A", protein="p.Arg389His", gene="GENE")
        assert result["变异提及"] is True
        assert "pathogenic" in result["致病性"].lower()


# ═══════════════════════════════════════════════════════════════════════════════
# INTEGRATION TEST — against expected JSON (informational)
# ═══════════════════════════════════════════════════════════════════════════════

class TestExpectedOutput:
    """
    Informational comparison against data/pdf/variant_results.json.
    These are NOT strict pass/fail — they print the diff for review.
    """

    def test_33374015_vs_expected(self, pdf_33374015, expected_results):
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from server import extract_variant_info

        expected = expected_results.get("PMID:33374015.pdf", {})

        # Build cdna/protein from matched keywords in expected
        cdna, protein = "", ""
        for kw in expected.get("匹配关键词", []):
            if not cdna and kw.startswith("c."): cdna = kw
            if not protein and kw.startswith("p."): protein = kw

        result = extract_variant_info(
            str(pdf_33374015), cdna=cdna, protein=protein, gene=expected.get("基因", "")
        )

        # Print a readable diff
        print("\n" + "=" * 60)
        print("PMID:33374015.pdf — actual vs expected")
        print("=" * 60)
        for field, exp_val in expected.items():
            if exp_val in ("", None): continue
            act_val = result.get(field)
            match = "✓" if exp_val == act_val else "✗"
            print(f"  {match} {field}:")
            print(f"      expected: {str(exp_val)[:80]}")
            print(f"      actual:   {str(act_val)[:80]}")

    def test_34426522_vs_expected(self, pdf_34426522, expected_results):
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from server import extract_variant_info

        expected = expected_results.get("PMID:34426522.pdf", {})

        result = extract_variant_info(
            str(pdf_34426522), cdna="", protein="", gene=expected.get("基因", "")
        )

        print("\n" + "=" * 60)
        print("PMID:34426522.pdf — actual vs expected")
        print("=" * 60)
        for field, exp_val in expected.items():
            if exp_val in ("", None): continue
            act_val = result.get(field)
            match = "✓" if exp_val == act_val else "✗"
            print(f"  {match} {field}:")
            print(f"      expected: {str(exp_val)[:80]}")
            print(f"      actual:   {str(act_val)[:80]}")
