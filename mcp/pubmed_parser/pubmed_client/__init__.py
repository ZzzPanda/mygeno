"""
PubMed Variant Extractor package - split from pubmed_client.py
"""

# Import all public functions/classes from all modules
from .constants import (
    PHENOTYPE_MAP,
    DISEASE_SUBTYPE_MAP,
    LAB_FINDINGS_MAP,
    TOPIC_KEYWORD_MAP,
    MAX_RETRIES,
    DAILY_SITE_LIMIT,
    PAUSE_AFTER_N,
    EXTRA_PAUSE_SECS,
    MIN_DELAY,
    MAX_DELAY,
    _get_state_file,
    load_site_counts,
    save_site_counts,
    check_and_increment,
    random_sleep,
    periodic_long_pause,
)

from .xml_parser import (
    AA_3TO1,
    AA_1TO3,
    KNOWN_VARIANT_NAMES,
)

from .network import (
    safe_request,
    fetch_ncbi_abstract,
    fetch_europe_pmc_fulltext,
    fetch_pmc_fulltext,
    fetch_europe_pmc_text,
    fetch_pmc_html,
    fetch_ncbi_gene_info,
)

from .xml_parser import (
    _empty_result,
    parse_ncbi_xml,
    parse_pmc_europe_xml,
    _extract_tables_from_xml,
    parse_pmc_html,
    _extract_tables_from_html,
    split_sentences,
)

from .variants import (
    build_variant_keywords,
    expand_protein_keywords,
    _build_descriptive_matches,
    check_variant_match,
    find_variant_sentences,
    infer_variant_type,
)

from .extractors import (
    extract_pathogenicity,
    extract_zygosity,
    extract_inheritance,
    extract_patient_phenotypes,
    extract_clinical_details,
    extract_co_variants,
    extract_phase_evidence,
    _extract_parental_variant,
    _extract_parental_details,
    extract_patient_count,
    extract_variant_features,
    extract_info_for_variant,
    _parse_patient_sentence,
    _parse_table_row_for_patient,
    _extract_pathogenicity_from_table,
    _extract_table_features,
)

from .summarizer import (
    generate_summary_paragraph,
    _generate_literature_summary,
    generate_one_sentence_summary,
    _generate_study_background,
    _generate_excel_csv,
)

from .cli import main

# Public API
__all__ = [
    # Constants
    "PHENOTYPE_MAP",
    "DISEASE_SUBTYPE_MAP",
    "LAB_FINDINGS_MAP",
    "TOPIC_KEYWORD_MAP",
    "AA_3TO1",
    "AA_1TO3",
    "KNOWN_VARIANT_NAMES",
    "MAX_RETRIES",
    "DAILY_SITE_LIMIT",
    "PAUSE_AFTER_N",
    "EXTRA_PAUSE_SECS",
    "MIN_DELAY",
    "MAX_DELAY",
    # Network functions
    "fetch_ncbi_abstract",
    "fetch_europe_pmc_fulltext",
    "fetch_pmc_fulltext",
    "fetch_europe_pmc_text",
    "fetch_pmc_html",
    "fetch_ncbi_gene_info",
    "safe_request",
    # Site count/delay functions
    "random_sleep",
    "periodic_long_pause",
    "check_and_increment",
    "load_site_counts",
    "save_site_counts",
    # XML/HTML parsing
    "split_sentences",
    "parse_ncbi_xml",
    "parse_pmc_europe_xml",
    "_extract_tables_from_xml",
    "parse_pmc_html",
    "_extract_tables_from_html",
    # Variant functions
    "build_variant_keywords",
    "expand_protein_keywords",
    "check_variant_match",
    "find_variant_sentences",
    "infer_variant_type",
    # Extraction functions
    "extract_pathogenicity",
    "extract_zygosity",
    "extract_inheritance",
    "extract_patient_phenotypes",
    "extract_clinical_details",
    "extract_co_variants",
    "extract_phase_evidence",
    "extract_patient_count",
    "extract_variant_features",
    "extract_info_for_variant",
    # Summary generation
    "generate_summary_paragraph",
    "generate_one_sentence_summary",
    "_generate_excel_csv",
    # Main
    "main",
    # Internal helpers (included as requested)
    "_empty_result",
    "_get_state_file",
    "_build_descriptive_matches",
    "_extract_parental_variant",
    "_extract_parental_details",
    "_parse_patient_sentence",
    "_parse_table_row_for_patient",
    "_extract_pathogenicity_from_table",
    "_extract_table_features",
    "_generate_literature_summary",
    "_generate_study_background",
]