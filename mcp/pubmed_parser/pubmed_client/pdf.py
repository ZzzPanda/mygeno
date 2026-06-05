"""
PDF text and table extraction using pymupdf.
"""

import pymupdf


def extract_text_from_pdf(pdf_path: str) -> str:
    """
    Extract full text from PDF using pymupdf.

    Args:
        pdf_path: Path to the PDF file.

    Returns:
        Full text content of the PDF.
    """
    text_parts = []
    try:
        with pymupdf.open(pdf_path) as doc:
            for page in doc:
                page_text = page.get_text()
                if page_text:
                    text_parts.append(page_text)
    except Exception:
        return ""

    return "\n".join(text_parts)


def extract_tables_from_pdf(pdf_path: str) -> list[dict]:
    """
    Extract tables from PDF using pymupdf find_tables() method.

    Args:
        pdf_path: Path to the PDF file.

    Returns:
        List of dicts with keys: "caption", "rows" (list of lists).
        Each table dict contains the table caption and rows as lists of cell strings.
    """
    tables = []
    try:
        with pymupdf.open(pdf_path) as doc:
            for page_num, page in enumerate(doc):
                # Use find_tables() to detect tables on this page
                text_page = page.get_text("text")
                tables_found = page.find_tables()
                if tables_found and tables_found.tables:
                    for table in tables_found.tables:
                        # Extract table bbox
                        bbox = table.bbox
                        # Get caption from above the table (look for text near table)
                        caption = ""
                        # Try to extract caption from nearby text
                        # The table object has a header property
                        header = table.header
                        # Build rows from table data
                        rows = []
                        for row in table.extract():
                            # row is a list of cell contents
                            cleaned_row = []
                            for cell in row:
                                if cell is None:
                                    cleaned_row.append("")
                                else:
                                    cleaned_row.append(str(cell).strip())
                            if any(c for c in cleaned_row):  # Skip fully empty rows
                                rows.append(cleaned_row)

                        if rows:
                            tables.append({
                                "id": f"page_{page_num + 1}",
                                "caption": caption,
                                "rows": rows,
                            })

                # Also try extracting tables using the legacy method if no tables found
                if not tables or not tables_found.tables:
                    # Fallback: try extracting tables from blocks
                    pass

    except Exception as e:
        raise RuntimeError(f"Failed to extract tables from PDF: {e}")

    return tables