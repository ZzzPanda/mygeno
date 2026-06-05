"""
Variant keyword building and matching functions.
"""

import re

from .xml_parser import AA_3TO1, AA_1TO3, KNOWN_VARIANT_NAMES, split_sentences


def _build_descriptive_matches(target_cdna, target_protein):
    """
    为缺失、插入、错义等变异类型构建描述性短语的正则模式列表。
    例如 "deletion of a serine at amino acid 189", "3bp deletion at codon 565"
    返回 list of (regex_pattern, description) tuples。
    """
    patterns = []
    cdna_clean = (target_cdna or "").replace("c.", "").replace(" ", "")
    protein_clean = (target_protein or "").replace("p.", "").replace(" ", "")

    # 提取 cDNA 位置
    cdna_pos = ""
    m = re.search(r'(\d+)', cdna_clean)
    if m:
        cdna_pos = m.group(1)

    # 提取蛋白位置
    protein_pos = ""
    protein_aa = ""
    if protein_clean:
        m = re.match(r'([A-Za-z]{1,3})(\d+)', protein_clean)
        if m:
            protein_aa = m.group(1)
            protein_pos = m.group(2)
            # 扩展三字母
            if len(protein_aa) == 1:
                protein_aa3 = AA_1TO3.get(protein_aa, protein_aa)
            else:
                protein_aa3 = protein_aa

    # 缺失 (del)
    if "del" in cdna_clean or "del" in protein_clean:
        if protein_pos:
            patterns.append((
                rf"deletion\s+of\s+(?:a\s+)?{protein_aa3}\s*(?:residue\s*)?at\s+(?:amino\s*acid\s*|position\s*|residue\s*)?"
                rf"(?:codon\s*)?(?:Ser)?\s*{protein_pos}",
                f"Ser{protein_pos}del"
            ))
            # 也匹配简单数字形式
            patterns.append((
                rf"deletion\s+of\s+[A-Za-z]+\s+at\s+(?:amino\s*acid\s*|position\s*|residue\s*)?"
                rf"\s*{protein_pos}",
                f"p.{protein_pos}del"
            ))
        if cdna_pos:
            patterns.append((
                rf"(?:3\s*bp|3bp|three\s*bp|in-frame)\s+deletion\s+.*?{cdna_pos}",
                f"c.{cdna_pos}del"
            ))

    # 无义 (nonsense/stop)
    if any(t in cdna_clean.lower() for t in ["ter", "stop", "*", "x"]) or \
       any(t in protein_clean for t in ["Ter", "Stop", "*", "X"]):
        if protein_pos:
            patterns.append((
                rf"(?:nonsense|stop)\s+(?:mutation|variant)\s+at\s+(?:amino\s*acid\s*|position\s*)?\s*{protein_pos}",
                f"p.{protein_pos}*"
            ))

    # 错义 (missense)
    if ">" in cdna_clean or (protein_clean and not any(t in protein_clean for t in ["del", "Ter", "Stop", "*", "X", "fs"])):
        if protein_pos and protein_aa:
            patterns.append((
                rf"(?:amino\s*acid\s*)?substitution\s+at\s+(?:amino\s*acid\s*|position\s*|residue\s*)?\s*{protein_pos}",
                f"p.{protein_pos}"
            ))

    return patterns


def expand_protein_keywords(protein_str):
    """
    给定蛋白变异字符串（如 p.Arg389His 或 p.R389H），
    生成所有可能的关键词变体，用于在文献中匹配。

    返回 list[str]，包含所有应搜索的关键词。
    """
    if not protein_str:
        return []

    keywords = set()
    keywords.add(protein_str)  # 原始输入

    # 去除 "p." 前缀
    bare = protein_str
    if bare.startswith("p."):
        bare = bare[2:]
    keywords.add(bare)

    # 解析三字母格式: Arg389His / Arg389Ter / Arg382Stop
    m3 = re.match(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|Stop|\*|X)', bare)
    if m3:
        ref3, pos, alt3 = m3.groups()
        ref1 = AA_3TO1.get(ref3, ref3[0])
        alt1 = AA_3TO1.get(alt3, alt3[0] if alt3 not in ("Ter", "Stop", "*", "X") else "*")

        # 单字母格式: R389H
        keywords.add(f"{ref1}{pos}{alt1}")
        # p. 前缀单字母: p.R389H
        keywords.add(f"p.{ref1}{pos}{alt1}")
        # 无p.三字母: Arg389His
        keywords.add(f"{ref3}{pos}{alt3}")
        # 带空格: R 389 H
        keywords.add(f"{ref1} {pos} {alt1}")
        # 星号/终止密码子变体
        if alt3 in ("Ter", "Stop", "*", "X"):
            keywords.add(f"{ref1}{pos}*")
            keywords.add(f"{ref1}{pos}Ter")
            keywords.add(f"{ref1}{pos}X")
            keywords.add(f"{ref1}{pos}Stop")
            keywords.add(f"p.{ref1}{pos}*")
            keywords.add(f"p.{ref1}{pos}Ter")
            keywords.add(f"p.{ref1}{pos}X")
            keywords.add(f"p.{ref1}{pos}Stop")
            keywords.add(f"{ref3}{pos}*")
            keywords.add(f"{ref3}{pos}X")
            keywords.add(f"{ref3}{pos}Stop")

    # 解析单字母格式: R389H
    m1 = re.match(r'([A-Z*])(\d+)([A-Z*])', bare)
    if m1:
        ref1, pos, alt1 = m1.groups()
        ref3 = AA_1TO3.get(ref1, ref1)
        alt3 = AA_1TO3.get(alt1, alt1)
        # 三字母格式
        if len(ref3) == 3 and len(alt3) == 3:
            keywords.add(f"{ref3}{pos}{alt3}")
            keywords.add(f"p.{ref3}{pos}{alt3}")

    return [k for k in keywords if k]


def build_variant_keywords(cdna_str, protein_str, transcript_str=""):
    """
    构建完整的变异搜索关键词列表。

    策略：
    1. 精确 cDNA: c.1166G>A
    2. 精简 cDNA: 1166G>A, c1166G>A
    3. cDNA 数字: 1166（用于转录本版本差异匹配）
    4. 蛋白所有变体: p.Arg389His, p.R389H, R389H, Arg389His 等
    5. 转录本+变异组合（如果有转录本信息）

    返回 dict: {"exact": [...], "fuzzy": [...], "protein": [...]}
    """
    exact = set()
    fuzzy = set()
    protein_keywords = set()

    # --- cDNA 关键词 ---
    if cdna_str:
        exact.add(cdna_str)
        # 去除 c. 前缀: 1166G>A
        bare_cdna = cdna_str.replace("c.", "")
        fuzzy.add(bare_cdna)
        # 去掉空格: c1166G>A
        no_dot = cdna_str.replace("c.", "c")
        fuzzy.add(no_dot)
        # 提取纯数字: 1166（用于跨转录本匹配）
        num_m = re.search(r'c\.(\d+)', cdna_str)
        if num_m:
            fuzzy.add(num_m.group(1))

    # --- 蛋白关键词 ---
    if protein_str:
        protein_keywords.update(expand_protein_keywords(protein_str))

    # --- 描述性短语模式 ---
    descriptive_patterns = _build_descriptive_matches(cdna_str, protein_str)

    # --- 历史命名 ---
    # 基于基因名查找已知历史命名
    historical_names = []
    # 从 cdna_str 或 protein_str 推断基因名（调用方传入 gene）
    # 这里我们直接检查所有已知基因
    for gene_symbol, variants in KNOWN_VARIANT_NAMES.items():
        for hist_name, hgvs in variants.items():
            # 检查 HGVS 是否与当前目标变异匹配
            if (cdna_str and hgvs.get("cdna", "").replace(" ", "") == cdna_str.replace(" ", "")) or \
               (protein_str and hgvs.get("protein", "").replace(" ", "") == protein_str.replace(" ", "")):
                historical_names.append(hist_name)

    return {
        "exact": sorted(exact),
        "fuzzy": sorted(fuzzy),
        "protein": sorted(protein_keywords),
        "descriptive": descriptive_patterns,
        "historical": sorted(set(historical_names)),
        "all": sorted(exact | fuzzy | protein_keywords),
    }


def check_variant_match(text, keywords):
    """
    检查文本中是否包含目标变异的任何关键词变体。
    返回 (matched: bool, matched_keywords: list[str])
    """
    text_lower = text.lower()
    matched = []
    for kw in keywords["all"]:
        if kw.lower() in text_lower:
            matched.append(kw)
    # 描述性短语
    for pattern, desc in keywords.get("descriptive", []):
        if re.search(pattern, text, re.IGNORECASE):
            matched.append(desc)
    # 历史命名
    for hist_name in keywords.get("historical", []):
        if hist_name.lower() in text_lower:
            matched.append(hist_name)
    return len(matched) > 0, matched


def find_variant_sentences(full_text, keywords):
    """
    找出所有提及目标变异的句子。

    策略：
    1. 先匹配精确关键词（cDNA精确匹配、蛋白变体）
    2. 如果精确匹配不够，尝试"数字+氨基酸位置"的模糊匹配（转录本版本差异）
    3. 匹配描述性短语（如 "deletion of serine at position 189"）
    4. 匹配历史命名（如 "G6PD Tsukui"）
    """
    sentences = split_sentences(full_text)
    variant_sentences = []
    seen = set()
    all_matched = []

    for s in sentences:
        s_lower = s.lower()
        matched_here = []

        # 精确匹配：cDNA 关键词
        for kw in keywords["exact"]:
            if kw.lower() in s_lower:
                matched_here.append(kw)
        for kw in keywords["fuzzy"]:
            # fuzzy cDNA 需要确认附近有 c. 或基因名上下文
            if kw.lower() in s_lower:
                if re.search(r'\bc[\.<]', s, re.IGNORECASE) or re.search(r'[A-ZTCG]\s*>\s*[A-ZTCG*]', s):
                    matched_here.append(kw)

        # 蛋白关键词匹配
        for kw in keywords["protein"]:
            if kw.lower() in s_lower:
                matched_here.append(kw)

        # 描述性短语匹配
        for pattern, desc in keywords.get("descriptive", []):
            if re.search(pattern, s, re.IGNORECASE):
                matched_here.append(desc)

        # 历史命名匹配
        for hist_name in keywords.get("historical", []):
            if hist_name.lower() in s_lower:
                matched_here.append(hist_name)

        if matched_here:
            if s not in seen:
                seen.add(s)
                variant_sentences.append(s)
                all_matched.extend(matched_here)

    # 模糊匹配（转录本版本差异）：如果精确匹配不够，尝试用蛋白位置 + cDNA数字搜索
    if len(variant_sentences) < 3:
        # 提取 cDNA 中的数字部分
        cdna_nums = set()
        for kw in keywords["fuzzy"]:
            num_m = re.match(r'^(\d+)', kw)
            if num_m:
                cdna_nums.add(num_m.group(1))

        for num in cdna_nums:
            # 搜索 "c." + 附近数字 + 基因上下文
            pattern = rf'c\.\s*{num}\s*[\w><\-\+]+'
            for s in sentences:
                if s not in seen and re.search(pattern, s):
                    seen.add(s)
                    variant_sentences.append(s)
                    all_matched.append(f"c.{num}*")

    return variant_sentences, list(set(all_matched))


def infer_variant_type(cdna, protein, sentences):
    """推断变异类型。"""
    types_found = []
    if cdna:
        if re.search(r'[A-ZTCG*]\s*>\s*[A-ZTCG*]', cdna):
            if protein and re.search(r'(Ter|\*|X\b)', protein):
                types_found.append("无义突变 (nonsense)")
            elif protein and "fs" in protein:
                types_found.append("移码突变 (frameshift)")
            elif protein and "Met1" in protein:
                types_found.append("起始密码子丢失 (start loss)")
            else:
                types_found.append("错义突变 (missense)")
        if re.search(r'del', cdna, re.IGNORECASE):
            if "fs" in (protein or ""):
                types_found.append("移码突变 (frameshift)")
            elif re.search(r'_\d+del', cdna) and '>[A-Z]' not in cdna:
                types_found.append("缺失 (deletion)")
        if re.search(r'ins', cdna, re.IGNORECASE):
            types_found.append("插入 (insertion)")
        if re.search(r'dup', cdna, re.IGNORECASE) and not re.search(r'del|ins|>', cdna):
            types_found.append("重复 (duplication)")
        if re.search(r'[-+]\d+', cdna) and '>' in cdna:
            types_found.append("剪接位点突变 (splice site)")
        if re.search(r'=', cdna):
            types_found.append("同义突变 (silent)")

    if not types_found:
        context = " ".join(sentences).lower()
        type_patterns = [
            ("无义突变 (nonsense)", r'nonsense|stop.*(gain|codon)|premature.*stop'),
            ("剪接位点突变 (splice site)", r'splice.*(site|donor|acceptor)|splicing.*defect'),
            ("移码突变 (frameshift)", r'frameshift|fs\*|fs\.|p\.[A-Z][a-z]{2}\d+fs'),
            ("起始密码子丢失 (start loss)", r'start.*(loss|lost)|initiation.*codon|met1[>\s]'),
            ("同义突变 (silent)", r'silent|synonymous'),
        ]
        for vtype, pattern in type_patterns:
            if re.search(pattern, context):
                types_found.append(vtype)

    return ", ".join(types_found) if types_found else "未指明"