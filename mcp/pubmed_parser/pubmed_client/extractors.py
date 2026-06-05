"""
Extraction functions for variant information from text and tables.
"""

import re

from .constants import PHENOTYPE_MAP, DISEASE_SUBTYPE_MAP, LAB_FINDINGS_MAP
from .xml_parser import _empty_result
from .variants import find_variant_sentences, infer_variant_type


def extract_pathogenicity(sentences, full_text_lower):
    """提取致病性评级。"""
    vs_lower = " ".join(sentences).lower()
    path_terms = []
    pathogenicity_map = [
        ("致病 (pathogenic)", r'\bpathogenic\b'),
        ("可能致病 (likely pathogenic)", r'likely\s+pathogenic'),
        ("良性 (benign)", r'\bbenign\b'),
        ("可能良性 (likely benign)", r'likely\s+benign'),
        ("意义不明 (VUS)", r'variant\s+of\s+uncertain\s+significance|\bvus\b'),
        ("有害 (damaging)", r'\bdamaging\b'),
        ("有害 (deleterious)", r'\bdeleterious\b'),
        ("致病 (disease-causing)", r'disease.?causing'),
    ]
    for term, pattern in pathogenicity_map:
        if re.search(pattern, vs_lower) and term not in path_terms:
            path_terms.append(term)
    if not path_terms:
        for term, pattern in pathogenicity_map:
            m = re.search(pattern, vs_lower)
            if m:
                pos = m.start()
                context = vs_lower[max(0, pos - 300):pos + 300]
                if re.search(r'\bc[\.<]|p\.', context):
                    path_terms.append(term)
    return ", ".join(path_terms) if path_terms else "未指明"


def extract_zygosity(sentences, full_text_lower, target_cdna="", target_protein="", tables=None, keywords=None):
    """提取遗传方式（合子状态）。v8: 添加邻近度评分，防止不同变异的合子状态误匹配。

    策略：
    1. 仅在同时包含目标变异的句子中匹配合子关键词
    2. 如果句子提到其他变异（不同的 c./p.），降低权重
    3. 表格行中的 Het/Hom 标记优先
    4. 如果存在 "in trans with" 证据 → 排除纯合
    """
    zygosity_keywords = [
        ("复合杂合 (compound heterozygous)", r'compound\s*heterozyg'),
        ("纯合 (homozygous)", r'homozyg'),
        ("杂合 (heterozygous)", r'heterozyg'),
        ("双等位基因 (biallelic)", r'biallelic'),
        ("单等位基因 (monoallelic)", r'monoallelic'),
        ("半合子 (hemizygous)", r'hemizyg'),
    ]

    # 构建目标变异标识（用于区分"这个变异"和"其他变异"）
    target_cdna_short = (target_cdna or "").replace("c.", "").replace(" ", "").lower()
    target_protein_short = (target_protein or "").replace("p.", "").replace(" ", "").lower()

    def _has_target_variant_in_sentence(s):
        """检查句子是否实际提到目标变异（而不仅是其他变异）。"""
        s_lower = s.lower()
        s_clean = s_lower.replace(" ", "")
        # 检查 cDNA
        if target_cdna_short and target_cdna_short in s_clean:
            return True
        # 检查蛋白关键词
        if target_protein_short and target_protein_short in s_clean:
            return True
        # 检查扩展关键词
        if keywords:
            for kw in keywords.get("exact", []):
                if kw.lower() in s_lower:
                    return True
            for kw in keywords.get("fuzzy", []):
                if kw.lower() in s_lower:
                    # 确认附近有 cDNA 上下文
                    if re.search(r'\bc[\.<]', s, re.IGNORECASE) or re.search(r'[A-ZTCG]\s*>\s*[A-ZTCG*]', s):
                        return True
        return False

    def _count_other_variants_in_sentence(s):
        """计数句子中不同于目标变异的其他变异数量。"""
        s_clean = s.replace(" ", "").lower()
        count = 0
        # 找所有 c.X 格式变异
        other_cdnas = re.findall(r'c\.([\d\w_*><+=\-\+]+)', s_clean)
        for c in other_cdnas:
            c_norm = c.lower()
            if target_cdna_short and c_norm != target_cdna_short:
                count += 1
        # 找所有 p.X 格式变异
        other_prots = re.findall(r'p\.([\w\*?]+)', s_clean)
        for p in other_prots:
            p_norm = p.lower()
            if target_protein_short and p_norm != target_protein_short:
                count += 1
        return count

    # ---- 第1层：句子级别邻近匹配（最高权重） ----
    scored_matches = []
    for s in sentences:
        s_lower = s.lower()
        has_target = _has_target_variant_in_sentence(s)
        other_count = _count_other_variants_in_sentence(s)

        for z_name, z_pattern in zygosity_keywords:
            m = re.search(z_pattern, s_lower)
            if m:
                # 基础分
                score = 10
                # 句子包含目标变异 → +100
                if has_target:
                    score += 100
                # 句子有其他变异 → 扣分
                score -= other_count * 20
                # "in trans" 存在 → 这是复合杂合，不是纯合
                if re.search(r'in\s+trans\b', s_lower) and z_name.startswith("纯合"):
                    score -= 50  # 大幅降权
                if re.search(r'in\s+trans\b', s_lower) and z_name.startswith("复合杂合"):
                    score += 50  # 大幅加权
                # 如果同时有"纯合"和"in trans"，改为大幅加权复合杂合（而非降权纯合）
                if re.search(r'in\s+trans\b', s_lower) and z_name.startswith("纯合"):
                    score += 50  # in trans 否定纯合，应提升复合杂合优先级
                scored_matches.append((z_name, score, s, m.start()))
                break  # 每句只取第一个匹配的合子关键词

    # ---- 第2层：表格级别匹配（高权重） ----
    table_zygosity = None
    if tables:
        for table_info in tables:
            for ri, row in enumerate(table_info["rows"]):
                row_text = "\t".join(str(cell or "") for cell in row)
                row_lower = row_text.lower()
                # 该行是否包含目标变异（仅当keywords可用时检查）
                has_target = False
                if keywords:
                    for kw in keywords.get("all", []):
                        if kw.lower() in row_lower:
                            has_target = True
                            break
                    # 如果keywords存在但没有匹配，则跳过该行
                    if not has_target:
                        continue
                # 从该行提取合子状态
                for cell in row:
                    cell_text = str(cell or "").strip().lower()
                    if cell_text in ("het", "heterozygous", "heterozygote"):
                        table_zygosity = "杂合 (heterozygous)"
                    elif cell_text in ("hom", "homozygous", "homozygote"):
                        table_zygosity = "纯合 (homozygous)"
                    elif cell_text in ("compound het", "compound heterozygous", "compound heterozygote"):
                        table_zygosity = "复合杂合 (compound heterozygous)"

                # v10: 同表行多变异检测 —— 同一行有两个不同蛋白变异 → 复合杂合
                if not table_zygosity:
                    # 收集该行所有蛋白变异
                    row_prot_variants = set()
                    prot_re_1letter = re.compile(r'p\.\s*([A-Z]\d+[A-Z\*])')
                    prot_re_3letter = re.compile(r'p\.\s*([A-Z][a-z]{2}\d+[A-Z][a-z]{2})')
                    for cell in row:
                        cell_str = str(cell or "").strip().replace(" ", "")
                        for m in prot_re_1letter.finditer(cell_str):
                            row_prot_variants.add(m.group(1))
                        for m in prot_re_3letter.finditer(cell_str):
                            row_prot_variants.add(m.group(1))
                    # 如果有≥2个不同的蛋白变异，则为复合杂合
                    if len(row_prot_variants) >= 2:
                        table_zygosity = "复合杂合 (compound heterozygous)"

                # v11: 跨行变异补充检测 —— 下一行仅有变异（无患者ID）→ 复合杂合
                # 处理每行一个等位基因的表格格式（如 c.3385C>T 在一行，c.858+2T>A 在下一行）
                if not table_zygosity:
                    rows_list = table_info.get("rows", [])
                    if ri + 1 < len(rows_list):
                        next_row = rows_list[ri + 1]
                        next_first_cell = str(next_row[0] or "") if next_row else ""
                        next_row_text = "\t".join(str(cell or "") for cell in next_row)
                        # 下一行不包含患者 ID（无人口学信息），仅有变异
                        has_patient_id_next = bool(re.search(r'^\d+$', next_first_cell)) or \
                            bool(re.search(r'(?:F\d+|P\d+|OX\d+|R\d+)', next_row_text))
                        is_variant_row = bool(re.search(r'^(c\.|p\.|g\.|n\.)', next_first_cell)) or \
                            (not has_patient_id_next and re.search(r'(c\.[\w\.\-\+>]+|p\.\w+\d+\w+)', next_row_text))
                        if is_variant_row and not has_patient_id_next:
                            # 收集下一行的变异（蛋白 + cDNA，因为可能是剪接位点如c.858+2T>A）
                            next_prot_variants = set()
                            next_cdna_variants = set()
                            cdna_re = re.compile(r'(c\.[\w\.\-\+>]+)')
                            for cell in next_row:
                                cell_str = str(cell or "").strip().replace(" ", "")
                                for m in prot_re_1letter.finditer(cell_str):
                                    next_prot_variants.add(m.group(1))
                                for m in prot_re_3letter.finditer(cell_str):
                                    next_prot_variants.add(m.group(1))
                                for m in cdna_re.finditer(cell_str):
                                    next_cdna_variants.add(m.group(1))
                            # 检查下一行是否有不同于目标cDNA的变异
                            target_cdna_clean = (target_cdna or "").replace(" ", "")
                            has_different_allele = False
                            for cv in next_cdna_variants:
                                if cv.replace(" ", "") != target_cdna_clean:
                                    has_different_allele = True
                                    break
                            combined = row_prot_variants | next_prot_variants
                            if len(combined) >= 2 or has_different_allele:
                                table_zygosity = "复合杂合 (compound heterozygous)"

    # ---- 第3层：全文级别匹配（最低权重，有邻近限制） ----
    full_text_matches = []
    if not scored_matches:
        for z_name, z_pattern in zygosity_keywords:
            for m in re.finditer(z_pattern, full_text_lower):
                pos = m.start()
                # 检查前后 300 字符范围内是否有 cDNA 上下文
                context = full_text_lower[max(0, pos - 300):pos + 300]
                nearby_cdna = re.search(r'\bc[\.<\s]', context)
                nearby_protein = re.search(r'\bp\.', context)
                if nearby_cdna or nearby_protein:
                    # 再检查是否提到了目标变异
                    if target_cdna_short and target_cdna_short in context.replace(" ", ""):
                        score = 5
                    elif target_protein_short and target_protein_short in context.replace(" ", ""):
                        score = 5
                    else:
                        score = 1  # 附近有变异但不确定是目标变异
                    full_text_matches.append((z_name, score))
                    break

    # ---- 综合判定 ----
    zygosity_found = []

    # 表格标记优先（最精确）
    if table_zygosity:
        zygosity_found.append(table_zygosity)

    # 句子级别结果（按分数排序）
    seen_names = set(zygosity_found)
    for z_name, score, s, pos in sorted(scored_matches, key=lambda x: -x[1]):
        if z_name not in seen_names and score > 0:
            zygosity_found.append(z_name)
            seen_names.add(z_name)

    # 全文级别结果（仅在没有句子级别匹配时）
    if not zygosity_found:
        for z_name, score in sorted(full_text_matches, key=lambda x: -x[1]):
            if z_name not in seen_names:
                zygosity_found.append(z_name)
                seen_names.add(z_name)

    # "in trans" 修正：如果有 trans 证据但结果包含纯合，转为复合杂合
    has_trans_evidence = False
    for s in sentences:
        if re.search(r'in\s+trans\b', s, re.IGNORECASE):
            has_trans_evidence = True
            break
    if has_trans_evidence and any("纯合" in z for z in zygosity_found):
        zygosity_found = ["复合杂合 (compound heterozygous)" if "纯合" in z else z for z in zygosity_found]
        if "复合杂合 (compound heterozygous)" not in zygosity_found:
            zygosity_found.insert(0, "复合杂合 (compound heterozygous)")

    return ", ".join(zygosity_found) if zygosity_found else "未指明"


def extract_inheritance(sentences, full_text_lower):
    """提取遗传模式。"""
    vs_lower = " ".join(sentences).lower()
    inheritance_map = {
        "常染色体隐性遗传": r'autosomal\s+recessive',
        "常染色体显性遗传": r'autosomal\s+dominant',
        "X连锁遗传": r'x.?linked',
        "新发突变 (de novo)": r'de\s+novo',
    }
    inheritance_found = []
    for inh, pattern in inheritance_map.items():
        if re.search(pattern, vs_lower):
            inheritance_found.append(inh)
    if not inheritance_found:
        for inh, pattern in inheritance_map.items():
            if re.search(pattern, full_text_lower):
                inheritance_found.append(inh)
    return ", ".join(inheritance_found) if inheritance_found else "未指明"


def extract_patient_phenotypes(sentences, tables=None, keywords=None):
    """从变异相关句子和表格中提取表型。"""
    combined = " ".join(sentences)
    combined_lower = combined.lower()
    found = set()

    for kw, zh in PHENOTYPE_MAP.items():
        if kw.lower() in combined_lower:
            found.add(zh)

    # 从表格中提取（如果表格行包含变异关键词）
    if tables and keywords:
        for table_info in tables:
            for row in table_info["rows"]:
                row_text = "\t".join(str(cell or "") for cell in row)
                row_lower = row_text.lower()
                has_variant = any(kw.lower() in row_lower for kw in keywords["all"])
                if has_variant:
                    for pkw, zh in PHENOTYPE_MAP.items():
                        if pkw.lower() in row_lower:
                            found.add(zh)

    return sorted(found)


def extract_clinical_details(sentences, tables=None, keywords=None):
    """从变异相关句子中提取深层次临床细节。
    返回 dict 包含: disease_subtype, onset_age, lab_findings,
                  progression, variant_frequency, exon_domain
    """
    combined = " ".join(sentences)
    combined_lower = combined.lower()
    details = {}

    # 1. 疾病亚型
    subtypes = []
    for en_pattern, zh_label in DISEASE_SUBTYPE_MAP.items():
        if en_pattern in combined_lower and zh_label not in subtypes:
            subtypes.append(zh_label)
    if subtypes:
        details["disease_subtypes"] = subtypes

    # 2. 发病年龄
    onset_patterns = [
        r'(?:onset|presented|diagnosed|first\s+symptoms?|first\s+neurological|manifested)\s+at\s+(?:age\s+)?(\d+[\-–]\d*\s*(?:months?|years?|yrs?|月|岁))',
        r'(?:onset|presented)\s+(?:at|between)\s+(\d+[\-–]\d*\s*(?:months?|years?|yrs?))',
        r'(?:age\s+at\s+(?:onset|diagnosis|presentation))\s*(?:was|is|:)?\s*(\d+[\-–]?\d*\s*(?:months?|years?|yrs?))',
        r'(\d+[\-–]\d+)\s*(?:years?|yrs?)\s+(?:old|age)\s+(?:at\s+)?(?:onset|presentation)',
        r'(?:onset|presented)\s+at\s+(\d+)\s*(?:years?|yrs?|months?)',
    ]
    for pat in onset_patterns:
        m = re.search(pat, combined_lower)
        if m:
            age_raw = m.group(1).strip()
            age_cn = age_raw.replace('years', '岁').replace('year', '岁').replace('yrs', '岁').replace('months', '个月').replace('month', '个月')
            details["onset_age"] = f"首发年龄{age_cn}"
            break

    # 3. 实验室/生化发现
    lab_found = []
    for en_label, zh_label in LAB_FINDINGS_MAP.items():
        if en_label in combined_lower:
            ctx_pat = re.compile(
                r'[^.]*' + re.escape(en_label) + r'[^.]*(?:elevated|increased|reduced|decreased|absent|normal|positive|negative|abnormal|classic|variant|mild|severe|moderate|\d+\.?\d*\s*(?:ng|ug|mg|pg)\s*/?\s*[mud]?[lL])[^.]*\.',
                re.IGNORECASE
            )
            ctx_m = ctx_pat.search(combined)
            if ctx_m:
                ctx_text = ctx_m.group(0).strip()
                # Classify lab result direction
                if re.search(r'(?:elevated|increased|high)', ctx_text, re.IGNORECASE):
                    lab_found.append(f"{zh_label}水平升高")
                elif re.search(r'(?:reduced|decreased|low|absent)', ctx_text, re.IGNORECASE):
                    lab_found.append(f"{zh_label}水平降低")
                elif re.search(r'(?:classic|variant)', ctx_text, re.IGNORECASE):
                    direction = "经典型" if "classic" in ctx_text.lower() else "变异型"
                    lab_found.append(f"{zh_label}判定为{direction}")
                else:
                    lab_found.append(f"{zh_label}异常")
            else:
                lab_found.append(zh_label)
    if lab_found:
        details["lab_findings"] = lab_found[:4]

    # 4. 疾病进展
    prog_patterns = [
        (r'(progressive).{0,30}(deterioration|course|decline|neurological)', "进行性加重"),
        (r'(rapidly)\s*(progressive|progressing)', "快速进展"),
        (r'(slowly)\s*(progressive|progressing)', "缓慢进展"),
        (r'(stable)\s*(course|disease)', "病程稳定"),
    ]
    for pat, label in prog_patterns:
        if re.search(pat, combined_lower):
            details["progression"] = label
            break

    # 5. 变异频率/新发性
    freq_patterns = [
        (r'private\s+mutation', "私有突变"),
        (r'novel\s+(?:mutation|variant|missense)', "新发突变"),
        (r'recurrent\s+mutation', "复发突变"),
        (r'found\s+in\s+only\s+(?:one|a single)\s+(?:patient|case|individual|family)', "未在队列中重复出现"),
        (r'not\s+(?:present|found|detected|identified)\s+in\s+(?:other|the\s+rest\s+of|control)', "未在对照或其他患者中检出"),
        (r'(?:not\s+included|absent\s+from).{0,40}(?:high.?frequency|common|recurrent|mutation\s+list)', "未包含在高频突变列表中"),
        (r'(?:minor allele frequency|MAF).{0,20}(?:below|less than|<)\s*0\.0+1', "人群频率极低（MAF<0.01）"),
        (r'(?:first|initially)\s+(?:described|reported|identified)', "首次报道"),
    ]
    freq_items = []
    for pat, label in freq_patterns:
        if re.search(pat, combined_lower):
            if label not in freq_items:
                freq_items.append(label)
    if freq_items:
        details["frequency_info"] = freq_items

    # 6. 外显子/结构域
    exon_m = re.search(r'(?:exon\s*)(\d{1,2}[A-Za-z]?)', combined, re.IGNORECASE)
    if not exon_m:
        # "E9" 格式 — 仅匹配1-2位数字的缩写，且后跟非字母数字字符（避免误匹配 p.E451K 等蛋白变化）
        exon_m = re.search(r'\bE(\d{1,2})\b', combined)
    if exon_m:
        details["exon"] = exon_m.group(1)
    domain_m = re.search(
        r'(?:sterol.sensing|SSD|helical|cytoplasmic|lumenal|NTD|MLD|CTD|cysteine.?rich|luminal|transmembrane|LumenC?|TM\d*)\s*(domain|loop|lumen|region|helix)?',
        combined, re.IGNORECASE
    )
    if domain_m:
        domain_raw = domain_m.group(0).strip()
        # 过滤掉基因名误匹配
        domain_lower = domain_raw.lower()
        if domain_lower not in ('npc1', 'abca4', 'npc2'):
            domain_cn = (domain_raw
                .replace('domain', '结构域').replace('Domain', '结构域')
                .replace('loop', '环').replace('Loop', '环')
                .replace('lumen', '腔').replace('Lumen', '腔')
                .replace('helix', '螺旋').replace('region', '区域')
                .replace('LumenC', '腔C').replace('luminal', '腔内')
                .replace('transmembrane', '跨膜').replace('cytoplasmic', '胞质')
                .replace('cysteine-rich', '半胱氨酸富集').replace('cysteine rich', '半胱氨酸富集')
                .replace('sterol-sensing', '固醇感应'))
            details["domain"] = domain_cn

    return details


def extract_co_variants(sentences, tables=None, target_cdna="", target_protein="", variant_keywords=None):
    """从变异相关句子/表格中提取与目标变异配对的共存变异。v8: 增加 in trans with 句式 + 表格同行 ID 配对。

    支持格式：
    1. (;) 配对: p.[(G607R(;)R2040Q)] / c.5461-10T>C(;)p.(R2040Q)
    2. in trans with: "c.1761-2A>G was detected in trans with c.5882G>A"
    3. compound heterozygosity with: "found in compound heterozygosity with c.X"
    4. 表格同行配对: 同一患者 ID 的两行
    返回列表，每个元素为 dict: {"cdna": ..., "蛋白变异": ...}。"""
    # v8: 标准化 Unicode 空白字符（thin space 等）为普通空格，确保正则匹配不中断
    def _norm_ws(t):
        return re.sub(r'[ -‏]', ' ', t)

    combined = " ".join(sentences)
    combined = _norm_ws(combined)

    # 也对 sentences 做标准化（用于后续句子处理）
    sentences_norm = [_norm_ws(s) for s in sentences]

    # 1. 构建 protein → cDNA 映射表
    protein_to_cdna = {}
    mapping_re = re.compile(
        r'(?:exon/ivs\s*\d+\s+)?(c\.[\w\.\-\+>]+)\s+'
        r'((?:p\.)?\(?\s*[A-Za-z]{0,3}\d+\s*[A-Za-z*;]+\)?)',
        re.IGNORECASE
    )
    # 扩大映射来源：不仅是 variant_sentences，也包括所有 sentences
    all_text = " ".join(sentences_norm)  # v8: 使用标准化后的句子
    # 也从表格构建映射
    if tables:
        for t in tables:
            for row in t.get("rows", []):
                row_text = "\t".join(str(c or "") for c in row)
                all_text += "\t" + _norm_ws(row_text)

    for m in mapping_re.finditer(all_text):
        cdna = re.sub(r'\s+', '', m.group(1))
        protein = re.sub(r'\s+', '', m.group(2))
        if protein.startswith('p.'):
            protein_to_cdna[protein] = cdna
        else:
            protein_to_cdna[f"p.{protein}"] = cdna

    for p_full, c_val in list(protein_to_cdna.items()):
        short = re.sub(r'^p\.\(?|\)?$', '', p_full)
        protein_to_cdna[short] = c_val

    # 2. 提取 (;) 配对的共存变异
    pair_re = re.compile(
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)\s*\(;\)\s*'
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)',
        re.IGNORECASE
    )

    # v8: "in trans with" / "compound heterozygosity with" 语句模式
    trans_with_re = re.compile(
        r'(?:detected|found|identified|was|is|were|are)\s+in\s+trans\s+with\s+'
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)',
        re.IGNORECASE
    )
    # 反向: "c.X was detected in trans with [target]"
    trans_with_target_re = re.compile(
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)\s+'
        r'(?:was|is|were|are)\s+(?:detected|found|identified)\s+in\s+trans\s+with',
        re.IGNORECASE
    )
    compound_het_with_re = re.compile(
        r'(?:found|detected|identified|was|is|were|are)\s+in\s+compound\s+heterozygosity\s+with\s+'
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)',
        re.IGNORECASE
    )
    compound_het_state_re = re.compile(
        r'(?:in\s+a\s+)?compound\s+heterozygous\s+state\s+with\s+'
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)',
        re.IGNORECASE
    )

    # v9: "compound heterozygous for X and Y" / "compound heterozygous for X/Y" 句式
    compound_het_for_and_re = re.compile(
        r'compound\s+heterozygous\s+for\s+'
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)'
        r'\s+(?:and|/)\s+'
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)',
        re.IGNORECASE
    )
    # "X and Y were in compound heterozygous state" / "X and Y in compound heterozygosity"
    two_var_compound_re = re.compile(
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)'
        r'\s+(?:and|/)\s+'
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)'
        r'\s+(?:were|was|are|is)\s+(?:found\s+)?(?:in\s+)?compound\s+heterozyg',
        re.IGNORECASE
    )
    # HGVS allele format: c.[1761-2A>G];[5512C>T] or c.[1761-2A>G(+)5512C>T]
    hgvs_allele_re = re.compile(
        r'c\.\[([\w\.\-\+>]+)\]\s*[;,+]\s*\[([\w\.\-\+>]+)\]',
        re.IGNORECASE
    )
    # "with c.X on the other/second/opposite allele"
    other_allele_re = re.compile(
        r'(?:with|and)\s+'
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)'
        r'\s+on\s+the\s+(?:other|second|opposite)\s+allele',
        re.IGNORECASE
    )
    # "carried X and Y mutations" / "carrying X and Y"
    carried_and_re = re.compile(
        r'(?:carried|carrying|harbored|harboring|had)\s+'
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)'
        r'\s+(?:and|/)\s+'
        r'(c\.[\w\.\-\+>]+|\(?(?:p\.\s*)?[\[\(]?\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>+\-\.\w]*(?:\s*[\]\)])?)',
        re.IGNORECASE
    )

    # v9: in-cis 复杂等位基因 (complex allele) 配对模式
    # [p.D252N;p.S369X] — HGVS 方括号分号格式
    cis_bracket_re = re.compile(
        r'\[([^\];]*?p\.\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>\+\-\.\w]*)\s*;\s*'
        r'(p\.\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>\+\-\.\w]*[^\]]*?)\]',
        re.IGNORECASE
    )
    # in-cis mutations: 两个变异在同一句中，由 and 或 ; 连接
    # "two in-cis mutations p.D252N and p.S369X"
    in_cis_and_re = re.compile(
        r'(?:two|both|2)\s+(?:in[\s\-]cis\s+)?(?:missense|nonsense|frameshift|splice|point\s+)?mutations?\s+'
        r'(p\.\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>\+\-\.\w]*)'
        r'\s+(?:and|,)\s+'
        r'(p\.\s*[A-Za-z]{1,3}\d+\s*[A-Za-z*>\+\-\.\w]*)',
        re.IGNORECASE
    )
    # allele carrying two in-cis mutations [...]
    in_cis_allele_re = re.compile(
        r'(?:allele|one\s+allele)\s+(?:carrying|carried|with|had|containing)\s+'
        r'(?:two|both)\s+in[\s\-]cis\s+mutations?\s+\[([^\]]+)\]',
        re.IGNORECASE
    )

    target_protein_clean = (target_protein or "").replace("p.", "").replace("(", "").replace(")", "").strip()
    target_key_parts = []
    if target_protein_clean:
        target_key_parts.append(target_protein_clean)
    if variant_keywords:
        for kw in variant_keywords.get("protein", []):
            kw_clean = kw.replace("p.", "").replace("(", "").replace(")", "").replace(" ", "").strip()
            if kw_clean and kw_clean not in target_key_parts:
                target_key_parts.append(kw_clean)
        for kw in variant_keywords.get("fuzzy", []):
            kw_clean = kw.strip()
            if kw_clean and kw_clean.isdigit() and kw_clean not in target_key_parts:
                target_key_parts.append(kw_clean)

    def _is_target(v):
        """判断一个变异片段是否是目标变异。"""
        vc = re.sub(r'\s+', '', v.strip('()[]p. '))
        tc = (target_cdna or "").replace(" ", "")
        if vc.startswith('c.') and vc.replace("c.", "") == tc.replace("c.", ""):
            return True
        cp = vc.replace('p.', '').replace('(', '').replace(')', '').replace('[', '').replace(']', '')
        for kp in target_key_parts:
            if kp.replace(" ", "").lower() == cp.lower():
                return True
        if cp and any(kp.replace(" ", "") in cp for kp in target_key_parts if len(kp) >= 3):
            return True
        return False

    co_variants = []
    seen_keys = set()

    def _extract_from_text(text):
        """从文本中提取 (;) 配对的共存变异。"""
        results = []
        for m in pair_re.finditer(text):
            left = m.group(1)
            right = m.group(2)
            if _is_target(left) and not _is_target(right):
                co = right
            elif _is_target(right) and not _is_target(left):
                co = left
            else:
                continue
            results.append(_to_co_dict(co))
        return results

    def _extract_trans_with(text):
        """从文本提取 'in trans with' 语句中的共存变异。"""
        results = []
        for m in trans_with_re.finditer(text):
            co = m.group(1)
            if not _is_target(co):
                results.append(_to_co_dict(co))
        for m in compound_het_with_re.finditer(text):
            co = m.group(1)
            if not _is_target(co):
                results.append(_to_co_dict(co))
        for m in compound_het_state_re.finditer(text):
            co = m.group(1)
            if not _is_target(co):
                results.append(_to_co_dict(co))
        # 反向: "c.X was detected in trans with [target]" → c.X 是共存变异
        for m in trans_with_target_re.finditer(text):
            prefix_part = m.group(1)
            # 检查此句后面是否提到了目标变异
            rest_of_sentence = text[m.end():]
            has_target_nearby = False
            for kp in target_key_parts:
                if kp.lower() in rest_of_sentence.lower().replace(" ", ""):
                    has_target_nearby = True
                    break
            if target_cdna:
                cdna_short = target_cdna.replace("c.", "").replace(" ", "").lower()
                if cdna_short in rest_of_sentence.lower().replace(" ", ""):
                    has_target_nearby = True
            if has_target_nearby and not _is_target(prefix_part):
                results.append(_to_co_dict(prefix_part))
        return results

    def _extract_in_cis_with(text):
        """v9: 从文本提取 'in-cis' 复杂等位基因中的共存变异。
        支持格式:
        - [p.D252N;p.S369X] — HGVS 方括号分号 (allele notation)
        - two in-cis mutations p.D252N and p.S369X
        - allele carrying two in-cis mutations [p.D252N;p.S369X]
        """
        results = []
        # 1. 方括号 [p.X;p.Y] 格式
        for m in cis_bracket_re.finditer(text):
            left = m.group(1).strip()
            right = m.group(2).strip()
            if _is_target(left) and not _is_target(right):
                results.append(_to_co_dict(right))
            elif _is_target(right) and not _is_target(left):
                results.append(_to_co_dict(left))
        # 2. "two in-cis mutations X and Y" 格式
        for m in in_cis_and_re.finditer(text):
            left = m.group(1).strip()
            right = m.group(2).strip()
            if _is_target(left) and not _is_target(right):
                results.append(_to_co_dict(right))
            elif _is_target(right) and not _is_target(left):
                results.append(_to_co_dict(left))
        # 3. "allele carrying two in-cis mutations [...]" → 用分号拆分方括号内容
        for m in in_cis_allele_re.finditer(text):
            content = m.group(1)
            # 拆分分号分隔的变异
            parts = re.split(r'\s*;\s*', content)
            for part in parts:
                part = part.strip().rstrip(',')
                if part and not _is_target(part):
                    # 提取 p. 或 c. 变异
                    var_m = re.search(r'(c\.[\w\.\-\+>]+|p\.[A-Za-z]{1,3}\d+[A-Za-z*>\+\-\.\w]*)', part)
                    if var_m:
                        results.append(_to_co_dict(var_m.group(1)))
        return results

    def _to_co_dict(co_str):
        co_clean = re.sub(r'\s+', '', co_str.strip('()[] '))
        if co_clean.startswith('c.'):
            return {"cdna": co_clean, "蛋白变异": None}
        protein_norm = co_clean
        if not protein_norm.startswith('p.'):
            protein_norm = f"p.{protein_norm}"
        protein_norm = protein_norm.replace('(', '').replace(')', '').replace('[', '').replace(']', '')
        cdna_mapped = protein_to_cdna.get(protein_norm)
        if not cdna_mapped:
            short = protein_norm.replace('p.', '')
            cdna_mapped = protein_to_cdna.get(short)
        return {"cdna": cdna_mapped, "蛋白变异": protein_norm if not cdna_mapped else None}

    # 3. 从句子中提取
    for s in sentences_norm:
        s_lower = s.lower()
        has_target = False
        if target_cdna and target_cdna.replace(" ", "").lower() in s_lower.replace(" ", ""):
            has_target = True
        elif target_key_parts:
            for kp in target_key_parts:
                if kp.lower() in s_lower:
                    has_target = True
                    break
        if not has_target:
            continue

        # (;) 格式
        for co_d in _extract_from_text(s):
            key = co_d.get("cdna") or co_d.get("蛋白变异")
            if key and key not in seen_keys:
                seen_keys.add(key)
                co_variants.append(co_d)

        # v8: "in trans with" 格式
        for co_d in _extract_trans_with(s):
            key = co_d.get("cdna") or co_d.get("蛋白变异")
            if key and key not in seen_keys:
                seen_keys.add(key)
                co_variants.append(co_d)

        # v9: "in-cis" 复杂等位基因格式
        for co_d in _extract_in_cis_with(s):
            key = co_d.get("cdna") or co_d.get("蛋白变异")
            if key and key not in seen_keys:
                seen_keys.add(key)
                co_variants.append(co_d)

    # 4. 从表格中提取（包括同患者 ID 配对）
    if tables:
        # 4a. 先在每一行中查找 (;) 或 in trans with
        for table_info in tables:
            for row in table_info["rows"]:
                row_text = "\t".join(str(cell or "") for cell in row)
                row_lower = row_text.lower()
                has_target = False
                for kp in target_key_parts:
                    if kp.lower() in row_lower:
                        has_target = True
                        break
                if target_cdna:
                    cdna_short = target_cdna.replace("c.", "").replace(" ", "")
                    if cdna_short.lower() in row_lower:
                        has_target = True
                if not has_target:
                    continue

                for co_d in _extract_from_text(row_text):
                    key = co_d.get("cdna") or co_d.get("蛋白变异")
                    if key and key not in seen_keys:
                        seen_keys.add(key)
                        co_variants.append(co_d)

                # v9: also check table rows for in-trans / in-cis patterns
                for co_d in _extract_trans_with(row_text):
                    key = co_d.get("cdna") or co_d.get("蛋白变异")
                    if key and key not in seen_keys:
                        seen_keys.add(key)
                        co_variants.append(co_d)
                for co_d in _extract_in_cis_with(row_text):
                    key = co_d.get("cdna") or co_d.get("蛋白变异")
                    if key and key not in seen_keys:
                        seen_keys.add(key)
                        co_variants.append(co_d)

        # 4b. v8: 表格同患者 ID 配对 —— 相邻行共享相同患者/家系 ID，一行有目标变异，另一行有另一个变异
        for table_info in tables:
            rows = table_info.get("rows", [])
            if len(rows) < 2:
                continue
            # 找到表头，识别患者 ID 列和变异列
            header = rows[0] if rows else []
            patient_id_cols = []
            cdna_cols = []
            protein_cols = []
            for ci, cell in enumerate(header):
                cell_text = str(cell or "").strip().lower()
                # 清理 [H] 前缀
                cell_text = cell_text.replace("[h]", "")
                if any(kw in cell_text for kw in ["proband", "patient", "family", "pedigree", "id", "subject"]):
                    patient_id_cols.append(ci)
                if any(kw in cell_text for kw in ["nucleotide", "cdna", "dna change", "cdna change", "mutation"]):
                    cdna_cols.append(ci)
                if any(kw in cell_text for kw in ["amino acid", "protein", "aa change", "protein change"]):
                    protein_cols.append(ci)

            if not patient_id_cols:
                # 无正式 ID 列时，使用每行第一个非空单元格作为患者 ID
                patient_id_cols = [0]

            # 按患者 ID 分组行
            patient_groups = {}
            last_pid = None  # v11: 跟踪上一个有效患者 ID，用于处理跨行变异
            for ri, row in enumerate(rows[1:], 1):  # 跳过表头
                row_text = "\t".join(str(cell or "") for cell in row)
                # 获取患者 ID
                pid = None
                for pci in patient_id_cols:
                    if pci < len(row):
                        cell_val = str(row[pci] or "").strip()
                        if cell_val:
                            # v11: 如果"患者ID"列的值是变异格式（c./p.开头），
                            # 说明这是跨行变异续行（如 c.858+2T>A），不应作为患者ID
                            if not re.search(r'^(c\.|p\.|g\.|n\.)', cell_val):
                                pid = cell_val
                                break
                if not pid:
                    # 尝试从行中提取 ID 模式
                    id_match = re.search(r'(?:F\d+|P\d+|OX\d+|R\d+|[A-Z]+\d+)', row_text)
                    if id_match:
                        pid = id_match.group(0)
                if not pid:
                    # v11: 检查是否为变异续行（无患者ID，仅有变异如 c.858+2T>A）
                    # 合并到上一个患者组，确保多行变异（每行一个等位基因）被正确配对
                    first_cell = str(row[0] or "").strip() if row else ""
                    is_variant_continuation = bool(re.search(r'^(c\.|p\.|g\.|n\.)', first_cell))
                    if is_variant_continuation and last_pid and last_pid in patient_groups:
                        pid = last_pid
                    else:
                        pid = f"_row_{ri}"

                if pid not in patient_groups:
                    patient_groups[pid] = []
                patient_groups[pid].append((ri, row, row_text))
                last_pid = pid

            # 对于每个患者，如果一组行包含目标变异，提取该组中的其他变异
            for pid, group_rows in patient_groups.items():
                if len(group_rows) < 2:
                    continue
                has_target = False
                target_row_idx = -1
                for gi, (ri, row, row_text) in enumerate(group_rows):
                    row_lower = row_text.lower()
                    for kp in target_key_parts:
                        if kp.lower() in row_lower:
                            has_target = True
                            target_row_idx = gi
                            break
                    if target_cdna:
                        cdna_short = target_cdna.replace("c.", "").replace(" ", "")
                        if cdna_short.lower() in row_lower:
                            has_target = True
                            target_row_idx = gi
                            break
                if not has_target:
                    continue

                # 提取同组其他行中的变异
                for gi, (ri, row, row_text) in enumerate(group_rows):
                    if gi == target_row_idx:
                        continue
                    # 从该行提取变异
                    for cci in cdna_cols:
                        if cci < len(row):
                            cell_val = str(row[cci] or "").strip()
                            cdna_m = re.search(r'(c\.[\w\.\-\+>]+)', cell_val.replace(" ", ""))
                            if cdna_m and not _is_target(cdna_m.group(1)):
                                cdna_val = cdna_m.group(1)
                                if cdna_val not in seen_keys:
                                    seen_keys.add(cdna_val)
                                    co_variants.append({"cdna": cdna_val, "蛋白变异": None})
                                    break
                    for pci in protein_cols:
                        if pci < len(row):
                            cell_val = str(row[pci] or "").strip()
                            prot_m = re.search(r'(p\.(?:[A-Z][a-z]{2}|[A-Z])\d+(?:[A-Z][a-z]{2}|[A-Z\*])(?:Ter|fs|\*)?)', cell_val.replace(" ", ""))
                            if prot_m:
                                prot_val = prot_m.group(1)
                                key = prot_val
                                if key not in seen_keys:
                                    # 查找映射
                                    cdna_map = protein_to_cdna.get(prot_val)
                                    seen_keys.add(key)
                                    co_variants.append({"cdna": cdna_map, "蛋白变异": None if cdna_map else prot_val})
                                    break

                    # fallback: 从行文本提取 c./p. (仅在列提取未找到时)
                    col_extracted = False
                    if cdna_cols:
                        max_col = max(cdna_cols)
                        for cci in cdna_cols:
                            if cci < len(row) and row[cci] and not _is_target(str(row[cci]).strip()):
                                col_extracted = True
                                break
                    if not col_extracted:
                        cdna_fallback = re.search(r'(c\.[\w\.\-\+>]+)', row_text.replace(" ", ""))
                        if cdna_fallback and not _is_target(cdna_fallback.group(1)):
                            cdna_val = cdna_fallback.group(1)
                            if cdna_val not in seen_keys:
                                seen_keys.add(cdna_val)
                                co_variants.append({"cdna": cdna_val, "蛋白变异": None})

    # 4c. v10: 同表行多变异提取 —— 同一行中存在两个不同变异（不同列）
    # 处理如 p.G1961E  p.R1129C 出现在同一患者行的相邻列
    for table_info in tables:
        for row in table_info.get("rows", []):
            row_text = "	".join(str(cell or "") for cell in row)
            row_lower = row_text.lower()
            # 检查该行是否包含目标变异
            has_target = False
            for kp in target_key_parts:
                if kp.lower() in row_lower:
                    has_target = True
                    break
            if target_cdna:
                cdna_short = target_cdna.replace("c.", "").replace(" ", "")
                if cdna_short.lower() in row_lower.replace(" ", ""):
                    has_target = True
            if not has_target:
                continue

            # 扫描该行所有单元格，找其他变异
            # 蛋白变异正则（同时支持三字母和单字母格式）
            prot_re_3letter = re.compile(r'(p\.\s*[A-Z][a-z]{2}\d+[A-Z][a-z]{2}(?:\*|Ter|fs)?)')
            prot_re_1letter = re.compile(r'(p\.\s*[A-Z]\d+[A-Z\*])')
            cdna_re_cell = re.compile(r'(c\.[\w\.\-\+>]+)')

            for cell in row:
                cell_str = str(cell or "").strip()
                if not cell_str:
                    continue
                cell_clean = cell_str.replace(" ", "")

                # 尝试三字母格式
                for prot_m in prot_re_3letter.finditer(cell_clean):
                    prot_val = prot_m.group(1)
                    if not _is_target(prot_val):
                        key = prot_val
                        if key not in seen_keys:
                            seen_keys.add(key)
                            cdna_map = protein_to_cdna.get(prot_val)
                            if not cdna_map:
                                short = prot_val.replace("p.", "")
                                cdna_map = protein_to_cdna.get(short)
                            co_variants.append({"cdna": cdna_map, "蛋白变异": None if cdna_map else prot_val})

                # 尝试单字母格式
                for prot_m in prot_re_1letter.finditer(cell_clean):
                    prot_val = prot_m.group(1)
                    if not _is_target(prot_val):
                        key = prot_val
                        if key not in seen_keys:
                            seen_keys.add(key)
                            cdna_map = protein_to_cdna.get(prot_val)
                            if not cdna_map:
                                short = prot_val.replace("p.", "")
                                cdna_map = protein_to_cdna.get(short)
                            co_variants.append({"cdna": cdna_map, "蛋白变异": None if cdna_map else prot_val})

                # 尝试 cDNA 格式
                for cdna_m in cdna_re_cell.finditer(cell_clean):
                    cdna_val = cdna_m.group(1)
                    if not _is_target(cdna_val):
                        key = cdna_val
                        if key not in seen_keys:
                            seen_keys.add(key)
                            co_variants.append({"cdna": cdna_val, "蛋白变异": None})

    return co_variants


def extract_phase_evidence(sentences, co_variants=None, zygosity=""):
    """
    全面提取正反式（cis/trans）配置和等位基因相位（phase）信息。

    从正文和表格内容中搜索以下证据类型：
    1. 亲本检测证据（最高置信度）
    2. 反式 (trans) 配置证据
    3. 顺式 (cis) 配置证据
    4. 相位未知/不确定的证据
    5. 亲本源具体信息

    Returns dict:
        - trans_evidence: bool
        - cis_evidence: bool
        - phase_status: confirmed_in_trans / confirmed_in_cis / presumed_in_trans / presumed_in_cis / phase_not_determined / not_assessed / not_applicable
        - phase_confidence: confirmed / presumed / unknown / not_applicable
        - parental_testing: bool
        - maternal_variant: str | None
        - paternal_variant: str | None
        - phase_evidence_sentences: list of {type, label, sentence}
        - phase_detail: 中文详细描述
    """
    combined = " ".join(sentences)
    combined_lower = combined.lower()

    result = {
        "trans_evidence": False,
        "cis_evidence": False,
        "phase_status": "not_assessed",
        "phase_confidence": "not_applicable",
        "parental_testing": False,
        "maternal_variant": None,
        "paternal_variant": None,
        "phase_evidence_sentences": [],
        "phase_detail": "",
    }

    evidence_sentences = []

    # ── 1. 反式 (trans) 配置证据 ──
    trans_patterns = [
        (r'in\s+trans\b', "in trans"),
        (r'on\s+opposite\s+allele', "on opposite alleles"),
        (r'trans\s+configuration', "trans configuration"),
        (r'on\s+different\s+allele', "on different alleles"),
        (r'biallelic\s+in\s+trans', "biallelic in trans"),
        (r'compound\s+heterozygous.*confirmed', "compound heterozygous confirmed"),
        (r'confirmed.*compound\s+heterozyg', "confirmed compound heterozygous"),
    ]
    for pattern, label in trans_patterns:
        for s in sentences:
            if re.search(pattern, s, re.IGNORECASE):
                result["trans_evidence"] = True
                evidence_sentences.append({"type": "trans", "label": label, "sentence": s.strip()[:400]})
                break
        if result["trans_evidence"]:
            break

    # ── 2. 顺式 (cis) 配置证据 ──
    cis_patterns = [
        (r'in[\s\-]cis\b', "in cis"),
        (r'on\s+(?:the\s+)?same\s+allele', "on the same allele"),
        (r'cis\s+configuration', "cis configuration"),
        (r'complex\s+allele', "complex allele"),
        (r'double\s+mutant\s+allele', "double mutant allele"),
        (r'both\s+variants\s+on\s+(?:the\s+)?same', "both variants on same"),
    ]
    for pattern, label in cis_patterns:
        for s in sentences:
            if re.search(pattern, s, re.IGNORECASE):
                result["cis_evidence"] = True
                evidence_sentences.append({"type": "cis", "label": label, "sentence": s.strip()[:400]})
                break
        if result["cis_evidence"]:
            break

    # ── 3. 亲本检测证据 ──
    parental_patterns = [
        (r'inherited\s+from\s+(?:the\s+)?mother', "inherited from mother"),
        (r'inherited\s+from\s+(?:the\s+)?father', "inherited from father"),
        (r'maternally\s+inherited', "maternally inherited"),
        (r'paternally\s+inherited', "paternally inherited"),
        (r'mother\s+(?:was|is)\s+(?:a\s+)?(?:carrier|heterozygous|homozygous)', "mother carrier"),
        (r'father\s+(?:was|is)\s+(?:a\s+)?(?:carrier|heterozygous|homozygous)', "father carrier"),
        (r'parental\s+(?:testing|analysis|study|segregation|origin)', "parental testing"),
        (r'segregation\s+analysis', "segregation analysis"),
        (r'co[\-\s]?segregat', "co-segregation"),
        (r'trio\s+(?:analysis|sequencing|WES|WGS)', "trio analysis"),
        (r'parents?\s+(?:were|are)\s+(?:tested|analyzed|genotyped|sequenced)', "parents tested"),
    ]
    for pattern, label in parental_patterns:
        for s in sentences:
            m = re.search(pattern, s, re.IGNORECASE)
            if m:
                result["parental_testing"] = True
                evidence_sentences.append({"type": "parental", "label": label, "sentence": s.strip()[:400]})
                # 提取亲本源
                s_lower = s.lower()
                if 'mother' in s_lower or 'maternal' in s_lower:
                    # 尝试提取母源变异
                    if result["maternal_variant"] is None:
                        result["maternal_variant"] = _extract_parental_variant(s, "maternal")
                if 'father' in s_lower or 'paternal' in s_lower:
                    if result["paternal_variant"] is None:
                        result["paternal_variant"] = _extract_parental_variant(s, "paternal")
                break
        if result["parental_testing"]:
            # 继续搜索更多亲本信息
            for s in sentences:
                _extract_parental_details(s, result)

    # ── 4. 相位未知/不确定证据 ──
    phase_unknown_patterns = [
        (r'phase\s+(?:not\s+determined|unknown|uncertain|could\s+not\s+be\s+determined)', "phase not determined"),
        (r'(?:cis|trans)\s+not\s+determined', "cis/trans not determined"),
        (r'allelic\s+phase\s+(?:unknown|not\s+determined)', "allelic phase unknown"),
        (r'phase\s+(?:was|is|remains)\s+(?:unable\s+to\s+be\s+)?(?:determined|unknown)', "phase unknown"),
    ]
    phase_unknown_found = False
    for pattern, label in phase_unknown_patterns:
        for s in sentences:
            if re.search(pattern, s, re.IGNORECASE):
                phase_unknown_found = True
                evidence_sentences.append({"type": "phase_unknown", "label": label, "sentence": s.strip()[:400]})
                break
        if phase_unknown_found:
            break

    # ── 5. 新发突变 (de novo) 证据 ──
    de_novo_found = False
    for s in sentences:
        if re.search(r'de\s+novo\b', s, re.IGNORECASE):
            de_novo_found = True
            evidence_sentences.append({"type": "de_novo", "label": "de novo", "sentence": s.strip()[:400]})
            break

    # ── 6. 相位状态判定 ──
    has_co_variants = co_variants and len(co_variants) > 0
    is_homozygous = "纯合" in zygosity or "homozyg" in zygosity.lower()

    if result["trans_evidence"]:
        if result["parental_testing"]:
            result["phase_status"] = "confirmed_in_trans"
            result["phase_confidence"] = "confirmed"
            result["phase_detail"] = "经亲本检测确认位于反式位置（不同等位基因）"
        else:
            result["phase_status"] = "confirmed_in_trans"
            result["phase_confidence"] = "confirmed"
            result["phase_detail"] = "文献明确指出处于反式位置"
    elif result["cis_evidence"]:
        result["phase_status"] = "confirmed_in_cis"
        result["phase_confidence"] = "confirmed"
        result["phase_detail"] = "文献明确指出处于顺式位置（同一等位基因，复杂等位基因）"
    elif phase_unknown_found:
        result["phase_status"] = "phase_not_determined"
        result["phase_confidence"] = "unknown"
        result["phase_detail"] = "文献明确指出相位未确定"
    elif is_homozygous and not has_co_variants:
        result["phase_status"] = "not_applicable"
        result["phase_confidence"] = "not_applicable"
        result["phase_detail"] = "纯合变异，不存在相位问题"
    elif has_co_variants:
        # 复合杂合，隐性遗传默认推定反式
        # v11: 如果同时标记为纯合但有共存变异，以共存变异为准（纯合判定可能来自不完整行分析）
        result["phase_status"] = "presumed_in_trans"
        result["phase_confidence"] = "presumed"
        result["phase_detail"] = "复合杂合变异，推定为反式位置（未经亲本验证）"
    elif not has_co_variants:
        result["phase_status"] = "not_applicable"
        result["phase_confidence"] = "not_applicable"
        result["phase_detail"] = "单变异，不存在相位问题"
    else:
        result["phase_status"] = "not_assessed"
        result["phase_confidence"] = "unknown"
        result["phase_detail"] = "文献未涉及相位信息"

    # 如果 de novo，覆盖部分判定
    if de_novo_found:
        result["phase_detail"] += "；该变异为新发突变 (de novo)"

    result["phase_evidence_sentences"] = evidence_sentences

    # ── 兼容旧接口：反式确认 ──
    old_trans_confirmed = result["trans_evidence"]

    return result, old_trans_confirmed


def _extract_parental_variant(sentence, parent_type):
    """从句子中提取亲本携带的具体变异。"""
    # 匹配模式如: "mother carried c.123A>G" / "paternal allele: c.456C>T"
    patterns = [
        rf'{parent_type}.*?(c\.[\d\w_*><+=delinsup\-\+]+)',
        rf'{parent_type}.*?(p\.(?:[A-Z][a-z]{{2}}|[A-Z])\d+(?:[A-Z][a-z]{{2}}|[A-Z\*])(?:Ter|fs|\*)?)',
    ]
    for pat in patterns:
        m = re.search(pat, sentence, re.IGNORECASE)
        if m:
            return m.group(1)
    return None


def _extract_parental_details(sentence, result):
    """从句子中补充亲本详细信息。"""
    s_lower = sentence.lower()

    # 检测母源携带状态
    mother_patterns = [
        (r'mother.*?(?:c\.|p\.)', 'mother_carrier'),
        (r'maternal\s+allele.*?(?:c\.|p\.)', 'maternal_carrier'),
    ]
    for pat, label in mother_patterns:
        if re.search(pat, s_lower, re.IGNORECASE):
            if result["maternal_variant"] is None:
                m = re.search(r'(c\.[\d\w_*><+=delinsup\-\+]+)', sentence)
                if m:
                    result["maternal_variant"] = m.group(1)
            result["parental_testing"] = True
            break

    # 检测父源携带状态
    father_patterns = [
        (r'father.*?(?:c\.|p\.)', 'father_carrier'),
        (r'paternal\s+allele.*?(?:c\.|p\.)', 'paternal_carrier'),
    ]
    for pat, label in father_patterns:
        if re.search(pat, s_lower, re.IGNORECASE):
            if result["paternal_variant"] is None:
                m = re.search(r'(c\.[\d\w_*><+=delinsup\-\+]+)', sentence)
                if m:
                    result["paternal_variant"] = m.group(1)
            result["parental_testing"] = True
            break


def extract_patient_count(sentences, tables=None, keywords=None, target_cdna=""):
    """统计携带目标变异的患者数量。v8: 改进表格计数，按患者 ID 或家系 ID 去重。

    策略：
    1. 句子中提取数字（如 "5 patients"）
    2. 表格中按唯一患者 ID 计数（更准确）
    """
    combined = " ".join(sentences)
    combined_lower = combined.lower()
    count = 0

    # 从句子中提取数字
    count_patterns = [
        r'(\d+)\s*(patients?|subjects?|individuals?|cases?|families?|probands?)',
        r'(\d+)\s*heterozygous',
        r'(\d+)\s*homozygous',
        r'(\d+)\s*carriers?',
        r'(\d+)\s*affected',
        r'n\s*=\s*(\d+)',
    ]
    for pattern in count_patterns:
        m = re.search(pattern, combined_lower)
        if m:
            count = max(count, int(m.group(1)))

    # v8: 从表格中按唯一患者/家系 ID 计数
    table_count = 0
    if tables and keywords:
        target_cdna_short = (target_cdna or "").replace("c.", "").replace(" ", "").lower()
        for table_info in tables:
            rows = table_info.get("rows", [])
            if len(rows) < 2:
                continue

            # 找到患者 ID 列
            header = rows[0] if rows else []
            patient_id_col = 0  # 默认第一列
            for ci, cell in enumerate(header):
                cell_text = str(cell or "").strip().lower().replace("[h]", "")
                if any(kw in cell_text for kw in ["proband", "patient", "family", "pedigree", "id", "subject"]):
                    patient_id_col = ci
                    break

            # 收集含有目标变异的唯一患者/家系 ID
            variant_patients = set()
            for row in rows[1:]:  # 跳过表头
                row_text = "\t".join(str(cell or "") for cell in row)
                row_lower = row_text.lower()
                has_variant = False
                for kw in keywords.get("all", []):
                    if kw.lower() in row_lower:
                        has_variant = True
                        break
                if target_cdna_short and target_cdna_short in row_lower.replace(" ", ""):
                    has_variant = True

                if has_variant:
                    # 提取患者 ID
                    pid = None
                    if patient_id_col < len(row):
                        pid = str(row[patient_id_col] or "").strip()
                    if not pid:
                        # 尝试从行中提取 ID 模式
                        id_match = re.search(r'(?:F\d+|P\d+|OX\d+|R\d+|[A-Z]+[:\-]?\d+)', row_text)
                        if id_match:
                            pid = id_match.group(0)
                    if not pid:
                        # fallback: 使用行文本 hash
                        pid = f"_row_{hash(row_text) % 10000}"
                    variant_patients.add(pid)

            if variant_patients:
                # 每个唯一的家系/患者 ID 计为 1
                # 但如果同一患者有两行（复合杂合），也只计 1
                table_count = max(table_count, len(variant_patients))

    return max(count, table_count, 1) if (count > 0 or table_count > 0) else 0


def extract_variant_features(sentences, cdna, protein):
    """提取变异特征（CpG位点、NMD、人群频率、新型变异等）。"""
    combined = " ".join(sentences)
    combined_lower = combined.lower()
    features = {}

    # 新型/首次报道
    if re.search(r'(?:novel|first.{0,20}(?:mutation|report|describe|identif))', combined_lower):
        features["是否为新型变异"] = True

    # CpG 二核苷酸位点
    m = re.search(r'[^.\n]{0,150}CpG\s*(?:di)?nucleotide[^.\n]{0,150}', combined, re.IGNORECASE)
    if m:
        features["CpG位点"] = True
        features["CpG原文"] = m.group(0).strip()

    # 人群频率
    m = re.search(r'(?:\d+)\s+(?:control|allele|chromosome|individual).{0,50}(?:not\s+found|not\s+present|absent)', combined, re.IGNORECASE)
    if m:
        features["人群频率信息"] = True
        features["人群频率原文"] = m.group(0).strip()

    m2 = re.search(r'(?:gnomAD|ExAC|1000\s+Genomes|dbSNP).{0,100}(?:absent|not\s+found|frequency)', combined, re.IGNORECASE)
    if m2:
        features["人群频率信息"] = True
        features["人群频率原文"] = m2.group(0).strip()

    # NMD
    if re.search(r'nonsense.{0,10}mediated\s+(?:mRNA\s+)?decay|NMD|premature\s+stop\s+codon.{0,100}(?:decay|degrad|NMD)', combined, re.IGNORECASE):
        features["NMD相关信息"] = True
        m = re.search(r'[^.\n]{0,300}(?:NMD|nonsense.{0,10}mediated|premature.{0,50}decay)[^.\n]{0,300}', combined, re.IGNORECASE | re.DOTALL)
        if m:
            features["NMD原文"] = m.group(0).strip()

    # 转染/功能实验
    if re.search(r'(?:transfect|expression|functional\s+assay|enzyme\s+activit|complementation).{0,80}(?:wild.?type|WT|mutant|plasmid|vector)', combined, re.IGNORECASE):
        features["转染/功能实验"] = True

    # 生物信息学预测
    pred_tools = {}
    for tool, pat in [("SIFT", r'SIFT[^.]{0,50}(?:deleterious|damaging|tolerated)'),
                       ("PolyPhen", r'PolyPhen[^.]{0,50}(?:damaging|benign)'),
                       ("CADD", r'CADD[^.]{0,30}score[^.]{0,20}\d'),
                       ("MutationTaster", r'MutationTaster[^.]{0,50}(?:disease.?causing|benign)')]:
        m = re.search(pat, combined, re.IGNORECASE)
        if m:
            pred_tools[tool] = m.group(0).strip()
    if pred_tools:
        features["生物信息学预测"] = pred_tools

    return features


def _extract_pathogenicity_from_table(table_variant_rows):
    """从表格行中推断致病性分类（识别如 DP, HLP, LP, VUS 等分类标记）。"""
    classification_map = {
        "dp": "致病 (pathogenic)",
        "hlp": "可能致病 (likely pathogenic)",
        "lp": "可能致病 (likely pathogenic)",
        "hyp": "可能致病 (likely pathogenic)",
        "i": "意义不明 (VUS)",
        "ln": "可能良性 (likely benign)",
        "vus": "意义不明 (VUS)",
    }
    for row in table_variant_rows:
        for cell in row:
            cell_lower = str(cell or "").strip().lower()
            if cell_lower in classification_map:
                return classification_map[cell_lower]
    return "未指明"


def _parse_table_row_for_patient(row, target_cdna=None, target_protein=None):
    """从表格行中解析患者信息。根据常见表格列结构推断字段。"""
    detail = {}
    # 解析常见表格列：Ped ID, Pt ID, Gene, Exon, Codon, cDNA, Protein, Classification, ...
    clean_cells = [str(c or "").strip() for c in row]
    clean_cells = [c for c in clean_cells if c]  # 移除空值

    if len(clean_cells) >= 2:
        # 第一个非空的通常是 pedigree/family ID
        if re.match(r'^[MP]\d+', clean_cells[0]):
            detail["pedigree_id"] = clean_cells[0]
            # 第二个可能是 patient ID
            if len(clean_cells) >= 2 and re.match(r'^(OX|R)\d+', clean_cells[1]):
                detail["patient_id"] = clean_cells[1]

    # 找到包含 cDNA 和 protein 的列
    for cell in clean_cells:
        cell_clean = cell.replace(" ", "")
        if re.match(r'c\.\d+[A-ZTCG*><delinsdup]+', cell_clean):
            detail["variant"] = cell
        if re.match(r'p\.(?:[A-Z][a-z]{2}|[A-Z])\d+(?:[A-Z][a-z]{2}|[A-Z\*])(?:Ter|fs|\*)?', cell):
            detail["protein_change"] = cell

    if target_cdna:
        detail.setdefault("variant", target_cdna)
    if target_protein:
        detail.setdefault("protein_change", target_protein)

    if len(detail) > 1:
        return detail
    return None


def _extract_table_features(features, table_variant_rows):
    """从表格行中提取额外特征（交叉引用次数等）。"""
    for row in table_variant_rows:
        for cell in row:
            cell_str = str(cell or "").strip()
            # 检测交叉引用次数标记，如 "4×" 或 "23×" 或 "1×"
            m = re.match(r'^(\d+)[×xX]', cell_str)
            if m:
                features["交叉引用次数"] = m.group(1)
                return


def extract_info_for_variant(result, target_gene, target_cdna, target_protein, keywords, pdf_result=None):
    """仅针对目标变异提取所有信息。同时搜索正文和表格。v9: 支持PDF搜索结果交叉引用。"""
    full_text = result.get("全文", "")
    # v8: Unicode 空白字符标准化 + 变异符号间距修复
    full_text = re.sub(r'[-‏]', ' ', full_text)
    full_text = re.sub(r'\bc\.\s*([\d\w\.\-\+]+)\s*>\s*([\d\w\.\-\+]+)', r'c.\1>\2', full_text)
    full_text = re.sub(r'\bp\.\s*\(?([A-Za-z]{1,3})\s*(\d+)\s*([A-Za-z\*]+)\)?',
               lambda m: f"p.{m.group(1)}{m.group(2)}{m.group(3)}", full_text)
    # 同时标准化表格数据
    tables = result.get("tables", [])
    if tables:
        for t in tables:
            for row in t.get("rows", []):
                for i in range(len(row)):
                    if row[i]:
                        row[i] = re.sub(r'[-‏　]', ' ', str(row[i]))
    full_text_lower = full_text.lower()

    result["基因"] = target_gene.upper() if target_gene else ""
    result["cDNA变异"] = target_cdna
    result["蛋白变异"] = target_protein

    # 搜索变异相关句子（正文）
    variant_sentences, matched_kws = find_variant_sentences(full_text, keywords)

    # 同时搜索表格中的变异
    table_sentences = []
    table_matched_kws = []
    table_variant_rows = []
    if tables:
        for table_info in tables:
            for row in table_info["rows"]:
                row_text = "\t".join(str(cell or "") for cell in row)
                row_text = re.sub(r"[-‏ 　]", " ", row_text)  # v8: 去除 Unicode 空白
                row_lower = row_text.lower()
                for kw in keywords["all"]:
                    if kw.lower() in row_lower:
                        table_sentences.append(row_text)
                        table_matched_kws.append(kw)
                        table_variant_rows.append(row)
                        break

    # 合并正文和表格的匹配结果
    all_sentences = variant_sentences + table_sentences
    all_matched_kws = list(set(matched_kws + table_matched_kws))

    result["相关句子"] = all_sentences[:10]
    result["匹配关键词"] = all_matched_kws

    if not all_sentences:
        # v9: 交叉引用PDF搜索结果
        if pdf_result and pdf_result.get("变异提及") and (pdf_result.get("正文匹配句") or pdf_result.get("表格匹配行")):
            print(f"  [交叉引用] 在线API未检测到变异，但PDF搜索已确认提及")
            pdf_sentences = list(pdf_result.get("正文匹配句", []))
            for row in pdf_result.get("表格匹配行", []):
                if isinstance(row, list):
                    pdf_sentences.append("\t".join(str(c or "") for c in row))
                elif isinstance(row, str):
                    pdf_sentences.append(row)
            all_sentences = [s for s in pdf_sentences if isinstance(s, str) and s.strip()]
            variant_sentences = all_sentences
            table_sentences = []  # PDF表格行已合并到all_sentences中
            result["全文来源"] = result.get("全文来源", "") + "+pdf_crossref"
            if pdf_result.get("匹配关键词"):
                all_matched_kws = list(set(pdf_result["匹配关键词"]))
                result["匹配关键词"] = all_matched_kws
            # 将PDF文本也合并到全文以支持后续提取
            pdf_text = " ".join(all_sentences)
            if pdf_text:
                full_text = full_text + "\n[PDF全文]\n" + pdf_text
                full_text_lower = full_text.lower()
            # 不在此处return，继续执行下方提取逻辑
        else:
            result["变异提及"] = False
            result["变异类型"] = "不适用"
            result["致病性"] = "不适用"
            result["合子状态"] = "不适用"
            result["临床表型"] = "不适用"
            result["功能验证"] = "不适用"
            return result

    result["变异提及"] = True

    # 变异类型
    result["变异类型"] = infer_variant_type(target_cdna, target_protein, variant_sentences)

    # 致病性（同时从表格和文字提取）
    result["致病性"] = extract_pathogenicity(variant_sentences, full_text_lower)
    if result["致病性"] == "未指明" and table_variant_rows:
        result["致病性"] = _extract_pathogenicity_from_table(table_variant_rows)

    # 遗传方式（v8: 添加邻近度参数）
    result["合子状态"] = extract_zygosity(variant_sentences, full_text_lower, target_cdna, target_protein, tables, keywords)

    # 遗传模式
    result["遗传模式"] = extract_inheritance(variant_sentences, full_text_lower)

    # 临床表型
    phenotypes = extract_patient_phenotypes(variant_sentences, tables, keywords)
    result["临床表型"] = "、".join(phenotypes) if phenotypes else "未指明"

    # v9: 深层临床细节（疾病亚型、发病年龄、实验室发现、变异频率等）
    clinical_details = extract_clinical_details(variant_sentences, tables, keywords)
    # 同时从表格行中提取临床细节
    if tables:
        for t in tables:
            for row in t.get("rows", []):
                row_text = "\t".join(str(cell or "") for cell in row)
                has_variant = any(kw.lower() in row_text.lower() for kw in keywords["all"])
                if has_variant:
                    row_details = extract_clinical_details([row_text])
                    for key, val in row_details.items():
                        if isinstance(val, list) and key in clinical_details and isinstance(clinical_details.get(key), list):
                            for item in val:
                                if item not in clinical_details[key]:
                                    clinical_details[key].append(item)
                        elif not clinical_details.get(key):
                            clinical_details[key] = val
    result["临床详情"] = clinical_details

    # 患者详情
    patient_details = []
    for s in variant_sentences:
        detail = _parse_patient_sentence(s, target_cdna, target_protein)
        if detail:
            patient_details.append(detail)
    # 也从表格行提取患者信息
    for row in table_variant_rows:
        detail = _parse_table_row_for_patient(row, target_cdna, target_protein)
        if detail:
            patient_details.append(detail)
    result["患者详情"] = patient_details

    # 共存变异
    result["共存变异"] = extract_co_variants(variant_sentences, tables, target_cdna, target_protein, keywords)

    # ── v7: 正反式/相位提取（替代旧 extract_trans_evidence）──
    phase_result, trans_confirmed = extract_phase_evidence(
        all_sentences, result.get("共存变异", []), result.get("合子状态", "")
    )
    result["反式确认"] = trans_confirmed
    result["顺式确认"] = phase_result["cis_evidence"]
    result["相位状态"] = phase_result["phase_status"]
    result["相位置信度"] = phase_result["phase_confidence"]
    result["亲本检测"] = phase_result["parental_testing"]
    result["母源变异"] = phase_result["maternal_variant"]
    result["父源变异"] = phase_result["paternal_variant"]
    result["相位证据"] = phase_result["phase_evidence_sentences"]
    result["相位详情"] = phase_result["phase_detail"]

    # 从相位证据中追加相关句子（避免遗漏关键亲本/相位描述）
    for ev in phase_result.get("phase_evidence_sentences", []):
        ev_sentence = ev.get("sentence", "")
        if ev_sentence and ev_sentence not in all_sentences:
            all_sentences.append(ev_sentence)

    # 患者数量（优先从表格计数）
    result["患者数量"] = extract_patient_count(variant_sentences, tables, keywords, target_cdna)

    # 变异特征
    result["变异特征"] = extract_variant_features(variant_sentences, target_cdna, target_protein)
    # 表格中的交叉引用次数等信息
    if table_variant_rows:
        _extract_table_features(result["变异特征"], table_variant_rows)

    # 功能验证
    func_found = []
    func_details = []
    func_keywords_map = {
        "体外表达实验": r'in\s+vitro\s+(?:experiment|study|assay|expression)',
        "酶活性检测": r'enzyme\s+activit|residual\s+activit|catalytic\s+activit',
        "蛋白表达分析": r'protein\s+expres|western\s+blot|immunoblot',
        "mRNA分析": r'mrna\s+(?:level|expression|analysis)|transcript.*level',
        "剪接分析": r'splic(e|ing)\s+(?:assay|analysis|minigene)',
        "结构建模": r'structural\s+model|homology\s+model|3d\s+model|pymol|swiss.?model',
        "功能实验": r'functional\s+(?:assay|study|test|characterization)',
        "共分离分析": r'segregation\s+(?:analysis|study)|co.?segregat',
    }
    for func_name, func_pattern in func_keywords_map.items():
        for s in variant_sentences:
            if re.search(func_pattern, s, re.IGNORECASE):
                if func_name not in func_found:
                    func_found.append(func_name)
                    func_details.append(s.strip())
                break
    result["功能验证"] = "、".join(func_found) if func_found else "未进行"
    result["功能验证详情"] = func_details[:5]

    # 表格信息摘要
    if tables:
        result["表格数量"] = len(tables)
        result["表格摘要"] = [
            {"id": t["id"], "caption": t["caption"], "rows": len(t["rows"])}
            for t in tables
        ]

    return result


def _parse_patient_sentence(sentence, target_cdna=None, target_protein=None):
    """从单个句子中尝试解析患者个体信息。"""
    detail = {}

    patient_id_match = re.search(
        r'(?:patient|case|proband|family|病例|患者)\s*(\d+|[IVX]+)',
        sentence, re.IGNORECASE
    )
    if patient_id_match:
        detail["patient_id"] = patient_id_match.group(1)

    if re.search(r'homozyg', sentence, re.IGNORECASE):
        detail["zygosity"] = "homozygous"
    elif re.search(r'compound\s*heterozyg', sentence, re.IGNORECASE):
        detail["zygosity"] = "compound heterozygous"
    elif re.search(r'heterozyg', sentence, re.IGNORECASE):
        detail["zygosity"] = "heterozygous"

    phenotype_words = [
        "developmental delay", "intellectual disability", "growth retardation",
        "hypotonia", "seizure", "encephalopathy", "feeding difficulty",
        "speech delay", "motor delay", "dysmorphic", "microcephaly",
        "hepatomegaly", "jaundice", "lethargy", "vomiting",
        "poor feeding", "failure to thrive", "metabolic crisis",
        "xanthoma", "hypercholesterolemia", "sitosterolemia",
        "atherosclerosis", "thrombocytopenia", "anemia",
    ]
    phenos = [p for p in phenotype_words if p in sentence.lower()]
    if phenos:
        detail["phenotype"] = ", ".join(phenos)

    age_match = re.search(r'(\d+\s*(?:months?|years?|days?|weeks?|old|月|岁))', sentence, re.IGNORECASE)
    if age_match:
        detail["age"] = age_match.group(1)

    sex_match = re.search(r'\b(male|female|boy|girl|son|daughter|男|女)\b', sentence, re.IGNORECASE)
    if sex_match:
        detail["sex"] = sex_match.group(1).capitalize()

    if target_cdna:
        detail["variant"] = target_cdna
    if target_protein:
        detail["protein_change"] = target_protein

    if len(detail) > 1:
        return detail
    return None