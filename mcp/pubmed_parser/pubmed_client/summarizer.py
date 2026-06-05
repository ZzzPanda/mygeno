"""
Summary generation functions for PubMed variant extraction results.
"""

import re
import os
import csv

from .constants import TOPIC_KEYWORD_MAP


def generate_summary_paragraph(result):
    """根据遗传模式和文献实际数据生成详细中文总结段落（v9增强版）。
    整合变异基本信息+患者详情+临床细节+文献来源。
    """
    zygosity = result.get("合子状态", "")
    cdna = result.get("cDNA变异", "") or ""
    protein = result.get("蛋白变异", "") or ""
    gene = result.get("基因", "")
    var_char = result.get("变异特征", {})
    patient_count = result.get("患者数量", 0)
    phenotype = result.get("临床表型", "")
    func_valid = result.get("功能验证", "")
    co_variants = result.get("共存变异", [])
    patient_details = result.get("患者详情", [])
    clinical_details = result.get("临床详情", {})
    variant_type = result.get("变异类型", "")

    # 相位信息
    phase_status = result.get("相位状态", "")
    parental_testing = result.get("亲本检测", False)
    maternal_var = result.get("母源变异")
    paternal_var = result.get("父源变异")

    variant_display = f"{cdna}（{protein}）" if cdna and protein else (cdna or protein or "该变异")

    # 蛋白缩写形式 p.D501Y
    protein_short = protein.replace("p.", "") if protein else ""
    if protein_short:
        variant_display_short = f"{cdna}（p.{protein_short}）"
    else:
        variant_display_short = variant_display

    # 构建参考文献字符串
    ref_parts = []
    if result.get("作者"):
        authors = result["作者"]
        if len(authors) > 3:
            ref_parts.append(f"{authors[0]} 等")
        else:
            ref_parts.append(", ".join(authors))
    if result.get("期刊"):
        ref_parts.append(result["期刊"])
    if result.get("发表年份"):
        ref_parts.append(result["发表年份"])
    pmid = result.get("PMID", "")
    ref_str = ". ".join(ref_parts) if ref_parts else "待补充"

    is_recessive = any(k in zygosity for k in ["复合杂合", "纯合", "compound heterozygous", "homozygous", "双等位基因"])

    parts = []

    # ── 1. 变异身份识别 ──
    identity_parts = [variant_display_short]
    if variant_type and variant_type != "未指明":
        var_type_cn = variant_type.replace("missense", "错义突变").replace("nonsense", "无义突变").replace("frameshift", "移码突变").replace("splicing", "剪接突变")
        identity_parts.append(f"是{gene}基因的{var_type_cn}")
    else:
        identity_parts.append(f"是{gene}基因突变")

    if clinical_details.get("exon"):
        exon_info = f"位于第{clinical_details['exon']}号外显子"
        if clinical_details.get("domain"):
            exon_info += f"（{clinical_details['domain']}）"
        identity_parts.append(exon_info)

    if clinical_details.get("disease_subtypes"):
        parts.append(f"{'，'.join(identity_parts)}；患者临床分型为{'/'.join(clinical_details['disease_subtypes'][:3])}")
    else:
        parts.append(f"{'，'.join(identity_parts)}")

    # ── 2. 患者/队列背景 ──
    patient_info = []
    if patient_details:
        for pd_ in patient_details[:3]:
            if pd_.get("patient_id"):
                patient_info.append(f"患者{pd_['patient_id']}")
                break
    if patient_count > 0 and not patient_info:
        patient_info.append(f"共{patient_count}例患者")
    if patient_info:
        parts.append(f"该变异出现于{'，'.join(patient_info)}")

    # ── 3. 合子状态 + 共存变异 ──
    if is_recessive:
        if "compound" in zygosity.lower() or "复合杂合" in zygosity:
            zyg_label = "呈复合杂合状态"
            if co_variants:
                co_strs = []
                for cv in co_variants[:3]:
                    if isinstance(cv, dict):
                        c = cv.get("cdna") or ""
                        p = cv.get("蛋白变异") or ""
                        if c and p:
                            co_strs.append(f"{c}（{p}）")
                        else:
                            co_strs.append(c or p or str(cv))
                    else:
                        co_strs.append(str(cv))
                zyg_label += f"，与{'/'.join(co_strs)}突变共存"
            parts.append(zyg_label)
        elif "homozygous" in zygosity.lower() or "纯合" in zygosity:
            parts.append("呈纯合状态")
        else:
            parts.append(f"合子状态：{zygosity}")
    else:
        if zygosity and zygosity != "未指明":
            parts.append(f"合子状态：{zygosity}")

    # ── 4. 相位信息 ──
    if phase_status == "confirmed_in_trans":
        if parental_testing:
            trans_detail = "经亲本检测确认处于反式位置（分别来自父母双方）"
            if maternal_var:
                trans_detail += f"，其中{variant_display}遗传自母亲"
            if paternal_var:
                trans_detail += f"，共分离变异遗传自父亲"
            parts.append(trans_detail)
        else:
            parts.append("文献明确两个变异处于反式位置")
    elif phase_status == "presumed_in_trans":
        parts.append("推定为反式位置（未经亲本验证）")
    elif phase_status == "confirmed_in_cis":
        if co_variants:
            co_strs = []
            for cv in co_variants[:3]:
                if isinstance(cv, dict):
                    c = cv.get("cdna") or ""
                    p = cv.get("蛋白变异") or ""
                    if c and p:
                        co_strs.append(f"{c}（{p}）")
                    else:
                        co_strs.append(c or p or str(cv))
                else:
                    co_strs.append(str(cv))
            parts.append(f"处于顺式位置（同一等位基因），与{'/'.join(co_strs)}构成复杂等位基因")
        else:
            parts.append("处于顺式位置（同一等位基因），构成复杂等位基因")

    # ── 5. 发病年龄 + 临床表型 ──
    if clinical_details.get("onset_age"):
        parts.append(clinical_details["onset_age"])

    if phenotype and phenotype != "未指明":
        pheno_list = phenotype.split("、") if "、" in phenotype else [phenotype]
        pheno_display = "、".join(pheno_list[:6])
        parts.append(f"临床表现为{pheno_display}")

    if clinical_details.get("progression"):
        parts.append(clinical_details["progression"])

    # ── 6. 实验室发现 ──
    if clinical_details.get("lab_findings"):
        lab_items = clinical_details["lab_findings"]
        if lab_items:
            parts.append("；".join(lab_items[:3]))

    # ── 7. 变异频率/新发性 ──
    if clinical_details.get("frequency_info"):
        freq_items = clinical_details["frequency_info"]
        parts.append("；".join(freq_items[:3]))
    elif var_char.get("是否为新型变异"):
        parts.append("该变异为本文首次报道的新型突变")
    if var_char.get("人群频率信息"):
        parts.append("在正常对照人群中未检出")

    # ── 8. 功能验证 ──
    if func_valid and func_valid != "未进行":
        parts.append(f"功能实验（{func_valid}）提示可能影响蛋白功能")
    if var_char.get("转染/功能实验"):
        parts.append("转染/功能实验支持该突变的致病性")

    # ── 9. NMD相关 ──
    if var_char.get("NMD相关信息"):
        nmd_text = var_char.get("NMD原文", "").lower()
        if any(neg in nmd_text for neg in ["not undergo", "did not undergo", "no evidence", "no nmd"]):
            parts.append("该突变虽产生提前终止密码子，但未引发无义介导的mRNA降解")
        else:
            parts.append("该突变产生提前终止密码子，可能引发无义介导的mRNA降解")

    if var_char.get("生物信息学预测"):
        pred_tools = ", ".join(var_char["生物信息学预测"].keys())
        parts.append(f"生物信息学预测（{pred_tools}）提示该变异可能有害")

    # ── 组装 ──
    paragraph = "。".join(parts) + f"。参考文献：{ref_str}。[PMID:{pmid}]" if parts else f"未获取到足够信息。参考文献：{ref_str}。[PMID:{pmid}]"

    return paragraph


def _generate_literature_summary(result):
    """未提及目标变异时，从文献标题/摘要/全文生成 >200 字纯中文文献简介。

    使用预置的 TOPIC_KEYWORD_MAP 进行主题词匹配，确保每次运行结果可重复。
    不输出英文原文，全部使用中文关键词和结构化描述。
    """
    title = result.get("标题", "").strip()
    abstract = result.get("摘要", "").strip()
    fulltext = result.get("全文", "").strip()
    journal = result.get("期刊", "").strip()
    year = result.get("发表年份", "").strip()
    authors = result.get("作者", [])
    mesh_terms = result.get("MeSH术语", [])

    text = (title + " " + abstract + " " + fulltext[:5000])
    text_lower = text.lower()
    title_lower = title.lower()

    # 1. 主题关键词匹配
    matched_topics = []
    seen = set()
    for en, zh in TOPIC_KEYWORD_MAP.items():
        if en in title_lower and zh not in seen:
            matched_topics.append(zh)
            seen.add(zh)
    if len(matched_topics) < 6:
        for en, zh in TOPIC_KEYWORD_MAP.items():
            if en in text_lower and zh not in seen:
                matched_topics.append(zh)
                seen.add(zh)
                if len(matched_topics) >= 6:
                    break
    for term in mesh_terms[:5]:
        term_lower = term.lower()
        for en, zh in TOPIC_KEYWORD_MAP.items():
            if en in term_lower and zh not in seen:
                matched_topics.append(zh)
                seen.add(zh)
                if len(matched_topics) >= 8:
                    break

    # 2. 判断研究类型和研究对象
    study_type = ""
    if any(w in text_lower for w in ["case report", "case series"]):
        study_type = "病例报告/病例系列"
    elif any(w in text_lower for w in ["meta-analysis", "systematic review"]):
        study_type = "Meta分析/系统综述"
    elif any(w in text_lower for w in ["carrier screening", "preconception"]):
        study_type = "携带者筛查"
    elif any(w in text_lower for w in ["prenatal", "fetal", "fetus", "antenatal"]):
        study_type = "产前诊断"
    elif any(w in text_lower for w in ["newborn screening"]):
        study_type = "新生儿筛查"
    elif any(w in text_lower for w in ["cohort"]):
        study_type = "队列研究"
    elif any(w in text_lower for w in ["whole-exome", "exome sequencing", "whole-genome", "genome sequencing"]):
        study_type = "基于高通量测序的基因检测研究"
    elif any(w in text_lower for w in ["gene panel", "panel testing", "targeted sequencing"]):
        study_type = "基因包检测研究"
    elif any(w in text_lower for w in ["genotype-phenotype"]):
        study_type = "基因型-表型关联分析"
    elif any(w in text_lower for w in ["functional"]):
        study_type = "功能学研究"
    elif any(w in text_lower for w in ["population", "carrier frequency", "genetic prevalence", "genetic architecture"]):
        study_type = "人群遗传学调查"
    elif any(w in text_lower for w in ["review", "overview"]):
        study_type = "综述"
    else:
        study_type = "遗传学/临床研究"

    # 研究对象推断
    population = ""
    for kw, label in [("chinese", "中国"), ("japanese", "日本"), ("korean", "韩国"),
                       ("thai", "泰国"), ("turkish", "土耳其"), ("iranian", "伊朗"),
                       ("dutch", "荷兰"), ("german", "德国"), ("french", "法国"),
                       ("spanish", "西班牙"), ("italian", "意大利"), ("british", "英国"),
                       ("polish", "波兰"), ("danish", "丹麦"), ("russian", "俄罗斯"),
                       ("indian", "印度"), ("african", "非洲")]:
        if kw in text_lower:
            population = label
            break

    # 样本量推断
    sample_size = ""
    for m in re.finditer(r'(\d[\d,]*)\s*(patients?|subjects?|individuals?|cases?|fetuses?|families?|exomes?|genomes?|samples?)', text_lower):
        num_str = m.group(1).replace(",", "")
        try:
            n = int(num_str)
            if 5 <= n <= 1000000:
                label_map = {"patient": "例患者", "patients": "例患者",
                             "subject": "例受试者", "subjects": "例受试者",
                             "individual": "例个体", "individuals": "例个体",
                             "case": "个病例", "cases": "个病例",
                             "fetus": "例胎儿", "fetuses": "例胎儿",
                             "family": "个家系", "families": "个家系",
                             "exome": "例外显子组", "exomes": "例外显子组",
                             "genome": "例基因组", "genomes": "例基因组",
                             "sample": "例样本", "samples": "例样本"}
                label = label_map.get(m.group(2), "例")
                sample_size = f"共纳入{n}{label}"
                break
        except ValueError:
            continue

    # 3. 组装纯中文简介
    parts = []

    # 研究主题
    if matched_topics:
        parts.append("该文献为" + "、".join(matched_topics[:6]) + "相关研究")

    # 研究类型
    if study_type:
        parts.append(f"研究类型为{study_type}")
    if population:
        parts.append(f"研究对象为{population}人群")
    if sample_size:
        parts.append(sample_size)

    # 期刊信息
    author_str = ""
    if authors:
        first_author = authors[0].split()[-1] if authors else ""
        if first_author and len(first_author) > 1:
            et_al = "等" if len(authors) > 1 else ""
            author_str = f"第一作者为{first_author}{et_al}"
    if journal and year:
        parts.append(f"发表于{year}年《{journal}》")
        if author_str:
            parts.append(author_str)
    elif year:
        parts.append(f"发表于{year}年")
        if author_str:
            parts.append(author_str)
    elif author_str:
        parts.append(author_str)

    # MeSH 术语
    if mesh_terms:
        mesh_cn = []
        for term in mesh_terms[:5]:
            term_lower = term.lower()
            for en, zh in TOPIC_KEYWORD_MAP.items():
                if en in term_lower and zh not in mesh_cn and zh not in seen:
                    mesh_cn.append(zh)
                    break
            else:
                mesh_cn.append(term)
        if mesh_cn:
            parts.append(f"MeSH关键词包括{'、'.join(mesh_cn[:5])}")

    # 标题关键信息（从标题中提取中文关键词描述，不使用英文原文）
    if title:
        title_keywords = []
        title_lower2 = title.lower()
        # 提取标题中的基因名
        for gene_candidate in re.findall(r'\b([A-Z][A-Z0-9]{2,}(?:\s*and\s*[A-Z][A-Z0-9]{2,})?)\b', title):
            if gene_candidate.lower() not in ("the", "and", "for", "not", "but", "are", "was", "were", "both",
                                               "type", "without", "with", "from", "that", "this", "role", "novel",
                                               "method", "case", "gene", "study", "analysis", "short", "rib",
                                               "major", "renal", "retinal", "common", "cause", "exome"):
                if gene_candidate not in title_keywords:
                    title_keywords.append(gene_candidate)
        if title_keywords:
            parts.append(f"文献标题涉及基因{'、'.join(title_keywords[:4])}")

    # 研究内容描述（基于主题词推断，不使用英文原文）
    content_desc_parts = []
    if matched_topics:
        content_desc_parts.append(f"围绕{'、'.join(matched_topics[:4])}展开")
    if study_type == "携带者筛查":
        content_desc_parts.append("分析相关致病变异在人群中的携带频率与分布特征")
    elif study_type == "产前诊断":
        content_desc_parts.append("通过影像学和分子检测手段对胎儿进行产前评估与遗传学诊断")
    elif study_type == "基因包检测研究":
        content_desc_parts.append("采用靶向基因包对相关疾病基因进行系统性检测与变异分析")
    elif study_type in ("基于高通量测序的基因检测研究", "基因检测/测序研究"):
        content_desc_parts.append("利用高通量测序技术系统鉴定致病基因变异并分析基因型-表型关联")
    elif study_type == "队列研究":
        content_desc_parts.append("系统收集病例队列进行遗传学分析和临床特征总结")
    elif study_type == "病例报告/病例系列":
        content_desc_parts.append("报道新发或罕见病例的临床特征及遗传学发现")
    elif study_type == "人群遗传学调查":
        content_desc_parts.append("在特定人群中进行遗传流行病学调查和致病变异谱分析")
    elif study_type == "功能学研究":
        content_desc_parts.append("通过体内外实验探讨相关基因和变异的分子功能与致病机制")
    elif study_type == "Meta分析/系统综述":
        content_desc_parts.append("系统检索并整合已有文献数据，进行定量或定性综合分析")
    elif study_type == "综述":
        content_desc_parts.append("对相关领域的研究进展、致病机制和临床管理进行系统回顾与总结")
    else:
        content_desc_parts.append("通过遗传学方法研究疾病的分子基础、变异谱和临床特征")

    if content_desc_parts:
        parts.append("，".join(content_desc_parts))

    summary = "。".join(parts) + "。"

    # 确保 >200 字符：不足时从标题关键词和主题词补充
    if len(summary) < 200:
        extra_parts = []
        # 补充更多主题关键词
        extra_topics = [t for t in matched_topics if t not in summary][:5]
        if extra_topics:
            extra_parts.append(f"还涉及{'、'.join(extra_topics)}等领域")
        # 期刊信息
        if journal and journal not in summary:
            extra_parts.append(f"发表于《{journal}》")
        if extra_parts:
            summary = summary[:-1] + "。" + "。".join(extra_parts) + "。"

    # 仍不足时用标题中文翻译描述补充
    if len(summary) < 200:
        # 用已匹配的主题词反推文献内容方向
        if matched_topics:
            more = f"本文献对于理解{'、'.join(matched_topics[:3])}的遗传基础、分子机制及临床管理具有参考价值"
        else:
            more = "本文献为该领域的遗传学研究和临床实践提供了有价值的参考数据"
        summary = summary[:-1] + "。" + more + "。"

    return summary


def generate_one_sentence_summary(result):
    """生成一句话概括。未提及变异时委托给文献简介生成器（>200 字）。"""
    mentioned = result.get("变异提及", False)
    if not mentioned:
        return _generate_literature_summary(result)

    gene = result.get("基因", "").strip()
    cdna = result.get("cDNA变异", "").strip()
    phenotype = result.get("临床表型", "")
    zygosity = result.get("合子状态", "")

    disease = ""
    if phenotype and phenotype != "未指明":
        pheno_items = [p.strip() for p in phenotype.split("、") if p.strip()]
        disease = pheno_items[0] if pheno_items else ""

    parts = []
    if gene:
        parts.append(gene)
    if cdna:
        parts.append(cdna)

    extra = ""
    if zygosity and zygosity != "未指明":
        extra = f"（{zygosity}）"

    if disease:
        return f"{disease}相关研究，{' '.join(parts)}{extra}。"
    return f"{' '.join(parts)}{extra}相关遗传学分析。" if parts else "无法生成概括"


def _generate_study_background(result):
    """生成文献背景描述（是什么研究），所有文献通用。

    基于标题、摘要、MeSH术语和TOPIC_KEYWORD_MAP生成约50-100字中文描述。
    预置逻辑确保每次运行结果可重复。
    """
    title = result.get("标题", "").strip()
    abstract = result.get("摘要", "").strip()
    fulltext = result.get("全文", "").strip()
    journal = result.get("期刊", "").strip()
    year = result.get("发表年份", "").strip()
    mesh_terms = result.get("MeSH术语", [])
    authors = result.get("作者", [])

    text = (title + " " + abstract[:2000] + " " + fulltext[:2000])
    text_lower = text.lower()
    title_lower = title.lower()

    # 提取主题关键词
    matched_topics = []
    seen = set()
    for en, zh in TOPIC_KEYWORD_MAP.items():
        if en in title_lower and zh not in seen:
            matched_topics.append(zh)
            seen.add(zh)
    if len(matched_topics) < 4:
        for en, zh in TOPIC_KEYWORD_MAP.items():
            if en in text_lower and zh not in seen:
                matched_topics.append(zh)
                seen.add(zh)
                if len(matched_topics) >= 4:
                    break

    # 判断研究类型
    study_type = ""
    if any(w in text_lower for w in ["case report", "case series"]):
        study_type = "病例报告/病例系列"
    elif any(w in text_lower for w in ["cohort", "cohort study"]):
        study_type = "队列研究"
    elif any(w in text_lower for w in ["meta-analysis", "systematic review"]):
        study_type = "Meta分析/系统综述"
    elif any(w in text_lower for w in ["carrier screening", "preconception"]):
        study_type = "携带者筛查研究"
    elif any(w in text_lower for w in ["prenatal", "fetal", "fetus", "antenatal"]):
        study_type = "产前诊断研究"
    elif any(w in text_lower for w in ["newborn screening"]):
        study_type = "新生儿筛查研究"
    elif any(w in text_lower for w in ["whole-exome", "exome sequencing", "whole-genome", "genome sequencing", "gene panel", "panel testing", "next-generation"]):
        study_type = "基因检测/测序研究"
    elif any(w in text_lower for w in ["genotype-phenotype", "genotype phenotype"]):
        study_type = "基因型-表型关联研究"
    elif any(w in text_lower for w in ["functional", "functional analysis", "functional characterization"]):
        study_type = "功能学研究"
    elif any(w in text_lower for w in ["population", "genetic architecture", "carrier frequency", "genetic prevalence"]):
        study_type = "人群遗传学研究"
    elif any(w in text_lower for w in ["mutation", "mutational", "variant"]):
        study_type = "突变筛查/变异分析研究"
    elif any(w in text_lower for w in ["review", "overview"]):
        study_type = "综述"
    else:
        study_type = "遗传学/临床研究"

    # 组装背景
    parts = []
    if study_type:
        parts.append(study_type)
    if matched_topics:
        parts.append("涉及" + "、".join(matched_topics[:4]))
    if journal and year:
        parts.append(f"发表于{year}年《{journal}》")

    if not parts:
        return title[:100] if title else "文献信息不足"

    return "，".join(parts) + "。"


def _generate_excel_csv(results, output_path):
    """生成 Excel 兼容的 CSV 文件（UTF-8 BOM，零外部依赖，Excel/WPS 直接打开不乱码）。

    预置列：PMID、标题、是否提及此位点、患者数、致病性、关联合子状态、
            反式(trans)/顺式(cis)位点、患者临床表型、文献背景、总结
    """
    # CSV 字段转义：含逗号/换行/引号时用双引号包裹
    def esc(s):
        if s is None:
            return ""
        return str(s).replace("\n", " ").replace("\r", "")

    headers = [
        "PMID", "标题", "是否提及此位点", "患者数", "致病性",
        "关联合子状态", "反式(trans)/顺式(cis)位点", "患者临床表型",
        "文献背景(是什么研究)", "总结"
    ]

    rows = []
    for r in results:
        pmid = r.get("PMID", "")
        title = r.get("标题", "")
        mentioned = "是" if r.get("变异提及") else "否"
        patient_count = str(r.get("患者数量", 0) or 0)
        pathogenicity = r.get("致病性", "") or ""
        zygosity = r.get("合子状态", "") or ""

        # 反式(trans)/顺式(cis)位点 — v9: 显示 trans/cis 相位信息和共存变异详情
        phase_status = r.get("相位状态", "")
        trans_confirmed = r.get("反式确认", False)
        cis_confirmed = r.get("顺式确认", False)
        co_vars = r.get("共存变异", [])

        trans_parts = []

        # (a) 相位标签
        if phase_status == "confirmed_in_trans":
            trans_parts.append("确认反式")
        elif phase_status == "presumed_in_trans":
            trans_parts.append("推定反式")
        elif phase_status == "confirmed_in_cis":
            trans_parts.append("确认顺式")
        elif phase_status == "presumed_in_cis":
            trans_parts.append("推定顺式")
        elif trans_confirmed:
            trans_parts.append("确认反式")
        elif cis_confirmed:
            trans_parts.append("确认顺式")

        # (b) 共存变异 — 始终显示
        if co_vars:
            co_strs = []
            for cv in co_vars:
                if isinstance(cv, dict):
                    cdna_v = cv.get("cdna") or ""
                    prot_v = cv.get("蛋白变异") or ""
                    if cdna_v and prot_v:
                        prot_short = str(prot_v).replace("p.", "")
                        co_strs.append(f"{cdna_v}(p.{prot_short})")
                    elif cdna_v:
                        co_strs.append(cdna_v)
                    elif prot_v:
                        co_strs.append(str(prot_v))
                else:
                    co_strs.append(str(cv))
            if co_strs:
                trans_parts.append("; ".join(co_strs[:3]))

        if not trans_parts:
            if r.get("变异提及"):
                trans_parts.append("未确认")
            else:
                trans_parts.append("不适用")

        # (c) 亲本源
        if r.get("亲本检测"):
            maternal = r.get("母源变异") or ""
            paternal = r.get("父源变异") or ""
            if maternal:
                trans_parts.append(f"母源: {maternal}")
            if paternal:
                trans_parts.append(f"父源: {paternal}")

        trans_display = "; ".join(trans_parts) if trans_parts else "不适用"

        phenotype = r.get("临床表型", "") or ""
        if phenotype == "未指明":
            phenotype = ""

        background = _generate_study_background(r)

        summary = r.get("总结段落", "") or ""
        if not summary and not r.get("变异提及"):
            summary = r.get("一句话概括", "") or ""

        rows.append([
            esc(pmid),
            esc(title[:200]),
            esc(mentioned),
            esc(patient_count),
            esc(pathogenicity[:200]),
            esc(zygosity[:200]),
            esc(trans_display[:500]),
            esc(phenotype[:300]),
            esc(background[:500]),
            esc(summary[:2000]),
        ])

    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else ".", exist_ok=True)

    # UTF-8 BOM + CSV — Excel/WPS 原生支持，零依赖
    with open(output_path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)

    return output_path