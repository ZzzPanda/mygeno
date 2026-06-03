#!/usr/bin/env python3
"""
PubMed Variant Extractor v9
从PubMed/PMC文献中仅针对目标变异提取信息，生成标准化中文总结段落。

支持转录本版本差异识别（如 NM_022436.2 vs NM_022436.3 编号差异）。
v9: 修复 in-cis 识别（连字符匹配）、顺式复杂等位基因共存变异提取、Excel 列显示顺式标签。
纯Python标准库，零依赖。

Usage:
  python pubmed_extractor.py --pmids PMID1 PMID2 ... --gene GENE \
      --variant "c.NNN X>Y (p.AA)" [--transcript NM_xxx.x] [--output path.json]
"""

import argparse
import json
import time
import random
import urllib.request
import urllib.parse
import urllib.error
import xml.etree.ElementTree as ET
import sys
import re
import os
from datetime import date

# ── 编码 ──
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')
if sys.stderr.encoding != 'utf-8':
    sys.stderr.reconfigure(encoding='utf-8')

# ── 配置 ──
MAX_RETRIES = 2
DAILY_SITE_LIMIT = 500
PAUSE_AFTER_N = 15
EXTRA_PAUSE_SECS = 60
MIN_DELAY = 10
MAX_DELAY = 40

# ── 表型关键词映射 ──
PHENOTYPE_MAP = {
    "sitosterolemia": "谷固醇血症",
    "phytosterolemia": "植物固醇血症",
    "xanthoma": "黄色瘤",
    "hypercholesterolemia": "高胆固醇血症",
    "hyperlipidemia": "高脂血症",
    "atherosclerosis": "动脉粥样硬化",
    "coronary artery disease": "冠状动脉疾病",
    "myocardial infarction": "心肌梗死",
    "premature coronary": "早发冠心病",
    "cardiovascular disease": "心血管疾病",
    "stroke": "卒中",
    "hepatomegaly": "肝肿大",
    "splenomegaly": "脾肿大",
    "thrombocytopenia": "血小板减少",
    "anemia": "贫血",
    "stomatocyte": "口形红细胞",
    "arthritis": "关节炎",
    "growth retardation": "生长迟缓",
    "failure to thrive": "生长发育迟缓",
    "developmental delay": "发育迟缓",
    "seizure": "癫痫",
    "encephalopathy": "脑病",
    "hypotonia": "肌张力低下",
    "psychomotor retardation": "精神运动发育迟缓",
    "intellectual disability": "智力障碍",
    "mental retardation": "智力低下",
    "neurodegenerative": "神经退行性",
    "retinal dystrophy": "视网膜变性",
    "retinitis pigmentosa": "视网膜色素变性",
    "vision loss": "视力丧失",
    "visual impairment": "视力障碍",
    "blindness": "失明",
    "macular dystrophy": "黄斑变性",
    "macular atrophy": "黄斑萎缩",
    "chorioretinal atrophy": "脉络膜视网膜萎缩",
    "bull.s eye maculopathy": "牛眼样黄斑病变",
    "bull's eye": "牛眼样黄斑病变",
    "yellow-white flecks": "黄白色斑点",
    "beaten-bronze": "青铜样黄斑",
    "rpe atrophy": "RPE萎缩",
    "outer nuclear layer": "外核层变薄",
    "onl thinning": "外核层变薄",
    "erg group": "ERG分组异常",
    "electroretinogram": "视网膜电图异常",
    "flecks": "视网膜斑点",
    "chorioretinal": "脉络膜视网膜病变",
    "fundus flavimaculatus": "眼底黄色斑点",
    "bull's eye maculopathy": "牛眼样黄斑病变",
    "attenuation of retinal": "视网膜变薄",
    "atrophic macular": "黄斑萎缩",
    "central scotoma": "中心暗点",
    "color vision": "色觉异常",
    "nyctalopia": "夜盲",
    "photophobia": "畏光",
    "retinal flecks": "视网膜斑点",
    "pisciform": "鱼形斑点",
    "ataxia": "共济失调",
    "cerebellar ataxia": "小脑共济失调",
    "cerebellar atrophy": "小脑萎缩",
    "brain atrophy": "脑萎缩",
    "cognitive decline": "认知衰退",
    "epilepsy": "癫痫",
    "psychiatric": "精神症状",
    "liver failure": "肝衰竭",
    "cirrhosis": "肝硬化",
    "jaundice": "黄疸",
    "cholestasis": "胆汁淤积",
    "feeding difficulty": "喂养困难",
    "vomiting": "呕吐",
    "diarrhea": "腹泻",
    "lethargy": "嗜睡",
    "metabolic crisis": "代谢危象",
    "acidosis": "酸中毒",
    "coma": "昏迷",
    "fatigue": "乏力",
    "short stature": "身材矮小",
    "joint pain": "关节痛",
    "stunted growth": "生长迟缓",
}

# ── 疾病亚型关键词映射 ──
DISEASE_SUBTYPE_MAP = {
    "neonatal": "新生儿型", "perinatal": "围产期型",
    "early infantile": "早发型婴儿型", "late infantile": "迟发型婴儿型",
    "juvenile": "青少年型", "adolescent": "青少年型",
    "adult onset": "成人型", "adult-onset": "成人型",
    "classic": "经典型（CLASSIC）", "variant": "变异型（VARIANT）",
    "biochemical variant": "生化变异型",
    "mild": "轻型", "severe": "重型", "moderate": "中度",
    "intermediate": "中间型", "severe infantile": "重症婴儿型",
    "early infantile systemic lethal": "早发型婴儿系统性致死型",
}

# ── 实验室/生化检测关键词映射 ──
LAB_FINDINGS_MAP = {
    "filipin": "Filipin染色",
    "oxysterol": "氧化固醇",
    "cholestane": "胆甾烷三醇",
    "7-ketocholesterol": "7-酮胆固醇（7-KC）",
    "lysosphingomyelin": "溶血鞘磷脂",
    "chitotriosidase": "壳三糖酶",
    "sphingomyelinase": "鞘磷脂酶",
    "acid sphingomyelinase": "酸性鞘磷脂酶",
    "biomarker": "生物标志物",
    "bone marrow": "骨髓",
    "foam cell": "泡沫细胞",
    "sea blue": "海蓝组织细胞",
    "plasma": "血浆",
    "serum": "血清",
}

# ── 文献主题英文→中文关键词映射（用于未提及变异时生成 >200 字文献简介）──
TOPIC_KEYWORD_MAP = {
    # 骨骼 / 纤毛病
    "skeletal dysplasia": "骨骼发育不良",
    "skeletal ciliopath": "骨骼纤毛病",
    "short rib": "短肋",
    "short-rib": "短肋",
    "thoracic dystroph": "胸廓发育不良",
    "jeune": "Jeune综合征",
    "asphyxiating thoracic": "窒息性胸廓发育不良",
    "polydactyly": "多指(趾)",
    "skeletal disorder": "骨骼疾病",
    "osteogenesis": "成骨不全",
    "achondroplasia": "软骨发育不全",
    "bone dysplasia": "骨发育不良",
    "osteoporosis": "骨质疏松",
    "arthrogryposis": "关节挛缩",
    "scoliosis": "脊柱侧弯",
    # 纤毛病 / 肾脏
    "ciliopath": "纤毛病",
    "ciliary": "纤毛",
    "cilia": "纤毛",
    "intraflagellar": "鞭毛内转运",
    "nephronophthisis": "肾单位肾痨",
    "polycystic kidney": "多囊肾",
    "renal": "肾脏",
    "kidney": "肾脏",
    "cystic kidney": "囊性肾病",
    "bardet-biedl": "Bardet-Biedl综合征",
    "joubert": "Joubert综合征",
    "meckel": "Meckel综合征",
    "oral-facial-digital": "口面指综合征",
    "orofaciodigital": "口面指综合征",
    # DYNC2H1 / 动力蛋白
    "dync2h1": "DYNC2H1",
    "dynein": "动力蛋白",
    "cytoplasmic dynein": "胞质动力蛋白",
    "ift": "鞭毛内转运",
    "hedgehog": "Hedgehog信号",
    "sonic hedgehog": "SHH信号",
    # 视网膜 / 眼科
    "retinal degeneration": "视网膜变性",
    "retinal dystroph": "视网膜营养不良",
    "retinitis pigmentosa": "视网膜色素变性",
    "stargardt": "Stargardt病",
    "macular degeneration": "黄斑变性",
    "macular dystroph": "黄斑营养不良",
    "inherited retinal": "遗传性视网膜病",
    "retina": "视网膜",
    "photoreceptor": "光感受器",
    "blindness": "失明",
    "vision loss": "视力丧失",
    "visual impairment": "视力障碍",
    "cone-rod": "锥杆细胞",
    "rod-cone": "杆锥细胞",
    "ophthalmol": "眼科",
    "fundus": "眼底",
    "electroretinogram": "视网膜电图",
    "usher syndrome": "Usher综合征",
    "usher": "Usher综合征",
    # 神经 / 发育
    "neurodegenerative": "神经退行性",
    "ataxia": "共济失调",
    "cerebellar": "小脑",
    "intellectual disability": "智力障碍",
    "developmental delay": "发育迟缓",
    "seizure": "癫痫",
    "epilepsy": "癫痫",
    "encephalopathy": "脑病",
    "hypotonia": "肌张力低下",
    "microcephaly": "小头畸形",
    "autism": "自闭症",
    "parkinson": "帕金森",
    "alzheimer": "阿尔茨海默",
    "neuropathy": "神经病变",
    "spina bifida": "脊柱裂",
    "neural tube": "神经管",
    # 代谢 / 内分泌
    "metabolic disease": "代谢性疾病",
    "inborn error": "先天性代谢异常",
    "diabetes": "糖尿病",
    "mody": "MODY型糖尿病",
    "hyperinsulinism": "高胰岛素血症",
    "hypoglycemia": "低血糖",
    "obesity": "肥胖",
    "hypercholesterolemia": "高胆固醇血症",
    "hyperlipidemia": "高脂血症",
    "homocystinuria": "同型半胱氨酸尿症",
    "phenylketonuria": "苯丙酮尿症",
    "gaucher": "戈谢病",
    "pompe": "庞贝病",
    "fabry": "法布里病",
    "mucopolysaccharidosis": "黏多糖贮积症",
    "lysosomal": "溶酶体",
    "mitochondrial": "线粒体",
    # 心脏 / 血管
    "cardiovascular": "心血管",
    "cardiac": "心脏",
    "congenital heart": "先天性心脏",
    "cardiomyopathy": "心肌病",
    "arrhythmia": "心律失常",
    "long qt": "长QT综合征",
    "brugada": "Brugada综合征",
    "hypertension": "高血压",
    "thrombosis": "血栓",
    "stroke": "卒中",
    "aneurysm": "动脉瘤",
    # 血液
    "hemolytic": "溶血性",
    "anemia": "贫血",
    "thalassemia": "地中海贫血",
    "sickle cell": "镰状细胞",
    "hemophilia": "血友病",
    "thrombocytopenia": "血小板减少",
    "coagulation": "凝血",
    # 肝脏 / 消化
    "cholestasis": "胆汁淤积",
    "liver disease": "肝病",
    "cirrhosis": "肝硬化",
    "hepatic": "肝脏",
    "bile acid": "胆汁酸",
    "pfic": "进行性家族性肝内胆汁淤积",
    "jaundice": "黄疸",
    "hepatitis": "肝炎",
    # 皮肤
    "ichthyosis": "鱼鳞病",
    "epidermolysis": "大疱性表皮松解",
    "ectodermal dysplasia": "外胚层发育不良",
    "albinism": "白化病",
    # 免疫
    "immunodeficiency": "免疫缺陷",
    "autoimmune": "自身免疫",
    "scid": "重症联合免疫缺陷",
    # 肿瘤 / 癌症
    "hereditary cancer": "遗传性肿瘤",
    "cancer predisposition": "肿瘤易感",
    "breast cancer": "乳腺癌",
    "ovarian cancer": "卵巢癌",
    "colorectal cancer": "结直肠癌",
    "lynch syndrome": "Lynch综合征",
    "li-fraumeni": "Li-Fraumeni综合征",
    "retinoblastoma": "视网膜母细胞瘤",
    "neurofibromatosis": "神经纤维瘤病",
    "tuberous sclerosis": "结节性硬化症",
    "von hippel": "Von Hippel-Lindau综合征",
    # 产前 / 生殖
    "prenatal": "产前",
    "prenatal diagnosis": "产前诊断",
    "fetal": "胎儿",
    "fetus": "胎儿",
    "ultrasound": "超声",
    "antenatal": "产前",
    "carrier screening": "携带者筛查",
    "preconception": "孕前",
    "infertility": "不孕不育",
    "recurrent miscarriage": "复发性流产",
    "preimplantation": "植入前遗传学检测",
    # 方法学 / 通用
    "whole-exome": "全外显子组测序",
    "whole-genome": "全基因组测序",
    "exome sequencing": "外显子组测序",
    "genome sequencing": "基因组测序",
    "gene panel": "基因包",
    "panel testing": "基因包检测",
    "next-generation": "二代测序",
    "targeted sequencing": "靶向测序",
    "sanger": "Sanger测序",
    "genotype-phenotype": "基因型-表型",
    "genotype": "基因型",
    "phenotype": "表型",
    "mutational": "突变",
    "molecular diagnosis": "分子诊断",
    "genetic testing": "基因检测",
    "genetic diagnosis": "遗传学诊断",
    "functional characterization": "功能鉴定",
    "functional analysis": "功能分析",
    "functional": "功能",
    "diagnostic": "诊断",
    "prognosis": "预后",
    "natural history": "自然史",
    "newborn screening": "新生儿筛查",
    "population": "人群",
    "carrier frequency": "携带者频率",
    "genetic prevalence": "遗传患病率",
    "genetic architecture": "遗传架构",
    "association": "关联分析",
    "meta-analysis": "Meta分析",
    "systematic review": "系统综述",
    "cohort": "队列",
    "case report": "病例报告",
    "case series": "病例系列",
    "variant interpretation": "变异解读",
    "acmg": "ACMG",
    "pathogenic": "致病性",
    "novel variant": "新型变异",
    # 地域/人群
    "chinese": "中国",
    "china": "中国",
    "japanese": "日本",
    "korean": "韩国",
    "thai": "泰国",
    "indian": "印度",
    "turkish": "土耳其",
    "iranian": "伊朗",
    "dutch": "荷兰",
    "german": "德国",
    "french": "法国",
    "spanish": "西班牙",
    "italian": "意大利",
    "british": "英国",
    "polish": "波兰",
    "russian": "俄罗斯",
    "danish": "丹麦",
    "african": "非洲",
    "caucasian": "高加索",
    "hispanic": "西班牙裔",
    "ashkenazi": "阿什肯纳兹犹太",
}


# ── 每日站点访问计数 ──

def _get_state_file():
    today = date.today().isoformat()
    return os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        f".site_counts_{today}.json")


def load_site_counts():
    path = _get_state_file()
    if os.path.exists(path):
        try:
            with open(path, "r") as f:
                return json.load(f)
        except (json.JSONDecodeError, OSError):
            pass
    return {}


def save_site_counts(counts):
    try:
        with open(_get_state_file(), "w") as f:
            json.dump(counts, f)
    except OSError:
        pass


def check_and_increment(site_key):
    counts = load_site_counts()
    current = counts.get(site_key, 0)
    if current >= DAILY_SITE_LIMIT:
        return False, current
    counts[site_key] = current + 1
    save_site_counts(counts)
    return True, current + 1


# ── 延迟 ──

def random_sleep():
    secs = random.randint(MIN_DELAY, MAX_DELAY)
    print(f"  等待 {secs}s...")
    time.sleep(secs)


def periodic_long_pause(request_count):
    if request_count > 0 and request_count % PAUSE_AFTER_N == 0:
        print(f"  [已达 {PAUSE_AFTER_N} 次，额外暂停 {EXTRA_PAUSE_SECS}s]")
        time.sleep(EXTRA_PAUSE_SECS)


# ── 网络请求 ──

def safe_request(url, description="", retries=MAX_RETRIES):
    req = urllib.request.Request(url)
    req.add_header("User-Agent", "PubMedVariantExtractor/7.0")
    if description:
        print(f"  -> {description}")
    for attempt in range(retries + 1):
        try:
            if attempt > 0:
                wait = random.randint(20, 60) * attempt
                print(f"  重试 {attempt}/{retries}，等待 {wait}s...")
                time.sleep(wait)
            with urllib.request.urlopen(req, timeout=45) as response:
                return response.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(random.randint(30, 90))
                continue
        except Exception:
            pass
        if attempt < retries:
            time.sleep(random.randint(10, 40))
    return None


def fetch_ncbi_abstract(pmid):
    site = "eutils.ncbi.nlm.nih.gov"
    allowed, count = check_and_increment(site)
    if not allowed:
        print(f"\n  [警告] {site} 今日已访问 {count} 次（上限 {DAILY_SITE_LIMIT}）")
        return None
    params = {"db": "pubmed", "id": pmid, "rettype": "xml", "retmode": "text"}
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{urllib.parse.urlencode(params)}"
    return safe_request(url, f"NCBI abstract PMID:{pmid}")


def fetch_europe_pmc_fulltext(pmid):
    site = "www.ebi.ac.uk"
    allowed, count = check_and_increment(site)
    if not allowed:
        print(f"\n  [警告] {site} 今日已访问 {count} 次")
        return None
    check_url = (f"https://www.ebi.ac.uk/europepmc/webservices/rest/"
                 f"search?query=EXT_ID:{pmid}&format=json&resultType=core")
    data = safe_request(check_url, f"Europe PMC check PMID:{pmid}")
    if not data:
        return None
    try:
        result = json.loads(data)
        entries = result.get("resultList", {}).get("result", [])
        if not entries:
            return None
        entry = entries[0]
        if entry.get("isOpenAccess") == "Y" and entry.get("pmcid"):
            pmcid = entry["pmcid"]
            allowed2, _ = check_and_increment(site)
            if not allowed2:
                print(f"  [警告] {site} 今日已达上限")
                return None
            fulltext_url = (f"https://www.ebi.ac.uk/europepmc/webservices/rest/"
                           f"{pmcid}/fullTextXML")
            return safe_request(fulltext_url, f"Europe PMC XML {pmcid}")
    except (json.JSONDecodeError, KeyError):
        pass
    return None


def fetch_pmc_fulltext(pmid):
    """通过 PMC esearch 获取 PMC 全文 XML（通过 PMID 搜索对应 PMC 文章）"""
    site = "eutils.ncbi.nlm.nih.gov"
    allowed, count = check_and_increment(site)
    if not allowed:
        print(f"\n  [警告] {site} 今日已访问 {count} 次")
        return None
    # 使用 esearch 在 PMC 中搜索对应 PMID 的文章（比 elink 更准确）
    search_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                  f"?db=pmc&term={pmid}[pmid]&retmode=json")
    data = safe_request(search_url, f"PMC esearch PMID:{pmid}")
    if not data:
        return None
    try:
        result = json.loads(data)
        pmc_ids = result.get("esearchresult", {}).get("idlist", [])
        if not pmc_ids:
            return None
        xml_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                   f"?db=pmc&id={pmc_ids[0]}&rettype=xml&retmode=text")
        return safe_request(xml_url, f"PMC XML PMC{pmc_ids[0]}")
    except (json.JSONDecodeError, KeyError, IndexError):
        pass
    return None


# ── 转录本编号转换查询 ──

def fetch_ncbi_gene_info(transcript_id):
    """
    通过 NCBI Gene/RefSeq 接口查询转录本信息，
    返回该转录本对应的 CDS 起始位置等元数据，
    用于计算不同转录本版本间的 cDNA 编号偏移量。
    如果查询失败返回 None。
    """
    if not transcript_id:
        return None
    site = "eutils.ncbi.nlm.nih.gov"
    allowed, _ = check_and_increment(site)
    if not allowed:
        return None

    # 用 efetch 获取转录本的 FASTA + 元数据
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nucleotide&id={transcript_id}&rettype=fasta&retmode=text")
    data = safe_request(url, f"NCBI nucleotide {transcript_id}", retries=1)
    if not data:
        return None

    # 从 FASTA header 中提取 CDS 信息
    # Header 格式类似: >ref|NM_022436.3| Homo sapiens ATP binding cassette subfamily G member 5 (ABCG5), ...
    header_line = data.split('\n')[0] if data else ""
    info = {"transcript": transcript_id, "header": header_line}

    # 尝试用 esummary 获取更多结构化信息
    summary_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                   f"?db=nucleotide&id={transcript_id}&retmode=json")
    summary_data = safe_request(summary_url, f"NCBI summary {transcript_id}", retries=1)
    if summary_data:
        try:
            sresult = json.loads(summary_data)
            uid = list(sresult.get("result", {}).keys())
            if uid and uid[0] != "uids":
                summary_info = sresult["result"][uid[0]]
                info["title"] = summary_info.get("title", "")
                info["extra"] = summary_info.get("extra", "")
        except (json.JSONDecodeError, KeyError):
            pass

    return info


# ── 变异关键词扩展 ──

# 标准三字母 -> 单字母氨基酸映射
AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "ter": "*", "Xaa": "X", "X": "*",
    "*": "*",
}

AA_1TO3 = {v: k for k, v in AA_3TO1.items() if v != "X"}
AA_1TO3["*"] = "Ter"

# ── 已知变异历史命名 → HGVS 映射 ──
# 格式: { 基因符号: { 历史名称: {"cdna": "c.xxx", "protein": "p.xxx"} } }
KNOWN_VARIANT_NAMES = {
    "G6PD": {
        "Tsukui": {"cdna": "c.565_567del", "protein": "p.Ser189del"},
        "G6PD Tsukui": {"cdna": "c.565_567del", "protein": "p.Ser189del"},
    },
    # 可继续扩展其他基因的历史命名
}


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


# ── XML 解析 ──

def _empty_result():
    return {
        "PMID": "", "标题": "", "摘要": "", "全文": "",
        "全文来源": "abstract", "作者": [], "期刊": "",
        "发表年份": "", "MeSH术语": [], "一句话概括": "",
        "基因": "", "cDNA变异": "", "蛋白变异": "",
        "变异类型": "", "致病性": "", "合子状态": "",
        "患者详情": [], "临床表型": "",
        "遗传模式": "", "功能验证": "",
        "功能验证详情": [], "变异提及": False,
        "相关句子": [], "匹配关键词": [],
        "共存变异": [], "反式确认": False,
        # v7 新增：正反式/相位字段
        "顺式确认": False,
        "相位状态": "",           # confirmed_in_trans / confirmed_in_cis / presumed_in_trans / presumed_in_cis / phase_not_determined / not_assessed / not_applicable
        "相位置信度": "",         # confirmed / presumed / unknown / not_applicable
        "亲本检测": False,
        "母源变异": None,
        "父源变异": None,
        "相位证据": [],           # list of {type, label, sentence}
        "患者数量": 0, "变异特征": {},
        "表格数量": 0, "表格摘要": [],
        "总结段落": "", "临床详情": {},
    }


def parse_ncbi_xml(xml_text):
    result = _empty_result()
    if not xml_text:
        return None
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return None

    pmid_elem = root.find(".//PMID")
    if pmid_elem is not None:
        result["PMID"] = pmid_elem.text.strip()

    title_elem = root.find(".//ArticleTitle")
    if title_elem is not None:
        result["标题"] = "".join(title_elem.itertext()).strip()

    abstract_parts = []
    for elem in root.findall(".//AbstractText"):
        label = elem.get("Label", "")
        text = "".join(elem.itertext()).strip()
        abstract_parts.append(f"[{label}] {text}" if label else text)
    result["摘要"] = " ".join(abstract_parts)
    result["全文"] = result["摘要"]
    result["全文来源"] = "abstract"

    for author_elem in root.findall(".//Author"):
        ln = author_elem.find("LastName")
        fn = author_elem.find("ForeName")
        parts = []
        if ln is not None:
            parts.append(ln.text or "")
        if fn is not None:
            parts.append(fn.text or "")
        if parts:
            result["作者"].append(" ".join(parts))

    journal_elem = root.find(".//Journal/Title")
    if journal_elem is not None:
        result["期刊"] = journal_elem.text.strip()
    else:
        abbrev = root.find(".//Journal/ISOAbbreviation")
        if abbrev is not None:
            result["期刊"] = abbrev.text.strip()

    pub_date = root.find(".//PubDate/Year")
    if pub_date is not None:
        result["发表年份"] = pub_date.text.strip()
    else:
        medline = root.find(".//PubDate/MedlineDate")
        if medline is not None and medline.text:
            m = re.search(r'(\d{4})', medline.text)
            if m:
                result["发表年份"] = m.group(1)

    return result


def fetch_europe_pmc_text(pmid):
    """
    通过 Europe PMC 搜索 API 获取文章的摘要文本（作为 fallback）。
    当 NCBI 摘要和 PMC 全文都不可用时使用。
    返回 dict 包含标题、摘要、作者等信息；失败返回 None。
    """
    site = "www.ebi.ac.uk"
    allowed, count = check_and_increment(site)
    if not allowed:
        print(f"  [警告] {site} 今日已访问 {count} 次")
        return None
    url = (f"https://www.ebi.ac.uk/europepmc/webservices/rest/"
           f"search?query=EXT_ID:{pmid}&format=json&resultType=core")
    data = safe_request(url, f"Europe PMC text PMID:{pmid}", retries=1)
    if not data:
        return None
    try:
        result = json.loads(data)
        entries = result.get("resultList", {}).get("result", [])
        if not entries:
            return None
        entry = entries[0]
        pub_info = {
            "title": entry.get("title", ""),
            "abstract": entry.get("abstractText", ""),
            "journal": entry.get("journalTitle", ""),
            "year": str(entry.get("pubYear", "")),
            "authors": [],
            "pmid": pmid,
            "pmcid": entry.get("pmcid", ""),
            "isOpenAccess": entry.get("isOpenAccess", "N"),
        }
        author_str = entry.get("authorString", "")
        if author_str:
            pub_info["authors"] = [a.strip() for a in author_str.split(",") if a.strip()]
        return pub_info
    except (json.JSONDecodeError, KeyError):
        pass
    return None


def fetch_pmc_html(pmid):
    """
    从 PMC 网页获取文章的 HTML 全文（用于出版商限制 XML 时的 fallback）。
    PMC 网站对大部分文章提供免费 HTML 全文，即使 XML 下载受限。
    返回 HTML 文本；失败返回 None。
    """
    site = "www.ncbi.nlm.nih.gov"
    allowed, count = check_and_increment(site)
    if not allowed:
        print(f"  [警告] {site} 今日已访问 {count} 次")
        return None
    # 先通过 esearch 获取 PMCID
    search_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                  f"?db=pmc&term={pmid}[pmid]&retmode=json")
    data = safe_request(search_url, f"PMC esearch HTML PMID:{pmid}", retries=1)
    if not data:
        return None
    try:
        result = json.loads(data)
        pmc_ids = result.get("esearchresult", {}).get("idlist", [])
        if not pmc_ids:
            return None
        pmcid = pmc_ids[0]
    except (json.JSONDecodeError, KeyError, IndexError):
        return None

    # 获取 PMC HTML 页面
    html_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
    allowed2, _ = check_and_increment(site)
    if not allowed2:
        return None
    html_data = safe_request(html_url, f"PMC HTML {pmcid}", retries=1)
    return html_data


def parse_pmc_html(html_text, result):
    """从 PMC HTML 页面提取正文文本和表格数据。"""
    if not html_text or not result:
        return result

    # HTML 实体解码
    html_text = html_text.replace("&gt;", ">").replace("&lt;", "<").replace("&amp;", "&")
    html_text = html_text.replace("&quot;", "\"").replace("&#x0002C;", ",")
    html_text = html_text.replace("&nbsp;", " ")

    text_parts = []

    # 1. 提取 abstract
    abs_match = re.search(r'<div[^>]*class="[^"]*abstract[^"]*"[^>]*>(.*?)</div>', html_text, re.DOTALL | re.IGNORECASE)
    if not abs_match:
        abs_match = re.search(r'<section[^>]*id="[^"]*abstract[^"]*"[^>]*>(.*?)</section>', html_text, re.DOTALL | re.IGNORECASE)
    if abs_match:
        abs_text = re.sub(r'<[^>]+>', ' ', abs_match.group(1))
        abs_text = re.sub(r'\s+', ' ', abs_text).strip()
        if abs_text:
            text_parts.append(abs_text)

    # 2. 从 <section class="body main-article-body"> 提取正文（PMC 标准结构）
    body_match = re.search(
        r'<section[^>]*class="[^"]*\bbody\b[^"]*"[^>]*>(.*?)</section>',
        html_text, re.DOTALL | re.IGNORECASE
    )
    if body_match:
        # 排除表格部分（表格单独提取）
        body_html = re.sub(r'<table[^>]*>.*?</table>', '', body_match.group(1), flags=re.DOTALL | re.IGNORECASE)
        # 提取所有 <p> 标签文本
        paras = re.findall(r'<p[^>]*>(.*?)</p>', body_html, re.DOTALL)
        for p in paras:
            p_text = re.sub(r'<[^>]+>', ' ', p)
            p_text = re.sub(r'\s+', ' ', p_text).strip()
            if len(p_text) > 50:
                text_parts.append(p_text)

    # 3. 通用 fallback：从所有 <p> 标签提取（排除表头/页脚杂讯）
    if not text_parts:
        all_paras = re.findall(r'<p[^>]*>(.*?)</p>', html_text, re.DOTALL)
        for p in all_paras:
            p_text = re.sub(r'<[^>]+>', ' ', p)
            p_text = re.sub(r'\s+', ' ', p_text).strip()
            # 排除 PMC 界面杂讯
            if len(p_text) > 80 and not re.match(
                r'(An official website|The \.gov|Federal government|The site is secure|'
                r'Access keys|NCBI Homepage|MyNCBI|PubMed|PMC|Follow|Share|Connect|'
                r'Disclaimer|Copyright|FOIA|Privacy|NLM|National Library|NIH|HHS|'
                r'Vulnerability Disclosure|Accessibility|Careers|Nondiscrimination)',
                p_text, re.IGNORECASE
            ):
                text_parts.append(p_text)

    if text_parts:
        full_text = " ".join(text_parts)
        if len(full_text) > len(result.get("全文", "")):
            result["全文"] = full_text
            result["全文来源"] = "pmc_html"

    # 提取表格数据
    tables = _extract_tables_from_html(html_text)
    if tables:
        result["tables"] = tables

    return result


def _extract_tables_from_html(html_text):
    """从 HTML 文本提取表格（正则方式，零依赖）。"""
    tables = []
    # 找所有 table 元素
    table_matches = list(re.finditer(r'<table[^>]*>(.*?)</table>', html_text, re.DOTALL | re.IGNORECASE))
    for table_match in table_matches:
        table_content = table_match.group(1)
        # 查找 caption（可能在 table 之前或内部）
        caption = ""
        cap_match = re.search(r'<caption[^>]*>(.*?)</caption>', table_content, re.DOTALL | re.IGNORECASE)
        if cap_match:
            caption = re.sub(r'<[^>]+>', ' ', cap_match.group(1)).strip()
            caption = re.sub(r'\s+', ' ', caption)

        rows = []
        tr_matches = re.finditer(r'<tr[^>]*>(.*?)</tr>', table_content, re.DOTALL | re.IGNORECASE)
        for tr in tr_matches:
            tr_content = tr.group(1)
            cells = []
            for cell_match in re.finditer(r'<(?:td|th)[^>]*>(.*?)</(?:td|th)>', tr_content, re.DOTALL | re.IGNORECASE):
                cell_text = re.sub(r'<[^>]+>', ' ', cell_match.group(1))
                cell_text = re.sub(r'\s+', ' ', cell_text).strip()
                if re.match(r'<(?:td|th)[^>]*>', tr_content[:50]):
                    # 标记 header
                    cell_text = f"[H]{cell_text}"
                cells.append(cell_text)
            if cells:
                rows.append(cells)

        if rows:
            tables.append({
                "id": "",
                "caption": caption,
                "rows": rows,
                "raw_text": "\n".join(
                    "\t".join(c for c in row) for row in rows
                ),
            })
    return tables


def parse_pmc_europe_xml(xml_text, result):
    if not xml_text or not result:
        return result

    # 检测出版商限制（如 ASN 出版社会在 XML 中标注限制）
    if "does not allow downloading of the full text in XML form" in xml_text:
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            return result
        # 提取 front 部分的文本（标题、摘要等）
        front_texts = []
        for elem in root.findall(".//front//*"):
            if elem.text and len(elem.text.strip()) > 20:
                front_texts.append(elem.text.strip())
        if front_texts:
            result["全文"] = " ".join(front_texts)
            result["全文来源"] = "pmc_restricted"
        return result

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return result

    # 提取全文文本（优先 body，同时提取摘要）
    text_parts = []

    # 1. 先提取 PMC XML 中的摘要（front 部分）
    for abs_elem in root.findall(".//abstract/abstract-text"):
        txt = " ".join(abs_elem.itertext()).strip()
        if txt:
            text_parts.append(txt)
    if not text_parts:
        for abs_elem in root.findall(".//AbstractText"):
            txt = "".join(abs_elem.itertext()).strip()
            if txt:
                text_parts.append(txt)

    # 2. 提取 body 文本
    body_texts = []
    body = root.find(".//body")
    if body is not None:
        for elem in body.iter():
            if elem.text:
                body_texts.append(elem.text.strip())
            if elem.tail:
                body_texts.append(elem.tail.strip())

    # 3. 如果 body 不存在或为空，尝试从 <sec> 元素提取（非JATS标准结构）
    if not body_texts:
        for sec in root.findall(".//sec"):
            for elem in sec.iter():
                if elem.text and len(elem.text.strip()) > 10:
                    body_texts.append(elem.text.strip())
                if elem.tail and len(elem.tail.strip()) > 10:
                    body_texts.append(elem.tail.strip())
    if not body_texts:
        for elem in root.iter():
            if elem.tag not in ("article", "front", "body", "back"):
                if elem.text:
                    body_texts.append(elem.text.strip())
                if elem.tail:
                    body_texts.append(elem.tail.strip())

    # 合并：摘要 + body
    if body_texts:
        text_parts.extend(body_texts)

    # 如果 PMC 完全没有文本，保留已有的摘要
    if not text_parts:
        return result

    full_text = " ".join(filter(None, text_parts))
    if full_text:
        result["全文"] = full_text
        result["全文来源"] = "pmc_fulltext"

    # 提取表格数据
    tables = _extract_tables_from_xml(root)
    if tables:
        result["tables"] = tables

    return result


def _extract_tables_from_xml(root):
    """从 XML 全文提取表格数据（行列表）。支持 table-wrap (JATS) 和 table 元素。"""
    tables = []

    def _process_table_element(table_elem, caption=""):
        """处理单个表格元素（可能是 <table> 或 <table-wrap> 内的 <table>）"""
        rows = []
        for tr in table_elem.findall(".//tr"):
            cells = []
            for cell in tr.findall("th"):
                txt = " ".join(t.strip() for t in cell.itertext() if t.strip())
                cells.append(f"[H]{txt}")
            for cell in tr.findall("td"):
                txt = " ".join(t.strip() for t in cell.itertext() if t.strip())
                cells.append(txt)
            if cells:
                rows.append(cells)
        return rows

    # JATS 格式: <table-wrap> 包裹 <table>，caption 在 table-wrap 层
    for tw in root.findall(".//table-wrap"):
        tw_id = tw.get("id", "")
        caption = ""
        caption_el = tw.find("caption")
        if caption_el is not None:
            caption = " ".join(
                t.strip() for t in caption_el.itertext() if t.strip()
            )
        # 在 table-wrap 内找 <table>
        table_el = tw.find("table")
        if table_el is not None:
            rows = _process_table_element(table_el, caption)
        else:
            # table-wrap 内可能直接有 <tr> (某些格式)
            rows = _process_table_element(tw, caption)
        if rows:
            tables.append({
                "id": tw_id,
                "caption": caption,
                "rows": rows,
                "raw_text": "\n".join(
                    "\t".join(c for c in row) for row in rows
                ),
            })

    # 直接的 <table> 元素（不在 table-wrap 内的，如 HTML 格式）
    processed_tables = set()
    for tw in root.findall(".//table-wrap"):
        table_el = tw.find("table")
        if table_el is not None:
            processed_tables.add(table_el)

    for table_elem in root.findall(".//table"):
        if table_elem in processed_tables:
            continue  # 已通过 table-wrap 处理
        caption_el = table_elem.find(".//caption")
        caption = ""
        if caption_el is not None:
            caption = " ".join(
                t.strip() for t in caption_el.itertext() if t.strip()
            )
        rows = _process_table_element(table_elem, caption)
        if rows:
            tables.append({
                "id": table_elem.get("id", ""),
                "caption": caption,
                "rows": rows,
                "raw_text": "\n".join(
                    "\t".join(c for c in row) for row in rows
                ),
            })

    return tables


# ── 句子拆分 ──

def split_sentences(text):
    """按句子拆分，保护 c.NNN 格式不被截断。"""
    protected = re.sub(r'\bc\.(?=\d)', 'c<DOT>', text)
    raw_sentences = re.split(r'(?<=[.!?])\s+', protected)
    return [s.replace('c<DOT>', 'c.').strip() for s in raw_sentences if len(s.strip()) >= 15]


# ── 变异信息提取（核心逻辑） ──

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


def extract_pathogenicity(sentences, full_text_lower):
    """提取致病性评级。"""
    vs_lower = " ".join(sentences).lower()
    path_terms = []
    pathogenicity_map = [
        ("致病 (pathogenic)", r'\bpathogenic\b'),
        ("可能致病 (likely pathogenic)", r'likely\s+pathogenic'),
        ("良性 (benign)", r'\bbenign\b'),
        ("可能良性 (likely benign)", r'likely\s+benign'),
        ("意义不明 (VUS)", r'variant\s+of\s+uncertain\s+significance|\bVUS\b'),
        ("有害 (damaging)", r'\bdamaging\b'),
        ("有害 (deleterious)", r'\bdeleterious\b'),
        ("致病 (disease-causing)", r'disease.?causing'),
    ]
    for term, pattern in pathogenicity_map:
        if re.search(pattern, vs_lower) and term not in path_terms:
            path_terms.append(term)
    if not path_terms:
        for term, pattern in pathogenicity_map:
            m = re.search(pattern, full_text_lower)
            if m:
                pos = m.start()
                context = full_text_lower[max(0, pos - 300):pos + 300]
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
                scored_matches.append((z_name, score, s, m.start()))
                break  # 每句只取第一个匹配的合子关键词

    # ---- 第2层：表格级别匹配（高权重） ----
    table_zygosity = None
    if tables and keywords:
        for table_info in tables:
            for ri, row in enumerate(table_info["rows"]):
                row_text = "\t".join(str(cell or "") for cell in row)
                row_lower = row_text.lower()
                # 该行是否包含目标变异
                has_target = False
                for kw in keywords.get("all", []):
                    if kw.lower() in row_lower:
                        has_target = True
                        break
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
                        next_first_cell = str(next_row[0] or "").strip() if next_row else ""
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
        return re.sub(r'[ -\u200F   　]', ' ', t)

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
    full_text = re.sub(r'[ -\u200F  　]', ' ', full_text)
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
                        row[i] = re.sub(r'[ -\u200F  　]', ' ', str(row[i]))
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
                row_text = re.sub(r"[ -‏  　]", " ", row_text)  # v8: 去除 Unicode 空白
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


# ── 总结段落生成 ──

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


# ── 主流程 ──

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
    import csv

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

def main():
    parser = argparse.ArgumentParser(
        description="从PubMed/PMC文献中仅针对目标变异提取信息"
    )
    parser.add_argument("--pmids", nargs="+", required=True, help="PubMed ID列表")
    parser.add_argument("--gene", default=None, help="目标基因名称")
    parser.add_argument("--variant", default=None,
                        help='目标变异，如 "c.1166G>A (p.Arg389His)"')
    parser.add_argument("--transcript", default=None,
                        help="转录本ID，如 NM_022436.3")
    parser.add_argument("--output", default="pubmed_variant_results.json",
                        help="输出JSON文件路径")
    parser.add_argument("--excel-dir", default=r"D:\claude_code\project1\文献提取结果",
                        help="Excel汇总表输出目录，默认 D:\\claude_code\\project1\\文献提取结果")
    parser.add_argument("--pdf-results", default=None,
                        help="PDF搜索结果的JSON文件路径（用于交叉引用，补充在线API遗漏的变异检出）")
    args = parser.parse_args()

    # 解析目标变异
    target_cdna = ""
    target_protein = ""
    if args.variant:
        cdna_m = re.search(r'c\.([\d\w_*><+=\-\+]+)', args.variant)
        if cdna_m:
            target_cdna = f"c.{cdna_m.group(1)}"
        prot_m = re.search(r'p\.\s*([\w?*]+)', args.variant)
        if prot_m:
            target_protein = f"p.{prot_m.group(1)}"

    if not args.gene:
        print("提示: 建议提供 --gene 参数以精确定位目标基因")

    # 去重
    seen = set()
    unique_pmids = []
    for pmid in args.pmids:
        if pmid not in seen:
            seen.add(pmid)
            unique_pmids.append(pmid)
    if len(unique_pmids) < len(args.pmids):
        print(f"去重: {len(args.pmids)} -> {len(unique_pmids)} 篇")

    # 构建变异关键词
    keywords = build_variant_keywords(target_cdna, target_protein, args.transcript or "")
    print(f"\n目标变异: {target_cdna} {target_protein}")
    print(f"搜索关键词 ({len(keywords['all'])} 个):")
    print(f"  精确: {keywords['exact']}")
    print(f"  模糊: {keywords['fuzzy']}")
    print(f"  蛋白: {keywords['protein']}")

    # 转录本信息查询
    if args.transcript:
        print(f"\n查询转录本信息: {args.transcript}")
        tx_info = fetch_ncbi_gene_info(args.transcript)
        if tx_info:
            print(f"  转录本信息: {tx_info.get('title', '')}")
        time.sleep(5)

    # 站点访问计数状态
    counts = load_site_counts()
    print(f"\n今日 ({date.today()}) 站点访问计数:")
    for site_key, count in counts.items():
        status = "[已达上限]" if count >= DAILY_SITE_LIMIT else "[OK]"
        print(f"   {site_key}: {count}/{DAILY_SITE_LIMIT} {status}")

    # v9: 加载PDF搜索结果用于交叉引用
    pdf_data = {}
    if args.pdf_results and os.path.exists(args.pdf_results):
        with open(args.pdf_results, 'r', encoding='utf-8') as f:
            pdf_list = json.load(f)
            for pr in pdf_list:
                pid = str(pr.get("PMID", "")).strip()
                # 回退: 从文件名提取 PMID (如 PubMed22326530.pdf → 22326530)
                if not pid:
                    fname = pr.get("文件", "")
                    m = re.search(r'(\d{7,9})', fname)
                    if m:
                        pid = m.group(1)
                if pid:
                    # 如果同一个 PMID 有多个条目（文件名不同），保留提及变异的
                    if pid in pdf_data and not pr.get("变异提及"):
                        continue
                    pdf_data[pid] = pr
        print(f"\nPDF交叉引用: 已加载 {len(pdf_data)} 条PDF搜索结果")
    elif args.pdf_results:
        print(f"\n[警告] PDF结果文件不存在: {args.pdf_results}")

    results = []
    request_count = 0
    stopped_early = False

    for i, pmid in enumerate(unique_pmids):
        print(f"\n{'='*80}")
        print(f"PMID: {pmid}  ({i+1}/{len(unique_pmids)})")
        print(f"{'='*80}")

        # Step 1: 获取NCBI摘要
        xml_text = fetch_ncbi_abstract(pmid)
        if xml_text is None and load_site_counts().get("eutils.ncbi.nlm.nih.gov", 0) >= DAILY_SITE_LIMIT:
            stopped_early = True
            break
        request_count += 1
        periodic_long_pause(request_count)

        result = None
        if xml_text:
            result = parse_ncbi_xml(xml_text)

        # NCBI 摘要获取失败时，回退到 Europe PMC
        if not result:
            print(f"  NCBI 摘要获取失败，尝试 Europe PMC...")
            epmc_info = fetch_europe_pmc_text(pmid)
            if epmc_info:
                result = _empty_result()
                result["PMID"] = pmid
                result["标题"] = epmc_info.get("title", "")
                result["摘要"] = epmc_info.get("abstract", "")
                result["全文"] = epmc_info.get("abstract", "")
                result["全文来源"] = "epmc_abstract"
                result["作者"] = epmc_info.get("authors", [])
                result["期刊"] = epmc_info.get("journal", "")
                result["发表年份"] = epmc_info.get("year", "")
                print(f"  已从 Europe PMC 获取摘要")
            else:
                print(f"  无法获取任何来源，跳过")
                continue

        # Step 2: 随机间隔
        random_sleep()
        request_count += 1
        periodic_long_pause(request_count)

        # Step 3: 尝试获取全文 (Europe PMC -> PMC)
        fulltext_xml = fetch_europe_pmc_fulltext(pmid)
        request_count += 1
        if not fulltext_xml:
            random_sleep()
            request_count += 1
            fulltext_xml = fetch_pmc_fulltext(pmid)

        if fulltext_xml:
            result = parse_pmc_europe_xml(fulltext_xml, result)
            print(f"  全文获取成功 (来源: {result['全文来源']})")
            if result.get("tables"):
                print(f"  提取到 {len(result['tables'])} 个表格")
            # 如果 PMC XML 受限，尝试从 PMC HTML 页面获取全文（网页版通常免费）
            if result.get("全文来源") == "pmc_restricted":
                print(f"  PMC XML 受出版商限制，尝试 PMC HTML 网页版...")
                html_data = fetch_pmc_html(pmid)
                if html_data:
                    result = parse_pmc_html(html_data, result)
                    print(f"  已获取 PMC HTML (来源: {result['全文来源']})")
                    if result.get("tables"):
                        print(f"  从 HTML 提取到 {len(result['tables'])} 个表格")
                else:
                    # 回退到 Europe PMC 摘要
                    print(f"  PMC HTML 获取失败，尝试 Europe PMC 摘要...")
                    epmc_text = fetch_europe_pmc_text(pmid)
                    if epmc_text and epmc_text.get("abstract"):
                        result["全文"] = epmc_text.get("abstract", "")
                        result["全文来源"] = "epmc_abstract"
                        print(f"  已获取 Europe PMC 摘要文本 ({len(result['全文'])} 字符)")
        else:
            print(f"  仅获取摘要")
            # NCBI 摘要成功但无全文，尝试补充 Europe PMC 摘要
            epmc_text = fetch_europe_pmc_text(pmid)
            if epmc_text and epmc_text.get("abstract"):
                result["全文"] = epmc_text.get("abstract", "")
                result["全文来源"] = "epmc_abstract"
                print(f"  已补充 Europe PMC 摘要文本 ({len(result['全文'])} 字符)")

        # Step 4: 仅针对目标变异提取信息
        pdf_match = pdf_data.get(str(pmid))
        result = extract_info_for_variant(result, args.gene, target_cdna, target_protein, keywords, pdf_match)

        # Step 5: 生成一句话概括和总结段落
        result["一句话概括"] = generate_one_sentence_summary(result)
        if result["变异提及"]:
            result["总结段落"] = generate_summary_paragraph(result)
        else:
            result["总结段落"] = result["一句话概括"]  # 未提及时使用文献简介作为总结

        # 打印结果
        print(f"\n  标题: {result['标题']}")
        print(f"  一句话概括: {result['一句话概括']}")
        print(f"  作者: {', '.join(result['作者'][:3])}{'...' if len(result['作者']) > 3 else ''}")
        print(f"  期刊: {result['期刊']} ({result['发表年份']})")
        print(f"  全文来源: {result['全文来源']}")

        if result["变异提及"]:
            print(f"\n  目标变异: {result['基因']} {result['cDNA变异']} {result['蛋白变异']}")
            print(f"    变异类型: {result['变异类型']}")
            print(f"    致病性:   {result['致病性']}")
            print(f"    遗传方式: {result['合子状态']}")
            print(f"    临床表型: {result['临床表型']}")
            print(f"    遗传模式: {result['遗传模式']}")
            print(f"    匹配关键词: {', '.join(result.get('匹配关键词', []))}")

            if result.get("共存变异"):
                co_var_display = []
                for cv in result['共存变异']:
                    if isinstance(cv, dict):
                        parts = [cv.get("cdna") or "", cv.get("蛋白变异") or ""]
                        co_var_display.append(" ".join(p for p in parts if p))
                    else:
                        co_var_display.append(str(cv))
                if co_var_display:
                    print(f"    共存变异: {', '.join(co_var_display)}")
            if result.get("反式确认"):
                print(f"    反式位置确认: 是")
            if result.get("顺式确认"):
                print(f"    顺式位置确认: 是")
            if result.get("相位状态"):
                print(f"    相位状态: {result['相位状态']} ({result.get('相位置信度', '')})")
                print(f"    相位详情: {result.get('相位详情', '')}")
            if result.get("亲本检测"):
                print(f"    亲本检测: 已进行")
                if result.get("母源变异"):
                    print(f"      母源变异: {result['母源变异']}")
                if result.get("父源变异"):
                    print(f"      父源变异: {result['父源变异']}")
            if result.get("相位证据"):
                print(f"    相位证据句数: {len(result['相位证据'])}")
            if result.get("患者数量"):
                print(f"    患者数量: {result['患者数量']}")

            if result['患者详情']:
                print(f"\n  患者详情:")
                for j, pd in enumerate(result['患者详情'], 1):
                    print(f"    患者{j}: {json.dumps(pd, ensure_ascii=False)}")

            if result.get("变异特征"):
                print(f"\n  变异特征:")
                for k, v in result["变异特征"].items():
                    print(f"    {k}: {v}")

            print(f"\n  功能验证: {result['功能验证']}")
            if result['功能验证详情']:
                for fd in result['功能验证详情']:
                    print(f"    - {fd[:150]}...")

            if result.get("相关句子"):
                print(f"\n  原文相关句:")
                for rs in result["相关句子"][:3]:
                    print(f"    > {rs[:200]}...")

            if result.get("总结段落"):
                print(f"\n  === 标准化总结段落 ===")
                print(f"  {result['总结段落']}")
                print(f"  === 结束 ===")
        else:
            print(f"\n  文献中未提及目标变异 {args.gene} {target_cdna}")
            if result.get("一句话概括"):
                print(f"\n  === 文献简介 ===")
                print(f"  {result['一句话概括']}")
                print(f"  === 结束 ===")

        results.append(result)

        if i < len(unique_pmids) - 1:
            random_sleep()
            request_count += 1
            periodic_long_pause(request_count)

    # 输出JSON
    output_path = os.path.join(os.getcwd(), args.output)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\n\nJSON output: {output_path}")
    print(f"共处理 {len(results)}/{len(args.pmids)} 篇文献")
    print(f"总请求数: {request_count}")

    # 生成 Excel 汇总表
    excel_dir = args.excel_dir
    safe_gene = re.sub(r'[\\/:*?"<>|]', '_', args.gene) if args.gene else "gene"
    safe_variant = re.sub(r'[\\/:*?"<>|]', '_', target_cdna) if target_cdna else "variant"
    excel_filename = f"{safe_gene}_{safe_variant}_文献汇总.csv"
    excel_path = os.path.join(excel_dir, excel_filename)
    try:
        _generate_excel_csv(results, excel_path)
        print(f"\nExcel 汇总表: {excel_path}")
    except Exception as e:
        print(f"\n[警告] Excel 汇总表生成失败: {e}")

    # 打印最终站点计数
    counts = load_site_counts()
    if counts:
        print(f"\n今日站点访问计数:")
        for site_key, count in counts.items():
            status = "[已达上限]" if count >= DAILY_SITE_LIMIT else "[OK]"
            print(f"   {site_key}: {count}/{DAILY_SITE_LIMIT} {status}")

    if stopped_early:
        print("\n因站点访问限制提前停止")

    # 汇总表格
    print(f"\n{'='*80}")
    print(f"汇总:")
    print(f"{'PMID':<10} {'提及':<6} {'变异类型':<16} {'致病性':<14} {'遗传方式':<14} {'患者数':<6}")
    print("-" * 80)
    for r in results:
        pmid = r.get("PMID", "?")[:10]
        mentioned = "是" if r.get("变异提及") else "否"
        vtype = (r.get("变异类型", "-") or "-")[:16]
        patho = (r.get("致病性", "-") or "-")[:14]
        zyg = (r.get("合子状态", "-") or "-")[:14]
        pcnt = r.get("患者数量", 0) or 0
        print(f"{pmid:<10} {mentioned:<6} {vtype:<16} {patho:<14} {zyg:<14} {pcnt:<6}")


if __name__ == "__main__":
    main()
