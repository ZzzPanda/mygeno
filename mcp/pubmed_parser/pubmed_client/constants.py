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