---
name: pubmed-extractor
description: Extract gene variant-specific information from PubMed articles with transcript version awareness. Uses NCBI E-utilities and Europe PMC APIs with stable, pre-written Python script. Mandatory extraction of cis/trans configuration and allelic phase. Outputs JSON + Excel summary to D:\claude_code\project1\文献提取结果.
---

# PubMed 变异文献提取器

从 PubMed/PMC 文献中**仅针对目标变异**提取相关信息，生成标准化中文总结段落和 Excel 汇总表。

## 核心原则

1. **使用预置 Python 脚本，无需临时编写代码。所有逻辑（关键词扩展、变异匹配、相位提取、文献简介生成、Excel 输出）均预写在脚本中，确保每次稳定运行并可重复结果。**
2. **必须提取文章中该位点的正反式（cis/trans）配置和等位基因相位（phase）信息。**
3. **未提及目标变异的文献，自动生成 >200 字的中文文献简介**（基于预置的 180+ 主题词映射表 + 标题/摘要/正文提取）。

脚本路径: `.claude/skills/pubmed-extractor/scripts/pubmed_extractor.py`
- 纯 Python 标准库，零外部依赖
- 所有逻辑预置在脚本中，不依赖 AI 临时生成代码
- Windows UTF-8 编码兼容
- 输入 PMID 自动去重
- **每次运行相同输入产生相同输出（确定性）**

## 使用方法

### 一行命令执行

```bash
python .claude/skills/pubmed-extractor/scripts/pubmed_extractor.py \
    --pmids PMID1 PMID2 ... \
    --gene GENE \
    --variant "c.NNN X>Y (p.AA)" \
    [--transcript NM_xxx.x] \
    [--output results.json] \
    [--excel-dir D:\path\to\output]
```

### 示例

```bash
# 基础用法（错义突变）
python .claude/skills/pubmed-extractor/scripts/pubmed_extractor.py \
    --pmids 36972931 --gene MFSD8 --variant "c.1444C>T (p.Arg482Ter)"

# 带转录本（推荐，可识别版本编号差异）
python .claude/skills/pubmed-extractor/scripts/pubmed_extractor.py \
    --pmids 20521169 28771437 --gene ABCG5 \
    --variant "c.1166G>A (p.Arg389His)" --transcript NM_022436.3

# 多 PMID 批量检索
python .claude/skills/pubmed-extractor/scripts/pubmed_extractor.py \
    --pmids 12842373 17048214 39019822 29801666 --gene PKD1 \
    --variant "c.1522T>C (p.Cys508Arg)" --transcript NM_001009944.3

# 指定 JSON 和 Excel 输出路径
python .claude/skills/pubmed-extractor/scripts/pubmed_extractor.py \
    --pmids 19177532 --gene MFSD8 --variant "c.1444C>T (p.Arg482Ter)" \
    --output custom_results.json \
    --excel-dir D:\claude_code\project1\文献提取结果
```

### 参数说明

| 参数 | 必需 | 说明 |
|------|------|------|
| `--pmids` | 是 | PubMed ID 列表，空格分隔，自动去重 |
| `--gene` | 是 | 目标基因名称 |
| `--variant` | 是 | 目标变异，格式如 `"c.1522T>C (p.Cys508Arg)"` |
| `--transcript` | 否 | 转录本 ID（推荐提供，可识别版本差异） |
| `--output` | 否 | 输出 JSON 路径，默认 `pubmed_variant_results.json` |
| `--excel-dir` | 否 | Excel 汇总表输出目录，默认 `D:\claude_code\project1\文献提取结果` |

## 工作流程（全自动）

| 步骤 | 操作 | 说明 |
|------|------|------|
| 1 | 关键词扩展 | 脚本自动将输入的变异信息扩展为 9+ 种搜索形式 |
| 2 | NCBI E-utilities | 获取 PubMed XML 摘要（标题、作者、期刊、MeSH 术语） |
| 3 | Europe PMC API | 获取开放获取全文 XML（含表格） |
| 4 | PMC API | 获取 PMC 全文 XML（含表格） |
| 5 | PMC HTML 网页 | 当 XML 受限时，自动抓取 PMC 网页版全文和表格 |
| 6 | 变异匹配 | 在正文和表格中搜索目标变异的所有关键词变体（双重搜索） |
| 7 | 信息提取 | 提取结构化信息（患者详情、致病性、遗传方式、功能验证等） |
| 8 | 正反式/相位提取 | 从文献中提取 cis/trans 配置、等位基因相位、亲本源信息 |
| 9 | 文献背景分析 | 基于 180+ 主题词映射表自动判断研究类型和领域 |
| 10 | 总结生成 | 提及变异→标准化中文总结段落；未提及→>200 字文献简介 |
| 11 | Excel 汇总表 | 生成 CSV 汇总表（UTF-8 BOM 编码，Excel/WPS 直接打开，中文不乱码） |
| 12 | JSON 详情 | 保存结构化 JSON（含所有字段和原始证据） |

## 脚本内置特性

- **速率限制**：随机间隔 10-40 秒 + 每站点每日 500 次上限（自动停止，计数持久化到 `scripts/.site_counts_{日期}.json`）
- **转录本版本感知**：通过 `--transcript` 参数自动查询转录本信息，识别不同版本间的 cDNA 编号差异
- **变异关键词自动扩展**：输入 `c.1522T>C (p.Cys508Arg)` 自动生成 9 个搜索关键词
- **全文表格解析**：从 JATS XML 和 PMC HTML 中提取表格行与患者数据
- **双重变异搜索**：正文句子匹配 + 表格行匹配，确保表格中的变异不被遗漏
- **正反式/相位提取**：自动检测 cis/trans 配置、亲本检测、等位基因相位分类
- **四级回退链**：NCBI 摘要 → Europe PMC 全文 → PMC 全文 → PMC HTML 网页 → Europe PMC 摘要
- **表型关键词映射**：内置 60+ 临床表型关键词的中英文映射
- **文献主题词映射**：内置 180+ 中英文主题词映射，覆盖骨骼纤毛、视网膜、神经、代谢、心脏、血液、肿瘤、产前、方法学等领域
- **文献背景自动分析**：基于标题/摘要/MeSH + 主题词映射，自动判断研究类型（病例报告/队列研究/携带者筛查/产前诊断/基因检测/功能学等）
- **Excel 汇总表**：UTF-8 BOM CSV 格式，Excel/WPS 直接打开，中文无乱码，零外部依赖

---

## Excel 汇总表输出

### 输出位置

默认路径：`D:\claude_code\project1\文献提取结果\{基因}_{变异}_文献汇总.csv`

可通过 `--excel-dir` 参数自定义输出目录。

### 表头列

| 列名 | 说明 |
|------|------|
| PMID | PubMed ID |
| 标题 | 文献标题（截取前 200 字符） |
| **是否提及此位点** | 是 / 否 |
| **患者数** | 携带目标变异的患者数量 |
| **致病性** | 致病性评级（致病/可能致病/良性/有害等） |
| **关联合子状态** | 纯合/杂合/复合杂合 |
| **反式(trans)位点** | 相位状态 + 共存变异 + 亲本源信息（如有） |
| **患者临床表型** | 携带目标变异患者的具体临床表型 |
| **文献背景(是什么研究)** | 自动分析的研究类型和领域（如"产前诊断研究，涉及骨骼发育不良、短肋、纤毛病，发表于2021年《Frontiers in genetics》"） |
| **总结** | 提及变异→标准化中文总结段落；未提及→>200 字文献简介 |

### 格式

- 输出格式为 **UTF-8 BOM CSV**（`.csv`），零外部依赖，Excel / WPS 直接打开，中文无乱码
- UTF-8 编码，中文无乱码
- 每次运行覆盖同名文件

---

## 文献简介（未提及变异时）

### 何时触发

当文献中**未搜索到目标变异**的任何关键词变体时，自动触发文献简介生成。

### 生成策略（预置逻辑，确保可重复）

1. **主题词匹配**：从标题中匹配 `TOPIC_KEYWORD_MAP`（180+ 条目）中的英文关键词，转换为中文主题标签
2. **全文补充**：标题匹配不足时，从摘要和正文前 5000 字符补充主题词（上限 6 个）
3. **MeSH 术语补充**：从 PubMed MeSH 术语中提取额外主题词（上限 8 个）
4. **研究类型判断**：基于预置规则自动判断（病例报告/队列研究/Meta分析/携带者筛查/产前诊断/基因检测/功能学/人群遗传学/突变筛查/综述）
5. **正文提取**：提取摘要或正文前若干句作为核心内容描述
6. **长度保证**：如不足 200 字符，自动从全文/摘要追加补充内容

### 简介示例

> 该文献为携带者筛查、孕前、基因包检测、外显子组测序、致病性、ACMG相关研究。发表于2019年《PLoS genetics》。核心内容：Pathogenic variants in autosomal recessive and X-linked recessive mendelian disorders were identified from 14,125 exomes. We optimized clinical exome design and parallel gene-testing for recessive genetic conditions in preconception carrier screening...

---

## 正反式与相位提取（核心功能）

### 为什么需要提取相位信息

对于复合杂合变异，必须确定两个变异是否位于同一等位基因（顺式/cis）还是不同等位基因（反式/trans），这直接影响致病性判读：

- **反式 (trans)**：两个变异在不同等位基因上 → 符合常染色体隐性遗传模式 → 支持致病性 (ACMG PM3)
- **顺式 (cis)**：两个变异在同一等位基因上 → 另一等位基因可能正常 → 可能为携带者状态 → 降低致病性权重 (ACMG BP2)
- **相位未确定**：无法判断 cis/trans → 致病性判读受限

### 相位分类体系

脚本自动将提取到的相位信息归入以下类别：

| 相位状态 | 英文标识 | 含义 | 置信度 |
|----------|----------|------|--------|
| **已确认反式** | `confirmed_in_trans` | 经亲本检测或分子手段确认位于不同等位基因 | 确认 (confirmed) |
| **已确认顺式** | `confirmed_in_cis` | 经亲本检测确认位于同一等位基因（复杂等位基因） | 确认 (confirmed) |
| **推定反式** | `presumed_in_trans` | 复合杂合，隐性遗传默认推定反式，但未经亲本验证 | 推定 (presumed) |
| **推定顺式** | `presumed_in_cis` | 两个变异在同一读长/同一克隆中检出 | 推定 (presumed) |
| **相位未确定** | `phase_not_determined` | 文献明确指出相位未知或无法确定 | 未知 (unknown) |
| **未评估** | `not_assessed` | 文献未涉及相位信息 | 无 |
| **不适用** | `not_applicable` | 纯合变异或仅涉及单一变异，无相位问题 | 无 |

### 提取策略（预置逻辑）

脚本在正文和表格中自动搜索以下证据类型：

**1. 亲本检测证据（最高置信度）**
```
"inherited from the mother" / "inherited from the father"
"maternally inherited" / "paternally inherited"
"mother was a carrier" / "father was a carrier"
"parental testing confirmed" / "parental segregation analysis"
"segregation analysis" / "co-segregation"
```

**2. 反式 (trans) 配置证据**
```
"in trans" / "on opposite alleles" / "trans configuration"
"on different alleles" / "biallelic in trans"
"compound heterozygous ... confirmed in trans"
```

**3. 顺式 (cis) 配置证据**
```
"in cis" / "on the same allele" / "cis configuration"
"complex allele" / "same parental allele" / "double mutant allele"
```

**4. 相位未知/不确定证据**
```
"phase not determined" / "phase unknown"
"cis or trans not determined" / "allelic phase could not be determined"
```

**5. 亲本源具体信息**
- 母源等位基因变异 vs. 父源等位基因变异
- 各亲本携带状态（携带者/患者/正常）
- 新发突变 (de novo) 排除亲本遗传

### 相位信息在总结段落中的体现

- **已确认反式**: "通过亲本检测确认，该变异与 c.XXX 处于反式位置（分别来自父母双方）"
- **推定反式**: "该变异与 c.XXX 组成复合杂合，推定为反式位置（未经亲本验证）"
- **已确认顺式**: "该变异与 c.XXX 处于顺式位置，构成复杂等位基因"
- **相位未确定**: "该变异与 c.XXX 的等位基因相位尚未确定"
- **亲本溯源**: "该变异遗传自母亲，c.XXX 遗传自父亲"

---

## Excel 反式(trans)位点列的内容规则

| 场景 | Excel 列显示内容 |
|------|-----------------|
| 确认反式 + 亲本验证 | `确认反式 (confirmed trans)；共存变异: c.152T>G (p.L51R)；母源: c.152T>G；父源: c.988C>T` |
| 确认反式，无亲本详情 | `确认反式 (confirmed trans)；共存变异: c.XXX` |
| 推定反式 | `推定反式 (presumed trans)；共存变异: c.XXX` |
| 纯合变异 | `不适用（纯合变异）` |
| 单变异/未评估 | `未确认` |
| 未提及变异 | `不适用` |

---

## 变异关键词扩展策略（预置逻辑）

输入 `c.1522T>C (p.Cys508Arg)` 时，脚本自动生成：

| 类型 | 生成的关键词 |
|------|-------------|
| 精确 cDNA | `c.1522T>C` |
| 精简 cDNA | `1522T>C`, `c1522T>C`, `1522` |
| 蛋白三字母 | `p.Cys508Arg`, `Cys508Arg` |
| 蛋白单字母 | `p.C508R`, `C508R` |
| 带空格变体 | `C 508 R` |
| 终止密码子变体 | `R508*`, `C508Ter`, `R508X` (针对无义突变额外生成) |
| 描述性短语 | "deletion of serine at position 189" 等（自动生成正则模式） |
| 历史命名 | 如 "G6PD Tsukui" → c.565_567del（基于内置映射表） |

## 限制场景与 fallback 行为

| 场景 | 来源标识 | 说明 |
|------|----------|------|
| 开放获取文章 | `pmc_fulltext` | 完整正文 + 表格 |
| 出版商限制 XML | `pmc_restricted` → `epmc_abstract` | 自动回退到 PMC HTML 网页版 |
| NCBI 摘要失败 | `epmc_abstract` | 自动回退到 Europe PMC 摘要 |
| 仅摘要可用 | `abstract` | 仅搜索摘要文本；文献简介质量受限 |

## 输出文件

### 1. JSON 详情文件（`--output`）

| 字段 | 类型 | 说明 |
|------|------|------|
| `PMID` | string | PubMed ID |
| `标题` | string | 文献标题 |
| `期刊` | string | 期刊名称 |
| `发表年份` | string | 发表年份 |
| `作者` | list | 作者列表 |
| `全文来源` | string | 文本来源 (`pmc_fulltext` / `pmc_html` / `epmc_abstract` / `abstract`) |
| `一句话概括` | string | 提及变异→变异相关概要；未提及→>200 字文献简介 |
| `变异提及` | bool | 是否提及目标变异 |
| `基因` | string | 目标基因名称 |
| `cDNA变异` | string | cDNA 改变（如 c.1522T>C） |
| `蛋白变异` | string | 蛋白质改变（如 p.Cys508Arg） |
| `匹配关键词` | list | 实际匹配到的关键词列表 |
| `变异类型` | string | 变异类型（错义/无义/移码/剪接等） |
| `致病性` | string | 致病性评级 |
| `合子状态` | string | 合子状态（纯合/杂合/复合杂合） |
| `遗传模式` | string | 遗传模式（常隐/常显/X连锁等） |
| `临床表型` | string | 临床表型 |
| `患者详情` | list | 携带目标变异的患者个体信息 |
| `患者数量` | int | 携带目标变异的患者数量 |
| `共存变异` | list | 共存变异（同一患者记录中的另一等位基因变异） |
| `相位状态` | string | 等位基因相位分类 |
| `相位置信度` | string | confirmed / presumed / unknown / not_applicable |
| `顺式确认` | bool | 是否有顺式（cis）位置证据 |
| `反式确认` | bool | 是否有反式（trans）位置证据 |
| `亲本检测` | bool | 是否进行了亲本验证实验 |
| `母源变异` | string/null | 母源等位基因的变异 |
| `父源变异` | string/null | 父源等位基因的变异 |
| `相位证据` | list | 支持相位判断的原文句子及证据类型 |
| `变异特征` | dict | 变异特征（CpG位点、NMD、人群频率、新型变异等） |
| `功能验证` | string | 功能验证方法 |
| `功能验证详情` | list | 功能验证详情（原文句子） |
| `相关句子` | list | 文献中直接提及该变异的所有句子 |
| `表格数量` | int | 提取到的表格数量 |
| `表格摘要` | list | 表格简要信息 |
| `总结段落` | string | 提及变异→标准化中文描述段落；未提及→>200 字文献简介 |

### 2. Excel 汇总表（`--excel-dir`）

默认保存至 `D:\claude_code\project1\文献提取结果\{基因}_{变异}_文献汇总.csv`

表头：PMID | 标题 | 是否提及此位点 | 患者数 | 致病性 | 关联合子状态 | 反式(trans)位点 | 患者临床表型 | 文献背景(是什么研究) | 总结

---

## 总结段落生成原则

**核心原则：仅输出文献中确实存在的信息，不杜撰任何内容。**

- 如果文献中没有反式位置证据，不生成相关句子
- 如果没有功能验证数据，不添加"尚未经过功能学验证"等表述
- 如果没有共分离数据，不生成相关句子
- 共存变异仅在目标变异同一行/同一段落中出现时才被提取
- **相位信息必须基于文献原文证据，不可臆测推定**
- 未提及变异时，总结段落为自动生成的 >200 字文献简介

### 隐性遗传（复合杂合/纯合）模板

```
共X例患者为该c.XXX变异与一个致病性或可能致病性c.XXX变异的复合杂合。
通过亲本检测确认，该变异与c.XXX处于反式位置（仅当有trans确认证据时生成）。
该变异与c.XXX的等位基因相位尚未经亲本验证（仅当相位未确定时生成）。
该变异遗传自母亲，c.XXX遗传自父亲（仅当有亲本溯源时生成）。
该先证者通过父母检测确认处于反式位置（与旧版兼容）。
该患者表型为xx、xx、xx等。
该变异已经过[实验名称]显示可能影响蛋白功能（仅当有功能验证时添加）。
参考文献：xxx
```

### 显性/半合子遗传模板

```
该变异已在X例先证者中报道。
该患者表型为xx、xx、xx。
该变异已经过[实验名称]显示可能影响蛋白功能（仅当有功能验证时添加）。
参考文献：xxx
```

---

## 转录本版本差异处理

### 问题描述

同一基因的不同转录本版本（如 NM_022436.2 vs NM_022436.3）由于 CDS 起始位置不同，会导致相同的氨基酸变异对应不同的 cDNA 编号。

| 基因 | 转录本 | cDNA 变异 | 蛋白变异 | 说明 |
|------|--------|----------|---------|------|
| ABCG5 | NM_022436.2 | c.1306G>A | p.Arg389His | 旧版转录本 |
| ABCG5 | NM_022436.3 | c.1166G>A | p.Arg389His | 当前 RefSeq |

### 脚本应对策略

1. **蛋白位置优先**：蛋白变异编号（如 p.Arg389His）在所有转录本版本中通常一致，脚本优先匹配蛋白关键词
2. **cDNA 数字模糊匹配**：当精确 cDNA 匹配不足时，搜索 cDNA 数字位置（如 `1522`）
3. **转录本信息查询**：通过 `--transcript` 参数可查询转录本元数据
4. **匹配关键词记录**：结果 JSON 中 `匹配关键词` 字段记录实际匹配到的关键词，便于追溯

### 建议

- **始终提供 `--transcript` 参数**，以确保使用正确的参考转录本
- 检查输出 JSON 中的 `匹配关键词` 字段，确认匹配到的具体关键词
- 注意无义突变在文献中可能使用多种命名（`Ter`, `*`, `X`, `Stop`），脚本已覆盖

## 注意事项

- 非开放获取文献仅能获取摘要（脚本自动标识来源，文献简介质量可能受限）
- 出版商限制 XML 时，脚本自动回退到 PMC HTML 网页版
- 如文献未提及目标变异，自动生成 >200 字文献简介（基于预置主题词映射 + 标题/摘要提取）
- 患者详情依赖文献是否包含个体患者描述
- **速率限制**计数保存在 `scripts/.site_counts_{日期}.json` 文件中，跨会话持久化
- **每次运行结果确定且可重复**（相同输入产生相同输出，包括文献简介内容）
- 脚本执行期间会打印详细进度和中间结果，便于监控
- **相位信息提取依赖文献原文**，如文献未做亲本验证则标记为"推定"而非"确认"
- Excel 汇总表采用 UTF-8 BOM CSV 格式（`.csv`），零外部依赖，Excel/WPS 直接打开，中文无乱码

## 维护与更新

脚本路径: `.claude/skills/pubmed-extractor/scripts/pubmed_extractor.py`

所有功能变更直接修改该脚本，不依赖 prompt engineering 或 AI 临时代码。修改后运行相同命令即可获得更新后的结果。

### 版本历史

| 版本 | 日期 | 变更 |
|------|------|------|
| v6 | 2026-04 | 转录本版本感知、PMC HTML fallback、描述性短语匹配 |
| v7 | 2026-05 | 正反式(cis/trans)/相位(phase)全面提取、亲本溯源、相位分类体系 |
| **v8** | 2026-05 | **文献背景自动分析（研究类型判断）、未提及变异时自动生成 >200 字文献简介（180+ 主题词映射表）、Excel 汇总表输出（UTF-8 BOM CSV 格式，零依赖）、文献简介和 Excel 输出均为预置逻辑确保可重复** |
| **v9** | 2026-05 | **修复 in-cis 识别：支持连字符 `in-cis` 匹配、提取顺式复杂等位基因共存变异（`[p.X;p.Y]` 格式）、Excel 列重命名为"反式(trans)/顺式(cis)位点"并显示顺式标签、总结段落包含顺式复杂等位基因详情** |