#!/usr/bin/env python3
"""
Sci-Hub PDF 下载器
根据 PMID 从 Sci-Hub 镜像站下载 PDF 全文。
流程：PMID → DOI (PubMed API) → Sci-Hub → PDF 下载
依次尝试三个镜像站：sci-hub.se → sci-hub.st → sci-hub.ru
"""

import argparse
import os
import random
import re
import sys
import time
from urllib.parse import urljoin

# Windows UTF-8 编码兼容
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')

try:
    import requests
except ImportError:
    print("错误: 需要安装 requests，请运行: pip install requests")
    sys.exit(1)

try:
    from bs4 import BeautifulSoup
except ImportError:
    print("错误: 需要安装 beautifulsoup4，请运行: pip install beautifulsoup4")
    sys.exit(1)

import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# ==================== 配置 ====================
OUTPUT_DIR = r"D:\claude_code\project1\sci"

MIRRORS = [
    "sci-hub.ru",
    "sci-hub.st",
    "sci-hub.se",
]

TIMEOUT = 30
DOWNLOAD_TIMEOUT = 120

# 人类行为模拟延迟（全部随机小数，禁止固定整数）
#  打开页面后等待 2~5 秒再操作
#  看到按钮到点击等待 0.8~2.5 秒
#  两次连续请求间隔 3~8 秒
#  一轮结束后休息 6~12 秒


# ==================== 工具函数 ====================

def create_session():
    session = requests.Session()
    session.verify = False
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125.0.0.0 Safari/537.36",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "Accept-Language": "en-US,en;q=0.9",
    })
    return session


def get_doi_from_pmid(pmid, session):
    """通过 PubMed E-utilities API 获取 PMID 对应的 DOI。"""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "pubmed",
        "id": pmid,
        "retmode": "json",
    }
    try:
        resp = session.get(url, params=params, timeout=TIMEOUT)
        resp.raise_for_status()
        # 打开页面后等待 2~5 秒再处理
        time.sleep(random.uniform(2, 5))
        data = resp.json()
        uid = str(pmid)
        article = data.get("result", {}).get(uid, {})
        doi = article.get("elocationid", "") or article.get("doi", "")
        if doi and doi.startswith("doi:"):
            doi = doi[4:].strip()
        if doi:
            return doi
    except Exception as e:
        print(f"    PubMed API 查询失败: {e}")
    return None


def find_pdf_url(soup, page_url):
    """从 Sci-Hub 文章页面中提取 PDF 下载链接。

    查找策略：
    1. <meta name="citation_pdf_url" content="...">
    2. <object type="application/pdf" data="...">
    3. <div class="download"><a href="...">
    4. <iframe src="...">
    5. <embed src="...">
    6. /storage/.../*.pdf 链接
    """
    # 策略1: meta citation_pdf_url
    meta = soup.find("meta", attrs={"name": "citation_pdf_url"})
    if meta and meta.get("content"):
        src = meta["content"]
        if src.startswith("//"):
            src = "https:" + src
        elif src.startswith("/"):
            src = urljoin(page_url, src)
        return src

    # 策略2: <object type="application/pdf">
    for obj in soup.find_all("object", attrs={"type": "application/pdf"}):
        data = obj.get("data", "")
        if data:
            if data.startswith("//"):
                data = "https:" + data
            elif data.startswith("/"):
                data = urljoin(page_url, data)
            return data

    # 策略3: <div class="download"> 中的链接
    download_div = soup.find("div", class_="download")
    if download_div:
        link = download_div.find("a", href=True)
        if link:
            href = link["href"]
            if href.startswith("//"):
                href = "https:" + href
            elif href.startswith("/"):
                href = urljoin(page_url, href)
            return href

    # 策略4: <iframe>
    for iframe in soup.find_all("iframe"):
        src = iframe.get("src", "")
        if src and not src.startswith("about:"):
            if src.startswith("//"):
                src = "https:" + src
            elif src.startswith("/"):
                src = urljoin(page_url, src)
            return src

    # 策略5: <embed>
    for embed in soup.find_all("embed"):
        src = embed.get("src", "")
        if src:
            if src.startswith("//"):
                src = "https:" + src
            return src

    # 策略6: <button onclick>
    for btn in soup.find_all("button"):
        onclick = btn.get("onclick", "")
        m = re.search(r"location\s*\.?\s*href\s*=\s*['\"]([^'\"]+)['\"]", onclick)
        if m:
            url = m.group(1)
            if url.startswith("//"):
                url = "https:" + url
            elif url.startswith("/"):
                url = urljoin(page_url, url)
            return url

    # 策略7: 页面文本中匹配 /storage/.../*.pdf
    page_text = str(soup)
    m = re.search(r'(/storage/[^\s"\']+\.pdf)', page_text, re.IGNORECASE)
    if m:
        return urljoin(page_url, m.group(1))

    # 策略8: 任意 .pdf 链接
    for a in soup.find_all("a", href=True):
        href = a["href"]
        if ".pdf" in href.lower():
            if href.startswith("//"):
                href = "https:" + href
            elif href.startswith("/"):
                href = urljoin(page_url, href)
            return href

    return None


def is_captcha_or_blocked(soup, html_text):
    """检测是否为真正的 CAPTCHA 拦截页面或屏蔽页面。

    Sci-Hub 正常文章页面也包含 altcha-widget（用于问题反馈），
    但同时会有 /storage/、citation_pdf_url 等 PDF 相关标记。
    只有在缺失这些 PDF 标记时才视为被拦截。
    """
    if len(html_text.strip()) < 200:
        return True

    # 如果有 PDF 相关标记，说明是正常的文章页面，不是拦截页
    if "/storage/" in html_text or "citation_pdf_url" in html_text:
        return False
    if "pdf2md" in html_text:
        return False
    if soup.find("object", attrs={"type": "application/pdf"}):
        return False

    # Cloudflare 拦截检测
    if soup.find("div", id="cf-wrapper") or soup.find("div", class_="cf-browser-verification"):
        return True

    # 纯 reCAPTCHA 拦截页（没有 PDF 标记的 CAPTCHA 页面）
    text_lower = html_text.lower()
    strong_captcha = ["recaptcha", "cf-captcha", "ddos protection"]
    for kw in strong_captcha:
        if kw in text_lower:
            return True

    return False


def try_mirror(mirror, doi, session):
    """从单个镜像获取 PDF。使用 DOI 访问 Sci-Hub。

    返回: (pdf_content, source_url) 或 (None, error_message)
    """
    protocols = ["http", "https"]

    for pi, proto in enumerate(protocols):
        if pi > 0:
            # 同镜像切换协议，短暂间隔
            time.sleep(random.uniform(1.0, 2.5))
        scihub_url = f"{proto}://{mirror}/{doi}"
        print(f"    尝试: {scihub_url}")

        try:
            resp = session.get(scihub_url, timeout=TIMEOUT, allow_redirects=True)
        except requests.exceptions.Timeout:
            print(f"      ✗ 超时")
            continue
        except requests.exceptions.ConnectionError:
            print(f"      ✗ 连接失败")
            continue
        except requests.exceptions.SSLError:
            continue
        except Exception as e:
            print(f"      ✗ 异常: {e}")
            continue

        if resp.status_code != 200:
            print(f"      ✗ HTTP {resp.status_code}")
            continue

        # 打开页面后等待 2~5 秒再操作（模拟人类阅读/反应时间）
        page_delay = random.uniform(2, 5)
        time.sleep(page_delay)

        html_text = resp.text
        soup = BeautifulSoup(html_text, "html.parser")

        if is_captcha_or_blocked(soup, html_text):
            print(f"      ✗ CAPTCHA/屏蔽页面 ({len(html_text)} 字节)")
            continue

        # 检查是否直接返回 PDF
        if resp.content[:5] == b"%PDF-":
            print(f"      ✓ 直接返回 PDF ({len(resp.content)/1024:.0f} KB)")
            return resp.content, scihub_url

        content_type = resp.headers.get("Content-Type", "")
        if "application/pdf" in content_type:
            return resp.content, scihub_url

        # 查找 PDF URL
        pdf_url = find_pdf_url(soup, resp.url)
        if not pdf_url:
            print(f"      ✗ 未找到 PDF 链接 ({len(html_text)} 字节)")
            continue

        print(f"      PDF: {pdf_url[:120]}...")

        # 看到按钮到点击间隔 0.8~2.5 秒
        time.sleep(random.uniform(0.8, 2.5))

        # 下载 PDF
        try:
            pdf_resp = session.get(
                pdf_url,
                timeout=DOWNLOAD_TIMEOUT,
                allow_redirects=True,
                headers={"Referer": scihub_url},
            )
        except requests.exceptions.Timeout:
            print(f"      ✗ PDF 下载超时")
            continue
        except Exception as e:
            print(f"      ✗ PDF 下载失败: {e}")
            continue

        if pdf_resp.status_code != 200:
            print(f"      ✗ PDF HTTP {pdf_resp.status_code}")
            continue

        # 验证是 PDF
        if pdf_resp.content[:5] == b"%PDF-" or "application/pdf" in pdf_resp.headers.get("Content-Type", ""):
            return pdf_resp.content, scihub_url

        print(f"      ✗ 返回内容非 PDF ({len(pdf_resp.content)} 字节)")

    return None, "所有协议均失败"


def validate_pdf(content):
    if not content or len(content) < 1000:
        return False, "文件太小"
    if content[:5] != b"%PDF-":
        return False, "非 PDF 格式"
    if b"%%EOF" not in content[-1024:]:
        return False, "PDF 不完整"
    return True, ""


def download_pmid(pmid, output_dir, session):
    pmid = str(pmid).strip()
    if not pmid.isdigit():
        print(f"错误: 无效的 PMID: {pmid}")
        return None

    output_path = os.path.join(output_dir, f"PubMed{pmid}.pdf")

    if os.path.exists(output_path):
        size_kb = os.path.getsize(output_path) / 1024
        print(f"⚠ 已存在: PubMed{pmid}.pdf ({size_kb:.0f} KB) → 跳过")
        return output_path

    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"PMID: {pmid}")

    # 步骤1: 获取 DOI
    print("  获取 DOI ...")
    doi = get_doi_from_pmid(pmid, session)
    if not doi:
        print("  ✗ 无法获取 DOI（PubMed API 无返回）")
        return None
    print(f"  DOI: {doi}")

    # 步骤2: 通过 DOI 从 Sci-Hub 下载
    for i, mirror in enumerate(MIRRORS, 1):
        if i > 1:
            # 切换镜像前休息 3~8 秒
            time.sleep(random.uniform(3, 8))
        print(f"\n  [{i}/{len(MIRRORS)}] 镜像: {mirror}")

        pdf_content, result = try_mirror(mirror, doi, session)

        if pdf_content:
            valid, msg = validate_pdf(pdf_content)
            if valid:
                with open(output_path, "wb") as f:
                    f.write(pdf_content)
                size_kb = len(pdf_content) / 1024
                print(f"\n  ✓ 成功! PubMed{pmid}.pdf ({size_kb:.0f} KB)")
                return output_path
            else:
                print(f"    ✗ 验证失败: {msg}")
        else:
            print(f"    ✗ {result}")

    print(f"\n  ✗ 所有镜像均无法下载")
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Sci-Hub PDF 下载器 — PMID → DOI → Sci-Hub → PDF",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例:\n  python scihub_downloader.py 20981092\n  python scihub_downloader.py 20981092 33090715",
    )
    parser.add_argument("pmids", nargs="*", help="PubMed ID，可多个")
    parser.add_argument("--pmids", dest="pmids_csv", help="逗号分隔 PMID 列表")
    parser.add_argument("--output", "-o", default=OUTPUT_DIR, help=f"输出目录（默认: {OUTPUT_DIR}）")

    args = parser.parse_args()

    pmids = list(args.pmids)
    if args.pmids_csv:
        pmids.extend(args.pmids_csv.split(","))
    pmids = list(dict.fromkeys(p.strip() for p in pmids if p.strip()))

    if not pmids:
        parser.print_help()
        print("\n错误: 请提供至少一个 PMID")
        sys.exit(1)

    print(f"PMID 数量: {len(pmids)}")
    print(f"Sci-Hub 镜像: {', '.join(MIRRORS)}")
    print(f"输出目录: {args.output}")

    session = create_session()

    # 预热：先访问 sci-hub.ru 首页获取 session cookie
    print("预热会话（访问 Sci-Hub 首页）...")
    try:
        session.get("http://sci-hub.ru/", timeout=15)
        time.sleep(random.uniform(2, 5))
        print("  预热完成\n")
    except Exception:
        pass

    success = 0
    failed = []
    for idx, pmid in enumerate(pmids):
        if idx > 0:
            # 一轮结束，休息 6~12 秒
            rest = random.uniform(6, 12)
            print(f"\n  ...休息 {rest:.1f}s 后继续下一轮...")
            time.sleep(rest)
        result = download_pmid(pmid, args.output, session)
        if result:
            success += 1
        else:
            failed.append(pmid)

    print(f"\n{'='*60}")
    print(f"完成: {success}/{len(pmids)} 成功")
    if failed:
        print(f"失败: {', '.join(failed)}")


if __name__ == "__main__":
    main()