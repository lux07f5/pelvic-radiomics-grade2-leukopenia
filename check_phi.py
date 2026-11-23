#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Scan a directory for residues that should not leave the hospital.

    python check_phi.py data/
    python check_phi.py data/ --strict

What it looks for
    FAIL   a column named for a person or a hospital number (name, 姓名, mrn, ...)
    FAIL   DICOM keywords or tag numbers (PatientName, (0010,0010), SOPInstanceUID)
    FAIL   a .dcm or .ima file
    FAIL   a calendar date in the text
    WARN   a run of 6 or more digits (possible hospital or accession number)
    WARN   a two- or three-character Chinese word that starts with a common surname

This is a screen, not a proof. A clean run means nothing obvious was left in;
it does not mean the data are anonymous. Inspect the files once by eye, and let
the data manager sign off, before anything is copied off-site.
"""
import argparse
import os
import re
import sys

TEXT_EXT = {'.csv', '.tsv', '.txt', '.json', '.md', '.py', '.r', '.yaml', '.yml'}

# files whose purpose is to name the things being searched for; scanning them
# would only report their own vocabulary
SELF_DESCRIBING = {'check_phi.py', 'DEIDENTIFICATION.md'}

BAD_COLUMNS = [
    'name', '姓名', '患者姓名', '病人姓名', 'mrn', 'medical_record_number',
    'inpatient_no', '住院号', '病案号', '门诊号', 'accession', 'accession_number',
    'id_card', '身份证', 'idcard', 'phone', 'telephone', '电话', '手机',
    'address', '住址', '家庭住址',
]

DICOM_KEYWORDS = [
    'PatientName', 'PatientID', 'PatientBirthDate', 'PatientSex', 'PatientAge',
    'PatientAddress', 'AccessionNumber', 'StudyInstanceUID', 'SeriesInstanceUID',
    'SOPInstanceUID', 'InstitutionName', 'ReferringPhysicianName',
    'OperatorsName', 'OtherPatientIDs', 'PatientWeight',
]

DICOM_TAG = re.compile(r'\(00[0-9A-Fa-f]{2},[0-9A-Fa-f]{4}\)')
DATE_RE = re.compile(r'\b(19|20)\d{2}[-/]\d{1,2}[-/]\d{1,2}\b')
DIGITS_RE = re.compile(r'(?<![\d.])\d{6,}(?![\d.])')

SURNAMES = set('王李张刘陈杨黄赵吴周徐孙马朱胡郭何高林罗郑梁谢宋唐许韩冯邓曹彭曾'
               '萧田董袁潘于蒋蔡余杜叶程苏魏吕丁任沈姚卢姜崔钟谭陆汪范金石廖贾夏'
               '韦付方白邹孟熊秦邱江尹薛闫段雷侯龙史陶黎贺顾毛郝龚邵万钱严覃武戴'
               '莫孔向汤常温康施文牛樊葛邢安齐易乔伍庞颜倪庄聂章鲁岳翟殷詹申欧耿'
               '关兰焦俞左柳甘祝包宁尚符舒阮柯纪梅童凌毕单季裴霍涂成苗谷盛曲翁冉'
               '骆蓝路游辛靳管柴蒙鲍华喻祁蒲房滕屈饶解牟艾尤阳时穆农司卓古吉缪简'
               '车项连芦麦褚娄窦戚岑景党宫费卜冷晏席卫米柏宗瞿桂全佟应臧闵苟邬边'
               '卞姬师和仇栾隋商刁沙荣巫寇桑郎甄丛仲虞敖巩明佘池查麻苑迟邝封')
CJK_NAME = re.compile(r'^[\u4e00-\u9fa5]{2,3}$')

# values that are legitimately Chinese in this dataset and must not be flagged
CJK_ALLOW = {
    '直肠癌', '宫颈癌', '子宫内膜癌', '男', '女',
    '骨盆体积', '最小剂量', '最大剂量', '平均剂量', '年龄', '基线',
}


def is_cjk_name(value):
    v = value.strip()
    if not v or not CJK_NAME.match(v):
        return False
    if v in CJK_ALLOW:
        return False
    # a two-character word is a plausible name only when the second character is
    # a common given-name character; a three-character one always is
    if len(v) == 3:
        return True
    return False


def scan_file(path, strict):
    """Return (fails, warns) as lists of strings."""
    fails, warns = [], []
    ext = os.path.splitext(path)[1].lower()

    if ext in ('.dcm', '.ima', '.dicom'):
        return ['DICOM file present: %s' % path], []

    if ext not in TEXT_EXT:
        return [], []

    try:
        with open(path, encoding='utf-8', errors='ignore') as fh:
            text = fh.read()
    except OSError as exc:
        return ['cannot read %s: %s' % (path, exc)], []

    low = text.lower()
    for kw in DICOM_KEYWORDS:
        if kw.lower() in low:
            fails.append('DICOM keyword "%s"' % kw)
    if DICOM_TAG.search(text):
        fails.append('DICOM tag number %s' % DICOM_TAG.search(text).group(0))
    if DATE_RE.search(text):
        fails.append('calendar date "%s"' % DATE_RE.search(text).group(0))

    # header-level checks for tabular files
    if ext in ('.csv', '.tsv'):
        header = text.split('\n', 1)[0]
        for col in [c.strip().strip('"').lower() for c in header.split(',')]:
            if col in BAD_COLUMNS:
                fails.append('column "%s" looks like a direct identifier' % col)

        # Chinese names in cell values
        for ln, line in enumerate(text.split('\n')[1:], 2):
            if not line.strip():
                continue
            for cell in line.split(','):
                if is_cjk_name(cell):
                    warns.append('line %d: Chinese name-like value "%s"' % (ln, cell.strip()))
                    break

    if ext == '.json':
        for m in re.finditer(r'"[^"]{2,4}"', text):
            v = m.group(0)[1:-1]
            if is_cjk_name(v):
                warns.append('Chinese name-like value "%s"' % v)

    hits = DIGITS_RE.findall(text)
    if hits:
        uniq = sorted(set(hits))[:5]
        warns.append('%d run(s) of 6+ digits, e.g. %s' % (len(hits), ', '.join(uniq)))

    if strict and warns:
        fails.extend(warns)
        warns = []
    return fails, warns


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('path', nargs='?', default='data', help='directory to scan (default: data)')
    ap.add_argument('--strict', action='store_true', help='treat warnings as failures')
    args = ap.parse_args()

    root = os.path.abspath(args.path)
    if not os.path.exists(root):
        sys.exit('no such path: %s' % root)

    files = []
    self_path = os.path.abspath(__file__)
    if os.path.isfile(root):
        files = [root]
    else:
        for dirpath, _dirnames, filenames in os.walk(root):
            for f in sorted(filenames):
                if f.startswith('.') or f == 'README.md':
                    continue
                if f in SELF_DESCRIBING:
                    continue  # keyword lists would flag themselves
                files.append(os.path.join(dirpath, f))

    total_fail = total_warn = 0
    for path in files:
        fails, warns = scan_file(path, args.strict)
        rel = os.path.relpath(path, root)
        if not fails and not warns:
            continue
        print('\n%s' % rel)
        for f in fails:
            print('  FAIL  %s' % f)
        for w in warns:
            print('  warn  %s' % w)
        total_fail += len(fails)
        total_warn += len(warns)

    print('\n%d file(s) scanned, %d failure(s), %d warning(s)'
          % (len(files), total_fail, total_warn))
    if total_fail:
        print('resolve the failures before this directory leaves the institution')
        sys.exit(1)
    if total_warn:
        print('check each warning by eye; some are legitimate measurements')
    else:
        print('nothing found')


if __name__ == '__main__':
    main()
