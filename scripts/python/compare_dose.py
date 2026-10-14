"""compare_dose.py - 通用化"两组 xlsx 数据对比"工具（去标识化版）

⚠️ DE-IDENTIFIED VERSION
原脚本含桌面绝对路径与 60 例患者姓名，已脱敏重构。
此脚本现在接收命令行参数，使用方式见下方。
"""
import sys
import os
from openpyxl import load_workbook
import json


def parse_xlsx(fp):
    """Parse an xlsx with one sheet containing patient data starting at row 3."""
    wb = load_workbook(fp, data_only=True, read_only=True)
    ws = wb[wb.sheetnames[0]]
    rows = list(ws.iter_rows(values_only=True))
    data = {}
    for r in rows[2:]:
        if r[0]:
            data[r[0]] = list(r[1:10])
    return data


def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_dose.py <group1_xlsx> <group2_xlsx> [old_json_for_comparison]")
        print()
        print("Both xlsx files must have format:")
        print("  Row 0: title (e.g. 'Group A')")
        print("  Row 1: column headers")
        print("  Row 2+: data (col 0 = patient ID, col 1-9 = 9 metric values)")
        print()
        print("Note: All patient identifiers should be replaced with P001-PNN codes")
        print("      before passing to this script (no PHI in input files).")
        sys.exit(1)

    new_data = {}
    for grp, fp in zip(['group1', 'group2'], sys.argv[1:3]):
        new_data[grp] = parse_xlsx(fp)

    print(f'Group 1: {len(new_data["group1"])} records')
    print(f'Group 2: {len(new_data["group2"])} records')

    # Optional comparison with old JSON
    if len(sys.argv) > 3:
        with open(sys.argv[3]) as f:
            old = json.load(f)
        # Patient list overlap analysis
        # (Adapt to your specific data structure as needed)
        print(f'\nLoaded old data from {sys.argv[3]}')


if __name__ == '__main__':
    main()
