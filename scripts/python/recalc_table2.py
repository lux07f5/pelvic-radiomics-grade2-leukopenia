"""recalc_table2.py - Wilcoxon rank-sum test for dose-volume metrics
by event status (60-patient cohort). DE-IDENTIFIED VERSION.

Usage:
    python recalc_table2.py <path_to_clinical_json> <path_to_dose_json>

Both JSON inputs must be in the format produced by the local pipeline.
The clinical JSON should be a list of patients with 'name' (de-identified)
and 'y' (0/1 outcome) fields. The dose JSON should have 'event' and
'nonevent' sublists with 'name' and 'vals' (9 metric values) fields.
"""
import json
import math
import statistics
import sys


def wilcoxon_p(a, b):
    """Two-sided Wilcoxon rank-sum (Mann-Whitney U) p-value, normal approximation."""
    na, nb = len(a), len(b)
    if na < 1 or nb < 1:
        return None
    allv = sorted([(v, 'a') for v in a] + [(v, 'b') for v in b])
    r = {v: i + 1 for i, (v, _) in enumerate(allv)}
    R1 = sum(r[v] for v in a)
    N = na + nb
    U1 = R1 - na * (na + 1) / 2
    U2 = na * nb - U1
    U = min(U1, U2)
    mu = na * nb / 2
    sigma = math.sqrt(na * nb * (N + 1) / 12)
    if sigma == 0:
        return 1.0
    z = (U - mu) / sigma
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


def main():
    if len(sys.argv) < 3:
        print("Usage: python recalc_table2.py <clinical.json> <dose.json>")
        sys.exit(1)

    with open(sys.argv[1]) as f:
        m = json.load(f)
    with open(sys.argv[2]) as f:
        dose = json.load(f)

    ev_names = [p['name'] for p in m if p['y'] == 1]
    ne_names = [p['name'] for p in m if p['y'] == 0]

    # vals: [0]Volume, [1]min, [2]max, [3]mean, [4]V5, [5]V10, [6]V20, [7]V35, [8]V40
    labels = ['Volume(cm3)', 'meandose(cGy)', 'V5(%)', 'V10(%)', 'V20(%)', 'V35(%)', 'V40(%)']
    idx_d = [0, 3, 4, 5, 6, 7, 8]
    ev_dose = [d for d in dose['event']['data'] if d['name'] in ev_names]
    ne_dose = [d for d in dose['nonevent']['data'] if d['name'] in ne_names]
    print(f'Event: n={len(ev_dose)}, Non-event: n={len(ne_dose)}')

    print(f'\n=== Table 2: Dose-volume metrics by event status ===')
    print(f'{"Metric":<16} {"Event":>14} {"Non-event":>14} {"P(Wilcox)":>10}')
    for li, lbl in enumerate(labels):
        i = idx_d[li]
        ev = [d['vals'][i] for d in ev_dose]
        ne = [d['vals'][i] for d in ne_dose]
        p = wilcoxon_p(ev, ne)
        sd_e = statistics.stdev(ev) if len(ev) > 1 else 0
        sd_n = statistics.stdev(ne) if len(ne) > 1 else 0
        print(f'{lbl:<16} {statistics.mean(ev):>7.1f}±{sd_e:>4.1f} '
              f'{statistics.mean(ne):>7.1f}±{sd_n:>4.1f} {p:>10.3f}')


if __name__ == '__main__':
    main()
