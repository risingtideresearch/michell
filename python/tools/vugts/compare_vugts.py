"""Compare this crate's 2-D coefficients and beam-wave loads with Vugts'
experiments and with SEAWAY's curves, as digitised from report 1213."""
import json, math, statistics

d = json.load(open('vugts_digitised.json'))
ours = {}
lines = open('ours.tsv').read().strip().split('\n')
head = lines[0].split('\t')
for ln in lines[1:]:
    c = ln.split('\t')
    if len(c) != len(head):
        continue
    row = dict(zip(head, c))
    ours.setdefault(row['case'], []).append({k: float(v) for k, v in row.items() if k != 'case'})

pages = {'circle': (9, 10), 'rect2': (11, 12), 'rect4': (13, 14), 'rect8': (15, 16)}


def place(box):
    x, y = box[0], box[1]
    col = 0 if x < 220 else (1 if x < 340 else 2)
    return col, y


def interp(case, key, w):
    rows = ours[case]
    for a, b in zip(rows, rows[1:]):
        if a["w'"] <= w <= b["w'"]:
            t = (w - a["w'"]) / (b["w'"] - a["w'"])
            va, vb = a[key], b[key]
            if key.startswith('ph'):
                dv = (vb - va + 180) % 360 - 180
                return va + t * dv
            return va + t * (vb - va)
    return None


def assign(case):
    hp, sp = pages[case]
    out = {}
    # heave page: rows by y
    subs = [s for s in d[str(hp)] if s['markers'] or s['curves']]
    ys = sorted(set(round(s['box'][1]) for s in subs))
    for s in subs:
        col, y = place(s['box'])
        row = ys.index(round(s['box'][1]))
        key = {(0, 0): 'a33', (0, 1): 'b33', (2, 0): 'F3', (2, 1): 'ph3'}.get((col, row))
        if key:
            out[key] = s
    subs = [s for s in d[str(sp)] if s['markers'] or s['curves']]
    ys = sorted(set(round(s['box'][1]) for s in subs))
    table = {(0, 0): 'a22', (0, 1): 'b22', (1, 0): 'a42', (1, 1): 'b42', (2, 0): 'F2', (2, 1): 'ph2',
             (0, 2): 'a44', (0, 3): 'b44', (1, 2): 'a24', (1, 3): 'b24', (2, 2): 'F4', (2, 3): 'ph4'}
    for s in subs:
        col, y = place(s['box'])
        row = min(range(len(ys)), key=lambda i: abs(ys[i] - s['box'][1]))
        key = table.get((col, row))
        if key and key not in out:
            out[key] = s
    return out


def stats(case, key, s):
    exp = [(m[1], m[2]) for m in s['markers'] if s['xr'][0] <= m[1] <= s['xr'][1]]
    errs = []
    for w, v in exp:
        o = interp(case, key, w)
        if o is None:
            continue
        if key.startswith('ph'):
            errs.append(abs((o - v + 180) % 360 - 180))
        else:
            errs.append(o - v)
    sea = []
    for curve in s['curves']:
        for w, v in curve[::2]:
            if w < 0.2:
                continue
            o = interp(case, key, w)
            if o is None:
                continue
            sea.append(abs((o - v + 180) % 360 - 180) if key.startswith('ph') else o - v)
    return exp, errs, sea


print(f"{'case':7} {'qty':4} {'n':>3} {'scale':>7} {'exp: median |err|':>18} {'(rel)':>6} {'SEAWAY: median |diff|':>22} {'(rel)':>6}")
for case in pages:
    for key, s in sorted(assign(case).items()):
        exp, errs, sea = stats(case, key, s)
        lo, hi = s['yr']
        scale = max(abs(v) for _, v in exp) if exp else (hi - lo)
        if key.startswith('ph'):
            scale = 180.0
        me = statistics.median([abs(e) for e in errs]) if errs else float('nan')
        ms = statistics.median([abs(e) for e in sea]) if sea else float('nan')
        print(f"{case:7} {key:4} {len(errs):3d} {scale:7.3f} {me:18.3f} {me/scale if scale else 0:6.0%} {ms:22.3f} {ms/scale if scale else 0:6.0%}")
