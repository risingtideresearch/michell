"""Digitise the vector plots of SEAWAY report 1213 (Vugts' cylinders).

For each page: find subplots from their tick labels, calibrate the axes,
and return experiment markers (small filled shapes) and SEAWAY curves (long
polylines) in data coordinates.
"""
import fitz, json, sys, re


def num(s):
    s = s.replace('−', '-').replace(',', '.')
    try:
        return float(s)
    except ValueError:
        return None


def subplots(pg):
    words = pg.get_text('words')
    labels = []
    for w in words:
        v = num(w[4])
        if v is not None and len(w[4]) <= 6:
            x0, y0, x1, y1 = w[:4]
            labels.append(dict(v=v, cx=0.5 * (x0 + x1), cy=0.5 * (y0 + y1), x1=x1, x0=x0, y1=y1))
    # y-tick columns: same right edge, ≥3 labels, values rising upward
    cols = []
    used = set()
    for i, a in enumerate(labels):
        if i in used:
            continue
        grp = [j for j, b in enumerate(labels) if abs(b['x1'] - a['x1']) < 1.0]
        # split by vertical gaps > 40
        grp.sort(key=lambda j: labels[j]['cy'])
        chunks, cur = [], [grp[0]]
        for j in grp[1:]:
            if labels[j]['cy'] - labels[cur[-1]]['cy'] < 40 and labels[j]['v'] < labels[cur[-1]]['v']:
                cur.append(j)
            else:
                chunks.append(cur)
                cur = [j]
        chunks.append(cur)
        for ch in chunks:
            if len(ch) >= 3:
                vals = [labels[j]['v'] for j in ch]  # top to bottom
                if all(vals[k] > vals[k + 1] for k in range(len(vals) - 1)):
                    cols.append(ch)
                    used.update(ch)
    rows = []
    for i, a in enumerate(labels):
        grp = [j for j, b in enumerate(labels) if abs(b['cy'] - a['cy']) < 1.0 and j not in used]
        grp.sort(key=lambda j: labels[j]['cx'])
        # split by horizontal gaps > 45
        chunks, cur = [], [grp[0]] if grp else []
        for j in grp[1:]:
            if labels[j]['cx'] - labels[cur[-1]]['cx'] < 45 and labels[j]['v'] > labels[cur[-1]]['v']:
                cur.append(j)
            else:
                chunks.append(cur)
                cur = [j]
        if cur:
            chunks.append(cur)
        for ch in chunks:
            if len(ch) >= 3 and ch not in rows:
                vals = [labels[j]['v'] for j in ch]
                if all(vals[k] < vals[k + 1] for k in range(len(vals) - 1)):
                    rows.append(ch)
    out = []
    for col in cols:
        cl = [labels[j] for j in col]
        right = max(l['x1'] for l in cl)
        bottom = max(l['cy'] for l in cl)
        best = None
        for row in rows:
            rl = [labels[j] for j in row]
            ry = rl[0]['cy']
            if 0 < ry - bottom < 25 and abs(rl[0]['cx'] - right) < 12:
                if best is None or ry < best[0]:
                    best = (ry, rl)
        if best is None:
            continue
        rl = best[1]
        fit = lambda pts: (lambda n, sx, sy, sxx, sxy: ((n * sxy - sx * sy) / (n * sxx - sx * sx), (sy * sxx - sx * sxy) / (n * sxx - sx * sx)))(
            len(pts), sum(p[0] for p in pts), sum(p[1] for p in pts), sum(p[0] ** 2 for p in pts), sum(p[0] * p[1] for p in pts))
        ax, bx = fit([(l['cx'], l['v']) for l in rl])
        ay, by = fit([(l['cy'], l['v']) for l in cl])
        box = (min(l['cx'] for l in rl) - 2, min(l['cy'] for l in cl) - 2, max(l['cx'] for l in rl) + 2, max(l['cy'] for l in cl) + 2)
        out.append(dict(box=box, ax=ax, bx=bx, ay=ay, by=by,
                        xr=(rl[0]['v'], rl[-1]['v']), yr=(cl[-1]['v'], cl[0]['v'])))
    return out


def extract(page_index, path='1213-ValidationSEAWAY.pdf'):
    d = fitz.open(path)
    pg = d[page_index]
    sps = subplots(pg)
    drawings = pg.get_drawings()
    res = []
    for sp in sps:
        x0, y0, x1, y1 = sp['box']
        to = lambda px, py: (sp['ax'] * px + sp['bx'], sp['ay'] * py + sp['by'])
        markers, curves = [], []
        for dr in drawings:
            r = dr['rect']
            cx, cy = 0.5 * (r.x0 + r.x1), 0.5 * (r.y0 + r.y1)
            inside = x0 <= cx <= x1 and y0 <= cy <= y1
            if not inside:
                continue
            kinds = set(it[0] for it in dr['items'])
            if dr['type'] in ('f', 'fs') and r.width < 8 and r.height < 8 and r.width > 1.5:
                shape = 'square' if kinds == {'re'} else 'circle'
                markers.append((shape, *to(cx, cy)))
            elif dr['type'] == 's' and len(dr['items']) >= 20:
                pts = []
                for it in dr['items']:
                    if it[0] == 'l':
                        pts.append(to(it[1].x, it[1].y))
                        pts.append(to(it[2].x, it[2].y))
                # de-duplicate consecutive
                clean = [pts[0]]
                for p_ in pts[1:]:
                    if abs(p_[0] - clean[-1][0]) > 1e-9 or abs(p_[1] - clean[-1][1]) > 1e-9:
                        clean.append(p_)
                curves.append(clean)
        res.append(dict(box=sp['box'], xr=sp['xr'], yr=sp['yr'], markers=markers, curves=curves))
    return res


if __name__ == '__main__':
    pages = [int(a) for a in sys.argv[1:]]
    allres = {}
    for p in pages:
        r = extract(p - 1)
        allres[p] = r
        for i, sp in enumerate(r):
            print(f'page {p} subplot {i}: box {tuple(round(v) for v in sp["box"])} x {sp["xr"]} y {sp["yr"]} markers {len(sp["markers"])} curves {len(sp["curves"])}')
    json.dump(allres, open('vugts_digitised.json', 'w'))
