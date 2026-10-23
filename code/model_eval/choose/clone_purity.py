import sys, json
sys.path.insert(0, '.')
import numpy as np, pandas as pd
from choose_clone import load_expr, call_clones
cells = open('export/cells.txt').read().strip().split('\n')
clone, st = call_clones(cells, 'cbc', 'cbc_libs.tsv')
m = pd.read_csv('export/meta.tsv', sep='\t', low_memory=False).set_index('cell').loc[cells]
m['clone'] = [clone.get(c) for c in cells]; m['lib'] = [int(c.split('_')[-1]) for c in cells]
h = m[m.clone.notna()]
g = h.groupby('clone')
size = g.size(); maj = g.gRNA.agg(lambda s: s.value_counts().iloc[0] / len(s)); nl = g.lib.nunique()
multi = size[(size >= 2) & (g.gRNA.nunique() > 1)].index
print('clones>=2:', int((size >= 2).sum()), 'multi-gRNA:', len(multi))
print('majority-gRNA share in multi-gRNA clones (cell-weighted):', float((maj[multi] * size[multi]).sum() / size[multi].sum()))
print('multi-gRNA clones by size:', pd.cut(size[multi], [1, 3, 10, 50, 1e6]).value_counts().sort_index().to_dict())
print('largest multi clones:', [(int(size[c]), int(g.get_group(c).gRNA.nunique()), round(float(maj[c]), 2), int(nl[c])) for c in size[multi].sort_values().index[-8:]])
print('clones spanning >1 library (size>=2):', float((nl[size >= 2] > 1).mean()))
