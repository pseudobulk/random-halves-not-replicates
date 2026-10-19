"""Which Tahoe-100M plates share drug conditions (only plates 6 and 14 are full replicates). Reads the sample metadata
under TAHOE_ROOT (config.env); prints the result. Usage: tahoe_plate_overlap.py"""
import ast, itertools
import pyarrow.parquet as pq
import config
sm = pq.read_table(config.tahoe_meta("sample_metadata.parquet")).to_pydict()
print(list(sm.keys()))
cond = {}
for s, d, dc, p in zip(sm["sample"], sm["drug"], sm["drugname_drugconc"], sm["plate"]):
    c = (d.strip(), str(ast.literal_eval(dc)[0][1]))
    if "DMSO" in c[0]:
        continue
    cond.setdefault(p, set()).add(c)
plates = sorted(cond, key=lambda x: int(x.replace("plate", "")))
print({p: len(cond[p]) for p in plates})
drugs = {p: {c[0] for c in cond[p]} for p in plates}
pairs = []
for a, b in itertools.combinations(plates, 2):
    nc, nd = len(cond[a] & cond[b]), len(drugs[a] & drugs[b])
    if nc or nd:
        pairs.append((a, b, nc, nd))
for x in sorted(pairs, key=lambda t: -t[2])[:15]:
    print(x)
alld = set().union(*drugs.values()); multi = [d for d in alld if sum(d in drugs[p] for p in plates) > 1]
print("drugs", len(alld), "on >1 plate", len(multi))
