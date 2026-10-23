"""Model pilot: score PIE's test predictions (pie eval save_predictions=true; checkpoint best_auprc) like the other
models. PIE predicts per (context, perturbation) a log2 fold change (lfc_pred) and an expression-level shift
(delta_p_pred) for every gene of its axis; both are mapped onto the pilot genes (genes missing from PIE's axis, e.g.
PuroR, are left out of the Pearson) and scored by pilot_common.score against truth A (same guide) and B (other guide).
Both readouts are reported. delta_p_pred is the same kind of quantity as the pilot effect (a shift in mean log1p
expression); lfc_pred is a pseudobulk log2 fold change, whose Pearson over all genes is dominated by low-expression
genes. Usage: pilot_pie_score.py <data_dir> <pred_parquet> <out_json>"""
import sys, json
import numpy as np
import pyarrow.parquet as pq
from pilot_common import design, effects, score, STATES


def main():
    d, pf, outp = sys.argv[1:4]
    des, res = effects(d); genes = res[STATES[1]]["genes"]
    t = pq.read_table(pf); meta = json.loads(t.schema.metadata[b"pie"]); pg = meta["genes"]["cd4"]
    col = {g: j for j, g in enumerate(pg)}; idx = np.array([col.get(g, -1) for g in genes]); mask = idx >= 0
    df = t.to_pandas(); out = {}
    for field in ("lfc_pred", "delta_p_pred"):
        pred = {}
        for _, r in df.iterrows():
            if r["context"] != "stim":
                continue
            v = np.asarray(r[field], np.float32); p = np.zeros(len(genes), np.float32); p[mask] = v[idx[mask]]
            pred[r["perturbation"]] = p
        out[field] = score(pred, des, res, genes, gene_mask=mask)
    json.dump(dict(model="PIE", checkpoint="best_auprc", genes_scored=int(mask.sum()), scores=out, rules=__doc__),
              open(outp, "w"), indent=1)
    print(int(mask.sum()), json.dumps(out))


if __name__ == "__main__":
    main()
