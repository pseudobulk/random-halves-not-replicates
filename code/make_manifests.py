"""Download manifest for the public data this analysis reads (README, "Data"). Each line:
url <TAB> destination path relative to DATA_ROOT <TAB> size in bytes <TAB> checksum type <TAB> checksum
('sha256' for Hugging Face LFS files, 'md5' for figshare, 'none' for S3 multipart uploads). pipelines/download.sh
fetches the files and checks them. Three inputs are not in the manifest and are fetched by hand (README): the Nadig et
al. Jurkat and HepG2 single-cell files (GEO GSE264667), the Human Protein Atlas cell-line table rna_celline.tsv and
the VCC 2025 H1 hESC training data (Arc Virtual Cell Atlas).
Usage: make_manifests.py <out.tsv>"""
import sys, json, re, urllib.request

HF = "https://huggingface.co"


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "manifest"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read(), r.headers


def hf_files(repo, path=""):
    """All files under path (recursive), following the tree API's Link pagination."""
    url = "%s/api/datasets/%s/tree/main/%s?recursive=true&expand=false" % (HF, repo, path)
    out = []
    while url:
        body, hdr = get(url)
        for f in json.loads(body):
            if f.get("type") == "file":
                lfs = f.get("lfs") or {}
                out.append((f["path"], f.get("size", 0), lfs.get("oid")))
        m = re.search(r'<([^>]+)>;\s*rel="next"', hdr.get("Link") or "")
        url = m.group(1) if m else None
    return out


def hf_rows(repo, local, prefix="", keep=None):
    rows = []
    for p, s, sha in hf_files(repo, prefix):
        if keep is None or p.split("/")[-1] in keep:
            rows.append(("%s/datasets/%s/resolve/main/%s" % (HF, repo, p), "%s/%s" % (local, p), s,
                         "sha256" if sha else "none", sha or "-"))
    return rows


def figshare_rows(article, names, local):
    body, _ = get("https://api.figshare.com/v2/articles/%d/files?page_size=100" % article)
    rows = []
    for f in json.loads(body):
        if f["name"] in names:
            rows.append((f["download_url"], "%s/%s" % (local, f["name"]), f["size"], "md5",
                         f.get("supplied_md5") or f.get("computed_md5")))
    return rows


def s3_rows(keys, local):
    base = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
    body, _ = get(base + "?list-type=2&prefix=marson2025_data/")
    sizes = dict(re.findall(r"<Key>([^<]+)</Key>.*?<Size>(\d+)</Size>", body.decode()))
    return [(base + k, "%s/%s" % (local, k.split("/")[-1]), int(sizes[k]), "none", "-") for k in keys]


def main():
    out = sys.argv[1]
    rows = []
    # CD4+ T-cell Perturb-seq (Zhu et al. 2026): single cells of donors 1-2, resting and stimulated 48 h; 4-donor pseudobulk
    rows += s3_rows(["marson2025_data/%s.assigned_guide.h5ad" % n for n in ("D1_Rest", "D1_Stim48hr", "D2_Rest", "D2_Stim48hr")], "cd4_sc")
    rows += s3_rows(["marson2025_data/GWCD4i.pseudobulk_merged.h5ad"], "cd4")
    # Tahoe-100M: all single-cell shards and the metadata tables; gene axis and DMSO means of the Arc PIE preprocessing
    rows += hf_rows("tahoebio/Tahoe-100M", "tahoe/Tahoe-100M", "data")
    rows += hf_rows("tahoebio/Tahoe-100M", "tahoe/Tahoe-100M", "metadata",
                    keep={"cell_line_metadata.parquet", "gene_metadata.parquet", "sample_metadata.parquet"})
    rows += hf_rows("arcinstitute/PIE_tahoe100m", "tahoe/PIE_tahoe100m", "preprocessed", keep={"meta.json", "ctrl_means.npy"})
    # Replogle et al. 2022 (figshare+ 20029387)
    rows += figshare_rows(20029387, {"K562_essential_raw_singlecell_01.h5ad", "rpe1_raw_singlecell_01.h5ad",
                                     "K562_gwps_raw_singlecell_01.h5ad"}, "replogle")
    with open(out, "w") as f:
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    print(len(rows), "files", round(sum(int(r[2]) for r in rows) / 1e9, 1), "GB")


if __name__ == "__main__":
    main()
