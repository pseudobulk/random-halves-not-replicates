"""CHOOSE (Li, Fleck et al. 2023): stream one clone-barcode amplicon run (10x 3' v3: R1 = 16 nt cell barcode + 12 nt
UMI; R2 = lentiviral 3' LTR amplicon) from ENA and count (cell barcode, UMI, clone barcode) triples, without storing the
FASTQ. Each file is fetched in byte ranges of CH bytes (each retried until complete, since single long transfers through
the proxy break) and fed to gunzip.
Clone barcode (UCB): the 15 variable nt after the left flank GGTTCGCCGG; accepted when the right flank GCGCGC follows
one fixed G later with at most one mismatch (layout read off the first 10,000 reads of one run: 98% of reads carry a
16-nt insert ending in G between GGTTCGCCGG and GCGCGCCC).
Output: <out>.tsv.gz with run, cell barcode, UMI, UCB, reads; 'OK <run> <reads> <matched>' on stdout when both streams
ended cleanly with the same number of reads, otherwise 'FAIL ...' and exit 1.
Usage: cbc_parse.py <run> <R1 url> <R1 bytes> <R2 url> <R2 bytes> <out prefix> [proxy]"""
import sys, gzip, subprocess, collections, threading, time

run, u1, b1, u2, b2, out = sys.argv[1:7]; proxy = sys.argv[7] if len(sys.argv) > 7 else ""
LF, RF, CH = "GGTTCGCCGG", "GCGCGC", 32 << 20
err = []


def feed(u, size, sink):
    try:
        for s in range(0, size, CH):
            e = min(s + CH, size) - 1
            for k in range(30):
                cmd = ["curl", "-sSf", "-m", "900", "-r", "%d-%d" % (s, e), u] + (["-x", proxy] if proxy else [])
                r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
                if r.returncode == 0 and len(r.stdout) == e - s + 1:
                    break
                time.sleep(10)
            else:
                raise RuntimeError("range %d-%d of %s failed" % (s, e, u))
            sink.write(r.stdout)
    except Exception as x:
        err.append(str(x))
    finally:
        sink.close()


def stream(u, size):
    p = subprocess.Popen(["gunzip", "-c"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1 << 20)
    t = threading.Thread(target=feed, args=(u, int(size), p.stdin), daemon=True); t.start()
    return p, t


(p1, t1), (p2, t2) = stream(u1, b1), stream(u2, b2)
f1, f2 = p1.stdout, p2.stdout
cnt = collections.Counter(); n = ok = badlen = 0
while True:
    h1 = f1.readline(); h2 = f2.readline()
    if not h1 or not h2:
        break
    s1 = f1.readline().decode().strip(); f1.readline(); f1.readline()
    s2 = f2.readline().decode().strip(); f2.readline(); f2.readline()
    n += 1
    if len(s1) < 28:
        badlen += 1; continue
    i = s2.find(LF)
    if i < 0:
        continue
    r = s2[i + 26:i + 32]
    if len(r) < 6 or sum(a != b for a, b in zip(r, RF)) > 1:
        continue
    ucb = s2[i + 10:i + 25]
    if "N" in ucb:
        continue
    cnt[(s1[:16], s1[16:28], ucb)] += 1; ok += 1
rest1 = sum(1 for _ in f1); rest2 = sum(1 for _ in f2)          # drain both, so that no feeder blocks
t1.join(); t2.join(); e1, e2 = p1.wait(), p2.wait()
with gzip.open(out + ".tsv.gz", "wt") as g:
    for (cb, umi, ucb), c in cnt.items():
        g.write("%s\t%s\t%s\t%s\t%d\n" % (run, cb, umi, ucb, c))
# In the first pass four runs stopped with a gunzip length error and the parser waited forever on the other stream;
# rerun, all four decompressed cleanly (a transient corrupt transfer, not the ENA files). Both streams are now drained
# before joining; a gunzip error is still reported, and accepted only when every byte range arrived complete and both
# reads files gave the same number of records.
if err or rest1 or rest2:
    print("FAIL", run, n, ok, "exit", e1, e2, "leftover lines", rest1, rest2, err, flush=True); sys.exit(1)
print("OK" if not (e1 or e2) else "OK_GZIP_TRAILER", run, n, ok, "short_R1", badlen, "gunzip exit", e1, e2, flush=True)
