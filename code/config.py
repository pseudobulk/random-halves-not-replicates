"""Paths shared by the scripts. All locations live in one file, config.env at the repository root (copy
config.example.env and edit it); pipelines/run_all.sh reads the same file. Environment variables of the same name
override it, and RHNR_CONFIG points to another file.
Only the Tahoe-100M location is read here, because a few scripts look up the Tahoe metadata tables themselves; every
other input and output is a command-line argument."""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _read(path):
    out = {}
    if os.path.exists(path):
        for line in open(path, encoding="utf-8"):
            line = line.split("#", 1)[0].strip()
            if "=" in line:
                k, v = line.split("=", 1)
                out[k.strip()] = v.strip().strip('"').strip("'")
    return out


def get(key):
    if os.environ.get(key):
        return os.environ[key]
    cfg = _read(os.environ.get("RHNR_CONFIG") or os.path.join(ROOT, "config.env"))
    if key == "TAHOE_ROOT" and key not in cfg:
        return os.path.join(get("DATA_ROOT"), "tahoe")
    if key not in cfg:
        raise KeyError("%s is not set: copy config.example.env to config.env and fill it in, or export %s" % (key, key))
    v = cfg[key]
    for k in sorted(cfg, key=len, reverse=True):   # ${DATA_ROOT}-style references to other keys, as in the shell
        x = os.environ.get(k) or cfg[k]
        v = v.replace("${%s}" % k, x).replace("$%s" % k, x)
    return v


def tahoe(*parts):
    """path under the Tahoe data directory (TAHOE_ROOT, default DATA_ROOT/tahoe; holds Tahoe-100M/ and PIE_tahoe100m/)"""
    return os.path.join(get("TAHOE_ROOT"), *parts)


def tahoe_meta(name):
    """one of the Tahoe-100M metadata tables, e.g. tahoe_meta("sample_metadata.parquet")"""
    return tahoe("Tahoe-100M", "metadata", name)
