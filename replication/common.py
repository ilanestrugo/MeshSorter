"""Shared helpers for the table scripts.

Every table in the manuscript is produced by running the standalone simulator
`meshsorter` once per cell and collecting the results.  Nothing here changes the
experiment; the parameters live in the table scripts themselves, written out in
full, so that a reader can see exactly what was run.
"""
import json, os, subprocess, sys, hashlib

HERE = os.path.dirname(os.path.abspath(__file__))
BIN  = os.environ.get("MESHSORTER_BIN", os.path.join(HERE, "meshsorter"))
# Worker threads per cell.  The result never depends on this; only the wall
# clock does.  Override with MESHSORTER_THREADS.
THREADS = int(os.environ.get("MESHSORTER_THREADS", os.cpu_count() or 2))


def require_binary():
    if not os.path.exists(BIN):
        sys.exit(f"{BIN} not found.  Build it first:\n"
                 f"    c++ -O2 -std=c++17 -pthread -o meshsorter meshsorter.cpp")


def run(n, m, dual=True, stagger=False, extra="turnaround", buffers=None,
        steps=4_000_000, reps=10, warmup=None, seed=20260901, loop=None,
        loops=None, dp=4, turn=4, df=4, width=2):
    """Run one configuration and return the parsed JSON result."""
    cmd = [BIN, "-n", str(n), "-m", str(m),
           "--dual" if dual else "--single",
           "-T", str(steps), "-R", str(reps), "--seed", str(seed),
           "--dp", str(dp), "--turn", str(turn), "--df", str(df),
           "--width", str(width),
           "-t", str(THREADS), "--json"]
    if loop     is not None: cmd += ["-L", str(loop)]
    if loops    is not None: cmd += ["--loops", ",".join(map(str, loops))]
    if stagger:              cmd += ["--stagger", "--extra", extra]
    if buffers  is not None: cmd += ["-b", buffers if isinstance(buffers, str)
                                     else ",".join(map(str, buffers))]
    if warmup   is not None: cmd += ["-w", str(warmup)]
    # Results are cached under results/cache so that rerunning a table, or
    # running two tables that share a panel, costs nothing.  The key covers
    # every argument that can change the answer.  Delete the directory, or set
    # MESHSORTER_NOCACHE=1, to force a fresh run.
    key = hashlib.sha1(" ".join(cmd[1:]).replace(f"-t {THREADS}", "").encode()).hexdigest()[:16]
    cache = os.path.join(HERE, "results", "cache", key + ".json")
    if not os.environ.get("MESHSORTER_NOCACHE") and os.path.exists(cache):
        with open(cache) as f:
            return json.load(f)
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"meshsorter failed on n={n} m={m}: {out.stderr.strip()}")
    res = json.loads(out.stdout)
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    with open(cache, "w") as f:
        json.dump(res, f)
    return res


def grid(ns, ms, **kw):
    """Run a full n by m grid, printing progress to stderr."""
    res = {}
    total = len(ns) * len(ms)
    for k, n in enumerate(ns):
        for l, m in enumerate(ms):
            res[(n, m)] = run(n, m, **kw)
            print(f"  n={n} m={m}  {res[(n,m)]['throughput']:.5f} "
                  f"+-{res[(n,m)]['halfwidth']:.5f}   "
                  f"[{k*len(ms)+l+1}/{total}]", file=sys.stderr, flush=True)
    return res


def max_halfwidth(*grids):
    return max(r["halfwidth"] for g in grids for r in g.values())


def latex_rows(ns, ms, value, gain=None):
    """One LaTeX body row per n.  `value` and `gain` are callables of (n, m)."""
    lines = []
    for n in ns:
        cells = []
        for m in ms:
            v = value(n, m)
            cells.append(f"{v:.4f}" if gain is None
                         else f"{v:.4f} ({100*(v/gain(n,m)-1):.1f})")
        lines.append(f"{n} & " + " & ".join(cells) + r"\\")
    return "\n".join(lines)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)
    print(f"wrote {path}", file=sys.stderr)
