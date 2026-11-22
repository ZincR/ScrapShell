"""Generate the stock: one deposit of flat scraps, in meters, shared by every model.

    python -m stock.generate                         # data/stock/deposit_s0.json, default parameters
    python -m stock.generate --seed 1 -o /tmp/deposit.json
    python -m stock.generate --config my_deposit.json     # override sizes, family parameters or limits

The deposit has n_scraps pieces, exactly a quarter from each family (stock/shapes.py). Sizes
(size = sqrt(area), in meters) are drawn from one distribution (stock/sizes.py), sorted, and dealt
out in consecutive groups of four, one per family in random order, so every family gets the same
size distribution. The order is then shuffled, so any run of ids is a balanced sample. Separate
random streams drive the sizes, the dealing and shuffling, and each family's shapes: changing one
family's parameters leaves the other families' scraps exactly as they were.

A config file (JSON) may set any of: n_scraps, sizes (a stock/sizes.py distribution), params
(keyword arguments per family) and limits (overrides of stock.shapes.Limits). The rest keep the
defaults in DEFAULTS. stock/place.py puts the deposit on a model.

File format (scrapshell-stock/2): top-level "name", "units" ("m"), "seed", "generator" (the
resolved configuration), "summary" and "scraps". Each scrap has "id", "family", "size" (m),
"polygon" (counter-clockwise vertices in m, not closed, with the incenter at the origin and the
long axis along x) and the descriptors of stock/describe.py (lengths in m, areas in m^2). Numbers
are rounded to 7 significant digits; the descriptors are computed on the rounded polygon.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from collections.abc import Mapping
from dataclasses import asdict
from pathlib import Path

import numpy as np

from stock import sizes
from stock.describe import describe
from stock.model_ref import ROOT
from stock.shapes import FAMILIES, Limits, check, generate_shape

FORMAT = "scrapshell-stock/2"
STOCK_DIR = ROOT / "data" / "stock"
DEFAULT_DEPOSIT = STOCK_DIR / "deposit_s0.json"
SIG_DIGITS = 7
DEFAULTS = {
    "n_scraps": 4000,
    "sizes": {"kind": "lognormal", "median": 0.5, "sigma": 0.3, "min": 0.15, "max": 1.5},
    "params": {},
    "limits": {},
}


def sig(x, digits: int = SIG_DIGITS):
    """Round floats (also inside lists, tuples, dicts and arrays) to significant digits."""
    if isinstance(x, (float, np.floating)):
        return float(f"{x:.{digits}g}")
    if isinstance(x, np.ndarray):
        return [sig(v, digits) for v in x.tolist()]
    if isinstance(x, (list, tuple)):
        return [sig(v, digits) for v in x]
    if isinstance(x, dict):
        return {k: sig(v, digits) for k, v in x.items()}
    return x


def resolve_config(config: Mapping | None = None) -> dict:
    """DEFAULTS overridden by config, checked; limits are filled in completely."""
    config = dict(config or {})
    unknown = set(config) - set(DEFAULTS)
    if unknown:
        raise ValueError(f"unknown config keys {sorted(unknown)}; allowed: {sorted(DEFAULTS)}")
    cfg = {**DEFAULTS, **config}
    n = cfg["n_scraps"]
    if not isinstance(n, int) or n <= 0 or n % len(FAMILIES):
        raise ValueError(f"n_scraps must be a positive multiple of {len(FAMILIES)} (one share per "
                         f"family), got {n!r}")
    sizes.validate(cfg["sizes"])
    unknown = set(cfg["params"]) - set(FAMILIES)
    if unknown:
        raise ValueError(f"unknown families {sorted(unknown)}; choose from {sorted(FAMILIES)}")
    return {"n_scraps": n, "families": list(FAMILIES), "sizes": cfg["sizes"],
            "params": {k: dict(v) for k, v in cfg["params"].items()},
            "limits": asdict(Limits(**cfg["limits"]))}


def generate_deposit(config: Mapping | None = None, seed: int = 0) -> dict:
    """The deposit document (see the module docstring)."""
    cfg = resolve_config(config)
    n, fams = cfg["n_scraps"], tuple(FAMILIES)
    streams = np.random.SeedSequence(seed).spawn(2 + len(fams))
    size_rng, deal_rng = (np.random.default_rng(s) for s in streams[:2])
    shape_rngs = {f: np.random.default_rng(s) for f, s in zip(fams, streams[2:])}

    size_sorted = np.sort(sizes.sample(size_rng, cfg["sizes"], n))
    # deal each consecutive group of len(fams) sizes out to the families, one each, in random order
    family_of = np.concatenate([deal_rng.permutation(len(fams)) for _ in range(n // len(fams))])
    position = deal_rng.permutation(n)              # id -> position in the sorted sizes

    limits = Limits(**cfg["limits"])
    stats = {f: Counter() for f in fams}
    scraps = []
    for i, pos in enumerate(position):
        f = fams[family_of[pos]]
        unit = generate_shape(shape_rngs[f], f, cfg["params"].get(f), limits, stats=stats[f])
        coords = np.array(sig(np.asarray(unit.exterior.coords)[:-1] * size_sorted[pos]))
        reason = check(coords, limits)
        if reason is not None:
            raise RuntimeError(f"scrap {i} became unusable after rounding: {reason}")
        d = describe(coords)
        scraps.append({"id": i, "family": f, "size": sig(math.sqrt(d.area)),
                       "polygon": coords.tolist(), **sig(d.to_dict())})

    size = np.array([s["size"] for s in scraps])
    summary = {
        "count": n,
        "totalArea": sig(sum(s["area"] for s in scraps)),
        "families": dict(Counter(s["family"] for s in scraps)),
        "size": sig({"min": size.min(), "p10": np.percentile(size, 10), "median": np.median(size),
                     "p90": np.percentile(size, 90), "max": size.max()}),
        "drawAcceptance": {f: sig(st["ok"] / max(sum(st.values()), 1)) for f, st in stats.items()},
    }
    return {"format": FORMAT, "name": f"deposit_s{seed}", "units": "m", "seed": seed,
            "generator": cfg, "summary": summary, "scraps": scraps}


def dumps_deposit(doc: Mapping) -> str:
    """JSON text with the header indented and one scrap per line (readable and diff-friendly)."""
    head = json.dumps({k: v for k, v in doc.items() if k != "scraps"}, indent=2)
    lines = ",\n    ".join(json.dumps(s, separators=(",", ":")) for s in doc["scraps"])
    return head[:-2] + ',\n  "scraps": [\n    ' + lines + "\n  ]\n}\n"


def load_deposit(path: Path = DEFAULT_DEPOSIT) -> dict:
    doc = json.loads(Path(path).read_text())
    if doc.get("format") != FORMAT:
        raise ValueError(f"{path} is not a {FORMAT} deposit (format {doc.get('format')!r})")
    return doc


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m stock.generate",
                                 description="Generate the deposit of flat scraps (in meters).")
    ap.add_argument("--seed", type=int, default=0, help="random seed (default 0)")
    ap.add_argument("--n-scraps", type=int, help=f"number of scraps (default {DEFAULTS['n_scraps']})")
    ap.add_argument("--config", type=Path, help="JSON file overriding n_scraps, sizes, params, limits")
    ap.add_argument("-o", "--output", type=Path, help="output file (default data/stock/deposit_s<seed>.json)")
    args = ap.parse_args(argv)

    config = json.loads(args.config.read_text()) if args.config else {}
    if args.n_scraps is not None:
        config["n_scraps"] = args.n_scraps
    try:
        doc = generate_deposit(config, args.seed)
    except ValueError as err:
        ap.error(str(err))
    out = args.output or STOCK_DIR / f"{doc['name']}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(dumps_deposit(doc))
    s = doc["summary"]
    print(f"{out}: {s['count']} scraps {s['families']}, total area {s['totalArea']:g} m^2, "
          f"size median {s['size']['median']:g} m (range {s['size']['min']:g}-{s['size']['max']:g} m)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
