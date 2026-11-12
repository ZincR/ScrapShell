"""Write the results table (from web/data/manifest.json) into README.md below <!-- results-table -->.

Usage: python scripts/results_table.py
"""
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
manifest = json.loads((ROOT / "web" / "data" / "manifest.json").read_text())
models = manifest["models"] if isinstance(manifest["models"], list) else [manifest["models"]]


def f(x, fmt):
    return "—" if x is None else format(x, fmt)


rows = [
    "| Model | Faces in | Doublings | LM iters | Singularities | Panels (hex + other) "
    "| Planarity max before → after | Planarity mean | Distance max / bbd |",
    "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
]
for m in models:
    k = m.get("metrics", {}) or {}
    stages = m["stages"] if isinstance(m["stages"], list) else [m["stages"]]
    data = json.loads((ROOT / "web" / "data" / m["file"]).read_text())
    doublings = data.get("params", {}).get("iter_opt")
    if max(stages) < 8:
        rows.append(f"| {m['label']} | {m['nf']} | {doublings} | {k.get('lmIterations', '—')} | — | — | "
                    f"stopped at stage {max(stages)} | — | — |")
        continue
    rows.append(
        f"| {m['label']} | {m['nf']} | {doublings} | {k['lmIterations']} | {k['nSing']} "
        f"| {k['nHex']} + {k['nOther']} | {f(k['planMaxBefore'], '.2f')}% → {f(k['planMax'], '.2f')}% "
        f"| {f(k['planMean'], '.3f')}% | {f(k['distMax'], '.4f')} |"
    )
table = "\n".join(rows)
note = (f"\n\nGenerated {manifest['generated']}; field solved with {manifest.get('solver', 'lsqnonlin')}. "
        "Paper (Fig. 24/25): 255–1449 faces, planarity max 0.3–1.0% (target 1%), mean 0.03–0.49%, "
        "distance 0.001–0.015 × bbd.\n")
readme = ROOT / "README.md"
text = readme.read_text()
marker = "<!-- results-table -->"
head, _, tail = text.partition(marker)
tail = re.sub(r"^.*?(?=\n## )", "", tail, count=1, flags=re.S)
readme.write_text(head + marker + "\n\n" + table + note + tail)
print(table)
