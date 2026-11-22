"""A target model's size and its baseline panels, read from the viewer data (web/data/<id>.json).

The meshes are in arbitrary units (bounding-box diagonal 2.0 on vase, 58.9 on sphere2con_s0).
stock/place.py uses the surface area and diagonal read here to give every model a physical size.
The baseline's planarized panels, where the baseline got that far, serve as a comparison: A_ref is
their median area.
"""
from __future__ import annotations

import base64
import json
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np
import shapely

ROOT = Path(__file__).resolve().parents[1]
WEB_DATA = ROOT / "web" / "data"

_DTYPES = {"f32": np.float32, "i32": np.int32, "u8": np.uint8}


def decode_array(packed: dict) -> np.ndarray:
    """Unpack a {dtype, shape, b64} array written by matlab/json_array.m (row-major)."""
    dtype = _DTYPES[packed["dtype"]]
    shape = packed["shape"] if isinstance(packed["shape"], list) else [packed["shape"]]
    raw = np.frombuffer(base64.b64decode(packed["b64"]), dtype=dtype).reshape(shape)
    return raw.astype(np.float64) if dtype is np.float32 else raw.astype(np.int64)


def flatten_panel(P: np.ndarray) -> np.ndarray:
    """Project a (nearly planar) 3D panel onto its best-fit plane: (k, 2), counter-clockwise."""
    c = P.mean(axis=0)
    _, _, W = np.linalg.svd(P - c)
    q = np.c_[(P - c) @ W[0], (P - c) @ W[1]]
    x, y = q[:, 0], q[:, 1]
    if np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)) < 0:   # clockwise: flip the frame
        q[:, 1] = -q[:, 1]
    return q


@dataclass(frozen=True)
class ModelRef:
    """One model's surface and its baseline panels, in the model's own units.

    panels is empty for a model the baseline did not planarize (car_r, 3holes_r).
    """
    model: str
    label: str
    source: str              # path of the viewer data, relative to the repo root
    surface_area: float
    bbd: float               # bounding-box diagonal of the input mesh
    panels: tuple            # flattened planarized panels, each (k, 2), counter-clockwise

    @cached_property
    def panel_areas(self) -> np.ndarray:
        return np.array([abs(shapely.Polygon(p).area) if len(p) >= 3 else 0.0 for p in self.panels])

    @cached_property
    def panel_circumradii(self) -> np.ndarray:
        """Radius of the smallest circle enclosing each panel."""
        return np.array([shapely.minimum_bounding_radius(shapely.MultiPoint(p)) for p in self.panels])

    @property
    def a_ref(self) -> float:
        """A_ref: the median baseline panel area. Raises ValueError if there are no panels."""
        if not self.panels:
            raise ValueError(f"model {self.model!r} has no planarized panels (the baseline stopped "
                             f"before stage 8)")
        return float(np.median(self.panel_areas))

    def to_dict(self) -> dict:
        d = {"source": self.source, "surfaceArea": self.surface_area, "bbd": self.bbd,
             "nPanels": len(self.panels)}
        if self.panels:
            d["panelAreaMedian"] = self.a_ref
            d["panelCircumradiusMedian"] = float(np.median(self.panel_circumradii))
        return d


def load_model_ref(model: str, data_dir: Path = WEB_DATA) -> ModelRef:
    """Read web/data/<model>.json. Raises ValueError if there is no such file."""
    path = Path(data_dir) / f"{model}.json"
    if not path.exists():
        known = sorted(p.stem for p in Path(data_dir).glob("*.json") if not p.stem.endswith(".meta")
                       and p.stem != "manifest")
        raise ValueError(f"no viewer data for model {model!r} at {path}; known models: {known}")
    data = json.loads(path.read_text())
    stages = data["stages"]

    X = decode_array(stages["input"]["V"])
    T = decode_array(stages["input"]["F"])
    tri = np.cross(X[T[:, 1]] - X[T[:, 0]], X[T[:, 2]] - X[T[:, 0]])
    surface_area = 0.5 * float(np.linalg.norm(tri, axis=1).sum())

    panels = ()
    if stages.get("planar", {}).get("frames"):
        V = decode_array(stages["planar"]["frames"][-1]["V"])
        deg = decode_array(stages["dual"]["D"]).ravel()
        idx = decode_array(stages["dual"]["idx"]).ravel()
        starts = np.concatenate([[0], np.cumsum(deg)])
        panels = tuple(flatten_panel(V[idx[a:b]]) for a, b in zip(starts[:-1], starts[1:]))

    try:
        source = str(path.resolve().relative_to(ROOT))
    except ValueError:
        source = str(path)
    return ModelRef(model=model, label=data.get("label", model), source=source,
                    surface_area=surface_area, bbd=float(data["bbd"]), panels=panels)
