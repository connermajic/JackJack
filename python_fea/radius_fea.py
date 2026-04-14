from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
from scipy.io import loadmat
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree


def _lazy_imports() -> Dict[str, Any]:
    mods: Dict[str, Any] = {}
    import skimage.filters
    import skimage.measure
    import skimage.morphology
    import nibabel as nib
    import tifffile
    import tetgen
    import meshio
    from skfem import Basis, MeshTet, asm
    from skfem.element import ElementTetP1, ElementVector
    from skfem.models.elasticity import lame_parameters, linear_elasticity

    mods["sk_filters"] = skimage.filters
    mods["sk_measure"] = skimage.measure
    mods["sk_morph"] = skimage.morphology
    mods["nib"] = nib
    mods["tifffile"] = tifffile
    mods["tetgen"] = tetgen
    mods["meshio"] = meshio
    mods["Basis"] = Basis
    mods["MeshTet"] = MeshTet
    mods["asm"] = asm
    mods["ElementTetP1"] = ElementTetP1
    mods["ElementVector"] = ElementVector
    mods["lame_parameters"] = lame_parameters
    mods["linear_elasticity"] = linear_elasticity
    return mods


@dataclass
class RadiusFEAParams:
    voxel_size_mm: Tuple[float, float, float] = (1.0, 1.0, 1.0)
    volume_origin_mm: Tuple[float, float, float] = (0.0, 0.0, 0.0)

    threshold: float | None = None
    invert_mask: bool = False
    morph_open_radius_vox: int = 1
    fill_holes: bool = True

    young_modulus_pa: float = 15e9
    poisson_ratio: float = 0.30

    axis_for_length: int = 2  # 0=x, 1=y, 2=z
    proximal_band_fraction: float = 0.03
    distal_band_fraction: float = 0.03

    force_total_n: Tuple[float, float, float] = (0.0, 0.0, -500.0)
    load_offset_mm: Tuple[float, float, float] = (2.0, 0.0, 0.0)
    load_patch_radius_mm: float = 3.0

    tet_max_volume_mm3: float | None = None
    write_mesh_vtu: bool = False


def load_converted_volume(path: str, mat_variable_name: str | None = None) -> np.ndarray:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)

    lower = p.name.lower()
    if lower.endswith(".nii") or lower.endswith(".nii.gz"):
        nib = _lazy_imports()["nib"]
        return np.asarray(nib.load(str(p)).get_fdata(), dtype=np.float64)

    if lower.endswith(".tif") or lower.endswith(".tiff"):
        tifffile = _lazy_imports()["tifffile"]
        arr = tifffile.imread(str(p))
        if arr.ndim == 2:
            arr = arr[:, :, None]
        return np.asarray(arr, dtype=np.float64)

    if lower.endswith(".mat"):
        d = loadmat(str(p))
        if mat_variable_name:
            if mat_variable_name not in d:
                raise KeyError(f"{mat_variable_name} not found in {path}")
            v = d[mat_variable_name]
        else:
            cands = [(k, v) for k, v in d.items() if isinstance(v, np.ndarray) and v.ndim == 3]
            if not cands:
                raise ValueError("No 3-D numeric array found in mat file.")
            cands.sort(key=lambda kv: kv[1].size, reverse=True)
            v = cands[0][1]
        return np.asarray(v, dtype=np.float64)

    raise ValueError(f"Unsupported input extension: {p.suffix}")


def segment_bone(volume: np.ndarray, params: RadiusFEAParams) -> Tuple[np.ndarray, float]:
    mods = _lazy_imports()
    sk_filters = mods["sk_filters"]
    sk_morph = mods["sk_morph"]

    v = -volume if params.invert_mask else volume
    if params.threshold is None:
        thr = float(sk_filters.threshold_otsu(v))
    else:
        thr = float(params.threshold)

    mask = v >= thr

    if params.fill_holes:
        mask = sk_morph.remove_small_holes(mask, area_threshold=64)

    if params.morph_open_radius_vox > 0:
        se = sk_morph.ball(params.morph_open_radius_vox)
        mask = sk_morph.binary_opening(mask, se)

    return mask.astype(bool), thr


def _project_offset_orthogonal_to_force(offset_mm: np.ndarray, force_n: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(force_n)
    if n <= 1e-12:
        raise ValueError("force_total_n must be non-zero")
    u = force_n / n
    return offset_mm - np.dot(offset_mm, u) * u


def mask_to_surface(mask: np.ndarray, params: RadiusFEAParams) -> Tuple[np.ndarray, np.ndarray]:
    mods = _lazy_imports()
    sk_measure = mods["sk_measure"]

    spacing = np.asarray(params.voxel_size_mm, dtype=np.float64)
    verts, faces, _, _ = sk_measure.marching_cubes(mask.astype(np.float32), level=0.5, spacing=spacing)
    verts = verts + np.asarray(params.volume_origin_mm)[None, :]
    return verts, faces.astype(np.int32)


def tetrahedralize_surface(surface_vertices_mm: np.ndarray, surface_faces: np.ndarray, params: RadiusFEAParams) -> Tuple[np.ndarray, np.ndarray]:
    tetgen = _lazy_imports()["tetgen"]
    tgen = tetgen.TetGen(surface_vertices_mm, surface_faces)
    kwargs: Dict[str, Any] = {"mindihedral": 20.0, "minratio": 1.5}
    if params.tet_max_volume_mm3 is not None:
        kwargs["maxvolume"] = float(params.tet_max_volume_mm3)
    nodes, tets = tgen.tetrahedralize(**kwargs)
    return np.asarray(nodes, dtype=np.float64), np.asarray(tets, dtype=np.int32)


def solve_homogeneous_fea(
    nodes_mm: np.ndarray,
    tets: np.ndarray,
    params: RadiusFEAParams,
) -> Dict[str, Any]:
    mods = _lazy_imports()
    MeshTet = mods["MeshTet"]
    Basis = mods["Basis"]
    asm = mods["asm"]
    ElementTetP1 = mods["ElementTetP1"]
    ElementVector = mods["ElementVector"]
    lame_parameters = mods["lame_parameters"]
    linear_elasticity = mods["linear_elasticity"]

    axis = int(params.axis_for_length)
    if axis not in (0, 1, 2):
        raise ValueError("axis_for_length must be 0, 1, or 2")

    coords = nodes_mm[:, axis]
    cmin = float(coords.min())
    cmax = float(coords.max())
    extent = cmax - cmin

    proximal_cut = cmin + params.proximal_band_fraction * extent
    distal_cut = cmax - params.distal_band_fraction * extent

    proximal_nodes = np.where(coords <= proximal_cut)[0]
    distal_nodes = np.where(coords >= distal_cut)[0]
    if proximal_nodes.size == 0 or distal_nodes.size == 0:
        raise RuntimeError("Could not detect proximal/distal node sets. Check axis and segmentation.")

    distal_centroid_mm = nodes_mm[distal_nodes].mean(axis=0)
    force_n = np.asarray(params.force_total_n, dtype=np.float64)
    offset_mm = _project_offset_orthogonal_to_force(np.asarray(params.load_offset_mm, dtype=np.float64), force_n)
    load_point_mm = distal_centroid_mm + offset_mm

    tree = cKDTree(nodes_mm[distal_nodes])
    in_patch = tree.query_ball_point(load_point_mm, r=float(params.load_patch_radius_mm))
    if len(in_patch) == 0:
        load_nodes = distal_nodes
    else:
        load_nodes = distal_nodes[np.asarray(in_patch, dtype=np.int32)]

    # skfem works in whatever units are used; switch to meters for SI consistency.
    nodes_m = nodes_mm * 1e-3
    mesh = MeshTet(nodes_m.T, tets.T)
    basis = Basis(mesh, ElementVector(ElementTetP1()))

    lam, mu = lame_parameters(params.young_modulus_pa, params.poisson_ratio)
    K = asm(linear_elasticity(lam, mu), basis)

    ndofs = K.shape[0]
    f = np.zeros(ndofs, dtype=np.float64)
    nodal_dofs = basis.nodal_dofs  # shape: (3, n_nodes)

    force_per_node = force_n / max(load_nodes.size, 1)
    for n in load_nodes:
        f[nodal_dofs[0, n]] += force_per_node[0]
        f[nodal_dofs[1, n]] += force_per_node[1]
        f[nodal_dofs[2, n]] += force_per_node[2]

    fixed_dofs = nodal_dofs[:, proximal_nodes].ravel(order="F")
    fixed_mask = np.zeros(ndofs, dtype=bool)
    fixed_mask[fixed_dofs] = True
    free = np.where(~fixed_mask)[0]

    Kff = K[free][:, free]
    ff = f[free]
    u = np.zeros(ndofs, dtype=np.float64)
    u[free] = spsolve(Kff, ff)

    ux = u[nodal_dofs[0]]
    uy = u[nodal_dofs[1]]
    uz = u[nodal_dofs[2]]
    umag = np.sqrt(ux**2 + uy**2 + uz**2)

    return {
        "load_point_target_mm": load_point_mm.tolist(),
        "distal_centroid_mm": distal_centroid_mm.tolist(),
        "force_total_n": force_n.tolist(),
        "num_load_nodes": int(load_nodes.size),
        "num_fixed_nodes": int(proximal_nodes.size),
        "max_displacement_m": float(umag.max()),
        "mean_displacement_m": float(umag.mean()),
        "node_displacements_m": np.vstack([ux, uy, uz]).T,
        "load_nodes": load_nodes,
        "fixed_nodes": proximal_nodes,
    }


def run_radius_fea_python(
    image_path: str,
    params: RadiusFEAParams,
    output_json: str | None = None,
    output_npz: str | None = None,
    mat_variable_name: str | None = None,
) -> Dict[str, Any]:
    vol = load_converted_volume(image_path, mat_variable_name=mat_variable_name)
    mask, threshold = segment_bone(vol, params)

    surface_v, surface_f = mask_to_surface(mask, params)
    nodes_mm, tets = tetrahedralize_surface(surface_v, surface_f, params)

    fea = solve_homogeneous_fea(nodes_mm, tets, params)

    result: Dict[str, Any] = {
        "input_path": image_path,
        "volume_shape": list(vol.shape),
        "threshold": float(threshold),
        "num_surface_vertices": int(surface_v.shape[0]),
        "num_surface_faces": int(surface_f.shape[0]),
        "num_tet_nodes": int(nodes_mm.shape[0]),
        "num_tets": int(tets.shape[0]),
        "params": asdict(params),
        "fea": {k: v for k, v in fea.items() if k not in ("node_displacements_m", "load_nodes", "fixed_nodes")},
    }

    if output_json:
        Path(output_json).write_text(json.dumps(result, indent=2), encoding="utf-8")

    if output_npz:
        np.savez_compressed(
            output_npz,
            nodes_mm=nodes_mm,
            tets=tets,
            displacements_m=fea["node_displacements_m"],
            load_nodes=fea["load_nodes"],
            fixed_nodes=fea["fixed_nodes"],
        )

    return result


def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Homogeneous distal-radius FEA from converted images (Python).")
    ap.add_argument("image_path", help="Input volume path (.mat, .nii/.nii.gz, .tif/.tiff)")
    ap.add_argument("--out-json", default="radius_fea_python_result.json")
    ap.add_argument("--out-npz", default="radius_fea_python_result.npz")
    ap.add_argument("--mat-variable", default=None)
    ap.add_argument("--force", nargs=3, type=float, default=[0.0, 0.0, -500.0], metavar=("FX", "FY", "FZ"))
    ap.add_argument("--offset", nargs=3, type=float, default=[2.0, 0.0, 0.0], metavar=("OX", "OY", "OZ"))
    ap.add_argument("--axis", type=int, default=2, choices=[0, 1, 2])
    ap.add_argument("--voxel", nargs=3, type=float, default=[1.0, 1.0, 1.0], metavar=("SX", "SY", "SZ"))
    ap.add_argument("--threshold", type=float, default=None)
    return ap.parse_args()


def main() -> None:
    args = _parse_args()
    params = RadiusFEAParams(
        voxel_size_mm=tuple(args.voxel),
        force_total_n=tuple(args.force),
        load_offset_mm=tuple(args.offset),
        axis_for_length=args.axis,
        threshold=args.threshold,
    )

    result = run_radius_fea_python(
        image_path=args.image_path,
        params=params,
        output_json=args.out_json,
        output_npz=args.out_npz,
        mat_variable_name=args.mat_variable,
    )

    p = result["fea"]["load_point_target_mm"]
    print(f"Load application target (mm): [{p[0]:.4f}, {p[1]:.4f}, {p[2]:.4f}]")
    print(f"Max displacement (m): {result['fea']['max_displacement_m']:.6e}")
    print(f"Wrote: {args.out_json} and {args.out_npz}")


if __name__ == "__main__":
    main()
