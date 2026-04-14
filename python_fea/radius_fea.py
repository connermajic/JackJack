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
    from skfem import Basis, MeshTet, MeshTri, asm
    from skfem.element import ElementTetP1, ElementTriP1, ElementVector
    from skfem.models.elasticity import lame_parameters, linear_elasticity
    import matplotlib.pyplot as plt

    mods["sk_filters"] = skimage.filters
    mods["sk_measure"] = skimage.measure
    mods["sk_morph"] = skimage.morphology
    mods["nib"] = nib
    mods["tifffile"] = tifffile
    mods["tetgen"] = tetgen
    mods["Basis"] = Basis
    mods["MeshTet"] = MeshTet
    mods["MeshTri"] = MeshTri
    mods["asm"] = asm
    mods["ElementTetP1"] = ElementTetP1
    mods["ElementTriP1"] = ElementTriP1
    mods["ElementVector"] = ElementVector
    mods["lame_parameters"] = lame_parameters
    mods["linear_elasticity"] = linear_elasticity
    mods["plt"] = plt
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

    analysis_mode: str = "3d"  # "3d" or "2d"
    slice_axis: int = 2
    slice_index: int | None = None
    mesh_step_2d: int = 3

    save_visualizations: bool = False
    viz_prefix: str = "radius_fea"


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
    thr = float(sk_filters.threshold_otsu(v)) if params.threshold is None else float(params.threshold)
    mask = v >= thr

    if params.fill_holes:
        if mask.ndim == 3:
            mask = sk_morph.remove_small_holes(mask, area_threshold=64)
        else:
            mask = sk_morph.remove_small_holes(mask, area_threshold=32)

    if params.morph_open_radius_vox > 0:
        se = sk_morph.ball(params.morph_open_radius_vox) if mask.ndim == 3 else sk_morph.disk(params.morph_open_radius_vox)
        mask = sk_morph.binary_opening(mask, se)

    return mask.astype(bool), thr


def _project_offset_orthogonal_to_force(offset: np.ndarray, force: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(force)
    if n <= 1e-12:
        raise ValueError("force_total_n must be non-zero")
    u = force / n
    return offset - np.dot(offset, u) * u


def _extract_2d_slice(volume: np.ndarray, mask: np.ndarray, params: RadiusFEAParams) -> Tuple[np.ndarray, np.ndarray, int]:
    axis = int(params.slice_axis)
    if axis not in (0, 1, 2):
        raise ValueError("slice_axis must be 0, 1, or 2")

    n = volume.shape[axis]
    idx = n // 2 if params.slice_index is None else int(np.clip(params.slice_index, 0, n - 1))

    if axis == 0:
        return volume[idx, :, :], mask[idx, :, :], idx
    if axis == 1:
        return volume[:, idx, :], mask[:, idx, :], idx
    return volume[:, :, idx], mask[:, :, idx], idx


def mask_to_surface(mask: np.ndarray, params: RadiusFEAParams) -> Tuple[np.ndarray, np.ndarray]:
    sk_measure = _lazy_imports()["sk_measure"]
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


def _build_2d_mesh_from_mask(mask2d: np.ndarray, dx_mm: float, dy_mm: float, step: int):
    MeshTri = _lazy_imports()["MeshTri"]

    h, w = mask2d.shape
    step = max(int(step), 1)

    xs_m = np.arange(0, w, step, dtype=np.float64) * dx_mm * 1e-3
    ys_m = np.arange(0, h, step, dtype=np.float64) * dy_mm * 1e-3
    if xs_m.size < 2 or ys_m.size < 2:
        raise ValueError("2D slice too small for meshing at chosen mesh_step_2d")

    mesh = MeshTri.init_tensor(xs_m, ys_m)
    tri = mesh.p[:, mesh.t]
    centroids = tri.mean(axis=1).T

    cx = np.clip(np.round((centroids[:, 0] * 1e3) / dx_mm).astype(int), 0, w - 1)
    cy = np.clip(np.round((centroids[:, 1] * 1e3) / dy_mm).astype(int), 0, h - 1)
    inside = mask2d[cy, cx]

    keep = np.where(inside)[0]
    if keep.size == 0:
        raise RuntimeError("No 2D elements remain after masking; adjust threshold or slice selection.")

    return mesh.with_elements(keep)


def solve_homogeneous_fea_3d(nodes_mm: np.ndarray, tets: np.ndarray, params: RadiusFEAParams) -> Dict[str, Any]:
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
    load_nodes = distal_nodes if len(in_patch) == 0 else distal_nodes[np.asarray(in_patch, dtype=np.int32)]

    nodes_m = nodes_mm * 1e-3
    mesh = MeshTet(nodes_m.T, tets.T)
    basis = Basis(mesh, ElementVector(ElementTetP1()))

    lam, mu = lame_parameters(params.young_modulus_pa, params.poisson_ratio)
    K = asm(linear_elasticity(lam, mu), basis)

    ndofs = K.shape[0]
    f = np.zeros(ndofs, dtype=np.float64)
    nodal_dofs = basis.nodal_dofs

    force_per_node = force_n / max(load_nodes.size, 1)
    for n in load_nodes:
        f[nodal_dofs[0, n]] += force_per_node[0]
        f[nodal_dofs[1, n]] += force_per_node[1]
        f[nodal_dofs[2, n]] += force_per_node[2]

    fixed_dofs = nodal_dofs[:, proximal_nodes].ravel(order="F")
    fixed_mask = np.zeros(ndofs, dtype=bool)
    fixed_mask[fixed_dofs] = True
    free = np.where(~fixed_mask)[0]

    u = np.zeros(ndofs, dtype=np.float64)
    u[free] = spsolve(K[free][:, free], f[free])

    ux = u[nodal_dofs[0]]
    uy = u[nodal_dofs[1]]
    uz = u[nodal_dofs[2]]
    umag = np.sqrt(ux**2 + uy**2 + uz**2)

    return {
        "mode": "3d",
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
        "nodes_mm": nodes_mm,
        "elements": tets,
    }


def solve_homogeneous_fea_2d(mask2d: np.ndarray, params: RadiusFEAParams, spacing_mm: Tuple[float, float]) -> Dict[str, Any]:
    mods = _lazy_imports()
    Basis = mods["Basis"]
    asm = mods["asm"]
    ElementTriP1 = mods["ElementTriP1"]
    ElementVector = mods["ElementVector"]
    lame_parameters = mods["lame_parameters"]
    linear_elasticity = mods["linear_elasticity"]

    dx_mm, dy_mm = spacing_mm
    mesh = _build_2d_mesh_from_mask(mask2d, dx_mm, dy_mm, params.mesh_step_2d)

    basis = Basis(mesh, ElementVector(ElementTriP1()))
    lam, mu = lame_parameters(params.young_modulus_pa, params.poisson_ratio)
    K = asm(linear_elasticity(lam, mu), basis)

    p = mesh.p.T
    axis2d = min(int(params.axis_for_length), 1)
    coords = p[:, axis2d] * 1e3
    cmin, cmax = float(coords.min()), float(coords.max())
    extent = cmax - cmin

    proximal_cut = cmin + params.proximal_band_fraction * extent
    distal_cut = cmax - params.distal_band_fraction * extent

    proximal_nodes = np.where(coords <= proximal_cut)[0]
    distal_nodes = np.where(coords >= distal_cut)[0]
    if proximal_nodes.size == 0 or distal_nodes.size == 0:
        raise RuntimeError("Could not detect proximal/distal node sets in 2D mesh.")

    distal_centroid_mm_2d = (p[distal_nodes].mean(axis=0) * 1e3)
    force2d = np.asarray(params.force_total_n[:2], dtype=np.float64)
    if np.linalg.norm(force2d) < 1e-12:
        force2d = np.array([0.0, -500.0], dtype=np.float64)
    offset2d = _project_offset_orthogonal_to_force(np.asarray(params.load_offset_mm[:2], dtype=np.float64), force2d)
    load_point_mm = distal_centroid_mm_2d + offset2d

    tree = cKDTree((p[distal_nodes] * 1e3))
    in_patch = tree.query_ball_point(load_point_mm, r=float(params.load_patch_radius_mm))
    load_nodes = distal_nodes if len(in_patch) == 0 else distal_nodes[np.asarray(in_patch, dtype=np.int32)]

    ndofs = K.shape[0]
    f = np.zeros(ndofs, dtype=np.float64)
    nodal_dofs = basis.nodal_dofs

    force_per_node = force2d / max(load_nodes.size, 1)
    for n in load_nodes:
        f[nodal_dofs[0, n]] += force_per_node[0]
        f[nodal_dofs[1, n]] += force_per_node[1]

    fixed_dofs = nodal_dofs[:, proximal_nodes].ravel(order="F")
    fixed_mask = np.zeros(ndofs, dtype=bool)
    fixed_mask[fixed_dofs] = True
    free = np.where(~fixed_mask)[0]

    u = np.zeros(ndofs, dtype=np.float64)
    u[free] = spsolve(K[free][:, free], f[free])

    ux = u[nodal_dofs[0]]
    uy = u[nodal_dofs[1]]
    umag = np.sqrt(ux**2 + uy**2)

    load_point_3d = [float(load_point_mm[0]), float(load_point_mm[1]), 0.0]
    return {
        "mode": "2d",
        "slice_shape": [int(mask2d.shape[0]), int(mask2d.shape[1])],
        "load_point_target_mm": load_point_3d,
        "distal_centroid_mm": [float(distal_centroid_mm_2d[0]), float(distal_centroid_mm_2d[1]), 0.0],
        "force_total_n": [float(force2d[0]), float(force2d[1]), 0.0],
        "num_load_nodes": int(load_nodes.size),
        "num_fixed_nodes": int(proximal_nodes.size),
        "max_displacement_m": float(umag.max()),
        "mean_displacement_m": float(umag.mean()),
        "node_displacements_m": np.vstack([ux, uy]).T,
        "load_nodes": load_nodes,
        "fixed_nodes": proximal_nodes,
        "nodes_mm": p * 1e3,
        "elements": mesh.t.T,
    }


def save_visualizations(
    params: RadiusFEAParams,
    volume: np.ndarray,
    mask: np.ndarray,
    fea: Dict[str, Any],
    slice_info: Dict[str, Any],
) -> Dict[str, str]:
    plt = _lazy_imports()["plt"]
    out: Dict[str, str] = {}

    if fea["mode"] == "2d":
        nodes = np.asarray(fea["nodes_mm"])
        tri = np.asarray(fea["elements"])
        disp = np.asarray(fea["node_displacements_m"])
        umag = np.sqrt((disp**2).sum(axis=1))

        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        ax[0].imshow(slice_info["slice_mask"], cmap="gray")
        lp = fea["load_point_target_mm"]
        ax[0].scatter([lp[0] / params.voxel_size_mm[0]], [lp[1] / params.voxel_size_mm[1]], c="r", s=50, label="Load target")
        ax[0].set_title("2D Slice Segmentation + Load Target")
        ax[0].legend(loc="lower right")

        tpc = ax[1].tripcolor(nodes[:, 0], nodes[:, 1], tri, umag, shading="flat")
        ax[1].set_title("2D Displacement Magnitude (m)")
        ax[1].set_aspect("equal")
        fig.colorbar(tpc, ax=ax[1])

        fig.tight_layout()
        out_path = f"{params.viz_prefix}_2d_summary.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        out["summary_plot"] = out_path
    else:
        nodes = np.asarray(fea["nodes_mm"])
        disp = np.asarray(fea["node_displacements_m"])
        umag = np.sqrt((disp**2).sum(axis=1))

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        sample = np.linspace(0, nodes.shape[0] - 1, min(6000, nodes.shape[0]), dtype=int)
        sc = ax.scatter(nodes[sample, 0], nodes[sample, 1], nodes[sample, 2], c=umag[sample], s=2, cmap="viridis")
        lp = np.asarray(fea["load_point_target_mm"])
        ax.scatter([lp[0]], [lp[1]], [lp[2]], c="r", s=80, marker="x")
        ax.set_title("3D Node Displacement Magnitude + Load Target")
        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")
        ax.set_zlabel("Z (mm)")
        fig.colorbar(sc, ax=ax, shrink=0.6, label="|u| (m)")
        out_path = f"{params.viz_prefix}_3d_displacement.png"
        fig.tight_layout()
        fig.savefig(out_path, dpi=160)
        plt.close(fig)
        out["displacement_plot"] = out_path

    return out


def run_radius_fea_python(
    image_path: str,
    params: RadiusFEAParams,
    output_json: str | None = None,
    output_npz: str | None = None,
    mat_variable_name: str | None = None,
) -> Dict[str, Any]:
    vol = load_converted_volume(image_path, mat_variable_name=mat_variable_name)
    mask, threshold = segment_bone(vol, params)

    slice_info: Dict[str, Any] = {}
    if params.analysis_mode.lower() == "2d":
        _, mask2d, idx = _extract_2d_slice(vol, mask, params)
        if params.slice_axis == 0:
            spacing2d = (params.voxel_size_mm[2], params.voxel_size_mm[1])
        elif params.slice_axis == 1:
            spacing2d = (params.voxel_size_mm[2], params.voxel_size_mm[0])
        else:
            spacing2d = (params.voxel_size_mm[0], params.voxel_size_mm[1])
        fea = solve_homogeneous_fea_2d(mask2d, params, spacing2d)
        slice_info = {"slice_axis": int(params.slice_axis), "slice_index": int(idx), "slice_mask": mask2d}
        surface_v = np.empty((0, 3), dtype=np.float64)
        surface_f = np.empty((0, 3), dtype=np.int32)
        nodes_mm = np.asarray(fea["nodes_mm"]) if fea["nodes_mm"].shape[1] == 2 else np.empty((0, 2))
        elems = np.asarray(fea["elements"])
    else:
        surface_v, surface_f = mask_to_surface(mask, params)
        nodes_mm, elems = tetrahedralize_surface(surface_v, surface_f, params)
        fea = solve_homogeneous_fea_3d(nodes_mm, elems, params)

    fea_public = {
        k: v
        for k, v in fea.items()
        if k not in ("node_displacements_m", "load_nodes", "fixed_nodes", "nodes_mm", "elements")
    }

    result: Dict[str, Any] = {
        "input_path": image_path,
        "analysis_mode": params.analysis_mode.lower(),
        "volume_shape": list(vol.shape),
        "threshold": float(threshold),
        "num_surface_vertices": int(surface_v.shape[0]),
        "num_surface_faces": int(surface_f.shape[0]),
        "num_nodes": int(np.asarray(fea["nodes_mm"]).shape[0]),
        "num_elements": int(np.asarray(fea["elements"]).shape[0]),
        "params": asdict(params),
        "slice": {k: v for k, v in slice_info.items() if k != "slice_mask"},
        "fea": fea_public,
    }

    plots = {}
    if params.save_visualizations:
        plots = save_visualizations(params, vol, mask, fea, slice_info)
        result["plots"] = plots

    if output_json:
        Path(output_json).write_text(json.dumps(result, indent=2), encoding="utf-8")

    if output_npz:
        np.savez_compressed(
            output_npz,
            nodes_mm=np.asarray(fea["nodes_mm"]),
            elements=np.asarray(fea["elements"]),
            displacements_m=np.asarray(fea["node_displacements_m"]),
            load_nodes=np.asarray(fea["load_nodes"]),
            fixed_nodes=np.asarray(fea["fixed_nodes"]),
        )

    return result


def _parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Homogeneous distal-radius FEA from converted images (Python).")
    ap.add_argument("image_path", help="Input volume path (.mat, .nii/.nii.gz, .tif/.tiff)")
    ap.add_argument("--out-json", default="radius_fea_python_result.json")
    ap.add_argument("--out-npz", default="radius_fea_python_result.npz")
    ap.add_argument("--mat-variable", default=None)
    ap.add_argument("--mode", default="3d", choices=["3d", "2d"], help="Analysis mode")
    ap.add_argument("--force", nargs=3, type=float, default=[0.0, 0.0, -500.0], metavar=("FX", "FY", "FZ"))
    ap.add_argument("--offset", nargs=3, type=float, default=[2.0, 0.0, 0.0], metavar=("OX", "OY", "OZ"))
    ap.add_argument("--axis", type=int, default=2, choices=[0, 1, 2])
    ap.add_argument("--slice-axis", type=int, default=2, choices=[0, 1, 2])
    ap.add_argument("--slice-index", type=int, default=None)
    ap.add_argument("--mesh-step-2d", type=int, default=3)
    ap.add_argument("--voxel", nargs=3, type=float, default=[1.0, 1.0, 1.0], metavar=("SX", "SY", "SZ"))
    ap.add_argument("--threshold", type=float, default=None)
    ap.add_argument("--save-plots", action="store_true")
    ap.add_argument("--viz-prefix", default="radius_fea")
    return ap.parse_args()


def main() -> None:
    args = _parse_args()
    params = RadiusFEAParams(
        voxel_size_mm=tuple(args.voxel),
        force_total_n=tuple(args.force),
        load_offset_mm=tuple(args.offset),
        axis_for_length=args.axis,
        threshold=args.threshold,
        analysis_mode=args.mode,
        slice_axis=args.slice_axis,
        slice_index=args.slice_index,
        mesh_step_2d=args.mesh_step_2d,
        save_visualizations=args.save_plots,
        viz_prefix=args.viz_prefix,
    )

    result = run_radius_fea_python(
        image_path=args.image_path,
        params=params,
        output_json=args.out_json,
        output_npz=args.out_npz,
        mat_variable_name=args.mat_variable,
    )

    p = result["fea"]["load_point_target_mm"]
    print(f"Mode: {result['analysis_mode']}")
    print(f"Load application target (mm): [{p[0]:.4f}, {p[1]:.4f}, {p[2]:.4f}]")
    print(f"Max displacement (m): {result['fea']['max_displacement_m']:.6e}")
    print(f"Wrote: {args.out_json} and {args.out_npz}")
    if "plots" in result:
        print(f"Plots: {result['plots']}")


if __name__ == "__main__":
    main()
