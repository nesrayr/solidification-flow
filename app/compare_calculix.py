import meshio
import numpy as np
from pathlib import Path

def compare(fenics_vtu, ccx_vtk):
    # FEniCS
    xf = meshio.read(fenics_vtu)
    try:
        u_f = xf.point_data["Uel"]
    except KeyError:
        raise KeyError(f"Поле 'Uel' не найдено в {fenics_vtu}. Есть только: {list(xf.point_data.keys())}")

    # CalculiX
    xc = meshio.read(ccx_vtk)
    pd = xc.point_data
    if "U" in pd:
        u_c = pd["U"]
    elif {"U1","U2","U3"}.issubset(pd):
        u_c = np.vstack((pd["U1"], pd["U2"], pd["U3"])).T
    else:
        raise KeyError(f"Поле смещения не найдено в {ccx_vtk}. Доступные поля: {list(pd.keys())}")

    # Считаем вектор ошибки и его норму
    diff = u_f - u_c
    e2   = np.linalg.norm(diff) / np.sqrt(diff.shape[0])
    einf = np.max(np.linalg.norm(diff, axis=1))
    print(f"L₂-ошибка: {e2:.3e}")
    print(f"L∞-ошибка: {einf:.3e}")

    # Генерируем error.vtu на той же сетке, что и CCX
    sim_dir    = ccx_vtk.parent
    mesh_xdmf  = sim_dir / "VTK" / "mesh_tet.xdmf"
    if not mesh_xdmf.exists():
        raise FileNotFoundError(f"Не нашёл сетку {mesh_xdmf}")
    mesh       = meshio.read(mesh_xdmf)
    points     = mesh.points
    cells_dict = mesh.cells_dict
    if "tetra" not in cells_dict:
        raise RuntimeError(f"В mesh_tet.xdmf нет тетра-элементов, есть: {list(cells_dict.keys())}")
    cells      = cells_dict["tetra"]

    error_vals = np.linalg.norm(diff, axis=1)

    error_mesh = meshio.Mesh(
        points=points,
        cells=[("tetra", cells)],
        point_data={"error": error_vals}
    )
    out_file = sim_dir / "error.vtu"
    meshio.write(str(out_file), error_mesh)
    print(f"Wrote error field to {out_file}")

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--fenics", required=True, help="Путь к Uel_*.vtu из FEniCS")
    p.add_argument("--ccx",    required=True, help="Путь к crosscheck.vtk из CalculiX")
    args = p.parse_args()
    compare(Path(args.fenics), Path(args.ccx))
