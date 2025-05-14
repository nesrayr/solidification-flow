import meshio
import numpy as np
from pathlib import Path

def compare(fenics_xdmf, ccx_vtu):
    # 1) FEniCS — читаем Uel
    xf = meshio.read(fenics_xdmf)
    try:
        u_f = xf.point_data["Uel"]
    except KeyError:
        raise KeyError(f"Поле 'Uel' не найдено в {fenics_xdmf}. Есть только: {list(xf.point_data.keys())}")

    # 2) CalculiX — ищем displacement
    xc = meshio.read(ccx_vtu)
    pd = xc.point_data
    # Вариант A: смещение под ключом "U"
    if "U" in pd:
        u_c = pd["U"]
    # Вариант B: компоненты под "U1","U2","U3"
    elif {"U1","U2","U3"}.issubset(pd):
        u_c = np.vstack((pd["U1"], pd["U2"], pd["U3"])).T
    else:
        raise KeyError(f"Поле смещения не найдено в {ccx_vtu}. Доступные поля: {list(pd.keys())}")

    # 3) Считаем ошибки
    diff = u_f - u_c
    e2   = np.linalg.norm(diff) / np.sqrt(diff.shape[0])
    einf = np.max(np.linalg.norm(diff, axis=1))
    print(f"L₂-ошибка: {e2:.3e}")
    print(f"L∞-ошибка: {einf:.3e}")

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--fenics", required=True, help="Путь к Uel_*.xdmf")
    p.add_argument("--ccx",    required=True, help="Путь к crosscheck.vtk/vtu")
    args = p.parse_args()
    compare(Path(args.fenics), Path(args.ccx))
