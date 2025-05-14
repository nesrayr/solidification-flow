import meshio
import numpy as np
from pathlib import Path

def compute_E_nu(lambda_, mu):
    E  = mu*(3*lambda_ + 2*mu)/(lambda_ + mu)
    nu = lambda_/(2*(lambda_ + mu))
    return E, nu

def export_inp(mesh_xdmf, temp_vtu, inp_path,
               lambda_=1e11, mu_=4e10, alpha_T=1e-5, T_ref=300.0,
               fix_tol=1e-8):

    # читаем сетку и температуру
    mesh   = meshio.read(mesh_xdmf)
    points = mesh.points
    cells  = mesh.cells_dict["tetra"]
    temp   = meshio.read(temp_vtu).point_data["T"]

    E, nu = compute_E_nu(lambda_, mu_)

    with open(inp_path, "w") as f:
        f.write("** cross-verification INP\n")
        f.write("*Heading\n")
        f.write("Calculated by export_to_calculix.py\n\n")

        # 1) Узлы
        f.write("*Node, nset=ALLNODES\n")
        for i, p in enumerate(points, start=1):
            f.write(f"{i}, {p[0]}, {p[1]}, {p[2]}\n")
        f.write("\n")

        # 2) Элементы
        f.write("*Element, type=C3D4, elset=ALL_ELEMENTS\n")
        for eid, elem in enumerate(cells, start=1):
            n = elem + 1
            f.write(f"{eid}, {n[0]}, {n[1]}, {n[2]}, {n[3]}\n")
        f.write("\n")

        # 3) Материал
        f.write("*Solid Section, elset=ALL_ELEMENTS, material=MAT1\n")
        f.write("*Material, name=MAT1\n")
        f.write("*Elastic\n")
        f.write(f"{E:.6e}, {nu:.6f}\n")
        f.write("*Expansion\n")
        f.write(f"{alpha_T}\n\n")

        # 4) Граничные условия
        f.write("*Nset, nset=FIXED_NODES\n")
        for i, p in enumerate(points, start=1):
            if abs(p[2]) < fix_tol:
                f.write(f"{i},\n")
        f.write("*Boundary\n")
        f.write("FIXED_NODES, 1, 3, 0.0\n\n")

        # 5) Начальная температура (реф.)
        f.write("*Initial Conditions, TYPE=TEMPERATURE\n")
        for i in range(1, len(points)+1):
            f.write(f"{i}, {T_ref:.6f}\n")
        f.write("\n")

        # 6) Статический шаг с тепловой нагрузкой
        f.write("*Step, name=ThermoElastic, nlgeom=NO\n")
        f.write("*Static\n\n")

        # 7) ВЫВОД В FRD
        f.write("*NODE FILE, OUTPUT=FRD\n")
        f.write("U\n\n")                # обязательно U, иначе расчёт в FRD не попадёт
        f.write("*EL FILE,   OUTPUT=FRD\n\n")

        # 8) Собственно поле температуры
        f.write("*Temperature\n")
        for i, Tval in enumerate(temp, start=1):
            f.write(f"{i}, {Tval:.6f}\n")
        f.write("*End Step\n")

    print(f"Wrote INP file: {inp_path}")

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--mesh",  required=True, help="VTK/mesh_tet.xdmf")
    p.add_argument("--temps", required=True, help="VTK/simulation_<t>/internal.vtu")
    p.add_argument("--out",   default="crosscheck.inp")
    args = p.parse_args()
    export_inp(Path(args.mesh), Path(args.temps), Path(args.out))
