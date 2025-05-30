import argparse
import json
import os
import re
import subprocess
import sys
import signal
from pathlib import Path

# -----------------------
# Функции обновления словарей OpenFOAM
# -----------------------

def update_dictionary_text(text, updates):
    """
    Обновляет числовые значения параметров в тексте словаря OpenFOAM.
    """
    for key, new_val in updates.items():
        pattern = r"(^\s*" + re.escape(key) + r"\s*(?:\[[^\]]*\])?\s+)([^\s;]+)(\s*;)"
        replacement = r"\g<1>" + str(new_val) + r"\g<3>"
        new_text, count = re.subn(pattern, replacement, text, flags=re.MULTILINE)
        if count > 0:
            text = new_text
            print(f"[INFO] Обновлён параметр '{key}' на {new_val} ({count} замен).")
        else:
            print(f"[WARNING] Параметр '{key}' не найден.")
    return text

def update_section_in_transport(text, section_name, updates):
    """
    Обновляет числовые значения в секции transportProperties.
    """
    pattern = r"(" + re.escape(section_name) + r"\s*\{\s*)([^}]*)(\s*\})"
    m = re.search(pattern, text, flags=re.MULTILINE)
    if not m:
        print(f"[WARNING] Раздел '{section_name}' не найден.")
        return text
    header, block, footer = m.groups()
    updated_block = update_dictionary_text(block, updates)
    return text.replace(header + block + footer, header + updated_block + footer)

# -----------------------
# Обёртка запуска команд в Docker
# -----------------------

def run_in_docker(command, mount_dir):
    """
    Выполняет команду внутри контейнера openfoam-docker
    """
    docker_cmd = ["openfoam-docker", f"-dir={mount_dir}", "-c", command]
    print(f"[INFO] Запуск: {command}")
    proc = subprocess.Popen(docker_cmd, start_new_session=True)
    ret = proc.wait()
    return ret

import os
import re
import meshio
import numpy as np
from dolfin import File as DOLFINFile, Mesh, XDMFFile, VectorFunctionSpace, FunctionSpace, Function, \
                  Constant, TrialFunction, TestFunction, Identity, div, sym, grad, inner, \
                  dot, dx, solve

def prepare_mesh(mesh_source_vtu, mesh_xdmf):
    """
    Читает mesh_source_vtu, разбивает все hexahedron → tetra вручную,
    пишет результат в mesh_xdmf, а затем загружает его в DOLFIN Mesh.
    """
    if not os.path.exists(mesh_xdmf):
        m = meshio.read(mesh_source_vtu)

        new_cells = []
        for cell_block in m.cells:
            if hasattr(cell_block, "type"):
                ctype, data = cell_block.type, cell_block.data
            else:
                ctype, data = cell_block

            if ctype == "hexahedron":
                # ручное разбиение каждого гекса в 5 тет
                # для каждой восьмёрки индексов делаем 5 строк в new_tets
                h = data
                # шаблон 5 тет:
                t0 = h[:, [0, 1, 3, 4]]
                t1 = h[:, [1, 2, 3, 6]]
                t2 = h[:, [1, 4, 5, 6]]
                t3 = h[:, [3, 4, 6, 7]]
                t4 = h[:, [1, 3, 4, 6]]
                new_tets = np.vstack([t0, t1, t2, t3, t4])
                new_cells.append(("tetra", new_tets))
            else:
                new_cells.append((ctype, data))

        m_tet = meshio.Mesh(points=m.points, cells=new_cells, point_data=m.point_data)
        meshio.write(mesh_xdmf, m_tet)

    mesh = Mesh()
    with XDMFFile(mesh_xdmf) as xdmf:
        xdmf.read(mesh)
    return mesh


def read_scalar_field(field_vtu, field_name, V):
    """Читает скалярное поле field_name из .vtu и переносит в FunctionSpace V."""
    m = meshio.read(field_vtu)
    data = m.point_data[field_name]
    f = Function(V)
    f.vector()[:] = data
    return f


def run_fenics(vtk_dir) -> float:
    dirs = [d for d in os.listdir(vtk_dir) if os.path.isdir(os.path.join(vtk_dir, d))]
    pat = re.compile(r"^simulation_(\d+(\.\d+)?)$")
    times = [(float(m.group(1)), d) for d in dirs if (m := pat.match(d))]
    if not times:
        print("[ERROR] Нет папок вида simulation_<time> в", vtk_dir)
        return 0

    times.sort(key=lambda x: x[0])
    tval, tdir = times[-1]
    print(f"[FEniCS] Работаем с шагом {tval} (папка {tdir})")

    sim_dir      = os.path.join(vtk_dir, tdir)
    internal_vtu = os.path.join(sim_dir, "internal.vtu")
    mesh_vtu     = os.path.join(vtk_dir, "mesh.vtu")
    mesh_xdmf    = os.path.join(vtk_dir, "mesh_tet.xdmf")

    # выбираем откуда брать сетку
    mesh_source = mesh_vtu if os.path.exists(mesh_vtu) else internal_vtu

    mesh = prepare_mesh(mesh_source, mesh_xdmf)

    dim = mesh.geometry().dim()
    V   = VectorFunctionSpace(mesh, "Lagrange", 1)
    S   = FunctionSpace(mesh, "Lagrange", 1)

    T      = read_scalar_field(internal_vtu, "T",           S)
    alpha1 = read_scalar_field(internal_vtu, "alpha.solid", S)

    alpha_T = Constant(1e-5)
    T_ref   = Constant(300.0)
    I       = Identity(dim)
    f       = -div(alpha_T*(T - T_ref)*I)

    u = TrialFunction(V)
    v = TestFunction(V)
    lambda_ = Constant(1e11)
    mu_     = Constant(4e10)
    def sigma(u):
        return lambda_*div(u)*I + 2*mu_*sym(grad(u))

    a = inner(sigma(u), sym(grad(v))) * dx
    L = dot(f, v) * dx

    Uel = Function(V)
    solve(a == L, Uel)

    out_dir = os.path.join(vtk_dir, "..", "FEM_results")
    os.makedirs(out_dir, exist_ok=True)

    step_i = int(tval)
    xdmf_out = os.path.join(out_dir, f"Uel_{step_i}.xdmf")

    with XDMFFile(xdmf_out) as xf:
        xf.write(Uel, tval)

    vtu_out = os.path.join(out_dir, f"Uel_{step_i}.vtu")

    coords = mesh.coordinates()
    cells  = mesh.cells()

    u_np = Uel.vector().get_local().reshape((-1, mesh.geometry().dim()))

    meshio.write(vtu_out, meshio.Mesh(
        points=coords,
        cells=[("tetra", cells)],
        point_data={"Uel": u_np}
    ))
    print(f"[FEniCS] Результаты сохранены: {xdmf_out} и {vtu_out}")
    return tval


# -----------------------
# Main
# -----------------------

def main():
    parser = argparse.ArgumentParser(description="Клиент OpenFOAM + FEniCS")
    parser.add_argument("--json", required=True, help="Путь к JSON с параметрами")
    parser.add_argument("--path", default="openfoam/solidification-flow/simulation/constant")
    parser.add_argument("--cross",  action="store_true", help="Включить cross-верификацию")
    args = parser.parse_args()

    if not os.path.isfile(args.json):
        print(f"[ERROR] {args.json} не найден")
        sys.exit(1)
    params = json.load(open(args.json))

    thermo_path = os.path.join(args.path, "thermoProperties")
    if os.path.isfile(thermo_path) and "thermoProperties" in params:
        text = open(thermo_path).read()
        text = update_dictionary_text(text, params["thermoProperties"])
        open(thermo_path, "w").write(text)
        print("[INFO] thermoProperties обновлён")

    transport_path = os.path.join(args.path, "transportProperties")
    if os.path.isfile(transport_path) and "transportProperties" in params:
        text = open(transport_path).read()
        for sec in ("solid", "liquid"):
            if sec in params["transportProperties"]:
                text = update_section_in_transport(text, sec, params["transportProperties"][sec])
        open(transport_path, "w").write(text)
        print("[INFO] transportProperties обновлён")

    project_root = os.path.abspath(os.path.join(args.path, "..", ".."))
    print(f"[INFO] Корень проекта: {project_root}")

    # Сборка и запуск OpenFOAM
    run_in_docker("cd solidificationFoam/solvers && ./Allwmake", project_root)
    run_in_docker("cd simulation && ./Allclean && ./Allrun", project_root)
    run_in_docker("cd simulation && reconstructPar && foamToVTK", project_root)

    # FEniCS
    vtk_dir = os.path.join(project_root, "simulation", "VTK")
    print("[INFO] Запуск FEniCS решения...")
    last_step = run_fenics(vtk_dir)

    # Кросс-верификация
    if args.cross:
      project_root = Path(args.path).parent.parent
      vtk_dir      = project_root / "simulation" / "VTK"
      mesh_xdmf = vtk_dir / "mesh_tet.xdmf"
      temp_vtu  = vtk_dir / f"simulation_{int(last_step)}" / "internal.vtu"
      inp_file  = project_root / "simulation" / "crosscheck.inp"

      subprocess.run([
          "python", "export_to_calculix.py",
          "--mesh", str(mesh_xdmf),
          "--temps", str(temp_vtu),
          "--out",   str(inp_file)
      ], check=True)

      # Запускаем CalculiX и конвертацию из .frd в .vtk внутри контейнера
      run_in_docker("cd simulation && ccx crosscheck", project_root)
      run_in_docker("cd simulation && ccx2paraview crosscheck.frd vtk", project_root)
      # Сравниваем полученные результаты
      vtk_file = project_root / "simulation" / "crosscheck.vtk"
      subprocess.run([
          "python", "compare_calculix.py",
          "--fenics", str(project_root / "simulation" / "FEM_results" / f"Uel_{int(last_step)}.vtu"),
          "--ccx",    str(vtk_file)
      ], check=True)
      sys.exit(0)

if __name__ == "__main__":
    main()
