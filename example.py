#!/bin/python3

from pathlib import Path
from pylammpsmpi import LammpsLibrary
import numpy as np

MPI_CORES = 1
OMP_THREADS = 4
OUT_DIR = Path('results')
INPUT_FILE = Path('fall700.input.data')
ZERO_LVL = 83.391
LATTICE = 5.43


def set_suffix():
    if OMP_THREADS <= 0:
        lmp.command("package gpu 0")
        lmp.command("suffix gpu")
    else:
        lmp.command(f"package omp {OMP_THREADS}")
        lmp.command("suffix omp")


def new_var(name, value):
    lmp.command(f"variable {name} delete")
    lmp.command(f"variable {name} equal {value}")


lmp = LammpsLibrary(cores=MPI_CORES)
lmp.command(f'log {OUT_DIR / "log.init"}')

set_suffix()

lmp.command("units metal")
lmp.command("dimension 3")
lmp.command("boundary p p m")
lmp.command("atom_style atomic")
lmp.command("atom_modify map yes")
lmp.command(f"read_data {INPUT_FILE}")

for i in range(1):
    run_num = i + 1
    lmp.command(
        f"lattice diamond {LATTICE} orient x 1 0 0 orient y 0 1 0 orient z 0 0 1"
    )

    new_var("step", 1e-3)
    new_var("C60_x", 0)
    new_var("C60_y", 0)
    new_var("C60_z", 15.3 * LATTICE + 20)
    new_var("C60_vel", -np.sqrt(8_000) * 5.174)
    new_var("box_width", 12)
    new_var("box_bottom", -16)
    new_var("Si_top", 15.3)
    new_var("temperature", 700)
    new_var("zero_lvl", ZERO_LVL)
    new_var("Si_lattice", LATTICE)

    lmp.command('molecule C60 "mol.C60"')
    lmp.command(
        "create_atoms 1 single ${C60_x} ${C60_y} ${C60_z} mol C60 1 units box"
    )

    lmp.command("variable Si_fixed equal 'v_box_bottom + 0.5'")
    lmp.command("region Si_fixed block -${box_width} ${box_width} -${box_width} ${box_width} "
                "${box_bottom} ${Si_fixed} units lattice")
    lmp.command("region clusters block -${box_width} ${box_width} -${box_width} ${box_width} 0 INF units lattice")

    lmp.command("pair_style tersoff/zbl")
    lmp.command("pair_coeff * * SiC.tersoff.zbl Si C")
    lmp.command("neighbor 3.0 bin")

    lmp.command("group C60 type 2")
    lmp.command("group Si type 1")
    lmp.command("group Si_fixed region Si_fixed")
    lmp.command("group nve subtract all Si_fixed")

    lmp.command('variable is_sputtered atom "z>v_zero_lvl"')

    lmp.command('fix nve nve nve')
    lmp.command('fix dt all dt/reset 1 $(v_step/10) ${step} 0.1')

    lmp.command(f'dump all all custom 20 {OUT_DIR / "all.dump"} id type x y z')
    lmp.command('velocity C60 set NULL NULL ${C60_vel} sum yes units box')

    lmp.run(200)

    lmp.command('group clusters variable is_sputtered')
    lmp.command('compute clusters clusters cluster/atom 3')
    lmp.command('compute mass clusters property/atom mass')

    lmp.command(
        f'dump clusters clusters custom 1 {OUT_DIR / "clusters.dump"} id x y z vx vy vz type c_clusters')
    lmp.command(f'dump final all custom 1 {OUT_DIR / "final.dump"} id x y z vx vy vz type c_clusters')
    lmp.run(0)

    atom_cluster = lmp.extract_compute("clusters", 1, 1)
    print("finished")
