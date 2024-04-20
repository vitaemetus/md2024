#!/bin/bash
#SBATCH --partition=GTX780
#SBATCH --job-name=mdcode_rdf
#SBATCH --comment="MD code run"

# lmp_mpi -in 3_hw/1_task/in.melt -v density 1.1 -v temp 7.5 > 3_hw/1_task/1_thermo_1.1.txt
lmp_mpi -in 3_hw/1_task/in.melt -v density 1.2 -v temp 15.0 > 3_hw/1_task/1_thermo_1.2.txt