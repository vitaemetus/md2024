#!/bin/bash
#SBATCH --partition=GTX780
#SBATCH --job-name=mdcode_rdf
#SBATCH --comment="MD code run"

./mdlj -ns 20000 -dt 0.001 -rho 0.7 -outf 200 -outdir ./2_hw/1_task/dump_0.7 -N 1000
./mdlj -ns 20000 -dt 0.001 -rho 0.001 -outf 200 -outdir ./2_hw/1_task/dump_0.001 -N 1000