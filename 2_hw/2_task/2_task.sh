#!/bin/bash
#SBATCH --partition=GTX780
#SBATCH --job-name=mdcode_rdf
#SBATCH --comment="MD code run"

./mdlj -ns 5000 -dt 0.001 -rho 0.7 -outf 1 -outdir ./2_hw/2_task/dump_0.7
./mdlj -ns 5000 -dt 0.001 -rho 0.001 -outf 1 -outdir ./2_hw/2_task/dump_0.001