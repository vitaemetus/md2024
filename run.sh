#!/bin/bash
#SBATCH --partition=GTX780
#SBATCH --job-name=mdcode_rdf
#SBATCH --comment="MD code run"
./mdlj -ns 5000 -outf 1 -thermof 5001 -dt 0.0001 -rho 1.2 -msd -outdir ./2_hw/3_task/dump_216_1.5