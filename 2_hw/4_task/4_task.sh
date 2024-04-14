#!/bin/bash
#SBATCH --partition=GTX780
#SBATCH --job-name=mdcode_rdf
#SBATCH --comment="MD code run"

./mdlj -ns 5000 -outf 1 -uf -thermof 100000 -dt 0.001 -rho 0.7 -T0 1.5 -rescale 1 -duration 5000 -outdir ./2_hw/4_task/dumpInit
./mdlj -ns 5000 -outf 1 -outdir ./2_hw/4_task/dump0.01 -uf -thermof 100000 -dt 0.01 -icf ./2_hw/4_task/dumpInit/5000.xyz 
./mdlj -ns 50000 -outf 10 -outdir  ./2_hw/4_task/dump0.001 -uf -thermof 100000 -dt 0.001 -icf ./2_hw/4_task/dumpInit/5000.xyz