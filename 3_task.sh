#!/bin/bash
#SBATCH --partition=GTX780
#SBATCH --job-name=mdcode_rdf
#SBATCH --comment="MD code run"

temp=(1.0 1.5 2.0)
for t in "${temp[@]}"; do
    ./mdlj -ns 10000 -outf 10000 -thermof 1 > ./2_hw/3_task/thermo$t -dt 0.001 -rho 0.7 -T0 $t -rescale 1 -duration 3000 -msd -msd_start 5000 -outdir ./2_hw/3_task/dump
done