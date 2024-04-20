#!/bin/bash
#SBATCH --partition=GTX780
#SBATCH --job-name=mdcode_rdf
#SBATCH --comment="MD code run"

gamma=(0.001 0.01 0.1 1.0 10.0 100.0)
for g in "${gamma[@]}"; do
    ./mdlj -ns 50000 -outf 50000 -thermof 1 > ./2_hw/5_task/thermo_lang/$g.txt -lang -gamma $g -msd -msd_start 10000 -seed 42 -Tl 2.0
done

tau=(0.001 0.01 0.1 1.0 10.0 100.0)
for t in "${tau[@]}"; do
    ./mdlj -ns 50000 -outf 50000 -thermof 1 > ./2_hw/5_task/thermo_ber/$t.txt -ber -tau $t -msd -msd_start 10000 -seed 42 -Tl 2.0
done