#!/bin/bash
#SBATCH --partition=GTX780
#SBATCH --job-name=mdcode_rdf
#SBATCH --comment="MD code run"

gamma=(0.001 0.01 0.1 1 10 100)
for g in "${gamma[@]}"; do
    ./mdlj -ns 10000 -outf 10000 -thermof 1 > ./2_hw/5_task/thermo_lang/$g.txt -lang -gamma $g -msd -msd_start 5000 -seed 42 -Tl 2.0
done

lambda=(0.001 0.01 0.1 1.0 10.0)
for l in "${lambda[@]}"; do
    ./mdlj -ns 10000 -outf 10000 -thermof 1 > ./2_hw/5_task/thermo_ber/$l.txt -ber -rescale $l -msd -msd_start 5000 -seed 42 -Tl 2.0
done

./mdlj -ns 10000 -outf 10000 -thermof 1 > ./2_hw/5_task/thermo_none/thermo.txt -msd -msd_start 5000 -seed 42
