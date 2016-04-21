  #!/bin/bash
for i in `seq 1 10`;
do
    bsub -M 10000 -n 12 -R "span[hosts=1]"  ./ihw_bonferroni_wasserman_normal_simulations_benchmark_E3.R $i
done

