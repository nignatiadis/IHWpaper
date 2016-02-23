  #!/bin/bash
for i in `seq 1 10`;
do
    bsub -M 10000 -n 8 -R "span[hosts=1]"  ./ihw_bonferroni_null_simulation_benchmark.R $i
done

