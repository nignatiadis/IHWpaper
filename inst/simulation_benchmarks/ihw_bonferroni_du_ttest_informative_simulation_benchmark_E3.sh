  #!/bin/bash
for i in `seq 1 20`;
do
    bsub -M 10000 -n 8 -R "span[hosts=1]"  ./ihw_bonferroni_du_ttest_informative_simulation_benchmark_E3.R $i
done

