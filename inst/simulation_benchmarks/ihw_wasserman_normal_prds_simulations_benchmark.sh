  #!/bin/bash
for i in `seq 1 60`;
do
    bsub -M 10000 -n 9 -R "span[hosts=1]"  ./ihw_wasserman_normal_prds_simulations_benchmark.R $i
done

