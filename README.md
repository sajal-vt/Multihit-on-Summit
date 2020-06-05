# Scaling out a Combinatorial Algorithm for Discovering Carcinogenic Gene Combinations #

We have designed a combinatorial algorithm to identify carcinogenic gene combinations for 31 cancer types from TCGA mutation data.
Since the runtime for this algorithm grows exponentially with the length of the combination

### Source Code ###

* compute_4hits.cu
* a compiler script (compile_script.sh) is included to compile the code to an executable
* ./compile_script.sh compute_4hits.cu compute_4hits.o

### Computing 4-hit combinations ###

* bsub ACC-100.lsf
* jsrun -n100 -a1 -c1 -g6 --bind=proportional-packed:1 --launch_distribution=packed stdbuf -o0 ./compute_4hits.o $data_dir/ACC.maf2dat.matrix.out.training $data_dir/manifest_normal_normal.txt.training.txt.geneSampleList 0.1 6 ACC-100-schedule.txt 100


### Data ###

* Tumor samples data: data/ACC.maf2dat.matrix.out.training, $data_dir/ACC.maf2dat.matrix.out.test
* Normal samples data: data/manifest_normal_normal.txt.training.txt.geneSampleList, data/manifest_normal_normal.txt.test.txt.geneSampleList

### Results ###

* Results.xlsx has 3 workbooks with following 3 results
* Classification performance
* 314 combinations
* Speedup on 100 Summit nodes