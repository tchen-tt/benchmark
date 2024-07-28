# Benchmark
identity transposon insert methods benchmark

# Workflow  
1. Generate simulated transposon insertion data
2. Generate genome sequencing data with different sequencing depths
3. Use BWA to align the data with different sequencing depths to the reference genome
4. Test different methods to identify transposon insertion sites

# Methods
1. [TEMP2](https://github.com/weng-lab/TEMP2)
2. [ngs_te_mapper2](https://github.com/bergmanlab/ngs_te_mapper2)
3. [RetroSeq](https://github.com/tk2/RetroSeq)
4. [TEi](https://github.com/tchen-tt/TEi)

# Reference  
\[1\] Keane T M, Wong K, Adams D J. RetroSeq: transposable element discovery from next-generation sequencing data[J]. Bioinformatics, 2013, 29(3): 389-390.  
\[2\] Han S, Basting P J, Dias G B, et al. Transposable element profiles reveal cell line identity and loss of heterozygosity in Drosophila cell culture[J]. Genetics, 2021, 219(2): iyab113.  
\[3\] Yu T, Huang X, Dou S, et al. A benchmark and an algorithm for detecting germline transposon insertions and measuring de novo transposon insertion frequencies[J]. Nucleic Acids Research, 2021, 49(8): e44-e44.  
