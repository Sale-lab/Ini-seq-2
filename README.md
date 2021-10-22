# Ini-seq-2
Ini-seq 2 origin caller
We devised a custom script in R (calling subroutines in Bedtools and Awk) to call origins based on the reads in the HL and LL fractions and attribute to them an efficiency score. Briefly, we started by counting sequencing reads in 100 bp windows for both HL and LL fractions. Read counts were normalised to the total number of reads in each sequencing library. We then kept only those windows with a log2(HL/LL) read ratio ≥ 2. The remaining windows that were ≤ 500 bp from each other were merged. Only domains ≥ 200bp were retained. These domains we now call islands. We recounted the sequencing tags for all islands and normalised the counts to total number of reads then calculated the efficiency score for each island (HL/(HL+LL)). For an island to be called an origin we set a cut off at ≥ 0.8 (corresponding to a log2  threshold ≥ 2). We finally computed the Z-score for each origin and used this to divide the origins into three equally sized groups, which we termed high, medium and low efficiency. 

