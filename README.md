# LILY
Detection of super-enhancers in cancer data    

---------------------------------------------------------------------------------------------------------------------------

LILY is a pipeline for detection of super-enhancers using H3K27ac ChIP-seq data, which includes explicit correction for copy number variation inherent to cancer samples. The pipeline is based on the ROSE algorithm originally developed by [the Young lab](http://younglab.wi.mit.edu/super_enhancer_code.html). 

---------------------------------------------------------------------------------------------------------------------------

**To cite:**

Heterogeneity of neuroblastoma cell identity revealed by transcriptional circuitries. Boeva et al., to be published.

**Documentation**

Here is step-by-step instuctions how to detect super-enhancers in cancer data with a correction for possible copy number abberations. I would advise you to add and use full paths to all files below. You may need to download some tools and R packages (e.g. the R package "rtracklayer").
1. Remove low quality reads (Q\<20 with `samtools view -q`) and remove duplicates from the BAM files ("samtools rmdup")  
`samtools view -u -q 20 $DATAIN/XXX.bam | samtools rmdup -s - $DATAOUT/$CL.K27ac.Q20.noDup.bam`  
`samtools view -u -q 20 $DATAIN/YYY.bam | samtools rmdup -s - $DATAOUT/$CL.Input.Q20.noDup.bam`  
2. Call peaks with HMCan: http://www.cbrc.kaust.edu.sa/hmcan/  
By default, the size of genomic regions to perform copy number normalization is 100Kb. In case you have 30M reads in the input, you can decrease this size to 10Kb if you wish. It is mostly important if you expect to have some amplifications in your tumor genomes. In this case, you will need to generate a GC content profile for such windows (use [this tool](http://www.cbrc.kaust.edu.sa/hmcan/GCCount.tar.gz)) or kindly ask Haitham Ashoor to do it for you.  
3. [Optional: bigWig files can be created automatically later by `runLILY.R`] 
Create bigWig files (.bw) from HMCan density (.wig) files:  
`file=$CL.K27ac.wig`   
`PATHTOUCSCTOOLS/UCSC_tools/wigToBigWig -clip $file hg19.fa.fai $file.bw`  
4. Run LILY; it will create $CL.K27ac.scores.bed where SE are sorted according to the strength:  
`cat PATHTO/runLILY.R | R --slave --args $CL.K27ac $OUTDIR 12500 2500 hg19_refseq.ucsc hg19.fa.fai`  
where `$CL.K27ac` is the path to the HMCan output for sample `$CL.K27ac`, `12500` is stitching distance, `2500` is distance around gene TSS to annotate promoters, `hg19_refseq.ucsc` is file with transcriptome information (which can be found at https://github.com/linlabbcm/rose2/tree/master/rose2/annotation) and `hg19.fa.fai` is a file with chromosome lengths to create bigwig files if they were not previousely created 
5. [Optional]  
You may wish to renormalize densities of all your K27ac samples together to be able to compare them afterwards.  
`cat PATHTO/renormalizeWig.R | R --slave --args HMCANoutputDir K27ac \>\> $MARK.log \>$MARK.log`  
Attention: This normalization does not take into account spike-in data.  

**Acknowledgements**

This work was supported by grants from the Institut National de la Sante et de la Recherche Medicale, the Institut Curie, the Ligue Nationale contre le Cancer (Equipe labellisee and CIT program).
