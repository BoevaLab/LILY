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
`samtools view -u -q 20 $DATAIN/XXX.bam | samtools rmdup -s - $DATAOUT/$CL.K27ac.Q20.noDup.bam  
samtools view -u -q 20 $DATAIN/YYY.bam | samtools rmdup -s - $DATAOUT/$CL.Input.Q20.noDup.bam`  
2. Call peaks with HMCan: http://www.cbrc.kaust.edu.sa/hmcan/  
By default, the size of genomic regions to perform copy number normalization is 100Kb. In case you have 30M reads in the input, you can decrease this size to 10Kb if you wish. It is mostly important if you expect to have some amplifications in your tumor genomes. In this case, you will need to generate a GC content profile for such windows (use [this tool](http://www.cbrc.kaust.edu.sa/hmcan/GCCount.tar.gz)) or kindly ask Haitham Ashoor to do it for you.  
3. Create GFF files with peaks for ROSE:  
`file=$INDIR/$CL.HK27ac _regions.bed  
SCORE=\`cat PATHTO/selectThreshold_for_HMCan_regions.R | R --slave --args $INDIR/$CL.HK27ac_peaks.narrowPeak`  
echo $SCORE  
perl PATHTO/bed_HMCan2gff.pl -f $file -v -minScore $SCORE -o $OUTDIR/$CL.HK27ac_regions.thresholded.gff`  
4. Index the resulting BAM files with `samtools index`  
5. Run standard ROSE on the indexed BAM files and $CL.HK27ac_regions.thresholded.gff  
`file=$CL.HK27ac_regions.thresholded.gff   
python PATHTOROSE/young_computation-rose-1a9bb86b5464/ROSE_main.py -g hg19 -i $file -r $BAMFOLDER/$CL.K27ac.Q20.noDup.bam -o $OUTDIR -t 2500 -c $BAMFOLDER/$CL.Input.Q20.noDup.bam`  
6. Create bigWig files (.bw) from HMCan density (.wig) files:  
`file=$CL.K27ac.wig  
PATHTOUCSCTOOLS/UCSC_tools/wigToBigWig -clip $file hg19.fa.fai $file.bw`  
7. [Optional]  
You may wish to renormalize densities of all your K27ac samples together to be able to compare them afterwards.  
`cat PATHTO/renormalizeWig.R | R --slave --args K27ac \>\> $MARK.log \>$MARK.log`  
Attention: This normalization does not take into account spike-in data.  
8. Run an alternative to ROSE using the .bw file (renormalized or not); it will create $CL.K27ac.scores.bed where SE are sorted:  
`Rscript PATHTO/getHMCanScore_BEDorGFF.R gff $CL\_12KB_STITCHED_TSS_DISTAL.gff $CL.K27ac.renorm.bw $CL.K27ac.scores.bed`  

**Acknowledgements**

This work was supported by grants from the Institut National de la Sante et de la Recherche Medicale, the Institut Curie, the Ligue Nationale contre le Cancer (Equipe labellisee and CIT program).
