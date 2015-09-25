#IOGA
====

##Iterative Organellar Genome Assembly

####Author: Rens Holmer

Update 18/08/2015

IOGA was published in The Biological Journal of the Linnean Society on 07/08/2015  
If you use it, please cite:  
> [Bakker et al. 2015][1], Herbarium genomics: plastome sequence assembly from a range of herbarium specimens using an Iterative Organelle Genome Assembly pipeline, Biol. J. Linnean Soc.

====

Typical runtime on 4 threads is ~20minutes.  
IOGA is now written in Python.  
IOGA now uses the BBmap suite to map reads and to do quality-filtering/adapter-trimming.  
IOGA now comes with a script to download and install dependencies: setup_IOGA.py

Dependencies: Python, BioPython, BBmap, SOAPdenovo2, SeqTK, SPAdes.py, ALE, Samtools 0.1.19, Picardtools

###INSTALL:

* run setup_IOGA.py to download dependencies, this creates IOGA_config.json
* run IOGA.py -h

###TODO: 
* BBmap outputs per contig coverage stats, need to use this in determining chloroplast inverted repeats
* A final step that blasts the assembly agains the input reference to filter out contigs with no hits at all is still required
* Random subsampling to counter excessive coverage is not implemented. If your sample has a lot of organellar reads, you probably want to reduce the number of reads to work with. This is generally the case, and it also speeds things up considerably, so you might want to do it anyway. Use ```seqtk sample [reads.fastq] 1000000 > [1million.reads.fastq]```to reduce excessive coverage.




[1]:http://onlinelibrary.wiley.com/doi/10.1111/bij.12642/abstract
