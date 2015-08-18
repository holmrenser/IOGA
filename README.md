#IOGA
====

##Iterative Organellar Genome Assembly

####Author: Rens Holmer

Update 18/08/2015

IOGA was published in The Biological Journal of the Linnean Society on 07/08/2015  
If you use it, please cite:  
> [Bakker et al. 2015][1], Herbarium genomics: plastome sequence assembly from a range of herbarium specimens using an Iterative Organelle Genome Assembly pipeline, Biol. J. Linnean Soc.

====

IOGA is now written in Python.  
IOGA now uses the BBmap suite to map reads and to do quality-filtering/adapter-trimming.  
IOGA now comes with a script to download dependencies: download_dependencies.sh, untried though.  

Dependencies: Python, BioPython, BBmap, SOAPdenovo2, SeqTK, SPAdes.py, ALE, Samtools 0.1.19, Picardtools

###INSTALL:

* run download_dependencies.sh
* add any added dependencies to the PATH (this is mentioned by the script)
* run IOGA.py -h

###TODO: 
* BBmap outputs per contig coverage stats, need to use this in determining chloroplast inverted repeats
* A final step that blasts the assembly agains the input reference to filter out contigs with no hits at all is still required
* Random subsampling to counter excessive coverage is not implemented, if your sample has a lot of organellar reads, use seqtk sample to reduce excessive coverage




[1]:http://onlinelibrary.wiley.com/doi/10.1111/bij.12642/abstract
