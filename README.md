# IOGA
## Iterative Organellar Genome Assembly

:warning: This project is no longer maintained, and no updates or new versions will be released :warning:


IOGA was used to assemble chloroplast genomes for a range of herbarium samples, and was published in The Biological Journal of the Linnean Society on 07/08/2015. This repository contains the code that was used for the paper, and mainly serves as documentation.

If you use IOGA, please cite:  
> [Bakker et al. 2015][1], Herbarium genomics: plastome sequence assembly from a range of herbarium specimens using an Iterative Organelle Genome Assembly pipeline, Biol. J. Linnean Soc.

---

- Typical runtime on 4 threads is ~20minutes.  
- Written in Python.  
- Uses the BBmap suite to map reads and to do quality-filtering/adapter-trimming. - Comes with a script to download and install dependencies: setup_IOGA.py

Dependencies: Python2, BioPython, BBmap, SOAPdenovo2, SeqTK, SPAdes.py, ALE, Samtools 0.1.19, Picardtools

### INSTALL:

* run setup_IOGA.py to download dependencies, this creates IOGA_config.json
* run IOGA.py -h

### NOTES: 
* BBmap outputs per contig coverage stats, this can be used to determine chloroplast inverted repeats
* A final step that blasts the assembly agains the input reference to filter out contigs with no hits at all is still required
* Random subsampling to counter excessive coverage is not implemented. If your sample has a lot of organellar reads, you probably want to reduce the number of reads to work with. This is generally the case, and it also speeds things up considerably, so you might want to do it anyway. Use ```seqtk sample [reads.fastq] 1000000 > [1million.reads.fastq]```to reduce excessive coverage.




[1]:http://onlinelibrary.wiley.com/doi/10.1111/bij.12642/abstract
