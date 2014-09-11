#!/bin/bash
#
#Author: Rens Holmer
#
#Dependencies:
#
##Bowtie2 [*]
##Seqtk [*]
##Trimmomatic [*]
##Quake [*]
##Jellyfish 1.1 [*]
##SOAPdenovo2 [*]
##MUMmer [ ]
##SPAdes3.0 [*]
##MIX [ ]
##ALE [*]
#
#wiggleplotter.py #custom python script to plot coverage.wig REQUIRES matplotlib 
#assemblyplotter.r #custom R script to plot assembly lengths REQUIRES seqinr library
#
#
#note: Dependencies have to be in $PATH, 
#
#

usage() { printf "
	Usage: $0  
	[-1 | --read1 <forward reads>] 
	[-2 | --read2 <reverse reads>]
	[-m | --machinename <sequencing machine ID in fastq files, default='HWI'>]	
	[-p | --prefix <output location and name, default='/curdir/prefix'>] 
	[-r | --reference <fasta formatted reference>] 
	[-s | --insertsize <expected insert size>] 
	[-t | --threads <number of threads>] 
	[--mix <use MIX.py to mix assemblies, default=FALSE>] 
	
	Orientation of inserts:
	[--fr <orientation of reads is forward/reverse (small insert size libraries)>]
	[--rf <orientation of reads is reverse/forward (mate pair libraries)>]
	[--ff <orientation of reads is forward/forward (no idea when this happens)>]\n\n" 1>&2; exit 1; }
	
	
TEMP=`getopt -o 1:2:m:p:r:s:t: --long read1:,read2:,machinename:,prefix:,reference:,insertsize:,thread:,mix,rf,fr,ff \
             -n 'IOGA_1.4' -- "$@"`

eval set -- "$TEMP"
		
READ1=
READ2=
DIRECTION="--fr"
MIX=false
PREFIX=
READ_PREFIX="HWI"
REFERENCE=
INSERTSIZE=
THREADS=16

		
while true; do
    case "$1" in
        -1 | --read1 ) READ1="$2" ; printf "READ1 set to $2\n" ; shift 2 ;;
        -2 | --read2 ) READ2="$2" ; printf "READ2 set to $2\n" ; shift 2 ;;
		-m | --machinename ) READ_PREFIX="$2" ; printf "READ_PREFIX set to $2\n" ; shift 2 ;;
		-p | --prefix ) PREFIX="$2" ; printf "PREFIX set to $2\n" ; shift 2 ;;
		-r | --reference ) REFERENCE="$2" ; printf "REFERENCE set to $2\n" ; shift 2 ;;
        -s | --insertsize ) INSERTSIZE="$2" ; printf "INSERTSIZE set to $2\n" ; shift 2 ;; 
		-t | --threads ) THREADS="$2" ; printf "THREADS set to $2\n" ; shift 2 ;;
		--fr ) DIRECTION="--fr" ; printf "DIRECTION set to fr\n" ; shift ;;
		--rf ) DIRECTION="--rf" ; printf "DIRECTION set to rf\n" ; shift ;;
		--ff ) DIRECTION="--ff" ; printf "DIRECTION set to ff\n" ; shift ;;
		--mix ) MIX=true ; printf "MIX set to true\n" ; shift ;;
        --) shift ; break ;;
		*)
			printf "Incorrect options provided! \n"
            usage
            ;;
    esac
done

printf "
	forward reads = ${READ1}
	reverse reads = ${READ2}
	reference file = ${REFERENCE}
	output to ${PREFIX}
	expected insert size = ${INSERTSIZE}
	machine name = ${READ_PREFIX} \n"

if 	[ -z "${READ1}" ] || \
	[ -z "${READ2}" ] || \
	[ -z "${PREFIX}" ] || \
	[ -z "${REFERENCE}" ] || \
	[ -z "${INSERTSIZE}" ]
then
	printf "ja"
    usage
fi



export MALLOC_CHECK_=3
	
if [ -e $PREFIX ]
then
	printf "$PREFIX already exists! Remove existing folder? [y/n] \n"
	read -n 1 -i "n" remove
	if [ $remove == "y" ]
	then
		rm -rf $PREFIX
		mkdir $PREFIX
	else
		printf "Exiting \n"
		exit
	fi
else
	mkdir $PREFIX
fi

FOLDERCOUNT=`echo $PREFIX | sed 's/\// /g' | wc -w`
COUNT=$(($FOLDERCOUNT + 1))
FILENAME=`echo $PREFIX | cut -d '/' -f $COUNT`
TEMP_FOLDER=`echo $PREFIX | sed 's/\/'$FILENAME'//g'`

if [ -e $READ1 ]
	then
	echo "forward reads found, extracting"
	if [ `file $READ1 | awk '{print $2}'` == gzip ]
		then
		gunzip -c $READ1 > $PREFIX/$FILENAME\_1.fastq
	else
		tar xvfzO $READ1 > $PREFIX/$FILENAME\_1.fastq
	fi
else
	echo "forward reads not found!"
	exit
fi

if [ -e $READ2 ]
	then
	echo "reverse reads found, extracting"
	if [ `file $READ2 | awk '{print $2}'` == gzip ]
		then
		gunzip -c $READ2 > $PREFIX/$FILENAME\_2.fastq
	else
		tar xvfzO $READ2 > $PREFIX/$FILENAME\_2.fastq
	fi
else
	echo "reverse reads not found!"
	exit
fi

if [ -e $REFERENCE ]
then
	REFERENCE=$REFERENCE
else
	printf "Reference not found! \n"
	exit
fi

OLDSIZE=0
NEWSIZE=1
i=1
CURDIR=$(pwd)

while [ $NEWSIZE -gt $OLDSIZE ]
do
	OLDSIZE=$NEWSIZE

	printf "\nround $i of Iterative Organellar Genome Assembly\n"
	if [ $i -gt 1 ]
		then
		printf "K$q was the best K-mer in the previous round\n\n"
	fi
	mkdir $PREFIX/$FILENAME.$i
	cd $PREFIX/$FILENAME.$i


	printf "Creating Bowtie2 index\n"
	bowtie2-build -q $REFERENCE $FILENAME.$i.bt2index
	printf "Mapping with Bowtie2\n"
	bowtie2 -q -p $THREADS --very-sensitive-local -x $FILENAME.$i.bt2index $DIRECTION -1 $PREFIX/$FILENAME\_1.fastq -2 $PREFIX/$FILENAME\_2.fastq -S $FILENAME.$i.sam
	printf "Filtering mapped reads\n"
	samtools view -Sh -F 4 $FILENAME.$i.sam > mappedreads.$FILENAME.$i.sam.temp

	if [ $i -gt 1 ]
		then
		#cat $PREFIX/$FILENAME.$((i-1))/mappedreads.$FILENAME.$((i-1)).sam mappedreads.$FILENAME.$i.sam | grep $READ_PREFIX | awk '{print $1}' | sort -u > $FILENAME.$i.names
		printf "Merging files\n"
		MergeSamFiles.jar INPUT= $PREFIX/$FILENAME.$((i-1))/mappedreads.$FILENAME.$((i-1)).sam INPUT= mappedreads.$FILENAME.$i.sam.temp OUTPUT= mappedreads.$FILENAME.$i.sam
		printf "Removing duplicates\n"
		MarkDuplicates.jar INPUT= mappedreads.$FILENAME.$i.sam OUTPUT= mappedreads.$FILENAME.$i.temp METRICS_FILE= $FILENAME.$i.duplicates.log REMOVE_DUPLICATES=true
		mv mappedreads.$FILENAME.$i.temp mappedreads.$FILENAME.$i.sam
		#rm mappedreads.$FILENAME.$i.sam.temp
	else
		#cat mappedreads.$FILENAME.$i.sam | grep $READ_PREFIX | awk '{print $1}' | sort -u > $FILENAME.$i.names
		printf "Removing duplicates\n"
		MarkDuplicates.jar INPUT= mappedreads.$FILENAME.$i.sam OUTPUT= mappedreads.$FILENAME.$i.temp METRICS_FILE= $FILENAME.$i.duplicates.log REMOVE_DUPLICATES=true
		mv mappedreads.$FILENAME.$i.temp mappedreads.$FILENAME.$i.sam
	fi

	#cat $PREFIX/$FILENAME.$i.names | sed -e 's/$/\/1/' > $PREFIX/$FILENAME.$i.names_1
	#cat $PREFIX/$FILENAME.$i.names | sed -e 's/$/\/2/' > $PREFIX/$FILENAME.$i.names_2
	printf "Writing fastq files\n"
	SamToFastq.jar INPUT= mappedreads.$FILENAME.$i.sam FASTQ= $FILENAME\_1_reads$i.fastq SECOND_END_FASTQ= $FILENAME\_2_reads$i.fastq

	#NEWSIZE=`ls -l $PREFIX/$FILENAME.$i/$FILENAME.$i.names | awk '{print $5}'`
	NEWSIZE=`wc -l $FILENAME\_2_reads$i.fastq`

	printf "NEWSIZE=$NEWSIZE \nOLDSIZE=$OLDSIZE\n"

	if [ $NEWSIZE -le $OLDSIZE ]
		then
		cd $CURDIR
		rm -rf $PREFIX/$FILENAME.$i
		lastround=$((i-1)) 
		printf "\nFINISHED after $lastround rounds\n"
		date
		break
	fi

	printf "Writing forward reads\n"
	seqtk subseq $PREFIX/$FILENAME\_1.fastq $FILENAME.$i.names > $FILENAME\_1_reads$i.fastq
	printf "Writing reverse reads\n"
	seqtk subseq $PREFIX/$FILENAME\_2.fastq $FILENAME.$i.names > $FILENAME\_2_reads$i.fastq

	rm $FILENAME.$i.sam

	#TODO: MAKE TRIMMOMATIC ACCESIBLE FROM ANYWHERE
	#TRIMMOMATIC_DIR="/home/biosys/Applications/Trimmomatic-0.32"
	TRIMMOMATIC_DIR=`which trimmomatic-0.32.jar | sed 's/\/trimmomatic-0\.32\.jar//g'`
	trimmomatic-0.32.jar PE -threads $THREADS $FILENAME\_1_reads$i.fastq $FILENAME\_2_reads$i.fastq $FILENAME.$i.FP.fastq $FILENAME.$i.FU.fastq $FILENAME.$i.RP.fastq $FILENAME.$i.RU.fastq ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/alladapter.fa:2:30:10 SLIDINGWINDOW:4:22 MINLEN:36

	#############################
	#checking Trimmomatic output#
	#############################

	if [ -s $PREFIX/$FILENAME.$i/$FILENAME.$i.FP.fastq ]
		then
		FP="$PREFIX/$FILENAME.$i/$FILENAME.$i.FP.fastq"
	else
		printf "FPsize = 0; NO PAIRED READS from Trimmomatic, EXITING!\n"
		date
		exit
	fi

	if [ -s $PREFIX/$FILENAME.$i/$FILENAME.$i.RP.fastq ]
		then
		RP="$PREFIX/$FILENAME.$i/$FILENAME.$i.RP.fastq"
	else
		printf "RPsize = 0; NO PAIRED READS from Trimmomatic, EXITING!\n"
		date
		exit
	fi

	if [ -s $PREFIX/$FILENAME.$i/$FILENAME.$i.FU.fastq ]
		then
		FU="$PREFIX/$FILENAME.$i/$FILENAME.$i.FU.fastq"
	else
		FU=""
	fi

	if [ -s $PREFIX/$FILENAME.$i/$FILENAME.$i.RU.fastq ]
		then
		RU="$PREFIX/$FILENAME.$i/$FILENAME.$i.RU.fastq"
	else
		RU=""
	fi

	##############
	#quake config#
	##############

	cat <<- EOF > quakeconfig.$FILENAME.$i
	$FP $RP
	$FU
	$RU
	EOF

	#######
	#quake#
	#######

	#TODO: CHECK QUAKE FUNCTIONALITY
	quake.py -f quakeconfig.$FILENAME.$i -k 9 -p $THREADS --log --hash_size 1000000

	rm quakeconfig*

	#######################
	#checking quake output#
	#######################

	if [ -s $PREFIX/$FILENAME.$i/$FILENAME.$i.FP.cor.fastq ]
		then
		q1="q1=$PREFIX/$FILENAME.$i/$FILENAME.$i.FP.cor.fastq"
		q2="q2=$PREFIX/$FILENAME.$i/$FILENAME.$i.RP.cor.fastq"
		if [ -s $PREFIX/$FILENAME.$i/$FILENAME.$i.FP.cor_single.fastq ]
			then
			q3="q=$PREFIX/$FILENAME.$i/$FILENAME.$i.FP.cor_single.fastq"
		else
			q3=""
		fi
		if [ -s $PREFIX/$FILENAME.$i/$FILENAME.$i.RP.cor_single.fastq ]
			then
			q4="q=$PREFIX/$FILENAME.$i/$FILENAME.$i.RP.cor_single.fastq"
		else
			q4=""
		fi
		if [ -s $PREFIX/$FILENAME.$i/$FILENAME.$i.FU.cor.fastq ]
			then 
			q5="q=$PREFIX/$FILENAME.$i/$FILENAME.$i.FU.cor.fastq"
		else
			q5=""
		fi
		if [ -s $PREFIX/$FILENAME.$i/$FILENAME.$i.RU.cor.fastq ]
			then
			q6="q=$PREFIX/$FILENAME.$i/$FILENAME.$i.RU.cor.fastq"
		else
			q6=""
		fi
	else
		printf "Quake failed, using trimmomatic output\n"
		if [ -s $PREFIX/$FILENAME.$i/$FILENAME.$i.FP.fastq ]
			then
			q1="q1=$PREFIX/$FILENAME.$i/$FILENAME.$i.FP.fastq"
			q2="q2=$PREFIX/$FILENAME.$i/$FILENAME.$i.RP.fastq"
		fi
		if [ -s $PREFIX/$FILENAME.$i/$FILENAME.$i.FU.fastq ]
			then
			q3="q=$PREFIX/$FILENAME.$i/$FILENAME.$i.FU.fastq"
		else
			q3=""
		fi
		if [ -s $PREFIX/$FILENAME.$i/$FILENAME.$i.RU.fastq ]
			then
			q4="q=$PREFIX/$FILENAME.$i/$FILENAME.$i.RU.fastq"
		else
			q4=""
		fi
		q5=""
		q6=""
	fi
	
	if [ $DIRECTION == '--rf' ]
	then
		reverse_seq='reverse_seq=1'
		pair_num_cutoff='pair_num_cutoff=5'
		map_len='map_len=35'
	else
		reverse_seq='reverse_seq=1'
		pair_num_cutoff='pair_num_cutoff=3'
		map_len='map_len=32'
	fi
	
	#############################
	#creating SOAPdenovo2 config#
	#############################

	cat <<- EOF > soapconfig.$FILENAME.$i
	max_rd_len=100
	[LIB]
	avg_ins=$INSERTSIZE
	$reverse_seq
	asm_flags=3
	rd_len_cutoff=100
	rank=0
	$pair_num_cutoff
	$map_len
	$q1
	$q2
	$q3
	$q4
	$q5
	$q6
	EOF

	######################
	#SOAPdenovo2 assembly#
	######################

	printf "Assembling with SOAPdenovo2\n"
	
	#TODO: EXPLORE SOAP PARAMETERS
	
	for j in {67..99..2}
	do
		date | awk '{print $4}'
		printf "K = $j\n"
		#SOAPdenovo-127mer all -s $PREFIX/$FILENAME.$i/soapconfig.$FILENAME.$i -o $PREFIX/$FILENAME.$i/$FILENAME.$i.$j -K $j -p $THREADS -a 100 -R -F -L 200 -c 0.2 1>$PREFIX/$FILENAME.$i/$FILENAME.$i.$j.err 2>$PREFIX/$FILENAME.$i/$FILENAME.$i.$j.log
		SOAPdenovo-127mer all -s $PREFIX/$FILENAME.$i/soapconfig.$FILENAME.$i -o $PREFIX/$FILENAME.$i/$FILENAME.$i.$j -K $j -p $THREADS 1>$PREFIX/$FILENAME.$i/$FILENAME.$i.$j.err 2>$PREFIX/$FILENAME.$i/$FILENAME.$i.$j.log
	done

	#############################################
	#figuring out the best assembly based on N50#
	#############################################

	cd $PREFIX/$FILENAME.$i
	
	q=`grep -Rin n50 *ics | grep :3[0-9]: | sort -k2n | tail -n 1 | awk 'BEGIN { FS = "." } ; { print $3 }'`

	INSERTSIZE=`grep -i insert  $FILENAME.$i.*.log | grep -i average | tail -n 1 | awk '{print $5}'` #update insert size

	cp $FILENAME.$i.$q.scafSeq bestassembly.$FILENAME.$i.$q.fa
	cp $FILENAME.$i.$q.scafStatistics bestassembly.$FILENAME.$i.$q.stats
	cp $FILENAME.$i.$q.log bestassembly.$FILENAME.$i.$q.log

	REFERENCE=$PREFIX/$FILENAME.$i/bestassembly.$FILENAME.$i.$q.fa

	########################
	#saving some statistics#
	########################
	#
	#grep -Rin n50 *ics | grep :3[0-9]: > stats.$PREFIX/$FILENAME.$i.n50.txt
	#grep -Rin longest_seq *ics | grep :[0-9]: > stats.$PREFIX/$FILENAME.$i.longest_seq.txt
	#grep -Rin size_includen *ics | grep :[0-9]: > stats.$PREFIX/$FILENAME.$i.size.txt
	#grep -Rin scaffold_num *ics > stats.$PREFIX/$FILENAME.$i.scaffold_num.txt
	#grep -Rin singleton_num *ics > stats.$PREFIX/$FILENAME.$i.singleton_num.txt

	#########################################
	#removing unnecessary SOAPdenovo2 output#
	#########################################

	tar cfvz reads.tgz *.fastq

	rm $FILENAME.$i*

	##########################################
	#removing all folders except the last two#
	##########################################

	if [ $i -gt 2 ]
		then
		echo -e "Removing folder $PREFIX/$FILENAME.$((i-2))"
		rm -rf $PREFIX/$FILENAME.$((i-2))
	fi

	#################
	#counter goes up#
	#################

	i=$((i+1))

	###################
	#finish while-loop#
	###################
done

if [ $lastround -gt 0 ]
	then
	
	cd $PREFIX/$FILENAME.$lastround
	
	mkdir SPAdes_standard
	mkdir SPAdes_K_33_55_77
	mkdir SPAdes_K_55_77_99
	mkdir SPAdes_K_89_99
	mkdir SPAdes_K_33_99
	mkdir SPAdes_K_97_99

	tar xvzf reads.tgz

	printf "Assembling final set of reads with SPAdes3.0\n"

	########################
	#checking file presence#
	########################

	if [ -s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.FP.cor.fastq ]
		then
		printf "Using Quake output\n"
		q1="--pe1-1 $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.FP.cor.fastq"
		q2="--pe1-2 $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.RP.cor.fastq"
		if [ -s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.FP.cor_single.fastq ]
			then
			q3="--pe1-s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.FP.cor_single.fastq"
		else
			q3=""
		fi
		if [ -s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.RP.cor_single.fastq ]
			then
			q4="--pe1-s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.RP.cor_single.fastq"
		else
			q4=""
		fi
		if [ -s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.FU.cor.fastq ]
			then 
			q5="--pe1-s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.FU.cor.fastq"
		else
			q5=""
		fi
		if [ -s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.RU.cor.fastq ]
			then
			q6="--pe1-s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.RU.cor.fastq"
		else
			q6=""
		fi
	else
		printf "Quake failed, using trimmomatic output\n"
		if [ -s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.FP.fastq ]
			then
			q1="--pe1-1 $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.FP.fastq"
			q2="--pe1-2 $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.RP.fastq"
		fi
		if [ -s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.FU.fastq ]
			then
			q3="--pe1-s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.FU.fastq"
		fi
		if [ -s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.RU.fastq ]
			then
			q4="--pe1-s $PREFIX/$FILENAME.$lastround/$FILENAME.$lastround.RU.fastq"
		fi
		q5=""
		q6=""
	fi

	echo -e "Running SPAdes standard"
	spades.py $q1 $q2 $q3 $q4 $q5 $q6 --careful -t $THREADS -m 60 -o SPAdes_standard &> SPAdes_standard/log
	echo -e "Running SPAdes K33,55,77"
	spades.py $q1 $q2 $q3 $q4 $q5 $q6 --careful -t $THREADS -m 60 -k 33,55,77 -o SPAdes_K_33_55_77 &> SPAdes_K_33_55_77/log
	echo -e "Running SPAdes K55,77,99"
	spades.py $q1 $q2 $q3 $q4 $q5 $q6 --careful -t $THREADS -m 60 -k 55,77,99 -o SPAdes_K_55_77_99 &> SPAdes_K_55_77_99/log
	echo -e "Running SPAdes K89,99"
	spades.py $q1 $q2 $q3 $q4 $q5 $q6 --careful -t $THREADS -m 60 -k 89,99 -o SPAdes_K_89_99 &> SPAdes_K_89_99/log
	echo -e "Running SPAdes K33,99"
	spades.py $q1 $q2 $q3 $q4 $q5 $q6 --careful -t $THREADS -m 60 -k 33,99 -o SPAdes_K_33_99 &> SPAdes_K_33_99/log
	echo -e "Running SPAdes K97,99"
	spades.py $q1 $q2 $q3 $q4 $q5 $q6 --careful -t $THREADS -m 60 -k 97,99 -o SPAdes_K_97_99 &> SPAdes_K_97_99/log


	mkdir assemblies
	cp -t $PREFIX/$FILENAME.$lastround/assemblies bestassembly*.fa

	for k in `ls -d SPAdes*`
	do
		cd $PREFIX/$FILENAME.$lastround/$k
		mv scaffolds.fasta $k.fa
		cp -t $PREFIX/$FILENAME.$lastround/assemblies $k.fa
	done

	cd $PREFIX/$FILENAME.$lastround/assemblies

	
	if [ $MIX == true ] #TOGGLE MIX
	then
		list=`find *.fa -maxdepth 1 -type f`
		echo -e $list > namelist
		for l in $list
		do
			echo -e $l
			for m in $list
			do
				echo -e $m
				preprocessing.py $l $m -o $l\_$m\_combined.fasta
				nucmer --maxmatch --nosimplify -p $l\_$m\_alignment $l\_$m\_combined.fasta $l\_$m\_combined.fasta
				show-coords -rcl $l\_$m\_alignment.delta > $l\_$m\_alignment.coords
				mkdir $l\_$m
				Mix.py -a $l\_$m\_alignment.coords -c $l\_$m\_combined.fasta -o $l\_$m
				rm *combined*
				rm *alignment*
			done
		done

		for n in `ls -d */ | sed 's/\///g'`
		do
			echo -e $n
			cd $PREFIX/$FILENAME.$lastround/assemblies/$n/Mix_results_A0_C500
			name=`echo -e $n | sed 's/\.fa//g' | sed 's/\//\.fa/g'`
			mv Mix_assembly.fasta $name.fa
			cp -t $PREFIX/$FILENAME.$lastround/assemblies $name.fa
		done
	fi
	
	cd $PREFIX/$FILENAME.$lastround/assemblies

	rm -rf */

	#TODO: BUILD PYTHON ASSEMBLYPLOTTER
	#Rscript --no-save /home/biosys/Applications/IOGA_1.4/assemblyplotter.r
	Rscript --no-save /mnt/nexenta/holme003/CODE/assemblyplotter.r

	cd $PREFIX/$FILENAME.$lastround

	forward=`find $(pwd) -name *FP* | grep -v cor`
	reverse=`find $(pwd) -name *RP* | grep -v cor`

	cd $PREFIX/$FILENAME.$lastround/assemblies

	#TODO: PLOT COVERAGE OF ALL ASSEMBLIES, USE SAMFILE USED FOR ALE
	#TODO: EXPLORE ALE BEHAVIOUR WITH DIFFERENT BOWTIE2 SETTINGS
		
	for o in `find *.fa -maxdepth 0 -type f`
	do
		mv $o $FILENAME.$o
		printf "$o \n Bowtie2 indexing \n"
		bowtie2-build -q $FILENAME.$o index.$FILENAME.$o
		printf "Bowtie2 mapping\n"
		#careful with reporting all and/or local alignments for ALE!
		#bowtie2 -q -a -p $THREADS --very-sensitive-local -x index.$n.$n -1 $forward -2 $reverse -S $n.$n.sam
		bowtie2 -q -p $THREADS --very-sensitive -x index.$FILENAME.$o -1 $forward -2 $reverse -S $FILENAME.$o.sam
		printf "ALE\n"
		ALE $FILENAME.$o.sam $FILENAME.$o $FILENAME.$o.ALE.txt &>> $FILENAME.ALE.log
		samtools view -bS -@ $THREADS $FILENAME.$o.sam > $FILENAME.$o.bam
		samtools sort -@ $THREADS $FILENAME.$o.bam $FILENAME.$o.sorted
		samtools depth $FILENAME.$o.sorted.bam > $FILENAME.$o.coverage.wig
		python /mnt/nexenta/holme003/CODE/wiggleplotter.py $FILENAME.$o.coverage.wig		#COVERAGEPLOTTER
		
		rm *wig
		rm *sam
		rm *bam
		rm *param*
		rm *index*
		rm -rf */

	done
	
	printf "#ALEscores for different assemblies, higher score is better\n#(Higher is closer to zero ;-))\n\n" > $FILENAME.ALESCORE
	grep -rin ale_score *ALE.txt | sort -k3nr >> $FILENAME.ALESCORE
	head -n 1 < $FILENAME.ALESCORE >> bestassemblies
	
	#rm *.txt
	
	cd $PREFIX/$FILENAME.$lastround

	rm *.fastq
	rm *.log
	rm soap*
	rm *.fa
	rm *.sam
	rm *.txt
	rm *.stats
	rm -rf SPA*/

	cd $PREFIX
	
	#TODO: CHECK IF THIS STEP IS NOT TOO GREEDY
	NUMBER_ROUNDS=`find $PREFIX/$FILENAME.* -maxdepth 0 -type d | wc -l` 
	if [ $NUMBER_ROUNDS -gt 1 ]
	then
		find $PREFIX/$FILENAME.* -maxdepth 0 -type d | sort | head -n 1 | xargs rm -rf
	fi
	
	if [ -e $READ1 ] && [ -e $PREFIX/$FILENAME\_1.fastq ]
	then
		printf "$READ1 exists \n removing $PREFIX/$FILENAME\_1.fastq \n\n"
		rm $PREFIX/$FILENAME\_1.fastq
	else
		printf "keeping $PREFIX/$FILENAME\_1.fastq \n\n"
	fi
	if [ -e $READ2 ] && [ -e $PREFIX/$FILENAME\_2.fastq ]
	then
		printf "$READ2 exists \n removing $PREFIX/$FILENAME\_2.fastq \n\n"
		rm $PREFIX/$FILENAME\_2.fastq
	else
		printf "keeping $PREFIX/$FILENAME\_2.fastq \n\n"
	fi
fi


date
