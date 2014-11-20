#!/bin/bash -f

INSTALLED=()

if [ -e exe ] 
	then 
	rm -rf exe
	mkdir exe
else
	mkdir exe
fi
pushd exe &> /dev/null

#BBMAP
if `type bbmap.sh >/dev/null 2>&1`
	then
	type bbmap.sh
	ADAPTERDIR="`type bbmap.sh | awk '{print $3}' | sed 's/bbmap.sh//g'`resources"
else
	wget http://downloads.sourceforge.net/project/bbmap/BBMap_33.89_java7.tar.gz &&
	tar -xvzf BBMap_33.89_java7.tar.gz &&
	rm BBMap_33.89_java7.tar.gz
	INSTALLED+=('exe/bbmap')
	ADAPTERDIR=${PWD}'/bbmap/resources'
fi
if [ -e $ADAPTERDIR/alladapters.fa.gz ]
	then
	rm $ADAPTERDIR/alladapters.fa.gz
fi
if [ -e $ADAPTERDIR/alladapters.fa ]
	then
	rm $ADAPTERDIR/alladapters.fa
fi
#for file in `ls "$ADAPTERDIR"/\*fa.gz`
for file in `find $ADAPTERDIR -name "*fa.gz"`
do
	gunzip $file
	unzipped_file=`echo $file | sed 's/.gz//g'`
	cat $unzipped_file >> $ADAPTERDIR/alladapters.fa
	gzip $unzipped_file
done
gzip $ADAPTERDIR/alladapters.fa
chmod +x $ADAPTERDIR/alladapters.fa.gz
ADAPTERDIR=`echo $ADAPTERDIR | sed 's/\//\\\\\//g'`
sed -i 's/REPLACE_ADAPTERDIR/'$ADAPTERDIR'\/alladapters.fa.gz/g' ../IOGA.py #> ../IOGA.py

#SOAPdenovo
if `type SOAPdenovo-127mer >/dev/null 2>&1`
	then
	type SOAPdenovo-127mer
else
	wget http://downloads.sourceforge.net/project/soapdenovo2/SOAPdenovo2/bin/r240/SOAPdenovo2-bin-LINUX-generic-r240.tgz &&
	tar -xvzf SOAPdenovo2-bin-LINUX-generic-r240.tgz &&
	rm SOAPdenovo2-bin-LINUX-generic-r240.tgz
	INSTALLED+=('exe/SOAPdenovo2-bin-LINUX-generic-r240')
fi

#SeqTK
if `type seqtk >/dev/null 2>&1`
	then
	type seqtk
else
	wget -O seqtk.zip https://github.com/lh3/seqtk/archive/master.zip &&
	unzip seqtk.zip &&
	pushd seqtk-master &> /dev/null &&
	make &&
	popd &> /dev/null &&
	rm seqtk.zip 
	INSTALLED+=('exe/seqtk-master')
fi

#Samtools
if `type samtools >/dev/null 2>&1`
	then
	type samtools
else
	wget http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2 &&
	tar -xvjf samtools-0.1.19.tar.bz2 &&
	pushd samtools-0.1.19 &> /dev/null &&
	make &&
	popd &> /dev/null &&
	rm samtools-0.1.19.tar.bz2
	INSTALLED+=('exe/samtools-0.1.19')
fi

#Picard tools
if `type picard.jar >/dev/null 2>&1`
	then
	type picard.jar
else
	wget https://github.com/broadinstitute/picard/releases/download/1.124/picard-tools-1.124.zip &&
	unzip picard-tools-1.124.zip &&
	pushd picard-tools-1.124 &> /dev/null &&
	chmod +x picard.jar &&
	popd &> /dev/null &&
	rm picard-tools-1.124.zip
	INSTALLED+=('exe/picard-tools-1.124')
fi




#ALE
if `type ALE >/dev/null 2>&1`
	then 
	type ALE
else
	wget -O ALE.zip https://github.com/sc932/ALE/archive/master.zip &&
	unzip ALE.zip &&
	pushd ALE-master &> /dev/null &&
	make &&
	popd &> /dev/null &&
	rm ALE.zip
	INSTALLED+=('exe/ALE-master/src')
fi


#SPAdes.py
if `type spades.py >/dev/null 2>&1`
	then
	type spades.py
else
	wget http://spades.bioinf.spbau.ru/release3.1.1/SPAdes-3.1.1-Linux.tar.gz &&
	tar -xvzf SPAdes-3.1.1-Linux.tar.gz &&
	rm SPAdes-3.1.1-Linux.tar.gz
	INSTALLED+=('exe/SPAdes-3.1.1-Linux/bin')
fi


if [ ${#INSTALLED[@]} -eq 0 ]
	then
	printf "All dependencies were already present\n"
else
	printf "\nAdd the following folders to your PATH:\n"
	for ((i=0;i<${#INSTALLED[@]};++i))
	do
		FILE=${INSTALLED[$i]}
		printf "${PWD}/$FILE\n"
	done
fi

popd &> /dev/null

sed -i 's/setup = 0/setup = 1/g' IOGA.py &&
printf '\nSetup successful\n'


