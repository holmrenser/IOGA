#!/usr/bin/python
import glob
import sys
import json
import wget
import os
import shutil
from subprocess import Popen,PIPE
import gzip
import contextlib

dependencies = {	'samtools':{'source':'http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2',
								'target':'IOGA_samtools.tbz',
								'dir':'samtools-0.1.19',
								'unpack':'tar -xvjf IOGA_samtools.tbz',
								'install':'make',
								'exe':'samtools'},
					'bbmap.sh':{'source':'http://downloads.sourceforge.net/project/bbmap/BBMap_33.89_java7.tar.gz',
								'target':'IOGA_bbmap.tgz',
								'dir':'bbmap',
								'unpack':'tar -xvzf IOGA_bbmap.tgz',
								'install':'',
								'exe':'bbmap.sh'},
					'SOAPdenovo-127mer':{	'source':'http://downloads.sourceforge.net/project/soapdenovo2/SOAPdenovo2/bin/r240/SOAPdenovo2-bin-LINUX-generic-r240.tgz',
											'target':'IOGA_SOAPdenovo.tgz',
											'dir':'SOAPdenovo2-bin-LINUX-generic-r240',
											'unpack':'tar -xvzf IOGA_SOAPdenovo.tgz',
											'install':'',
											'exe':'SOAPdenovo-127mer'},
					'seqtk':{	'source':'https://github.com/lh3/seqtk/archive/master.zip',
								'target':'IOGA_seqtk.zip',
								'dir':'seqtk-master',
								'unpack':'unzip IOGA_seqtk.zip',
								'install':'make',
								'exe':'seqtk'},
					'picard.jar':{	'source':'https://github.com/broadinstitute/picard/releases/download/1.124/picard-tools-1.124.zip',
									'target':'IOGA_picard.zip',
									'dir':'picard-tools-1.124',
									'unpack':'unzip IOGA_picard.zip',
									'install':'chmod +x picard.jar',
									'exe':'picard.jar'},
					'ALE':{	'source':'https://github.com/sc932/ALE/archive/master.zip',
							'target':'IOGA_ALE.zip',
							'dir':'ALE-master',
							'unpack':'unzip IOGA_ALE.zip',
							'install':'make',
							'exe':'src/ALE'},
					'spades.py':{	'source':'http://spades.bioinf.spbau.ru/release3.1.1/SPAdes-3.1.1-Linux.tar.gz',
									'target':'IOGA_SPAdes.tgz',
									'dir':'SPAdes-3.1.1-Linux',
									'unpack':'tar -xvzf IOGA_SPAdes.tgz',
									'install':'',
									'exe':'bin/spades.py'}
					}

def main():
	"""
	Download dependencies for IOGA.py
	"""
	config = {}
	try:
		os.mkdir('exe')
	except OSError:
		shutil.rmtree('exe')
		os.mkdir('exe')
	original = os.getcwd()
	os.chdir('exe')
	curdir = os.getcwd()
	for dep in dependencies:
		print dep
		source = dependencies[dep]['source']
		target = dependencies[dep]['target']
		directory = dependencies[dep]['dir']
		unpack = dependencies[dep]['unpack']
		install = dependencies[dep]['install']
		exe = os.path.abspath('{0}/{1}'.format(directory,dependencies[dep]['exe']))
		wget.download(source,out=target)
		p = Popen(unpack.split(),stdout=PIPE,stderr=PIPE)
		out,err = p.communicate()
		if err:
			print err
			quit()
		print
		if install:
			print install
			os.chdir(directory)
			p = Popen(install.split(),stdout=PIPE,stderr=PIPE)
			out,err = p.communicate()
			os.chdir(curdir)
		try:
			print 'trying {0}'.format(exe)
			if dep == 'picard.jar':
				e = 'java -jar {0}'.format(exe)
			else:
				e = exe
			p = Popen(e.split(),stdout=PIPE,stderr=PIPE)
			out,err = p.communicate()
			print 'succes'
		except OSError:
			print '{0} installation failed'.format(dep)
			quit()
		config[dep] = exe
		if dep == 'bbmap.sh':
			config['bbduk.sh'] = exe.strip('bbmap.sh') + 'bbduk.sh' 
	os.chdir(original)
	with contextlib.closing(gzip.open('exe/bbmap/alladapters.fa.gz','wb')) as outfile:
		for f in glob.glob('exe/bbmap/resources/*'):
			with contextlib.closing(gzip.open(f,'rb')) as infile:
				for line in infile:
					outfile.write(line)
	config['adapters'] = os.path.abspath('exe/bbmap/alladapters.fa.gz')

	with open('IOGA_config.json','w') as fh:
		json.dump(config,fh,indent=1)

if __name__ == '__main__':
	main()






