#!/usr/bin/python
"""
Author: Rens Holmer
Iterative Organellar Genome Assembly (IOGA)
r1
Dependencies: [bowtie2,bwa,bbmap],[trimmomatic,bbduk],seqtk,soapdenovo2,spades.py,ALE,BioPython
"""

#from __future__ import print_function
import matplotlib 
matplotlib.use('Agg')
from matplotlib import pyplot
from pylab import *

import argparse
import datetime
import gzip
import os
import shutil
import subprocess
import sys

from Bio import SeqIO

def check_dependencies():
	"""
	Check dependencies
	"""

def check_files():
	"""
	Check if all files specified are found
	"""

def rename_fasta(folder,prefix,ref):
	"""
	Give all fasta entries a unique name so SAM-headers don't get messed up
	"""
	sequences = SeqIO.parse(open(ref,'rU'),'fasta')
	count = 0
	temp_ref = folder + '/' + prefix + '.temp.fasta' 
	with open(temp_ref,'a') as outfile:
		for seq in sequences:
			count += 1
			seq.id = prefix + '.contig_' +str(count) 
			SeqIO.write(seq,outfile,'fasta')
	return temp_ref

def stats(fasta_file):
	"""
	Determine length, N50 and N50 index of assembly
	"""
	length = []
	N50 = 0
	N50_index = 0
	contig_no = 0
	sequences = SeqIO.parse(open(fasta_file,'rU'),'fasta')
	for seq in sequences:
		contig_no += 1
		length.append(len(seq.seq))
	length = sorted(length,reverse=True)
	for i in length:
		N50 += i
		N50_index += 1
		if N50 >= sum(length) / 2:
			N50 = i
			break
	return contig_no,sum(length),N50_index,N50

def run_bwa(folder,prefix,ref,forward,reverse,threads):
	"""
	Run bwa mem to map reads to reference, optional
	"""
	readgroup = '@RG\\tID:' + prefix + '\\tSM:' + prefix 
	samfile = folder + '/' + prefix + '.sam' 
	log = folder + '/' + prefix + '.bwa.log'  
	mappedreads = folder + '/' + prefix + '.mapped.sam' 
	print '['+prefix+']','BWA index'
	with open(os.devnull,'w') as fnull:
		subprocess.call(['bwa','index','-p',ref,ref],stderr=fnull)
	print '['+prefix+']','BWA mem'
	with open(samfile,'w') as output:
		with open(log,'w') as log:
			subprocess.call(['bwa','mem','-t',threads,'-R',readgroup,'-v','0',ref,forward,reverse],stdout=output,stderr=log)
	print '['+prefix+']','Filtering mapped reads'
	with open(os.devnull,'w') as fnull:
		with open(mappedreads,'w') as mappedreadsfile:
			subprocess.call(['samtools','view','-hS','-F','4',samfile],stdout=mappedreadsfile,stderr=fnull)
	os.remove(samfile)
	return os.path.abspath(mappedreads)

def run_bbmap(folder,prefix,ref,forward,reverse,threads):
	"""
	Run bbmap to map reads to reference, default, preferred
	"""
	readgroup = '@RG\\tID:' + prefix + '\\tSM:' + prefix 
	samfile = folder + '/' + prefix + '.sam' 
	log = folder + '/' + prefix + '.bbmap.log'  
	mappedreads = folder + '/' + prefix + '.mapped.sam' 
	covstats = folder + '/' + prefix + '.covstats.txt'
	basecov = folder + '/' + prefix + '.basecov.txt'
	print '['+prefix+']','BBmap'
	with open(os.devnull,'w') as fnull:
		subprocess.call(['bbmap.sh','ref='+ref,'in='+forward,'in2='+reverse,'threads='+threads,'outm='+mappedreads,'covstats='+covstats,'basecov='+basecov,'-Xmx10G'],stdout=fnull,stderr=fnull)
		#subprocess.call(['bbmap.sh','ref='+ref,'in='+forward,'in2='+reverse,'semiperfectmode=t','threads='+threads,'outm='+mappedreads,'-Xmx10G'],stdout=fnull,stderr=fnull)
	plot_coverage(basecov)
	return os.path.abspath(mappedreads)

def run_bowtie2(folder,prefix,ref,forward,reverse,threads):
	"""
	Run Bowtie2 to map reads to reference, not working yet
	"""
	readgroup = '@RG\\tID:' + prefix + '\\tSM:' + prefix 
	samfile = folder + '/' + prefix + '.sam' 
	log = folder + '/' + prefix + '.bbmap.log' 
	mappedreads = folder + '/' + prefix + '.mapped.sam' 
	print '['+prefix+']','Bowtie2 index'
	with open(os.devnull,'w') as fnull:
		subprocess.call(['bowtie2-build','-q','-c',ref,ref],stderr=fnull)

def plot_coverage(BBmap_coverage):
	y=[]
	contigbreaks=[]
	length=0
	with open(BBmap_coverage,'rU') as infile:
		for line in infile:
			if '#' not in line:
				length+=1
				line = line.split()
				y.append(int(line[-1]))
				if line[1] == '1':
					contigbreaks.append([length,'_'.join(line[0:-3])])

	with open('contigbreaks.txt','w') as testfile:
		for contig in contigbreaks:
			testfile.write(contig)

	cov_mean = mean(y)
	cov_stdev = std(y)
	ylim(0,1.5*max(y))
	ylim(0,1.5*max(y))
	for i in contigbreaks:
		pyplot.plot([i[0],i[0]],[0,1.5*max(y)],color='0.50',linestyle='-',linewidth=0.5)
	pyplot.plot(y,'k-')
	pyplot.plot([0,length],[cov_mean,cov_mean],'r-') #mean
	pyplot.plot([0,length],[cov_mean+cov_stdev,cov_mean+cov_stdev],'r--') #+stdev
	pyplot.plot([0,length],[cov_mean-cov_stdev,cov_mean-cov_stdev],'r--') #-stdev
	suptitle('Mean coverage = ' + str(cov_mean)+' SD = '+str(cov_stdev),fontsize=18)
	figure=gcf()
	figure.set_size_inches(20,10)
	savefig(BBmap_coverage+'.png',dpi=100)

def merge_mapped_reads(folder,prefix,current_sam,previous_sam=None):
	"""
	Merge SAM files from current and previous round, sort, remove duplicates
	"""
	with open(os.devnull,'w') as fnull:
		if previous_sam == None:
			os.rename(current_sam,folder + '/'+prefix+'.merged.sam')
		else:
			with open(folder + '/'+prefix+'.merged.sam','w') as temp:
				subprocess.call(['samtools','view','-HS',current_sam],stdout=temp,stderr=fnull)
				for line in subprocess.check_output(['samtools','view','-HS',previous_sam],stderr=fnull).split('\n'):
					if '@PG' not in line:
						if line.strip():
							temp.write(line+'\n')
			with open(folder + '/'+prefix+'.merged.sam','a') as temp:
				subprocess.call(['samtools','view','-S',current_sam],stdout=temp,stderr=fnull)
				subprocess.call(['samtools','view','-S',previous_sam],stdout=temp,stderr=fnull)
		print '['+prefix+']','Sorting SAM files'
		subprocess.call(['SortSam.jar','INPUT=',folder + '/merged.sam','OUTPUT=',folder + '/merged.sorted.sam','SORT_ORDER=','coordinate'],stdout=fnull,stderr=fnull)
		print '['+prefix+']','Removing duplicates'
		subprocess.call(['MarkDuplicates.jar','INPUT=',folder + '/merged.temp.sam','OUTPUT=',folder + '/' + prefix + '.merged.sam','METRICS_FILE=',folder + '/rmdup.temp.log','REMOVE_DUPLICATES=','true','ASSUME_SORTED=','true'],stdout=fnull,stderr=fnull)
	return os.path.abspath(folder + '/' + prefix + '.merged.sam')	

def extract_reads(folder,prefix,samfile,forward,reverse):
	"""
	Extract forward and reverse reads from samfile
	"""
	forward_out = os.path.abspath(folder + '/' + prefix + '.R1.fastq') 
	reverse_out = os.path.abspath(folder + '/' + prefix + '.R2.fastq')
	names = []
	print '['+prefix+']','Extracting reads'
	with open(os.devnull,'w') as fnull:
		for line in subprocess.check_output(['samtools','view','-S',samfile],stderr=fnull).split('\n'):
			if line.strip():
				line = line.split()
				names.append(line[0])
		names = set(names)
		with open(folder + '/' + prefix + '.names','w') as outfile:
			for name in names:
				outfile.write(name+'\n')
		print '['+prefix+']','Writing forward reads'
		with open(forward_out,'w') as outfile:
			subprocess.call(['seqtk','subseq',forward,folder + '/' + prefix + '.names'],stderr=fnull,stdout=outfile)
		print '['+prefix+']','Writing reverse reads'
		with open(reverse_out,'w') as outfile:
			subprocess.call(['seqtk','subseq',reverse,folder + '/' + prefix + '.names'],stderr=fnull,stdout=outfile)
	return forward_out,reverse_out

def run_bbduk(name,forward,reverse,threads):
	"""
	Use BBduk to quality-trim and remove adapters
	"""
	adapterdir='/mnt/nexenta/holme003/progs_nobackup/bbmap/resources/alladapters.fa.gz'
	FP = name + '.temp.FP.fastq'
	RP = name + '.temp.RP.fastq'
	print '['+name+']','BBduk'
	with open(os.devnull,'w') as fnull:
		subprocess.call(['bbduk.sh','ref='+adapterdir,'in='+forward,'in2='+reverse,'out='+FP,'out2='+RP,'threads='+threads,'k=25','ktrim=rl','qtrim=t','minlength=32','-Xmx10G'],stderr=fnull,stdout=fnull)
	return FP,RP

def run_trimmomatic(folder,prefix,forward,reverse,threads):
	"""
	Use trimmomatic to quality-trim and remove adapters
	"""
	adapterdir='/mnt/nexenta/holme003/progs_nobackup/Trimmomatic-0.32/adapters/'
	FP = os.path.abspath(folder + '/' + prefix + '.FP.fastq') 
	RP = os.path.abspath(folder + '/' + prefix + '.RP.fastq') 
	FU = os.path.abspath(folder + '/' + prefix + '.FU.fastq') 
	RU = os.path.abspath(folder + '/' + prefix + '.RU.fastq') 
	clip = 'ILLUMINACLIP:' + adapterdir + 'alladapter.fa:2:30:10' 
	trimmomatic_command   =['trimmomatic-0.32.jar','PE','-threads',threads,forward,reverse,FP,FU,RP,RU,clip,'SLIDINGWINDOW:4:22','MINLEN:36']
	print '['+prefix+']','Trimmomatic'
	with open(os.devnull,'w') as fnull:
		subprocess.call(trimmomatic_command,stdout=fnull,stderr=fnull)
	return FP,RP#,FU,RU

def run_quake(name,iteration,FP,FU,RP,RU,threads):
	quake_config = name + '.' + iteration + '.quake_config' 
	with open(quake_config,'w') as config_file:
		config_file.write(FP + ' ' + RP + '\n' + FU + '\n' + RU)
	quake_command = ['quake.py','-f',quake_config,'-k','31','-p',threads,'--hash_size','1000000']
	subprocess.call(quake_command)

def run_soapdenovo(folder,prefix,FP,RP,insertsize,threads):
	"""
	Assemble quality trimmed reads with SOAPdenovo, return assembly with largest scaffold N50
	"""
	n50 = []
	soap_config_file = folder + '/' + prefix + '.soap_config'
	with open(soap_config_file,'w') as config_file:
		config_file.write('max_read_len=100\n'+
			'[LIB]\n'+
			'avg_ins='+insertsize+'\n'+
			'reverse_seq=1\n'
			'asm_flags=3\n'
			'rd_len_cutoff=100\n'
			'rank=0\n'
			'pair_num_cutoff=3\n'
			'map_len=32\n'
			'q1='+FP+'\n'
			'q2='+RP)#+'\n'
			#'q='+FU+'\n'
			#'q='+RU)

	for k in xrange(33,99,4):
		k = str(k)
		print '['+prefix+']','Running SOAP, k =',k
		print '['+prefix+']',datetime.datetime.now().time().isoformat()
		soap_folder = folder + '/' + prefix + '.soap_' + k
		try:
			os.mkdir(soap_folder)
		except:
			print '['+prefix+']',soap_folder,'allready exists'
		with open(soap_folder + '/' + prefix + '.soap_' + k + '.log','w') as log:
			with open(soap_folder + '/' + prefix + '.soap_' + k + '.err','w') as err:
				subprocess.call(['SOAPdenovo-127mer','all','-s',soap_config_file,'-o',soap_folder + '/' + prefix + '.soap_' + k,'-K',k,'-p',threads,'-d','2'],stdout=log,stderr=err)
		try:
			with open(soap_folder + '/' + prefix + '.soap_' + k + '.scafStatistics','rU') as statfile:
				contig = 0
				for line in statfile:
					if "Contig" in line:
						contig = 1
					if "N50" in line and contig == 1:
						n50.append([k,line.split()[1],line.split()[2]])
						break
		except:
			print '['+prefix+']','SOAP k =',k,'failed'
	best=sorted(n50,key = lambda x: int(x[1]),reverse = True)[0][0]
	shutil.move(folder + '/' + prefix + '.soap_' + best + '/' + prefix + '.soap_' + best + '.scafSeq',folder + '/' + prefix + '.soap_' + best + '.ctg.fasta')
	for k in xrange(33,99,4):
		shutil.rmtree(folder + '/' + prefix + '.soap_' + str(k))
	return folder + '/' + prefix + '.soap_' + best + '.ctg.fasta',best

def IOGA_loop(name,ref,forward,reverse,insertsize,threads,maxrounds):
	"""
	Main loop, iterate between mapping and assembling until nothing new is found.
	"""
	newsize = 1
	oldsize = 0
	iteration = 0

	forward,reverse = run_bbduk(name,forward,reverse,threads)

	while True:
		iteration += 1
		if maxrounds > 0:
			if iteration == maxrounds:
				break
		if iteration > 2:
			shutil.rmtree(name + '.' + str(iteration-2))
		prefix = name + '.' + str(iteration)
		try:
			os.mkdir(prefix)
		except:
			print 'Folder',prefix,'exists'
		folder = os.path.abspath(prefix)		
		print 'Iteration',str(iteration)
		ref = rename_fasta(folder,prefix,ref)
		samfile = run_bbmap(folder,prefix,ref,forward,reverse,threads)
		#samfile = run_bwa(folder,prefix,ref,forward,reverse,threads)		
		if iteration > 1:
			rmdup_merged = merge_mapped_reads(folder,prefix,samfile,rmdup_merged)
		else:
			rmdup_merged = merge_mapped_reads(folder,prefix,samfile)		
		
		#forward_mapped,reverse_mapped = extract_reads(folder,prefix,rmdup_merged,forward,reverse)		
		FP,RP = extract_reads(folder,prefix,rmdup_merged,forward,reverse)

		newsize = os.path.getsize(FP)		
		print '['+prefix+']','Old size =',oldsize
		print '['+prefix+']','New size =',newsize		
		if newsize <= oldsize:
			break
		oldsize = newsize	

		best,k = run_soapdenovo(folder,prefix,FP,RP,insertsize,threads)
		contig_no,length,N50_index,N50 = stats(best)
		print '['+prefix+']','Best k =',k
		print '['+prefix+']','# contigs',contig_no
		print '['+prefix+']','Total length',length
		print '['+prefix+']','N50 index',N50_index
		print '['+prefix+']','N50',N50

		ref = best

	os.remove(forward)
	os.remove(reverse)

	return best,FP,RP,iteration

def run_spades(name,forward,reverse,threads):
	"""
	Assembly of final set of organellar reads with SPAdes.py, using several combinations of k.
	"""
	print 'Final assemblies'
	spades_assemblies=[]
	os.mkdir(name + '.final')
	k_list = [[33],[55],[77],[95],[33,55,77],[55,77,95],[89,99],[33,95],[33,37],[33,55,77,95],[57,61]]
	for k in k_list:
		n = [str(x) for x in k]
		folder = name + '.final/'+name+'.spades.'+'.'.join(n)
		os.mkdir(folder)
		print '['+name+']','SPAdes k =',','.join(n)
		with open(os.devnull,'w') as fnull:
			subprocess.call(['spades.py','-o',folder,'-1',forward,'-2',reverse,'-t',threads,'-k',','.join(n),'--careful'],stderr=fnull,stdout=fnull)
			try:
				shutil.move(folder+'/contigs.fasta',name+'.final/spades.k'+'.'.join(n)+'.contigs.fasta')
				spades_assemblies.append(os.path.abspath(name+'.final/spades.k'+'.'.join(n)+'.contigs.fasta'))
			except:
				print '['+name+']','Spades k',','.join(n),'contigs failed'
			try:
				shutil.move(folder+'/scaffolds.fasta',name+'.final/spades.k'+'.'.join(n)+'.scaffolds.fasta')
				spades_assemblies.append(os.path.abspath(name+'.final/spades.k'+'.'.join(n)+'.scaffolds.fasta'))
			except:
				print '['+name+']','Spades k',','.join(n),'scaffolds failed'
	return os.path.abspath(name + '.final'),spades_assemblies

def run_ALE(folder,prefix,samfile,assembly):
	"""
	Determine assembly likelyhood using ALE (Note: ALE favours assemblies that stack the chloroplast Inverted Repeat)
	"""
	with open(os.devnull,'w') as fnull:
		for line in subprocess.check_output(['ALE','--metagenome',samfile,assembly,folder+'/'+prefix+'.ALEOUT'],stderr=fnull).split('\n'):
			line = line.split(':')
			if 'Total ALE Score' in line:
				ALE_score = line[1]
	return [prefix,ALE_score]

def main(ref,name,forward,reverse,threads,insertsize,maxrounds):
	"""
	IOGA
	"""
	iteration = 0
	start = datetime.datetime.now().time().isoformat()

	out,FP,RP,final_iteration = IOGA_loop(name,ref,forward,reverse,insertsize,threads,maxrounds)

	final_folder,assemblies = run_spades(name,FP,RP,threads)
	shutil.move(out,final_folder+'/'+name+'.soap.ctg.fasta')
	assemblies.append(os.path.abspath(final_folder+'/'+name+'.soap.ctg.fasta'))
	ALE_score = []
	for assembly in assemblies:
		prefix = assembly.split('/')
		prefix = prefix[len(prefix)-1]
		prefix = name + '.' + prefix.split('.fasta')[0]
		samfile = run_bbmap(final_folder,prefix,assembly,FP,RP,threads)
		scores = run_ALE(final_folder,prefix,samfile,assembly) + list(stats(assembly))
		ALE_score.append(scores)

	ALE_score = sorted(ALE_score,key = lambda x: float(x[1]),reverse=True)

	with open(final_folder+'/'+name+'.statistics','w') as ALE_file:
		ALE_file.write('Filename ALE_score Contig_no Assembly_size N50_index N50\n')
		for score in ALE_score:
			ALE_file.write(score[0]+'\t'+score[1]+'\n')

	shutil.rmtree(name+'.'+str(final_iteration))
	shutil.rmtree(name+'.'+str(final_iteration-1))
	
	print '[IOGA]',ALE_score[0],' has the highest ALE score'
	print '[IOGA] Finished'
	print '[IOGA]',start
	print '[IOGA]',datetime.datetime.now().time().isoformat()


if __name__ == '__main__':
	#threads = "16"
	#insertsize = "230"
	#reference = sys.argv[1]
	#name = sys.argv[2]
	#forward = sys.argv[3]
	#reverse = sys.argv[4]

	if len(sys.argv) <= 1:
		print "IOGA requires input, try 'python IOGA.py -h'"
		quit()

		
	parser = argparse.ArgumentParser(description='IOGA')
	parser.add_argument('--reference','-r',help='reference file')
	parser.add_argument('--name','-n',help='sample name')
	parser.add_argument('--forward','-1',help='forward reads')
	parser.add_argument('--reverse','-2',help='reverse reads')
	parser.add_argument('--insertsize','-i',help='expected insertsize, default = 500',default='500')
	parser.add_argument('--threads','-t',help='number of threads, default = 1',default='1')
	parser.add_argument('--maxrounds','-m',type=int,help='maximum number of iterations, default = 0',default=0)

	args = parser.parse_args()


	main(args.reference,args.name,args.forward,args.reverse,args.threads,args.insertsize,args.maxrounds)



	

