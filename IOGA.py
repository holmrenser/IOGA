#!/usr/bin/python
"""
Author: Rens Holmer
Iterative Organellar Genome Assembly (IOGA)
r4
Dependencies: bbmap,bbduk,seqtk,soapdenovo2,spades.py,ALE,BioPython,picardtools,samtools
"""

import matplotlib 
matplotlib.use('Agg')
from matplotlib import pyplot
from pylab import *

import argparse
import datetime
import logging
import os
import re
import shutil
import subprocess
import sys

from Bio import SeqIO

import json

def check_dependencies():
	"""
	Check dependencies, currently done by download_dependencies.sh
	"""
	#horrible oneliner to get absolute path of IOGA folder
	IOGA_path = '/'.join(os.path.abspath(sys.argv[0]).split('/')[:-1]) 
	config_file = '{0}/IOGA_config.json'.format(IOGA_path)
	print config_file
	try:
		with open(config_file,'rU') as ch:
			config = json.load(ch)
	except IOError:
		print 'First run setup_IOGA.py'
		quit()
	missing = 0
	for dep in config:
		if not config[dep]:
			missing += 0
	if missing:
		print 'IOGA_config.json is malformatted'
		quit()
	return config


def main(ref,name,forward,reverse,threads,insertsize,maxrounds,verbosity):
	"""
	IOGA
	"""
	config = check_dependencies()
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
		temp_ref = '{0}/{1}.temp.fasta'.format(folder,prefix)
		with open(temp_ref,'a') as outfile:
			for seq in sequences:
				count += 1
				seq.id = '{0}.contig_{1}'.format(prefix,count)
				seq.name = ''
				seq.description = ''
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

	def run_bbmap(folder,prefix,ref,forward,reverse,threads):
		"""
		Run bbmap to map reads to reference, default, preferred
		"""
		readgroup = '@RG\\tID:{0}\\tSM:{0}\\tPG:IOGA Iterative Organellar Genome Assembly'.format(prefix)
		samfile = '{0}/{1}.sam'.format(folder,prefix)
		log = '{0}/{1}.bbmap.log'.format(folder,prefix)
		mappedreads = '{0}/{1}.mapped.sam'.format(folder,prefix)
		covstats = '{0}/{1}.covstats.txt'.format(folder,prefix)
		basecov = '{0}/{1}.basecov.txt'.format(folder,prefix)
		print '[{0}] BBmap'.format(prefix)
		with open(os.devnull,'w') as fnull:
			if verbosity:
				subprocess.call([config['bbmap.sh'],'ref='+ref,'in='+forward,'in2='+reverse,'threads='+threads,'outm='+mappedreads,'covstats='+covstats,'basecov='+basecov,'-Xmx10G','local=t','keepnames=t'])
			else:
				subprocess.call([config['bbmap.sh'],'ref='+ref,'in='+forward,'in2='+reverse,'threads='+threads,'outm='+mappedreads,'covstats='+covstats,'basecov='+basecov,'-Xmx10G','local=t','keepnames=t'],stdout=fnull,stderr=fnull)				
		plot_coverage(basecov)
		return os.path.abspath(mappedreads)

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
						contigbreaks.append([length,'_'.join(line[:-2])])
		cov_mean = mean(y)
		cov_stdev = std(y)

		pyplot.figure()
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
		pyplot.close()

	def merge_mapped_reads(folder,prefix,current_sam,previous_sam=None):
		"""
		Merge SAM files from current and previous round, sort, remove duplicates
		"""
		with open(os.devnull,'w') as fnull:
			if previous_sam == None:
				os.rename(current_sam,folder + '/'+prefix+'.merged.sam')
			else:
				with open(folder + '/'+prefix+'.merged.sam','w') as temp:
					subprocess.call([config['samtools'],'view','-HS',current_sam],stdout=temp,stderr=fnull)
					for line in subprocess.check_output([config['samtools'],'view','-HS',previous_sam],stderr=fnull).split('\n'):
						if '@PG' not in line:
							if line.strip():
								temp.write(line+'\n')
				with open(folder + '/'+prefix+'.merged.sam','a') as temp:
					subprocess.call([config['samtools'],'view','-S',current_sam],stdout=temp,stderr=fnull)
					subprocess.call([config['samtools'],'view','-S',previous_sam],stdout=temp,stderr=fnull)
			print '[{0}] Sorting SAM files'.format(prefix)
			subprocess.call(['java','-jar',config['picard.jar'],'SortSam','INPUT=',folder + '/merged.sam','OUTPUT=',folder + '/merged.sorted.sam','SORT_ORDER=','coordinate'],stdout=fnull,stderr=fnull)
			print '[{0}] Removing duplicates'.format(prefix)
			subprocess.call(['java','-jar',config['picard.jar'],'MarkDuplicates','INPUT=',folder + '/merged.temp.sam','OUTPUT=',folder + '/' + prefix + '.merged.sam','METRICS_FILE=',folder + '/rmdup.temp.log','REMOVE_DUPLICATES=','true','ASSUME_SORTED=','true'],stdout=fnull,stderr=fnull)
		return os.path.abspath('{0}/{1}.merged.sam'.format(folder,prefix))

	def extract_reads(folder,prefix,samfile,forward,reverse):
		"""
		Extract forward and reverse reads from samfile
		"""
		forward_out = os.path.abspath('{0}/{1}.R1.fastq'.format(folder,prefix)) 
		reverse_out = os.path.abspath('{0}/{1}.R2.fastq'.format(folder,prefix))
		names = []
		print '[{0}] Extracting reads'.format(prefix)
		with open(os.devnull,'w') as fnull:
			for line in subprocess.check_output([config['samtools'],'view','-S',samfile],stderr=fnull).split('\n'):
				if line.strip():
					line = line.split()
					name = re.sub('#$','',line[0]) #remove trailing '#'
					names.append(name)
			names = set(names)
			with open(folder + '/' + prefix + '.names','w') as outfile:
				for name in names:
					outfile.write(name+'\n')
			print '[{0}] Writing forward reads'.format(prefix)
			with open(forward_out,'w') as outfile:
				subprocess.call([config['seqtk'],'subseq',forward,folder + '/' + prefix + '.names'],stderr=fnull,stdout=outfile)
			print '[{0}] Writing reverse reads'.format(prefix)
			with open(reverse_out,'w') as outfile:
				subprocess.call([config['seqtk'],'subseq',reverse,folder + '/' + prefix + '.names'],stderr=fnull,stdout=outfile)
		return forward_out,reverse_out

	def run_bbduk(name,forward,reverse,threads):
		"""
		Use BBduk to quality-trim and remove adapters
		"""
		FP = '{0}.temp.filtered.R1.fastq'.format(name)
		RP = '{0}.temp.filtered.R2.fastq'.format(name) 
		print '[{0}] Quality trimming with BBduk'.format(name)
		with open(os.devnull,'w') as fnull:
			if verbosity:
				subprocess.call([config['bbduk.sh'],'ref='+config['adapters'],'in='+forward,'in2='+reverse,'out='+FP,'out2='+RP,'threads='+threads,'k=25','ktrim=rl','qtrim=t','minlength=32','-Xmx10G'])
			else:
				subprocess.call([config['bbduk.sh'],'ref='+config['adapters'],'in='+forward,'in2='+reverse,'out='+FP,'out2='+RP,'threads='+threads,'k=25','ktrim=rl','qtrim=t','minlength=32','-Xmx10G'],stderr=fnull,stdout=fnull)
		return FP,RP

	def run_soapdenovo(folder,prefix,FP,RP,insertsize,threads):
		"""
		Assemble quality trimmed reads with SOAPdenovo, return assembly with largest scaffold N50
		"""
		n50 = []
		soap_config_file = '{0}/{1}.soap_config'.format(folder,prefix)
		with open(soap_config_file,'w') as config_file:
			config_file.write('max_read_len=250\n'+
				'[LIB]\n'+
				'avg_ins='+insertsize+'\n'+
				'reverse_seq=0\n'
				'asm_flags=3\n'
				'rd_len_cutoff=250\n'
				'rank=0\n'
				'pair_num_cutoff=3\n'
				'map_len=32\n'
				'q1='+FP+'\n'
				'q2='+RP)

		for k in xrange(33,99,4):
			k = str(k)
			print '[{0}] Running SOAPdenovo2 k = {1} -- {2}'.format(prefix,k,datetime.datetime.now().time().isoformat())
			soap_folder = '{0}/{1}.soap_{2}'.format(folder,prefix,k)
			try:
				os.mkdir(soap_folder)
			except:
				print '[{0}] {1} allready exists, overwriting'.format(prefix,soap_folder)
			log = '{0}/{1}.soap_{2}.log'.format(soap_folder,prefix,k)
			err = '{0}/{1}.soap_{2}.err'.format(soap_folder,prefix,k)
			with open(log,'w') as log:
				with open(err,'w') as err:
					subprocess.call([config['SOAPdenovo-127mer'],'all','-s',soap_config_file,'-o',soap_folder + '/' + prefix + '.soap_' + k,'-K',k,'-p',threads,'-d','2'],stdout=log,stderr=err)
			try:
				statfile = '{0}/{1}.soap_{2}.scafStatistics'.format(soap_folder,prefix,k)
				with open(statfile,'rU') as statfile:
					contig = 0
					for line in statfile:
						if 'Contig' in line:
							contig = 1
						if 'N50' in line and contig == 1:
							n50.append([k,line.split()[1],line.split()[2]])
							break
			except:
				print '[{0}] SOAPdenovo2 k = {1} failed'.format(prefix,k)
		best=sorted(n50,key = lambda x: int(x[1]),reverse = True)[0][0]
		source = '{0}/{1}.soap_{2}/{1}.soap_{2}.scafSeq'.format(folder,prefix,best)
		target = '{0}/{1}.soap_{2}.ctg.fasta'.format(folder,prefix,best)
		shutil.move(source,target)
		for k in xrange(33,99,4):
			shutil.rmtree('{0}/{1}.soap_{2}'.format(folder,prefix,k))
		return target,best

	def run_spades(name,forward,reverse,threads):
		"""
		Assembly of final set of organellar reads with SPAdes.py, using several combinations of k.
		"""
		print 'Final assemblies'
		spades_assemblies=[]
		try:
			os.mkdir('{0}.final'.format(name))
		except:
			shutil.rmtree('{0}.final'.format(name))
			os.mkdir('{0}.final'.format(name))

		k_list = [[33],[55],[77],[95],[33,55,77],[55,77,95],[33,95],[33,55,77,95]]
		for k in k_list:
			n = [str(x) for x in k]
			folder = '{0}.final/{0}.spades.{1}'.format(name,'.'.join(n))
			os.mkdir(folder)
			print '[{0}] SPAdes k = {0}'.format('.'.join(n))
			with open(os.devnull,'w') as fnull:
				subprocess.call([config['spades.py'],'-o',folder,'-1',forward,'-2',reverse,'-t',threads,'-k',','.join(n),'--careful'],stderr=fnull,stdout=fnull)
				for i in 'contigs','scaffolds':		
					try:
						source = '{0}/{1}.fasta'.format(folder,i)
						target = '{0}.final/{0}.spades.k.{1}.{2}.fasta'.format(name,'.'.join(n),i)
						shutil.move(source,target)
						spades_assemblies.append(os.path.abspath(target))
					except:
						print '[{0}] SPAdes k{1} {2} failed'.format(name,'.'.join(n),i)
		return os.path.abspath(name + '.final'),spades_assemblies

	def run_ALE(folder,prefix,samfile,assembly):
		"""
		Determine assembly likelyhood using ALE (Note: ALE favours assemblies that stack the chloroplast Inverted Repeat)
		"""
		print '[{0}] ALE'.format(prefix)
		with open(os.devnull,'w') as fnull:
			for line in subprocess.check_output([config['ALE'],samfile,assembly,folder+'/'+prefix+'.ALEOUT'],stderr=fnull).split('\n'):
				line = line.split(':')
				if 'Total ALE Score' in line:
					ALE_score = line[1]
		return [prefix,ALE_score]


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
					iteration -= 1
					break
			if iteration > 2:
				shutil.rmtree(name + '.' + str(iteration-2))
			if len(name.split('/')) == 1:
				prefix = '{0}.{1}'.format(name,iteration)
				try:
					os.mkdir(prefix)
				except:
					print 'Folder {0} exists, overwriting'.format(os.path.abspath(prefix))
					shutil.rmtree(prefix)
					os.mkdir(prefix)
				folder = os.path.abspath(prefix)
			else:
				prefix = name.split('/')[-1] + '.' + str(iteration)
				folder = '/'.join(name.split('/')[:-1])+prefix
				try:
					os.mkdir(folder)
				except:
					print 'Folder {0} exists, overwriting'.format(os.path.abspath(folder))
					shutil.rmtree(folder)
					os.mkdir(folder)
		
			print 'Iteration',str(iteration)
			ref = rename_fasta(folder,prefix,ref)
			samfile = run_bbmap(folder,prefix,ref,forward,reverse,threads)
			if iteration > 1:
				rmdup_merged = merge_mapped_reads(folder,prefix,samfile,rmdup_merged)
			else:
				rmdup_merged = merge_mapped_reads(folder,prefix,samfile)		
			
			FP,RP = extract_reads(folder,prefix,rmdup_merged,forward,reverse)

			newsize = os.path.getsize(FP)		
			print '[{0}] Old size = {1}'.format(prefix,oldsize)
			print '[{0}] New size = {1}'.format(prefix,newsize)
			if newsize <= oldsize:
				break
			oldsize = newsize	

			best,k = run_soapdenovo(folder,prefix,FP,RP,insertsize,threads)
			contig_no,length,N50_index,N50 = stats(best)
			print '[{0}] Best k = {1}'.format(prefix,k)
			print '[{0}] #contigs = {1}'.format(prefix,contig_no)
			print '[{0}] Total length = {1}'.format(prefix,length)
			print '[{0}] N50 index = {1}'.format(prefix,N50_index)
			print '[{0}] N50 = {1}'.format(prefix,N50)

			ref = best

		os.remove(forward)
		os.remove(reverse)

		return best,FP,RP,iteration

	iteration = 0
	start = datetime.datetime.now().time().isoformat()

	source,FP,RP,final_iteration = IOGA_loop(name,ref,forward,reverse,insertsize,threads,maxrounds)

	final_folder,assemblies = run_spades(name,FP,RP,threads)

	target = '{0}/{1}.soap.ctg.fasta'.format(final_folder,name)
	shutil.move(source,target)
	assemblies.append(os.path.abspath(target))
	ALE_score = []
	for assembly in assemblies:
		prefix = assembly.split('/')
		prefix = prefix[len(prefix)-1].split('.fasta')[0]
		try:
			samfile = run_bbmap(final_folder,prefix,assembly,FP,RP,threads)
			scores = run_ALE(final_folder,prefix,samfile,assembly) + [str(x) for x in stats(assembly)]
			ALE_score.append(scores)
		except:
			print '[{0}] {1} not found'.format(prefix,assembly)

	ALE_score = sorted(ALE_score,key = lambda x: float(x[1]),reverse=True)
	ALE_file = '{0}/{1}.statistics'.format(final_folder,name)
	with open(ALE_file,'w') as outfile:
		outfile.write('Filename ALE_score Contig_no Assembly_size N50_index N50\n')
		for score in ALE_score:
			outfile.write('\t'.join(score)+'\n')

	shutil.rmtree('{0}.{1}'.format(name,final_iteration))
	shutil.rmtree('{0}.{1}'.format(name,final_iteration-1))
	
	print '[IOGA] final assembly statistics are in {}'.format(ALE_file)
 	print '[IOGA] Filename ALE_score Contig_no Assembly_size N50_index N50'
	print '[IOGA] {} {} {} {} {} {} has the highest ALE score'.format(*ALE_score[0])
	print '[IOGA] Started  {0}'.format(start)
	print '[IOGA] Finished {0}'.format(datetime.datetime.now().time().isoformat())


if __name__ == '__main__':
	check_dependencies()
	if len(sys.argv) <= 1:
		print "IOGA requires input, try 'python IOGA.py -h'"
		quit()
		
	parser = argparse.ArgumentParser(description='IOGA')
	
	parser.add_argument('--reference','-r', help = 'reference file')
	parser.add_argument('--name','-n', help = 'sample name, default = IOGA_RUN', default = 'IOGA_RUN')
	parser.add_argument('--forward','-1', help = 'forward reads')
	parser.add_argument('--reverse','-2', help = 'reverse reads')
	parser.add_argument('--insertsize','-i', help = 'expected insertsize, default = 500',  default = '500')
	parser.add_argument('--threads','-t', help = 'number of threads, default = 1', default = '1')
	parser.add_argument('--maxrounds','-m', help = 'maximum number of iterations, default = 0', type = int, default = 0)
	parser.add_argument('--verbose','-v', help = 'toggle verbosity', action = 'store_true', default = False)


	args = parser.parse_args()

	main(args.reference, args.name, args.forward, args.reverse, args.threads, args.insertsize, args.maxrounds, args.verbose)



	

