#!/usr/bin/env python

import os,sys

jars = {'SimSeq':'~/jars/SimSeq.jar',
		'SamToFastq':'~/jars/picard/SamToFastq.jar',
		'AddOrReplaceReadGroups':'~/jars/picard/AddOrReplaceReadGroups.jar',
		'CreateSequenceDictionary':'~/jars/picard/CreateSequenceDictionary.jar',
		'GenomeAnalysisTK':'~/jars/gatk/GenomeAnalysisTK.jar'}

java_ram = '1g'
example_dir = '~brant/src/20120723-SimSeq/examples'
insert_size,insert_stdev,r1_len,r2_len = 200,20,96,96

#bad idea--can't calc an error profile from non-reference individual!
#err = '~/code/AFS/err_profile.txt'

def sim_reads(source_reference,read_prefix,read_number,insert_size,insert_stdev,r1_len, \
				r2_len=None,r1_err=os.path.join(example_dir,'hiseq_mito_default_bwa_mapping_mq10_1.txt'), \
				r2_err=os.path.join(example_dir,'hiseq_mito_default_bwa_mapping_mq10_2.txt'),duplicate_probability=0.01):
	samf = '%s-simreads.sam' % read_prefix

	#SimSeq
	if r2_len is None:
		simcmd = 'java -Xmx%s -jar %s -1 %s --error %s --insert_size %s --insert_stdev %s --read_number %s --read_prefix %s --reference %s --duplicate_probability %s --out %s' % (java_ram,jars['SimSeq'],r1_len,r1_err,insert_size,insert_stdev,read_number,os.path.basename(read_prefix),source_reference,duplicate_probability,samf)
		sam2fqcmd = 'java -Xmx%s -jar %s INPUT=%s-simreads.sorted.bam FASTQ=%s_s_1_1_sequence.txt INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=SILENT' % (java_ram,jars['SamToFastq'],read_prefix,read_prefix)
	else:
		simcmd = 'java -Xmx%s -jar %s -1 %s -2 %s --error %s --error2 %s --insert_size %s --insert_stdev %s --read_number %s --read_prefix %s --reference %s --duplicate_probability %s --out %s' % (java_ram,jars['SimSeq'],r1_len,r2_len,r1_err,r2_err,insert_size,insert_stdev,read_number,os.path.basename(read_prefix),source_reference,duplicate_probability,samf)
		sam2fqcmd = 'java -Xmx%s -jar %s INPUT=%s-simreads.sorted.bam FASTQ=%s_s_1_1_sequence.txt SECOND_END_FASTQ=%s_s_1_2_sequence.txt INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=SILENT' % (java_ram,jars['SamToFastq'],read_prefix,read_prefix,read_prefix)

	#print >> sys.stderr, 'run simulation:\n\t%s' % cmd
	os.system(simcmd)

	#generate sam header
	header_dict = os.path.splitext(source_reference)[0]+'.dict'
	if not os.path.exists(header_dict):
		os.system('java -Xmx%s -jar %s R=%s O=%s' % (java_ram,jars['CreateSequenceDictionary'],source_reference,header_dict))

	os.system('cat %s %s | samtools view -bS -T %s - | samtools sort - %s-simreads.sorted' % (header_dict,samf,source_reference,read_prefix))

	os.system(sam2fqcmd)

	return '%s_s_1_1_sequence.txt' % read_prefix, r2_len and '%s_s_1_2_sequence.txt' % read_prefix or None

def bwa_aln(read_prefix,indiv,map_reference,r1_fq,r2_fq=None):

	os.system('bwa index %s' % map_reference)

	#generate sam header
	#header_dict = os.path.splitext(map_reference)[0]+'.dict'
	#if not os.path.exists(header_dict):
	#	os.system('java -Xmx%s -jar %s R=%s O=%s' % (java_ram,jars['CreateSequenceDictionary'],map_reference,header_dict))

	if r2_len is None:
		bwacmd = 'bwa aln %s %s > %s.sai; bwa samse %s %s.sai %s | samtools view -bS - | samtools sort - %s-aln.sorted' % (map_reference,r1_fq,r1_fq,map_reference,r1_fq,r1_fq,read_prefix)
	else:
		bwacmd = 'bwa aln %s %s > %s.sai; bwa aln %s %s > %s.sai; bwa sampe %s %s.sai %s.sai %s %s | samtools view -bS - | samtools sort - %s-aln.sorted' % (map_reference,r1_fq,r1_fq,map_reference,r2_fq,r2_fq,map_reference,r1_fq,r2_fq,r1_fq,r2_fq,read_prefix)

	os.system(bwacmd)

	os.system('java -Xmx%s -jar %s I=%s-aln.sorted.bam O=%s-aln.sorted.rg.bam ID=%s LB=lane1_%s PL=illumina PU=lane1 SM=%s' % (java_ram, jars['AddOrReplaceReadGroups'], read_prefix, read_prefix, read_prefix, read_prefix, indiv))
	os.system('samtools index %s-aln.sorted.rg.bam' % read_prefix)

	return '%s-aln.sorted.rg.bam' % read_prefix


if __name__ == "__main__":
	
	read_prefix, indiv, source_reference, map_reference, read_number, vcf = sys.argv[1:]

	r1_fq, r2_fq = sim_reads(source_reference,read_prefix,read_number,insert_size,insert_stdev,r1_len,r2_len) #,err,err)

	print >> sys.stderr,'\n',r1_fq,r2_fq,'created'
	
	bamf = bwa_aln(read_prefix,indiv,map_reference,r1_fq,r2_fq)
	
	#skip running GATK on each haplo
	#os.system('java -Xmx%s -jar %s -T UnifiedGenotyper -I %s -R %s -o %s -maxAlleles 1 -out_mode EMIT_ALL_SITES' % (java_ram, jars['GenomeAnalysisTK'], bamf, map_reference, vcf))
