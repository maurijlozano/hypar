#!/usr/bin/env python3
'''
HyPAR (Hypothetical Protein Annotation Reviser), tries to improve genome annotation by finding nearly related genomes, searching for SNP of the frameshift variant, and analyzing if those frameshifts generated hypothetical proteins in the genome of interest.
Version:1.0
Written by Mauricio J. Lozano
UNLP - CONICET - Instituto de Biotecnología y Biología Molecular (IBBM)
'''
VERSION=1.02
REF="None yet"
GITHUB="https://github.com/maurijlozano/HyPAR"

#modules
from Bio import Entrez
from Bio import SeqIO
from os import sys
import argparse
import re, os, subprocess, multiprocessing, shutil
import pandas as pd
import numpy as np
import requests
import gzip
import errno


#argument parsing
def parseArgs():
	'''
	Argument parsing is done.
	'''
	parser = argparse.ArgumentParser(description='HyPAR is a program aimed at improving the annotation of hypothetical proteins on bacterial genomes by searching single nucleotide insertions and deletions on a family of nearly related bacterial strains. Such SNI/D disrupt coding sequences producing small ORFs which might be annotated as Hypothetical proteins instead of larger pseudogenes.')
	parser.add_argument("-G", "--Genus",help="Query genome genus.",action="store", dest="queryGenus", required=True)
	parser.add_argument("-S", "--Specie",help="Query genome specie.",action="store", dest="querySpecies", required=True)
	parser.add_argument("-a", "--AccessionMode",help="Uses accession numbers for the query genome.", dest="acc", action='append', nargs='+', required=False)
	parser.add_argument("-f", "--Files",help="Uses local files for the query genome. Requires a genbank containing all chromosomes and plasmids, and a multifasta file with the sequences...", dest="files", action='append', nargs='+', required=False)
	parser.add_argument("-s", "--Strain",help="Query genome strain.",action="store", dest="queryStrain", required=True)
	parser.add_argument("-e", "--email",help="User email.",action="store", dest="yourEmail", required=True)
	parser.add_argument("-o", "--OutputDir",help="Output folder.",action="store", dest="output", required=False)
	parser.add_argument("-n", "--NumberOfStrains",help="Try to retrieve the genome of n strains to perform the analysis.",action="store", dest="nstrains", required=False)
	#parser.add_argument("-N", "--Nexus", help="Use for input tree in nexus format", action="store_true", dest="NEXUS", required=False)

	args = parser.parse_args()
	if not args.querySpecies or not args.queryGenus:
		sys.stdout.write("You must supply the query genome: -G 'genus' -S 'specie' parameters.\n")
		sys.exit()
	if not args.queryStrain:
		sys.stdout.write("You must supply the query genome strain: -s 'strain' parameter.\n")
		sys.exit()
	if not args.nstrains:
		args.nstrains = 5
	if not args.yourEmail:
		print("An email should be provided as required for NCBItools.\n")
	elif re.match('[^@]+@[^@]+\..{1,3}$',args.yourEmail):
		Entrez.email = args.yourEmail
	else:
		print("Incorrect email address. Email address should be: xxxxx@xxxx.xxx .\n")
		sys.exit()
	return args

#Print presentation

def printSoftName():
	print("\n\n\n")
	print("   ************************************************")
	print("   *****   Hypothetical Protein Annotation Reviser")
	print("   *****   Version: "+str(VERSION))
	print("   *****   Developed by Mauricio J. Lozano")
	print("   *****   github.com/maurijlozano")
	print("   ************************************************")
	print("   Please cite: "+REF)
	print("   Downloaded from: "+GITHUB)
	print("\n\n\n")

#Downloads files from NCBI ftp...
def downLoadNCBIftp(ftplink,outputDir,outFileName):
	httpsLink = re.sub("ftp://","https://",ftplink)
	fileName = outputDir+"/"+os.path.basename(httpsLink)
	try:
		with open(fileName, 'wb') as f:
			f.write(requests.get(httpsLink).content)
		with gzip.open(fileName, 'rb') as f:
			with open(outFileName,'wb') as f_out:
				f_out.write(f.read())
		os.remove(fileName)
	except OSError as e:
		if e.errno == errno.ENOENT:
			# handle file not found error.
			print("Something went wrong...")
			log.write("Something went wrong...\n")
		else:
			# Something else went wrong while trying to run `wget`
			print("Something went wrong, check internet connection...")
			log.write("Something went wrong, check internet connection...\n")
			raise
	pass

def getQuerySeqFromNCBIbyACC(acc):
	#Retrieve GB file for the query genome using accession numbers and NCBItools.
	#NCBItools
	queryID = queryGenus+"_"+querySpecies+"_"+queryStrain
	sequences = acc
	fileName="./"+outputDir+"/"+queryID+".gb"
	fileNameFasta="./"+outputDir+"/"+queryID+".fasta"
	if not os.path.isfile(fileName):
		f=open(fileName, "w+")
		for sequence in sequences:
			print("Downlading "+sequence+" query sequence in genebank format...")
			log.write("Downlading "+sequence+" query sequence in genebank format...\n")
			try:
				handle = Entrez.efetch(db="nucleotide", id=sequence, rettype="gbwithparts", retmode="text")
			except OSError as e:
				if e.errno == errno.ENOENT:
					# handle file not found error.
					print("Something went wrong... Is the accession number correct?")
					log.write("Something went wrong... Is the accession number correct?\n")
				else:
					# Something else went wrong while trying to run `wget`
					print("Something went wrong, check your accession number ...")
					log.write("Something went wrong, check your gaccession number ...\n")
				raise
			f.write(handle.read())
			handle.close()
		f.close()
	if not os.path.isfile(fileNameFasta):
		f=open(fileNameFasta, "w+")
		for sequence in sequences:
			print("Downlading "+sequence+" query sequence in fasta format...")
			log.write("Downlading "+sequence+" query sequence in fasta format...\n")
			try:
				handle = Entrez.efetch(db="nucleotide", id=sequence, rettype="fasta", retmode="text")
			except OSError as e:
				if e.errno == errno.ENOENT:
					# handle file not found error.
					print("Something went wrong... Is the accession number correct?")
					log.write("Something went wrong... Is the accession number correct?\n")
				else:
					# Something else went wrong while trying to run `wget`
					print("Something went wrong, check your accession number ...")
					log.write("Something went wrong, check your gaccession number ...\n")
				raise
			f.write(handle.read())
			handle.close()
		f.close()
	return queryID


def getQuerySeqFromNCBI(queryGenus,querySpecies,queryStrain):
	#Retrieve GB file for the query genome using NCBItools.
	#NCBItools
	queryID = '"'+queryGenus+' '+querySpecies+' '+queryStrain+'"'
	sequences = Entrez.read(Entrez.esearch(db="assembly", RetMax=10, term='"'+queryID+'"[Organism] AND (latest[filter] AND "complete genome"[filter] AND "full genome representation"[filter] AND (all[filter] AND all[filter] NOT anomalous[filter])) AND ("refseq has annotation"[Properties])', idtype="acc" ))
	if len(sequences['IdList']) == 0:
		sys.exit("Query sequence not found on NCBI, try using other strain...")
	summary=Entrez.read(Entrez.esummary(db="assembly", RetMax=10, id=sequences['IdList'][0]))
	ftpLink=summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
	if ftpLink != "":
		ftpLinkGB=ftpLink+"/"+os.path.basename(ftpLink)+"_genomic.gbff.gz"
		ftpLinkFNA=ftpLink+"/"+os.path.basename(ftpLink)+"_genomic.fna.gz"
	else:
		ftpLink=summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
		if ftpLink != "":
			ftpLinkGB=ftpLink+"/"+os.path.basename(ftpLink)+"_genomic.gbff.gz"
			ftpLinkFNA=ftpLink+"/"+os.path.basename(ftpLink)+"_genomic.fna.gz"
		else:
			print("Unable to find query sequence in genebank format...")
			sys.exit()
	queryID = queryGenus+"_"+querySpecies+"_"+queryStrain
	queryID = re.sub(' ','_',queryID)
	fileName="./"+outputDir+"/"+queryID+".gb"
	fileNameFasta="./"+outputDir+"/"+queryID+".fasta"
	if not os.path.isfile(fileName):
		# download ftp..
		print("Downlading query sequence in genebank format...")
		log.write("Downlading query sequence in genebank format...\n")
		downLoadNCBIftp(ftpLinkGB,outputDir,fileName)
	if not os.path.isfile(fileNameFasta):
		# download ftp..
		print("Downlading query sequence in fasta format...")
		log.write("Downlading query sequence in fasta format...\n")
		downLoadNCBIftp(ftpLinkFNA,outputDir,fileNameFasta)
	return queryID

#Extracts CDS from genbank file, given base name...
def genbankToFNA(gbfileID,outputDir):
	f=open(outputDir+"/"+gbfileID+".fna", "w+")
	f2=open(outputDir+"/"+gbfileID+".pos", "w+")
	for rec in SeqIO.parse(outputDir+"/"+gbfileID+".gb", "genbank"):
		if rec.features:
			for feature in rec.features:
				if feature.type == "CDS":
					feature_name=feature.qualifiers["locus_tag"][0]
					if feature.qualifiers.get("old_locus_tag"):
						feature_old_name=feature.qualifiers["old_locus_tag"][0]
					else:
						feature_old_name="ND"
					feature_seq=feature.location.extract(rec).seq
					feature_product=feature.qualifiers['product'][0]
					feature_pos=[str(feature.location.start),str(feature.location.end),str(feature.location.strand)]
					f.write(">" + feature_name + "\n" + str(feature_seq) + "\n")
					f2.write(feature_name + "\t" + feature_product + "\t" + feature_pos[0] + "\t" + feature_pos[1] + "\t" + feature_pos[2] + "\t" + feature_old_name + "\n" )
	f.close()
	f2.close()
	pass

def getSeqsFromNCBI(queryGenus,querySpecies,queryStrain,nstrains):
	#Retrieve fasta files for reference genomes nearly related to the query from genbank using NCBItools. This block retrieves and concatenates the genome refseq sequences for 5 strains if possible...
	# Create a list of all file paths for snippy.
	snippyRef=[]
	#NCBItools
	seqcount = 0
	Id = '"'+queryGenus+' '+querySpecies+'"'
	seqIds = Entrez.read(Entrez.esearch(db="assembly", RetMax=100, term='"'+Id+'"[Organism] AND (latest[filter] AND "complete genome"[filter] AND "refseq has annotation"[Properties]) NOT "'+queryStrain+'"', idtype="acc" ))
	if len(seqIds['IdList']) <= nstrains:
		nstrains = len(seqIds['IdList'])
		print("A total of " + str(len(seqIds['IdList'])) + " alternative " + Id + " strains where found.")
		log.write('A total of ' + str(len(seqIds['IdList'])) + ' alternative ' + Id + ' strains where found.\n')
		log.flush()
	if len(seqIds['IdList']) != 0:
		for ids in seqIds['IdList']:
			if 	seqcount < nstrains:
				#sequences = Entrez.read(Entrez.esearch(db="nucleotide", term='(txid'+ids+'[Organism]) NOT '+queryStrain+'[] NOT wgs[filter] AND refseq[filter]' ))
				summary=Entrez.read(Entrez.esummary(db="assembly", RetMax=10, id=ids))
				ftpLink=summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
				if ftpLink != "":
					ftpLinkGB=ftpLink+"/"+os.path.basename(ftpLink)+"_genomic.gbff.gz"
					ftpLinkFNA=ftpLink+"/"+os.path.basename(ftpLink)+"_genomic.fna.gz"
				else:
					ftpLink=summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
					if ftpLink != "":
						ftpLinkGB=ftpLink+"/"+os.path.basename(ftpLink)+"_genomic.gbff.gz"
						ftpLinkFNA=ftpLink+"/"+os.path.basename(ftpLink)+"_genomic.fna.gz"
					else:
						print("Unable to find reference sequence in genebank format...")
						exit()
				fileName="./"+outputDir+"/"+ids+".gb"
				fileNameFasta="./"+outputDir+"/"+ids+".fasta"
				if not os.path.isfile(fileName):
					# download ftp..
					print("Downlading "+ids+" sequence in genebank format...")
					log.write("Downlading "+ids+" reference sequence in genebank format...\n")
					downLoadNCBIftp(ftpLinkGB,outputDir,fileName)
				if not os.path.isfile(fileNameFasta):
					# download ftp..
					print("Downlading "+ids+" sequence in fasta format...")
					log.write("Downlading "+ids+" sequence in fasta format...\n")
					downLoadNCBIftp(ftpLinkFNA,outputDir,fileNameFasta)
				seqcount+=1
				snippyRef.append(ids)
		return snippyRef
	else:
		print("Unable to find reference sequences...")
		sys.exit()


#Run script
if __name__ == "__main__":
	#Presentation
	printSoftName()
	args = parseArgs()
	#
	#output folder creation
	queryGenus = args.queryGenus
	querySpecies = args.querySpecies
	queryStrain = args.queryStrain
	nstrains = int(args.nstrains)
	if args.output:
		outputDir=args.output
		if not os.path.exists(outputDir):
			os.mkdir(outputDir)
			print("Directory " , outputDir ,  " created.")
		else:
			print("Directory " , outputDir ,  " already exists.")
	else:
		outputDir=queryGenus+"_"+querySpecies
		if not os.path.exists(outputDir):
			os.mkdir(outputDir)
			print("Using default folder, "+outputDir+".")
		else:
			print("Directory " , outputDir ,  " already exists.")

	#Open log files
	log = open("./"+outputDir+'/HyPAR.log', 'a')
	log.write('HyPAR log file.\n')
	log.flush()
	#Get sequences for Snippy analysis...
	#Get query sequence by genus species strain parameters.
	print("Downloading query genome from NCBI...")
	log.write("Downloading query genome from NCBI...\n")
	if args.files:
		files = args.files
		for filewithPath in files[0]:
			if os.path.exists(filewithPath):
				file = os.path.basename(filewithPath)
				newPath="./"+outputDir+'/'+file
				shutil.copyfile(filewithPath,newPath)
		snippyQuery=os.path.splitext(os.path.basename(files[0][0]))[0]
	elif args.acc:
		acc=args.acc
		snippyQuery = getQuerySeqFromNCBIbyACC(acc[0])
	else:
		snippyQuery = getQuerySeqFromNCBI(queryGenus,querySpecies,queryStrain)
	#extract CDS sequences...
	print("Extracting CDS from reference genome genbank...")
	log.write("Extracting CDS from reference genome genbank...\n")
	genbankToFNA(snippyQuery,outputDir)
	#search NCBI for other strains for the same genus and specie
	print("Downloading genomes for other strains from NCBI...")
	log.write("Downloading genomes for other strains from NCBI...\n")
	snippyRef = getSeqsFromNCBI(queryGenus,querySpecies,queryStrain,nstrains)
	log.flush()
	#Run snippy to compare the reference against all Downloaded genomes
	f=open(outputDir+"/frameshiftVariantsGenes.fasta", "w+")
	f2=open(outputDir+"/frameshiftVariantsGenes.tab", "w+")
	print("Exporting CDS for genes with frameshift Variants...")
	log.write("Exporting CDS for genes with frameshift Variants...\n")
	for snippyR in snippyRef:
		try:
			ncpu=multiprocessing.cpu_count()
			#Query vs genomes -> searches for SNPs using all the reference genomes as REF and the query as Query for snippy.
			snippyQueryFile="./"+outputDir+"/"+snippyQuery+".fasta"
			snippyRefFile="./"+outputDir+"/"+snippyR+'.gb'
			snippyOutDir="./"+outputDir+"/"+snippyR
			if not os.path.isfile(snippyOutDir+"/snps.csv"):
				print("Running Snippy for "+ snippyRefFile +"...")
				log.write("Running Snippy for "+ snippyQueryFile +"...\n")
				subprocess.call(['snippy', '--cpus', str(int(ncpu/2)),'--outdir',snippyOutDir, '--ref', snippyRefFile,'--ctgs', snippyQueryFile,'--force'], stdout=log, stderr=log, shell=False)
			#loading tables and geting Locus_tag and variant rettype
			variants=pd.read_csv(snippyOutDir+'/snps.csv')
			frameshiftVariants=variants[variants['EFFECT'].str.contains("frameshift|stop_gained")==True].loc[:,["EFFECT","LOCUS_TAG"]]
			frameshiftVariantsGenes=pd.DataFrame({"locus_tag":frameshiftVariants["LOCUS_TAG"].unique()})
			#Extract frameshiftVariantsGenes and save to fasta format
			for rec in SeqIO.parse(snippyRefFile, "genbank"):
				if rec.features:
					for feature in rec.features:
						if feature.type == "CDS":
							if frameshiftVariantsGenes['locus_tag'].str.match(feature.qualifiers['locus_tag'][0]).any() and not 'hypothetical protein' in feature.qualifiers['product'][0]:
								feature_name=feature.qualifiers["locus_tag"][0]
								feature_seq=feature.location.extract(rec).seq
								feature_prod=feature.qualifiers['product'][0]
								f.write(">" + feature_name + "\n" + str(feature_seq) + "\n")
								f2.write(feature_name + "\t" + feature_prod + "\n")
		except OSError as e:
			if e.errno == errno.ENOENT:
				# handle file not found error.
				print("Something went wrong... Is snippy installed?")
				log.write("Something went wrong... Is snippy installed?\n")
			else:
				# Something else went wrong while trying to run `wget`
				print("Something went wrong, check your genbank files...")
				log.write("Something went wrong, check your genbank files...\n")
				raise
	f.close()
	f2.close()
	log.flush()
	#Blast frameshiftVariantsGenes.fasta against query genome
	try:
		ncpu=multiprocessing.cpu_count()
		dbIn="./"+outputDir+"/"+snippyQuery+".fna"
		dbInMap="./"+outputDir+"/"+snippyQuery+".fasta"
		dbOut="./"+outputDir+"/blastDB/"+snippyQuery+"DB"
		dbOutMap="./"+outputDir+"/blastDB/"+snippyQuery+"DBmap"
		dbOutDir="./"+outputDir+"/blastDB/"
		query=outputDir+"/frameshiftVariantsGenes.fasta"
		blastres=outputDir+"/frameshiftVariantsGenes.res"
		blastmap=outputDir+"/frameshiftVariantsGenes.map"
		blastmapQuery="./"+outputDir+"/"+snippyQuery+".pos"
		queryProduct="./"+outputDir+"/frameshiftVariantsGenes.tab"
		#if not os.path.isfile(outputDir+"/snps.csv"):
		print("Blasting against the Query genome...")
		log.write("Blasting against the Query genome...\n")
		subprocess.call(['makeblastdb', '-in',dbIn, '-out', dbOut, '-dbtype', 'nucl', '-hash_index' ], stdout=log, stderr=log, shell=False)
		subprocess.call(['makeblastdb', '-in',dbInMap, '-out', dbOutMap, '-dbtype', 'nucl', '-hash_index'], stdout=log, stderr=log, shell=False)
		subprocess.call(['blastn','-query', query, '-db', dbOut, '-out', blastres, '-evalue', '0.00000000001', '-num_threads', str(int(ncpu/2)), '-outfmt', '6 qseqid qaccver sseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs', '-max_target_seqs', '5'], stdout=log, stderr=log, shell=False)
		#This blast is to obtain position of candidate genes on the query genome, and for ploting...
		subprocess.call(['blastn','-query', query, '-db', dbOutMap, '-out', blastmap, '-evalue', '0.00000000001', '-num_threads', str(int(ncpu/2)), '-outfmt', '6 qseqid qaccver sseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs', '-max_target_seqs', '1'], stdout=log, stderr=log, shell=False)
		print("Processing blast results...")
		log.write("Processing blast results...\n")
		#Generate table of identified candidate genes...
		headerBlast=['qseqid', 'qaccver', 'sseqid', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs']
		blastResTab=pd.read_csv(blastres,sep='\t',names=headerBlast)
		blastQueryUnique=blastResTab.iloc[:,[2,0,14]]
		blastQueryUnique=blastQueryUnique[blastQueryUnique.iloc[:,2] < 100]
		blastResMapTab=pd.read_csv(blastmap,sep='\t',names=headerBlast)
		#Read product info from frameshiftVariantsGenes.tab...
		frameshiftVariantsGenesProduct=pd.read_csv(queryProduct,sep='\t',names=['qseqid','product'])
		#merge with position table...
		blastResMapTab=pd.merge(blastResMapTab, frameshiftVariantsGenesProduct, on='qseqid', how='left')
		#
		blastResMapTab=pd.DataFrame(blastResMapTab.groupby(['qseqid']).first())
		blastResMapTab=blastResMapTab.iloc[:,[1,9,10,13,14]]
		blastResMapTab.rename({'sseqid': "Query_Replicon_ID"}, axis='columns',inplace=True)
		blastResMapTab=blastResMapTab[blastResMapTab.iloc[:,3] > 90]
		#Reference genes position...
		blastResMapTabReference=pd.read_csv(blastmapQuery,sep='\t',names=['sseqid','product','start','end','strand','Old_Name'])
		#table left joins... there are some combinations...
		mergedTable=pd.merge(blastQueryUnique, blastResMapTabReference, on='sseqid', how='left')
		blastTableMerge=pd.merge(mergedTable, blastResMapTab, on='qseqid', how='left')
		blastTableMerge.columns=['Query_locus_tag' , 'Ref_locus_tag' , 'CDS_ref_cover' , 'Query_Product' , 'Start' , 'End' , 'Strand' ,'Old_Name', 'Query_Replicon_ID' , 'Query_Genome_coordinate_Start' , 'Query_Genome_coordinate_End' , 'Query_Genome_cover','Ref_Product']
		groupedTable=pd.DataFrame(blastTableMerge.groupby(['Query_locus_tag','Ref_locus_tag']).first())
		groupedTable=groupedTable.iloc[:,[1,9,0,2,3,4,5,6,7,8]]
		groupedTable.to_csv(outputDir+'/Candidate_CDS.csv')

	except OSError as e:
		if e.errno == errno.ENOENT:
			# handle file not found error.
			print("Something went wrong... Is Blast+ installed?")
			log.write("Something went wrong... Is Blast+ installed?\n")
		else:
			# Something else went wrong while trying to run `wget`
			print("Something went wrong...")
			log.write("Something went wrong...\n")
			raise
	print("Results are on "+ outputDir + " folder. Look for Candidate_CDS.csv file!")
	log.write("Results are on "+ outputDir + " folder. Look for Candidate_CDS.csv file!\n")
	log.flush()
	log.close()
