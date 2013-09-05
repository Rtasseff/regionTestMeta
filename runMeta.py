#######
# Python 2.7 command line funciton
# runs the meta analysis of multiple tests assuming a standard format 
# summeries (gets some stats on each test), merges (merges region lists together)
# and compares (stats for tests pairwsie comparision)
# comparisions only done for tests with at least 1 region at max(FDR_set) 
# merge only done for tests with at least 1 region at min(FDR_set) 
# 
#
# to run:
# python2.7 runMeta <input data file path> <output dir> <output test [1 = summary, 2 = 1 + merge, 3 = 2 + compare]> <ignore list, comma seperated list OR '*' for nothing> 
# creates a dir for the output, will not overwrite an existing dir 
# places output files in that dir
# created RAT 20130823
versionTag = 'runMeta_20130823'

import sys
import os
import numpy as np
import scipy.stats as stats

# used as sig cut off in all lists
FDR_merge = .1
# number of tests (or fraction) feature must appear to be merged
strict_merge = 2
# used to look at several values in the summary, max is the cutoff for comparisions
FDR_set = [.05, .1, .2]
# name of suammry file
summaryName = 'summary.tsv'
# name of top merge folder 
mergeDir = 'featureMerges'
# name of comparision folder
compareDir = 'testComparisions'



def fdr_bh(p,alpha=.05):
	"""Performs the Benjamini & Hochberg 1995
	multiple test correction for controlling
	the false discovery rate in familywise 
	analysis.  Tests must be independent or 
	positivly corrilated.
	p	original pvalues, np 1d array
	alpha 	threshold FDR, scalar float
	returns	h, regect or accept, 1d np bool array
	returus	p_adj, adjusted pvalues, 1d np array
	returns	pCrit, the critial p-value cut off
	"""
	m = len(p)
	if m>0:
		sortInd = np.argsort(p)
		pSort = p[sortInd]
		unsortInd = np.argsort(sortInd)
		pAdj = np.zeros(m)*np.nan
		gamma = (np.arange(m)+1)*(alpha/m)
		pTmp = m*pSort/(np.arange(m)+1)
		for i in range(m):
			pAdj[i] = np.min(pTmp[i:])

		pAdjUnsort = pAdj[unsortInd]
		rejSort = pSort<=gamma
		# find the largest value still under threshold
		# note they are sorted
		maxInd = np.sum(rejSort)-1
		if maxInd<0:
			pCrit = 0
		else:
			pCrit = pSort[maxInd]

		h = p<=pCrit
	else:
		pAdjUnsort = np.array([])
		h = np.array([])
		pCrit = np.nan
	
	return h,pAdjUnsort,pCrit


def NGD(nx,ny,nxy,N):
	lx = np.log(nx)
	ly = np.log(ny)
	lxy = np.log(nxy)
	lN = np.log(N)
	return ((np.max([lx,ly]) - lxy)/(lN-np.min([lx,ly])))

def enrich(n_pos,n_draw,total_pos,total_all):
	"""Standard enrichment test usign hypergeometric 
	distribution.
	ie I got n_pos red balls out of n_draw draws from
	a bag with total_pos red balls and total all balls of 
	any color, calculate the probability of drawing n_pos
	red balls at random.
	"""
	
	p = stats.hypergeom.sf(n_pos-1,total_all,total_pos,n_draw)
	return(p)

def region2gene_unq(regionList):
	n = len(regionList)
	geneList = np.zeros(n,dtype=regionList.dtype)
	for i in range(n):
		tmp = regionList[i]
		ind = tmp.find('_')
		geneList[i] = tmp[:ind]

	geneList = np.array(list(set(geneList)),dtype=str)
		
	return geneList

def region2gene_order(regionList):
	n = len(regionList)
	geneList = np.zeros(n,dtype=regionList.dtype)
	count = 0
	for i in range(n):
		tmp = regionList[i]
		ind = tmp.find('_')
		if not np.any(geneList == tmp[:ind]):
			geneList[count] = tmp[:ind]
			count = count+1

	geneList = geneList[:count]
		
	return geneList





def runSumm(data,header,outDir):
	n,m = data.shape
	# we will want to report which other analysis can be done
	# currently we do not want to do analysis on usless tests
	doComp = np.array(np.ones(m),dtype=bool)
	doMerg = np.array(np.ones(m),dtype=bool)

	print 'Running '+str(m)+' summary tests...'
	
	fout = open(outDir+'/'+summaryName,'w')
	# build the header
	fout.write('Test_ID')
	for fdr in FDR_set:
		fout.write('\t#Genes_FDR='+str(fdr))
	fout.write('\t1Perc_p\tmin_p')
	fout.write('\n')
	
	for i in range(m):
		#print 'running test '+str(i)+', '+header[i]
		fout.write(header[i])
		data[data[:,i]>1,i] = 1
		notNAN =  ~np.isnan(data[:,i])
		ind_1Perc = int(.01*np.sum(notNAN))
		# lets look at some FDR rates
		for fdr in FDR_set:
			h,_,_ = fdr_bh(data[notNAN,i],fdr)
			tmp = np.sum(h)
			fout.write('\t'+str(tmp))
			# check if test will be used elsewhere
			if fdr==np.max(FDR_set):
				if tmp==0:doComp[i]=False
			elif fdr==FDR_merge:
				if tmp==0:doMerg[i]=False
			
		
		tmp = np.sort(data[notNAN,i])[ind_1Perc]
		fout.write('\t'+str(tmp))
		tmp = np.min(data[notNAN,i])
		fout.write('\t'+str(tmp))
		fout.write('\n')
	fout.close()
	return doComp, doMerg

def runMerg(data,header,regionID,outDir):
	# create folder for output
	os.makedirs(outDir+'/'+mergeDir)
	n,m = data.shape
	print 'Running region merge on '+str(m)+' tests ...'
	# We are planning on doing a strict merge, a loose merge, and a rank aggrigation 
	# we are going to do the by phenotype and by member, with T -> NB
	# lets set those
	memb = []
	phen = []
	for i in range(m):
		bits = header[i].split('|')
		tmp = bits[2]
		if bits == 'T': tmp='NB'
		memb.append(tmp)
		phen.append(bits[1])
	memb = np.array(memb,dtype=str)
	phen = np.array(phen,dtype=str)

	membList = ['F','M','NB']
	nMemb = 3
	phenList = list(set(phen))
	nPhen = len(phenList)

	print 'Found '+str(nPhen)+' phenotypes to merge ...'

	curIndex = np.array(np.ones(m),dtype=bool)
	
	for curPhen in phenList:
		# reset index to all true
		curIndex_p = curIndex.copy()
		# make spot for the results
		curPhen_dir = curPhen.replace(':','_')
		os.makedirs(outDir+'/'+mergeDir+'/'+curPhen_dir)
		curIndex_p[phen!=curPhen]=False
		# now family members 
		for curMemb in membList:
			# reset to the original index for the phenotype
			curIndex_p_m = curIndex_p.copy()
			dest = outDir+'/'+mergeDir+'/'+curPhen_dir+'/'+curMemb
			os.makedirs(dest)
			curIndex_p_m[memb!=curMemb] = False
			
			# number fo tests in this feature merge
			nTests = np.sum(curIndex_p_m)
			# now we can grab the releveant data 
			curData = data[:,curIndex_p_m]

			if nTests>1:
				# merge time!
				# lets do the very simple list
				# combine all 
				h = np.zeros(n)
				rankScore = np.zeros(n)
				for i in range(nTests):
					# check the FDRs
					notNAN =  ~np.isnan(curData[:,i])
					h_tmp,_,_ = fdr_bh(curData[notNAN,i],FDR_merge)
					h[notNAN] = h[notNAN]+h_tmp
					# get the rank socre
					rankScore[notNAN] = rankScore[notNAN] + 1/(np.argsort(np.argsort(curData[notNAN,i]))+1.)

				# for the loose list 
				regionList = regionID[h>0]
				# save the regions
				np.savetxt(dest+'/regionList_loose.dat',regionList,fmt='%s')
				# get the unique genes 
				geneList = region2gene_unq(regionList)
				np.savetxt(dest+'/geneList_loose.dat',geneList,fmt='%s')
					
				# for the strict list 
				regionList = regionID[h>=strict_merge]
				# save the regions
				np.savetxt(dest+'/regionList_strict.dat',regionList,fmt='%s')
				# get the unique genes 
				geneList = region2gene_unq(regionList)
				np.savetxt(dest+'/geneList_strict.dat',geneList,fmt='%s')
				
				# order based on rank score	
				rankInd = np.argsort(rankScore)[::-1]
				regionList = regionID[rankInd]
				np.savetxt(dest+'/regionList_rankAgg.dat',regionList,fmt='%s')
				geneList = region2gene_order(regionList)
				np.savetxt(dest+'/geneList_rankAgg.dat',geneList,fmt='%s')		
			elif nTests==1:
				h = np.zeros(n)
				curData = curData[:,0]
				# check the FDRs
				notNAN =  ~np.isnan(curData)
				h_tmp,_,_ = fdr_bh(curData[notNAN],FDR_merge)
				h[notNAN] = h_tmp
				# for the loose list 
				regionList = regionID[h>0]
				# save the regions
				np.savetxt(dest+'/regionList_loose.dat',regionList,fmt='%s')
				# get the unique genes 
				geneList = region2gene_unq(regionList)
				np.savetxt(dest+'/geneList_loose.dat',geneList,fmt='%s')
	




def runComp(data,header,outDir):
	n = len(header)
	m = len(data)
	dest = outDir+'/'+compareDir
	os.makedirs(dest)
	comp_R = np.zeros((n,n)) + np.nan
	comp_NGD = np.zeros((n,n)) + np.nan
	comp_R_P = np.zeros((n,n)) + np.nan
	comp_NGD_P = np.zeros((n,n)) + np.nan

	FDR = np.max(FDR_set)

	for i in range(n):
		comp_R[i,i] = 1
		comp_NGD[i,i] = 0
		comp_R_P[i,i] = 0
		comp_NGD_P[i,i] = 0
		hi = np.zeros(m,dtype=int)
		h,_,_ = fdr_bh(data[:,i][~np.isnan(data[:,i])],FDR)
		hi[~np.isnan(data[:,i])] = h
		ni = np.sum(hi==1)
		for j in range(i):
			comp_R[i,j], comp_R_P[i,j] = stats.spearmanr(data[:,i],data[:,j])
			hj = np.zeros(m,dtype=int)
			h,_,_ = fdr_bh(data[:,j][~np.isnan(data[:,j])],FDR)
			hj[~np.isnan(data[:,j])] = h
			nj = np.sum(hj==1)
			nij = np.sum((hj+hi)==2)
			comp_NGD[i,j] = NGD(ni,nj,nij,m)
			comp_NGD_P[i,j] = enrich(nij,ni,nj,m)
			
			comp_NGD_P[j,i] = comp_NGD_P[i,j]
			comp_NGD[j,i] = comp_NGD[i,j]
			comp_R[j,i] = comp_R[i,j]
			comp_R_P[j,i] = comp_R_P[i,j]

	foutName = [dest+'/regionTest_comp_R',dest+'/regionTest_comp_R_P',dest+'/regionTest_comp_NGD',dest+'/regionTest_comp_NGD_P']
	allData = [comp_R,comp_R_P,comp_NGD,comp_NGD_P]

	for i in range(4):
		fout = open(foutName[i]+'.tsv','w')
		fout.write('test')
		for j in range(n):
			fout.write('\t'+header[j])
		fout.write('\n')
		for k in range(n):
			fout.write(header[k])
			for j in range(n):
				fout.write('\t'+str(allData[i][k,j]))
			fout.write('\n')
		fout.close()


def readParams(fpath):
	params = {}
	fin = open(fpath)
	for line in fin:
		if line[0]!='#':
			tmp = line.strip().split('=')
			params[tmp[0]] = tmp[1]

	return(params)

def main():
	global FDR_merge
	global strict_merge
	global FDR_set
	global summaryName
	global mergeDir
	global compareDir


	params_fpath = sys.argv[1] 	
	params = readParams(params_fpath)
	
	if params.has_key('FDR_merge'):
		FDR_merge = float(params['FDR_merge'])
	if params.has_key('strict_merge'):
		strict_merge = int(params['strict_merge'])
	if params.has_key('FDR_set'):
		FDR_set = np.array(params['FDR_set'].split(','),dtype=float)
	if params.has_key('summaryName'):
		summaryName = params['summaryName']
	if params.has_key('mergeDir'):
		mergeDir = params['mergeDir']
	if params.has_key('compareDir'):
		compareDir = params['compareDir']


		

	finPath = params['finPath']
	outDir = params['outDir']
	tests = int(params['tests'])
	if params.has_key('ignore'):
		ignore = params['ignore'].strip().split(',')
	else: ignore = ['*']
	if tests==1: print 'Doing Summary Only'
	elif tests==2: print 'Doing Summary and Merge'
	elif tests==3: print 'Doing Summary, Merge and Compare, may take some time'
	else: raise ValueError('sys arg 3 is wrong, must be 1, 2 or 3')
	
	
	os.makedirs(outDir)

	# write a short info file, readme	
	fout = open(outDir+'/README.txt','w')
	fout.write('This folder contains the meta analysis results for various region test.\nSoftware version tag = '+versionTag+'.\nInput data at '+finPath+'.\n')


 


	
	fin = np.loadtxt(finPath,dtype=str,delimiter='\t')
	data = np.array(fin[1:,1:],dtype=float)
	header = fin[0,1:]
	regionID = fin[1:,0]
	n = len(header)
	
	if ignore[0] != '*':
		m = len(ignore)
		keepInd = np.array(np.ones(n),dtype=bool)
		for i in range(n):
			for j in range(m):
				if header[i].find(ignore[j])>=0: keepInd[i] = False
		data = data[:,keepInd]
		header = header[keepInd]

	fout.write('Ignoring test_ID with any keywords = '+str(ignore)+'.\n')

	
	doComp, doMerg = runSumm(data,header,outDir)
	fout.write('Ran summary analysis:\n\tResults in '+summaryName+'.\n\tConsidered FDR values of '+str(FDR_set)+'.\n')
	if tests > 1:
		runMerg(data[:,doMerg],header[doMerg],regionID,outDir)
		fout.write('Ran region merge analysis:\n\tResults in '+mergeDir+'.\n\tConsidered tests with regions at FDR <= '+str(FDR_merge)+'.\n\tRank agg done with Borda method and geometric seris score.\n\tStric merge param = '+str(strict_merge)+'.\n')
	if tests > 2:
		runComp(data[:,doComp],header[doComp],outDir)
		fout.write('Ran comparision analysis:\n\tResults in '+compareDir+'.\n\tConsidered tests with regions at FDR <= '+str(max(FDR_set))+'.\n')

	fout.close()

 
	


if __name__ == '__main__':
	main()
