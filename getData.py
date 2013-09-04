import numpy as np

tag = '20130904'
workingDir = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/DF4_2/PPC'
outDir = workingDir+'/results/regionTests'
index = np.loadtxt(workingDir+'/index.tab',dtype=str)
names = index[:,0]
n = len(names)

finName = [ workingDir+'/results/regionTests/regionTest_data_20130830_v0.dat', workingDir+'/results/DF4_2_additional_20130819_hamScat.tsv','/titan/cancerregulome9/ITMI_PTB/users/mmiller/feature_matrix/2013_05_01_df4/2013_07_29_genomic_qc_output/pairwise_test_fm.tsv','/titan/cancerregulome9/ITMI_PTB/users/dgibbs/DF4.2/fbat/regions/fbat_region_combined.txt']




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






m = len(finName)
fin = []
for i in range(m):
	fin.append(open(finName[i]))


fout = open(outDir+'/regionTest_pValue_'+tag+'.dat','w')
fout.write('region_ID')
for i in range(m):
	line = fin[i].next()
	tmp = line.strip().split('\t')
	nLabel = len(tmp)
	for j in range(1,nLabel):
		fout.write('\t'+tmp[j])
fout.write('\n')

for i in range(n):
	fout.write(names[i])
	for j in range(m):
		line = fin[j].next()
		tmp = line.strip().split('\t')
		if tmp[0] != names[i]: raise ValueError('region ID out of order in '+finName[m]+', region '+names[i])
		for k in range(1,len(tmp)):
			fout.write('\t'+tmp[k])
	fout.write('\n')

fout.close()

for i in range(m):
	fin[i].close()

# Now Create the q-value matrix with additional information appended to front

data = np.loadtxt(outDir+'/regionTest_pValue_'+tag+'.dat',dtype=str,delimiter='\t')

colLables = data[1:,:]

data = np.array(data[1:,1:],dtype=float)

n,m = data.shape

for i in range(m):
	p = data[:,i]
	notNAN = ~np.isnan(p)
	q = np.zeros(n)+np.nan
	_,q_tmp,_ = fdr_bh(p[notNAN],.05)
	q[notNAN] = q_tmp
	data[:,i] = q

# now lets do the second round of output:





fout = open(outDir+'/regionTest_pValue_'+tag+'.dat','w')
fout.write('region_ID\tchr\tstart_pos\tstop_pos')
for i in range(m):
	fout.write('\t'+colLables[i])
fout.write('\n')

for i in range(n):
	fout.write(index[i,0]+'\t'+index[i,1]+'\t'+index[i,2]+'\t'+index[i,3])
	for j in range(m):
		tmp = '%05.4E' % (data[i,j])
		fout.write('\t'+tmp)
	fout.write('\n')

fout.close()


