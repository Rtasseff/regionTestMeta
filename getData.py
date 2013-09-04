import numpy as np

tag = '20130830'
workingDir = '/titan/cancerregulome9/ITMI_PTB/users/rtasseff/DF4/DF4_2/PPC'
outDir = workingDir+'/results/regionTests'
index = np.loadtxt(workingDir+'/index.tab',dtype=str)
names = index[:,0]
n = len(names)

finName = [ workingDir+'/results/regionTests/regionTest_data_20130830_v0.dat', workingDir+'/results/DF4_2_additional_20130819_hamScat.tsv','/titan/cancerregulome9/ITMI_PTB/users/mmiller/feature_matrix/2013_05_01_df4/2013_07_29_genomic_qc_output/pairwise_test_fm.tsv','/titan/cancerregulome9/ITMI_PTB/users/dgibbs/DF4.2/fbat/regions/fbat_region_combined.txt']

m = len(finName)
fin = []
for i in range(m):
	fin.append(open(finName[i]))


fout = open(outDir+'/regionTest_data_'+tag+'.dat','w')
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




