import pickle
import os
import numpy as np
import matplotlib.pyplot as pl

file_list = os.listdir('data')


spatial_data = [[] for j in range(7)]
polarization_data = [[] for j in range(7)]


for i in file_list:
	with open('data/{}/data_{}.pickle'.format(i,i),'rb') as file:

		data = pickle.load(file)

		for j in range(len(data[0])):

			spatial_data[j] += data[0][j]
			polarization_data[j] += data[1][j]

polarisation = 'OHPVMLR'
target = 'plots'

reso = 200
bin_range = 150*10**-6
num_of_bins = reso
bin_size = (bin_range)/num_of_bins

cls_intervals = [-(bin_range/2) + i*bin_size for i in range(num_of_bins)]

target = 'plots/{}_{}'.format(str(reso),str(bin_range*10**6))

os.system('mkdir {}'.format(target))

M = []

print ('...\n')

for i in range(7):

	tmp = []

	for j in range(7):

		#bin_dump
		pixel_array = np.zeros((len(cls_intervals),len(cls_intervals)))

		for k in range(len(spatial_data[i])):

			pixel_array[int((spatial_data[i][k][0]+(bin_range/2))/bin_size)][int((spatial_data[i][k][1]+(bin_range/2))/bin_size)] += abs(polarization_data[i][k][j])

		tmp.append(pixel_array)


	M.append(tmp)

print ('calculating Mueller matrix elements...\n')

mueller = np.zeros((4,4),dtype=list)

d = {'o':0,'h':1,'p':2,'v':3,'m':4,'l':5,'r':6}

mueller[0][0] = M[d['o']][d['o']]
mueller[0][1] = (M[d['h']][d['o']]-M[d['v']][d['o']])/2
mueller[0][2] = (M[d['p']][d['o']]-M[d['m']][d['o']])/2
mueller[0][3] = (M[d['l']][d['o']]-M[d['r']][d['o']])/2

mueller[1][0] = (M[d['o']][d['h']]-M[d['o']][d['v']])/2
mueller[1][1] = ((M[d['h']][d['h']]+M[d['v']][d['v']])/4) - ((M[d['h']][d['v']]+M[d['v']][d['h']])/4)
mueller[1][2] = ((M[d['p']][d['h']]+M[d['m']][d['v']])/4) - ((M[d['p']][d['v']]+M[d['m']][d['h']])/4)
mueller[1][3] = ((M[d['l']][d['h']]+M[d['r']][d['v']])/4) - ((M[d['l']][d['v']]+M[d['r']][d['h']])/4)

mueller[2][0] = (M[d['o']][d['p']]-M[d['o']][d['m']])/2
mueller[2][1] = ((M[d['h']][d['p']]+M[d['v']][d['m']])/4) - ((M[d['h']][d['m']]+M[d['v']][d['p']])/4)
mueller[2][2] = ((M[d['p']][d['p']]+M[d['m']][d['m']])/4) - ((M[d['p']][d['m']]+M[d['m']][d['p']])/4)
mueller[2][3] = ((M[d['l']][d['p']]+M[d['r']][d['m']])/4) - ((M[d['l']][d['m']]+M[d['r']][d['p']])/4)

mueller[3][0] = (M[d['o']][d['l']]-M[d['o']][d['r']])/2
mueller[3][1] = ((M[d['h']][d['l']]+M[d['v']][d['r']])/4) - ((M[d['h']][d['r']]+M[d['v']][d['l']])/4)
mueller[3][2] = ((M[d['p']][d['l']]+M[d['m']][d['r']])/4) - ((M[d['p']][d['r']]+M[d['m']][d['l']])/4)
mueller[3][3] = ((M[d['l']][d['l']]+M[d['r']][d['r']])/4) - ((M[d['l']][d['r']]+M[d['r']][d['l']])/4)


print ('plotting...')

for i in range(4):
	for j in range(4):

		pl.figure(figsize=(7,6))
		pl.pcolormesh(mueller[i][j],cmap='jet')
		pl.colorbar()
#		pl.show()

		pl.savefig('{}/mueller_{}{}'.format(target,str(i+1),str(j+1)))
		pl.close()

