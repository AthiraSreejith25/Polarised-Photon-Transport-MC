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
#os.system('mkdir {}'.format(target))
reso = 200



for i in range(7):
	for j in range(7):

		num_of_bins = reso
		bin_range = 150*10**-6
		bin_size = (bin_range)/num_of_bins

		cls_intervals = [-(bin_range/2) + i*bin_size for i in range(num_of_bins)]

		#bin_dump
		pixel_array = np.zeros((len(cls_intervals),len(cls_intervals)))

		for k in range(len(spatial_data[i])):

#			print (spatial_data[i][k][0])
#			print (pixel_array[int((spatial_data[i][k][0]+(bin_range/2))/bin_size)][int((spatial_data[i][k][1]+(bin_range/2))/bin_size)])
			pixel_array[int((spatial_data[i][k][0]+(bin_range/2))/bin_size)][int((spatial_data[i][k][1]+(bin_range/2))/bin_size)] += polarization_data[i][k][j]

		pl.figure(figsize=(7,6))
		pl.pcolormesh(pixel_array,cmap = 'jet')
		pl.colorbar()
#		pl.show()

		pl.savefig('{}/stokes_{}{}'.format(target,polarisation[i],polarisation[j]))
		pl.close()

