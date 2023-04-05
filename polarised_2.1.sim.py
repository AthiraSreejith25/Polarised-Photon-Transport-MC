# Monte Carlo simulation of polarised photon transport in scattering media

import pickle
from os import system as sys
import numpy as np
import matplotlib.pyplot as pl
import random
from datetime import datetime
from celluloid import Camera
from timeit import default_timer
import miepython

now = datetime.now()

target = './run/{}'.format(now.strftime('%H:%M:%S_%d-%m'))
sys('mkdir {}'.format(target))
start = default_timer()
time_taken = 0


polarisation = 'OHPVMLR'
random.seed(500)

n_photons = 500


def euler_matrix(k,sig):

	c, s, v = np.cos(sig), np.sin(sig), 1-np.cos(sig)
	x, y, z = k[0], k[1], k[2]

	return np.array([[x*x*v+c, y*x*v-z*s, z*x*v+y*s],[x*y*v+z*s, y*y*v+c, y*z*v-x*s],[x*z*v-y*s, y*z*v+x*s, z*z*v+c]])


class Medium:



	def __init__(self,scat,abso,g,refra):

		self.scat = scat #scattering coefficient ()
		self.abso = abso #absorption coefficient ()
		self.aniso = g #anisotropy ()
		self.refra = refra #refractive index ()

	def scat_phase_func(self,wavelength): #scattering phase function

		#re(m) refractive index and im(m) is extinction coeff
		m = complex(self.refra,-(self.abso*wavelength/(4*np.pi)))

#		x = 5.213 #size parameter

#		m = 1.55-0.1j
#		x = 18.621 #sphere diameter = 2 (microns), Re(n_sphere) = 1.59, Im(n_sphere) = -.1, Wavelength = 523nm, conc = .01
		x = 15.54 #sphere diameter = 2.02 (microns), Re(n_sphere) = 1.59, Im(n_sphere) = 0, Wavelength = 543nm, conc = .05

		theta = np.linspace(0,180,181) #this is the range of cosine of the scattering angles
		mu = np.cos(theta*np.pi/180)

		return miepython.mie_S1_S2(m,x,mu)

	def scatter(self):

		p = 0
		p_rand = 1

		while p_rand > p:

			alpha = np.arccos(2*random.random()-1)
			beta = random.random()*2*np.pi

			#S1(alpha) and S2(alpha) are used to calculate the elements of M(alpha)
			#s11 s12 s33 s34 must be real by definition

			S1, S2 = S1_S2[0][int(180*alpha/np.pi)], S1_S2[1][int(180*alpha/np.pi)]

			s11 = (abs(S2)**2 + abs(S1)**2)/2
			s12 = (abs(S2)**2 - abs(S1)**2)/2
			s33 = (S2 * S1.conjugate()).real
			s34 = (S2 * S1.conjugate()).imag


			p = s11*photon.stokes[0] + s12*(photon.stokes[1]*np.cos(2*beta)+photon.stokes[2]*np.sin(2*beta))
			p_rand = s11_0*photon.stokes[0] + s12_0*(photon.stokes[1]*np.cos(2*beta)+photon.stokes[2]*np.sin(2*beta))


		M = np.array([[s11,s12,0,0],[s12,s11,0,0],[0,0,s33,s34],[0,0,-s34,s33]])
		R = np.array([[1,0,0,0],[0,np.cos(2*beta),np.sin(2*beta),0],[0,-np.sin(2*beta),np.cos(2*beta),0],[0,0,0,1]])

		S_new = (np.matmul(M,np.matmul(R,photon.stokes)))

		#normalising the new stokes vector
		S_new[3] /= S_new[0]
		S_new[2] /= S_new[0]
		S_new[1] /= S_new[0]
		S_new[0] = 1

		v_new = np.dot(euler_matrix(photon.drxn,beta),photon.ref)
		u_new = np.dot(euler_matrix(v_new,np.pi*(alpha)/1800),photon.drxn)

		return (S_new,u_new,v_new)


class Photon:


	def __init__(self,medium,stokes,wavelength,core_diameter):

		self.medium = medium

		r = (random.uniform(-1,1)*core_diameter/2)
		theta = random.random()*2*np.pi

		self.dist = 0

		self.wavelength = wavelength

		#critical angle and NA (numerical aperture) are not used, can be included

		self.time_index = 0

		self.rec_path_x = [r*np.cos(theta)]
		self.rec_path_y = [r*np.sin(theta)]
		self.rec_path_z = [0]

		self.weight = 1
		self.pos_x = r*np.cos(theta)
		self.pos_y = r*np.sin(theta)
		self.pos_z = 0

		self.drxn = [0,0,1] #u
		self.speed = 3*(10**8)/self.medium.refra
		self.ref = [0,1,0] #v i	s this [0,1,0] for all u or only u = [0,0,1]??
		self.stokes = stokes #S
		self.threshold = 10**-3 #do we need to change it?


	def step(self):

		return (-np.log(random.random())/(self.medium.scat+self.medium.abso))/100

class Detector:

	stokes_list = [[1,0,0,0],[1,1,0,0],[1,0,1,0],[1,-1,0,0],[1,0,-1,0],[1,0,0,-1],[1,0,0,1]] #O H P V M L R

	def __init__(self,position,resolution):

		self.pos = position
		self.reso = resolution

		self.intensity_data = [[] for j in range(7)]
		self.spatial_data = [[] for j in range(7)]
		self.polarization_data = [[] for j in range(7)]

	def make_bins(self,spatial_data,num_of_bins,bin_range):

		bin_size = (bin_range)/num_of_bins

		return [-(bin_range/2) + i*bin_size for i in range(num_of_bins)]



#medium instantiated with given parameters
solution = Medium(11.88,0,0,1.33)


wavelength = 543*(10**9)

#scattering phase function for media for a given wavelength
S1_S2 = solution.scat_phase_func(wavelength)
s11_0 = (abs(S1_S2[1][0])**2 + abs(S1_S2[0][0])**2)/2
s12_0 = (abs(S1_S2[1][0])**2 - abs(S1_S2[0][0])**2)/2

#super list of all photon trajectories
photon_trajectories_x = [[] for i in range(7)]
photon_trajectories_y = [[] for i in range(7)]
photon_trajectories_z = [[] for i in range(7)]

times = [i*(10**-15) for i in range(101)] #list of times (0 to 101fs) #the  ~.3 microns per femtosecond

detector = Detector(0.5*10**-2,100)

for s in range(len(detector.stokes_list)): #the detector stokes_list will be the same of the sources, hence iterating through detector.stokes_list

	print ('------',polarisation[s]) #user update

	for i_photons in range(n_photons):

		print (i_photons/n_photons) #user update

		photon = Photon(solution,detector.stokes_list[s],wavelength,1200*10**-6) #core_diam 1.2mm

		while photon.weight > photon.threshold:

			length = photon.speed*times[photon.time_index]
			step = photon.step()


			if photon.dist + step >= length and photon.time_index + 1 < len(times):

				s1 = length - photon.dist

				photon.rec_path_x.append(photon.pos_x + s1*photon.drxn[0])
				photon.rec_path_y.append(photon.pos_y + s1*photon.drxn[1])
				photon.rec_path_z.append(photon.pos_z + s1*photon.drxn[2])

				photon.time_index += 1

			photon.weight -= (photon.weight*photon.medium.abso)/(photon.medium.scat+photon.medium.abso)

			photon.stokes, photon.drxn, photon.ref = photon.medium.scatter()


			if photon.pos_z + step*photon.drxn[2] >= detector.pos:

				#updating position of the photon
				photon.pos_x += photon.drxn[0]*(photon.pos_z - detector.pos)/photon.drxn[2]
				photon.pos_y += photon.drxn[1]*(photon.pos_z - detector.pos)/photon.drxn[2]
				photon.pos_z = detector.pos

#				print ('reached detector @ ', photon.pos_x, photon.pos_y)
				detector.spatial_data[s].append((photon.pos_x,photon.pos_y))
				detector.intensity_data[s].append(photon.weight)


				#I,Q,U,V --> O,H,P,V,M,L,R

				I, Q, U, V = photon.stokes[0],photon.stokes[1],photon.stokes[2],photon.stokes[3]

				O = photon.weight*I
				H = photon.weight*(Q+I)/2
				P = photon.weight*(U+I)/2
				V = photon.weight*(I-Q)/2
				M = photon.weight*(I-U)/2
				L = photon.weight*(V+I)/2
				R = photon.weight*(I-V)/2

				detector.polarization_data[s].append((O,H,P,V,M,L,R))


				photon.rec_path_x.append(photon.pos_x)
				photon.rec_path_y.append(photon.pos_y)
				photon.rec_path_z.append(photon.pos_z)

				photon.weight = 0 #Avada kadavra


			photon.pos_x += step*photon.drxn[0]
			photon.pos_y += step*photon.drxn[1]
			photon.pos_z += step*photon.drxn[2]


		photon_trajectories_x[s].append(photon.rec_path_x)
		photon_trajectories_y[s].append(photon.rec_path_y)
		photon_trajectories_z[s].append(photon.rec_path_z)



#pickling the data
with open('{}/data_{}.pickle'.format(target,now.strftime('%H:%M:%S_%d-%m')),'wb') as file:

	data = [detector.spatial_data,detector.polarization_data]
	pickle.dump(data,file)

print ('time taken: ',(default_timer()-start)/60,'min')

