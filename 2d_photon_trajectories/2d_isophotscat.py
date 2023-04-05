import numpy as np
import matplotlib.pyplot as pl
import random
from datetime import datetime
from celluloid import Camera
from timeit import default_timer

start = default_timer()

now = datetime.now()
random.seed(0)

n_photons = 10


class Medium:

	scale = 35

	def __init__(self,scat,abso,g,refra):

		self.scat = scat #scattering coefficient
		self.abso = abso #absorption coefficient
		self.aniso = g #anisotropy
		self.refra = refra #refractive index

	def rand_l_max(self):

		if self.abso == 0:

			return self.scale

		else:

			return -self.scale*(np.log(1-random.random()*self.abso))/self.abso


class Photon:


	def __init__(self,medium):

		self.medium = medium
		self.pos = [0,0]
		self.dist = 0
		self.time = 0
		self.path_x = [0]
		self.path_y = [0]
		self.l_max = medium.rand_l_max()
		self.next_drxn = [0,0]
		self.speed = 3*(10**8)/self.medium.refra

	def iso_rand_vec(self): #for isotropic medium

		phi = random.random()*2*np.pi

		cosphi = np.cos(phi)
		sinphi = np.sin(phi)

		return ((cosphi,sinphi))

	def step(self):

		return -np.log(random.random())/self.medium.scat

medium_1 = Medium(0.92,0.8,0,1.33) #medium instantiated with given parameters

photon_trajectories_x = []
photon_trajectories_y = []

for i_photons in range(n_photons):

	photon = Photon(medium_1) #photon instantiated in the given medium

	while photon.dist <= photon.l_max:

		photon.next_drxn = photon.iso_rand_vec()

		for coord in range(2):

			step = photon.step()
			photon.pos[coord] += step*photon.next_drxn[coord]

		photon.dist += step
		photon.path_x.append(photon.pos[0])
		photon.path_y.append(photon.pos[1])

	photon_trajectories_x.append(photon.path_x)
	photon_trajectories_y.append(photon.path_y)

sim_time = default_timer()

#data analysis and data visualisation

fig = pl.figure()
cam = Camera(fig)

longest = max([len(i) for i in photon_trajectories_x])

for i in range(longest):

	for j in range(len(photon_trajectories_x)):

		if i < len(photon_trajectories_x[j]):

			if i > 2:

				for alpha in range(1,4):

					pl.plot(photon_trajectories_x[j][i-4+alpha:i-2+alpha],photon_trajectories_y[j][i-4+alpha:i-2+alpha],'#1167b1{}'.format(str(hex(int(alpha*255/7)))[2:])) #blue

				pl.plot(photon_trajectories_x[j][i],photon_trajectories_y[j][i],color = '#03ac13{}'.format(str(hex(80))[2:]), marker = '.') #green

			else:

				pl.plot(photon_trajectories_x[j][:i+1],photon_trajectories_y[j][:i+1],'#1167b1{}'.format(str(hex(80))[2:])) #blue
				pl.plot(photon_trajectories_x[j][i],photon_trajectories_y[j][i],color = '#03ac13{}'.format(str(hex(80))[2:]), marker = '.') #green

		else:

			pl.plot(photon_trajectories_x[j],photon_trajectories_y[j],'#d21404{}'.format(str(hex(30))[2:])) #red

	cam.snap()
	print ('countdown:  {}'.format(longest-i))

print ('\nconverting individual frames to movie in mp4 format (use ffmpeg to convert to required format, if needed)\nplease wait...')
animation = cam.animate()
animation.save('{}.mp4'.format(now.strftime('%H-%M-%S_%d-%m')))


end = default_timer()

with open('./{}.info'.format(now.strftime('%H:%M:%S_%d-%m')),'w') as file:

	file.write("Photon trajectory animator -- shows animation of photons in a scattering and absorbing medium.\n\n")
	file.write("Program executed at {}\n\n".format(now.strftime('%H:%M:%S %d/%m/%Y')))

	file.write("Number of photons = {}\n".format(n_photons))
	file.write("Absorption coefficient = {}\n".format(medium_1.abso))
	file.write("Scattering coefficient = {}\n".format(medium_1.scat))
	file.write("Refractive index = {}\n".format(medium_1.refra))
	file.write("Anisotropy = {}\n\n".format(medium_1.aniso))

	file.write("Simulation runtime: {}\n".format(sim_time-start))
	file.write("Animation runtime: {}".format(end-sim_time))
