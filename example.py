#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy.linalg as lin
from RBF_module import rbf_module as rm

# Parameters
N = 1000
radius = 100e3
Pn = 1005
Pc = 950
A = 23.0
B = 1.5
rho_air = 1.15
OMEGA = 7.2722052166430395e-5
THETA_0 = 0.52359877559829882 
f = 2.0 * OMEGA * np.sin(THETA_0)
f = 0.0

# Evaluate profiles
#x = np.concatenate((np.linspace(-radius,-0.01,N),np.linspace(0.01,radius,N)),axis=0)
#radius = lambda x, y: np.abs(np.sqrt(x**2 + y**2)) * 1e-3
#p = Pc + (Pn - Pc) * np.exp(-A/(r)**B)
C = 1e1**2*A*B*(Pn-Pc)/(rho_air)

dist = lambda x_1, x_2: 6371*np.sqrt(  (np.sin(x_1[1]) - np.sin(x_2[1]))**2 + (np.cos(x_1[1])*np.sin(x_1[0]) - np.cos(x_2[1])*np.sin(x_2[0]))**2 + (np.cos(x_1[1])*np.cos(x_1[0]) - np.cos(x_2[1])*np.cos(x_2[0]))**2  )

direction = lambda x_1, x_2: np.array([-(x_2[1] - x_1[1]), x_2[0] - x_1[0] ])/lin.norm(np.array([x_2[0] - x_1[0], x_2[1] - x_1[1]])) if x_1 != x_2 else np.array([0,0])

#velocity function: takes distance from center in kilometers

v = lambda r: np.sqrt(C*np.exp(-A/r**B)/r**B + r**2 * f**2 / 4.0) - r * f / 2.0 if r != 0 else 0

wind_speed = lambda x_1, x_2: direction(x_1, x_2)*v(dist(x_1, x_2))

def test(epsilon = 1000000):
	#set up interpolating grid

	n_small = 22
	n_big = 200

	x = np.linspace(45.99, 46.01, n_small)
	y = np.linspace(19.99, 20.01, n_small)
	#eye center
	v_0 = (46, 20)

	values = np.array([wind_speed(v_0, np.array([i,j]))[0] for i in x for j in y])
	points = np.array([[i, j] for i in x for j in y])

	#set up interpolant
	rf_mat = rm.rbf_matrix(points, epsilon)
	coefficients = values.copy()
	rm.get_coefficients(rf_mat, coefficients)

	x_fine = np.linspace(45.99, 46.01, n_big)
	y_fine = np.linspace(19.99, 20.01, n_big)
	points_fine = np.array([[i, j] for i in x_fine for j in y_fine])
	eval_mat = rm.eval_over_array(points, coefficients, points_fine, epsilon)
	test_mat = np.array([wind_speed(v_0, np.array([i,j]))[0] for i in x_fine for j in y_fine])
	print np.max(np.abs(eval_mat - test_mat))

	eval_mat = eval_mat.reshape((n_big,n_big))
	test_mat = test_mat.reshape((n_big, n_big))
	plt.imshow(eval_mat, cmap = cm.coolwarm)
	plt.show()
	plt.imshow(test_mat, cmap = cm.coolwarm)
	plt.show()

	return (points_fine, test_mat, eval_mat, rf_mat)

#plt.figure(1)
#plt.plot(x,v)
#plt.title("Wind Velocity Profile")
#plt.xlabel('km')
#plt.ylabel('m/s')
#plt.axis([np.min(x),np.max(x),0.0,np.max(v)+5])
# plt.savefig('wind.pdf')

#plt.figure(2)
#plt.plot(x,p)
#plt.title("Pressure Profile")
#plt.xlabel('km')
#plt.ylabel('mb')
#plt.axis([np.min(x),np.max(x),Pc-5,Pn+5])
# plt.savefig('pressure.pdf')

#plt.show()
