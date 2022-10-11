import numpy as np
import matplotlib.pyplot as plt
import random as rnd
plt.style.use("ggplot")

"""Defining constants/values"""
q = 1
m = 1   

dt = 0.01
K = 1000

#defining initial values
x_0 = [0, 0]
v_0 = [5, 1]



""" Defining functions """
def rotation(t): #Gives the rotational matrix
    matrix = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    return matrix

def constant_field(pos_0, v_0): #Plots a particle in E(x) = E_0
    """ definerer verier """
    #creating empty arrays to fill with position and velocity
    position = np.zeros((K+1, 2))  
    velocity = np.zeros((K+1, 2))

    #defining initial values
    position[0] = np.array(pos_0) 
    velocity[0] = np.array(v_0)

    #defining electric and magnetic fields
    E = np.zeros((K+1)) #defining only E(x)
    B = 2 #defining as a constant

    #making a quick loop to make E constant
    for k in range(K):
        E[k] = 1
    #defining E(x) at t=0 (to E(x) = -x)
    #E[0] = -position[0,0]

    #defining theta 
    theta = - ((q * B * dt) / m)


    """Performing the half rotation backwards"""
    #defining the rotational matrix for the half rotation backwards
    theta_back = rotation(- (theta / 2))

    #defining the rotational matrix to be used in the for-loop
    rot_mat = rotation(theta)

    #performing the ACTUAL half rotation backwards
    v_min1 = np.matmul(theta_back, velocity[0])

    #performing the half deceleration
    velocity[0, 0] = v_min1[0] - ((q * dt * E[0]) / (2*m))
    velocity[0, 1] = v_min1[1]


    """Boris mover algorithm"""

    for k in range(K):
        half_acc = (q * dt * E[k]) / (2*m) #half acceleration stored as a temp variable
        v_min = np.array([velocity[k, 0] + half_acc, velocity[k, 1]]) #performing the half acceleration

        v_plus = np.matmul(rot_mat, v_min) #rotating

        velocity[k+1, 0] = v_plus[0] + half_acc #doing the second half acceleration in x-dir
        velocity[k+1, 1] = v_plus[1] #doing the second half acceleration in y-dir

        position[k+1, 0] = position[k, 0] + (velocity[k+1, 0] * dt) #updating x-pos
        position[k+1, 1] = position[k, 1] + (velocity[k+1, 1] * dt) #updatin y_pos
        
    plt.plot(position[:, 0], position[:, 1], label = "linear field")

def analytic(pos_0, v_0): #Plots analytic solution of ExB drifrt
    E = 1 #defining only E(x)
    B = 2 #defining as a constant (B = B_z)

    """ Defining values """
    #creating empty arrays to fill with position and velocity
    position = np.zeros((K+1, 2))  
    velocity = np.zeros((K+1, 2))

    #defining initial values
    position[0] = np.array(pos_0) 
    velocity[0] = np.array(v_0)

    omega = q*B/m #Defining Omega_c
    r_c = (m/(q*B))*np.sqrt((v_0[0]**2)+((v_0[1]+ (E/B))**2)) #Defining gyro-radius
    theta = np.arctan((-v_0[1] - (E/B))/v_0[0]) #Defining gyro-phase
    x_c0 = pos_0[0] - r_c*np.sin(theta) #defining guiding center x-pos at t=0
    y_c0 = pos_0[1] - r_c*np.cos(theta) #defining guiding center y-pos at t=0

    for k in range(K):
        position[k+1, 0] = x_c0 + r_c*np.sin(omega*k*dt + theta)  #updating x-pos
        position[k+1, 1] = y_c0 + r_c*np.cos(omega*k*dt + theta) - (E/B)*k*dt #updatin y_pos
        
    plt.plot(position[:, 0], position[:, 1], label = "Analytic ExB drift", ls='--')

def compare(N): #Compares analytic and numerical solution for N particles
    for n in range(N):
        a = rnd.randint(1, 30)
        b = rnd.randint(1, 30)
        constant_field([20*n, 3*n], [a, b])
        analytic([20*n, 3*n], [a, b])

'''Main function'''
if __name__ == "__main__": 
    oppg_4c = constant_field(x_0, v_0)
    plt.xlabel('x-position (m)')
    plt.ylabel('y-position (m)')
    plt.title('E = $E_x$, B = $B_z$')
    plt.axis("equal")
    plt.show()