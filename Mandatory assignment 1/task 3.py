''' This entire file is more or less useless, this is mereley a quick attempt to solve task 3, but it's not close to doing what its supposed to'''

import numpy as np
import matplotlib.pyplot as plt


"""Defining constants/values"""
q = 1
m = 1   

dt = 0.01
dx = 0.01
K = 1000 
J = 1000
#defining initial values
x_0 = [1, 1]
v_0 = [0, 1]
#starte def her som tar inn initial values?


""" definerer funksjoner """
#creating a function to give the rotational matrix
def rotation(t):
    matrix = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    return matrix

def linear_field(pos_0, v_0):
    """ Defining functions """
    #creating empty arrays to fill with position and velocity
    position = np.zeros((K+1, 2))  
    velocity = np.zeros((K+1, 2))

    #defining initial values
    position[0] = np.array(pos_0) 
    velocity[0] = np.array(v_0)

    #defining electric and magnetic fields
    E = np.zeros((K+1)) #defining only E(x)
    B = 0 #defining B as a constant

    #defining E(x) at t=0 (to E(x) = -x)
    x_grid = np.zeros(J)
    for j in range(J):
        x_grid[j] = j*dx

    '''My 'where' function'''
    def where(grid, pos):
        idx = np.abs(grid - pos).argmin()
    
        if pos >= grid[idx]:
            return np.array([grid[idx], grid[idx+1]])
        if pos < grid[idx]:
            return np.array([grid[idx], grid[idx-1]])

    x_range = where(x_grid, position[0, 0])
    x_min = x_range[0]
    x_plus = x_range[1]

    def E_field(x):
        return -(x-((J*dx)/2))

    E_min = E_field(x_min)
    E_plus = E_field(x_plus)
    E[0] = E_min*((x_plus-position[0, 0])/dx)+E_plus*((position[0, 0]-x_min)/dx)

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


    """
    Now the initial velocity is correct as it is a half step before the position so 
    remember that velocity[0] is half a step before position[0]
    """

    for k in range(K):
        half_acc = (q * dt * E[k]) / (2*m) #half acceleration stored as a temp variable
        v_min = np.array([velocity[k, 0] + half_acc, velocity[k, 1]]) #performing the half acceleration

        v_plus = np.matmul(rot_mat, v_min) #rotating

        velocity[k+1, 0] = v_plus[0] + half_acc #doing the second half acceleration in x-dir
        velocity[k+1, 1] = v_plus[1] #doing the second half acceleration in y-dir

        position[k+1, 0] = position[k, 0] + (velocity[k+1, 0] * dt) #updating x-pos
        position[k+1, 1] = position[k, 1] + (velocity[k+1, 1] * dt) #updatin y_pos


        '''My 'where' function'''
        def where(grid, pos):
            idx = np.abs(grid - pos).argmin()
    
            if pos >= grid[idx]:
                return np.array([grid[idx], grid[idx+1]])
            if pos < grid[idx]:
                return np.array([grid[idx], grid[idx-1]])

        x_range = where(x_grid, position[k+1, 0])
        x_min = x_range[0]
        x_plus = x_range[1]

        def E_field(x):
            return -(x-((J*dx)/2))

        E_min = E_field(x_min)
        E_plus = E_field(x_plus)
        E[k+1] = E_min*((x_plus-position[k+1, 0])/dx)+E_plus*((position[k+1, 0]-x_min)/dx)

    plt.plot(position[:, 0], position[:, 1], label = "linear field")
    # plt.show()


if __name__ == "__main__": 
    oppgave_a = linear_field(x_0, v_0)

    plt.title('Starting position: ({x}, {y}) Starting velocity: ({v_x}, {v_y})'
    .format(x = x_0[0], y = x_0[1], v_x = v_0[0], v_y = v_0[1]))
    plt.xlabel('x-position (m)')
    plt.ylabel('y-position (m)')
    plt.legend()
    plt.show()



