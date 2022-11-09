import numpy as n
import math as m
import matplotlib.pyplot as plt

'''defining constant values'''
q=1
J=32
rho = n.zeros(J)

def update_field(grid, N):
    L = 10 # length of domain in units of Debye length
    seed = 1234 # seed for the random number generator
    rng = n.random.default_rng(seed=seed)
    dx=L/J 
    particle_positions = rng.uniform(0, L, size=N) # creating N particles uniformly distributed

    for pos in particle_positions:
        # finding the closest grid points
        j_0 = m.floor(pos/dx)
        j_1 = j_0+1

        #Updating the charge density at the two closest grid points
        grid[j_0] += (q/dx)*((j_1*dx-pos)/dx)
        grid[j_1 % J] += (q/dx)*((pos-j_0*dx)/dx) 

    # Defining the theoretical charge density
    rho_theo = n.zeros(J)
    for i in range(J):
        rho_theo[i] = q*N/L
    
    plt.plot(n.arange(0, L, dx), grid,'o', label = 'charge density')
    plt.plot(n.arange(0,L, dx), rho_theo, '--', label = 'theoretical charge density')
    plt.xlabel('position')
    plt.ylabel('charge density')
    plt.legend()
    plt.show()    
    print(max(n.abs(grid-rho_theo)/rho_theo)) #printing the agreement between the two solutions
    return (n.abs(grid-rho_theo)/rho_theo)


if __name__ == "__main__": 
    update_field(rho, 10000)