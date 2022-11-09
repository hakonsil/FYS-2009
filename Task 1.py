import numpy as n
import matplotlib.pyplot as plt
import math as m
plt.style.use('ggplot')

'''defining constants'''
e0 = 8.854e-12 #Vacuum permittivity
L = 64 #Size of the grid
J = 2**8 #Points on the grid
dx = L/J #Spacing on the grid
grid_x = dx*n.arange(J) #Grid points

def task_c():
    grid_E = n.sin(2*n.pi*grid_x/L) + n.sin(6*n.pi*grid_x/L) #Electric field on the grid
    rho_i_e = e0*(2*n.pi*(3*n.cos((6*n.pi*grid_x)/L)+n.cos((2*n.pi*grid_x)/L)))/L #p_i - p_e (from Poisson's equation)

    rho_hat_i_e = n.fft.rfft(rho_i_e) #calculates fourier coefficients for p_i - p_e
    rho_hat_freq = n.fft.rfftfreq(J, dx) #calculates the frequency of the fourier coefficients

    Em_hat = n.zeros(len(rho_hat_freq), dtype=complex) #creating an array to fill with Em_hat
    Em_hat[1:] = -1.j*(1/(2*n.pi*e0*rho_hat_freq[1:]))*rho_hat_i_e[1:] #calculating the fourier coefficients for Em
        #the slicing is to ensure E_0 = 0
    
    Ex = n.fft.irfft(Em_hat) #calculating E(x) using inverse fourier transform on the fourier coefficients

    plt.plot(grid_x, grid_E, color = 'b', label = 'Analytical')
    plt.plot(grid_x, Ex, '--', color = 'r', label = 'Numerical')
    plt.ylabel('Electric field E(x)')
    plt.xlabel('x')
    plt.legend()
    plt.show()

def task_d(sigma):
    grid_E = (grid_x - L/2) * n.exp(-sigma*(grid_x - L/2)**2) #Electric field on the grid
    rho_i_e = e0* (n.exp(-sigma*(grid_x-L/2)**2)*(1 - 2*sigma*(grid_x-L/2)**2))#p_i - p_e (from Poisson's equation)

    rho_hat_i_e = n.fft.rfft(rho_i_e) #calculates fourier coefficients for p_i - p_e
    rho_hat_freq = n.fft.rfftfreq(J, dx) #calculates the frequency of the fourier coefficients

    Em_hat = n.zeros(len(rho_hat_freq), dtype=complex) #creating an array to fill with Em_hat
    Em_hat[1:] = -1.j*(1/(2*n.pi*e0*rho_hat_freq[1:]))*rho_hat_i_e[1:] #calculating the fourier coefficients for Em
        #the slicing is to ensure E_0 = 0
    
    Ex = n.fft.irfft(Em_hat) #calculating E(x) using inverse fourier transform on the fourier coefficients
    
    plt.plot(grid_x, Ex,color = 'r', label = 'Numerical, $sigma$={x}'.format(x=round(sigma, 4)))
    plt.plot(grid_x, grid_E, '--', color = 'b', label = 'Analytical')
    plt.ylabel('Electric field E(x)')
    plt.xlabel('x')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    task_c()
    task_d(5)
