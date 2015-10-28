# Function to compute energy loss

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import rv_continuous
from scipy.integrate import quad
import sys
import time


class EnergyLoss(object):

    '''
    So far for energy-independent dispersion relations m=0
    '''

    def __init__(self, m=0, n_particles=100, constant_loss=None,
                 glashow_unit=None):
        '''
        m: int determines the dispersion relation
        n_particles: intnumber of particles
        constant_loss: if None, the energy loss per decay
            is given by probability distribution defined
            below. If float, proportion of energy lost in each
            decay.
        glashow_unit : if None, the unit of energy in final_energy
        is 1 / (L * A), where L is the total distance, and A is the 
        energy independent factor in the decay rate. If float, 
        the unit of energy is the terminal energy of Cohen-Glashow with
        glashow_unit equal to k. 
        '''
        self.m = m
        self.n_particles = n_particles
        self.constant_loss = constant_loss
        self.glashow_unit = glashow_unit

    def decay_prob(self, E, x):
        '''
        E is the current energy in units of the terminal energy,
        x is the distance in units of the total flight path

        k=None means terminal energy is equal to 1 / (A * L),
        where A is the energy independent factor of the decay rate
        and L the total flight path. k = float means the terminal energy
        is the one defines by Cohen and Glashow: 1 / (A * (5+3m) * k * L)
        '''
        c = 1.0
    
        if self.glashow_unit:
            k = self.glashow_unit
            c = (5.0 * k)

        return np.exp(- E**5 * x / c)

    def final_energy(self, min_energy=0.1, max_energy=100, mode='log',
                     n_points=100, threshold=0):

        self.mode = mode
        self.n_points = n_points
        # initial energy in units of terminal energy
        if mode == 'log':
            self.initial_energies = np.logspace(np.log10(min_energy),
                                                np.log10(max_energy),
                                                self.n_points)
        else:
            self.initial_energies = np.linspace(min_energy, max_energy,
                                                self.n_points)
        self.mean_energy = np.zeros(n_points)
        self.std_energy = np.zeros(n_points)
        self.mean_decays = np.zeros(n_points)
        self.std_decays = np.zeros(n_points)

        # Initializa gamma distribution defined below
        gamma = gamma_gen(name='gamma', a=0, b=1)
        # set of particle energies that end with energies
        # higher than 1.5 times the terminal energy
        self.high = []

        for i, Ei in enumerate(self.initial_energies):

            if i % 10 == 0:
                print 'Initial energy = %.02f' % Ei
            # final energy in row 0, number decays row 1
            final = np.zeros((2, self.n_particles))
            for particle in xrange(self.n_particles):

                E = Ei
                n = 0
                x_2 = 0

                while x_2 < 1:

                    interval = np.min([1 / (10.0 * E**5), 1 - x_2])
                    possibilities = np.array([True, False])
                    weights = np.array([1 - self.decay_prob(E, interval),
                                        self.decay_prob(E, interval)])
                    decay = np.random.choice(possibilities, size=1, p=weights)
                    
                    if decay:
                        # l = 1 - k
                        if self.constant_loss:
                            l = 1 - self.constant_loss
                        else:
                            l = gamma.rvs()
                        E *= l
                        n += 1

                    if E <= threshold:
                        break

                    x_2 += interval

                final[0, particle] = E
                final[1, particle] = n

                if E > 1.5:

                    self.high.append((Ei, E))

            mean = final.mean(axis=1)
            std = final.std(axis=1)
            self.mean_energy[i] = mean[0]
            self.std_energy[i] = std[0]
            self.mean_decays[i] = mean[1]
            self.std_decays[i] = std[1]

    def decay_prob_prop(self, E, x):
        '''
        E is the current energy in units of the initial energy, 
        x is the distance in units of the initial mean free path
        '''

        return np.exp(- E**5 * x)

    def propagation(self, distance=1e+3, threshold=0):

        # total distance interms of initial free mean path
        self.distance = distance
        self.space_points = np.arange(0, distance, 0.1)

        # Initializa gamma distribution defined below
        gamma = gamma_gen(name='gamma', a=0, b=1)
        
        energy = np.zeros((self.n_particles, len(self.space_points)))
        for particle in xrange(self.n_particles):

            if particle % 10 == 0:
                print 'Particle number %d is flying' % particle

            # energy in units of initial energy
            E = 1
            x_0 = 0
            for i, x in enumerate(self.space_points):

                interval = x - x_0
                energy[particle, i] = E
                possibilities = np.array([True, False])
                weights = np.array([1 - self.decay_prob_prop(E, interval),
                                    self.decay_prob_prop(E, interval)])
                decay = np.random.choice(possibilities, size=1, p=weights)

                if decay:
                    # l = k - 1
                    if self.constant_loss:
                        l = 1 - self.constant_loss
                    else:
                        l = gamma.rvs()
                    E *= l

                if E <= threshold:
                    break

                x_0 = x

        self.current_energy = energy.mean(axis=0)
        self.std_current = energy.std(axis=0)


    def final_initial(self, std=False, n_decays=False):

        fig = plt.figure(figsize=(20,15))
        ax1 = plt.subplot(111)

        ax1.plot(self.initial_energies, self.mean_energy,
                 label='mean final energy within 5 sigmas',
                 color='blue')
        if std:
            ax1.plot(self.initial_energies,
                     self.mean_energy + 3*self.std_energy,
                     color='goldenrod', ls='--', label="mean + 3 std")

        confidence = 5.0 * self.std_energy / np.sqrt(self.n_particles)
        ax1.fill_between(self.initial_energies,
                         self.mean_energy - confidence,
                         self.mean_energy + confidence, alpha=0.2,
                         color='blue')

        initial_high = [i[0] for i in self.high]
        final_high = [i[1] for i in self.high]
        ax1.scatter(initial_high, final_high, color='green',
                    label='individual final energies > 1.5')

        ylabel = ''
        if n_decays:
            # ax2 = ax1.twinx()
            ax1.plot(self.initial_energies, self.mean_decays,
                     color='darkmagenta',
                     ls='--', label='mean number decays')
            ylabel += r'$ \, ,\,  \langle N \rangle $'
            # ax2.set_ylabel(r'$\langle N \rangle$')
            # ax2.legend(loc='best')
            # ax2.set_xscale('log')

        ax1.set_xscale('log')
        ax1.legend(loc='best')
        if self.glashow_unit:
            xlabel = '$E_i / E_t$'
            ylabel = r'$ E_f / E_t$' + ylabel
        else:
            xlabel = r'$E_i / E_L \, \left(E_L \equiv 1/ (A \times L), \,' + \
                r'\Gamma (E) = A\times E^{5} \right)$'
            ylabel = r'$ E_f / E_L$' + ylabel
        ax1.set_xlabel(xlabel)
        ax1.set_xlim((self.initial_energies[0], self.initial_energies[-1]))
        ax1.set_ylabel(ylabel)
        ax1.set_title('''
            %d particles per energy. %d different initial energies, 
            %s distribuited. Constant energy loss: %s ''' %
                      (self.n_particles, self.n_points, self.mode,
                       str(self.constant_loss)))

        # ax.set_xticks(np.arange(0, 10, 0.1))
        # ax.set_xticklabels(np.arange(0, 10, 0.1))
        # ax.xaxis.set_ticks(np.logspace(-1, 1, self.n_points/ 2))
        # plt.show((block=False))


    def final_propagation(self, std=False):

        fig, ax = plt.subplots(figsize=(15, 10), ncols=1)

        ax.plot(self.space_points, self.current_energy, color='b',
                label='mean energy within 5 sigmas')
        confidence = 5.0 * self.std_current / np.sqrt(self.n_particles)
        ax.fill_between(self.space_points, self.current_energy - confidence,
                        self.current_energy + confidence,
                        alpha=0.2, color='blue')

        if std:
            plt.plot(self.space_points,
                     self.current_energy + 3*self.std_current,
                     color='goldenrod', ls='--', label='mean + 3 std')

        ax.set_xscale('log')
        ax.set_xlim(xmin=0.0, xmax=self.distance)
        ax.set_ylim([0, 1])
        ax.set_xlabel(r'$x / \Gamma ^{-1} (E_i)$')
        ax.set_ylabel(r'$E(x) / E_i$')
        ax.legend(loc='best')
        ax.set_title('''
            %d particles. %d different spatial points, 
            Constant energy loss: %s ''' %
                     (self.n_particles, len(self.space_points),
                      str(self.constant_loss)))

        # plt.xlim(xmin=0.0, xmax=self.distance)

        plt.show()


class gamma_gen(rv_continuous):

    def _pdf(self, x):
        '''
        The fact that the distribution is zero for x<0 and x>1 is
        implemented through args a and b when initializating the class
        '''
        norm_factor = quad(lambda x: (1-x)**3 * (1 - x**3), 0, 1)[0]
        return (1-x)**3 * (1 - x**3) / norm_factor
