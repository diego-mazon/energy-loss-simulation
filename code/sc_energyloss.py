# script to run energyloss

from energyloss import EnergyLoss as EL
import sys


def main():
    if kind == 'constant':
        el_const = EL(
            n_particles=n_particles, constant_loss=True, glashow_unit=True)
        el_const.final_energy()
        print el_const.final_initial(std=True)

    if kind == 'fluctuating':
        el_var = EL(n_particles=n_particles)
        el_var.final_energy()
        el_var.final_initial()

if __name__ == '__main__':
    kind = sys.argv[1]
    n_particles = int(sys.argv[2])
    main()
