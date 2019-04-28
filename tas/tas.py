'''
This module generates thin-cell transmission spectra,
accounting for  cavity effects, transient atom dynamics following
depolarisation in atom-wall collisions, and atom-surface van der Waals
:math:`\propto 1/R^3` interactions.

Example:
    To generate simple thin cell spectra::

        from tas import *
        import matplotlib.pyplot as plt
        import numpy as np

        laserDetuning = np.linspace(-6500,-1500,60)  # (MHz)
        temperature = 228  # (degree Celsius)
        collisionalBroadening = 840  # (MHz)
        C3 = 2  # (kHz mum^3)
        cellLength = 80e-9  # (m)
        collisionalShift = 0  # (MHz)

        calc = ThinCellSpectra(nSurface=1.76)
        T = calc.getTransmission(laserDetuning,
                                   temperature,
                                   collisionalBroadening,
                                   collisionalShift,
                                   C3,
                                   cellLength)

        plt.figure()
        plt.plot(laserDetuning, T,
                 'b', label='Theoretical prediction')
        plt.xlabel('Detuning (GHz)')
        plt.ylabel('Transmission' )
        plt.legend(loc='lower right')
        plt.show()

'''

import numpy as np
#from pylab import *
from scipy.constants import Boltzmann as C_kb,\
                            atomic_mass,\
                            epsilon_0,\
                            hbar

class ThinCellSpectra():
    r"""
    Generates spectra for optically thin atom slabs in nano-cells.

    Includes atom-surface interactions, transient effects of atom dynamics
    following depolarisation at the walls of the cell, and cavity effects
    due to the cell walls.
    Neglects change of driving light power due to interaction with the
    atomic medium.

    Args:
        nSurface: (Optional) refractive index of the vapour cell surface
        wavelength: (Optional) transition wavelength (m)
        gamma0: (Optional) transition natural linewidth (MHz)
        atomMass: (Optional) mass of the atoms (atomic units)
        energyLevelsF: (Optional) array of offsets of energy levels,
            relative to center of gravity of HFS line (MHz)
        cg_coeff: (Optional) array of Clebsch–Gordan coefficients for the energy
            levels listed in `energyLevelsF`.
    Note:
        To include different transitions and/or atomic species
        change atom specific data with optional parameters of ThinCellSpectra
        during initialisation of the class instance.
        Optional parameters are all set by default to conditions in the
        experiment (Cs D1, :math:`F=4 \rightarrow F=3,4` ).

    """

    def __init__(self, nSurface=1.76,
                 wavelength=894.59295986e-9,
                 gamma0=4.5612,
                 atomMass=132.905451931,
                 energyLevelsF=[-3510.916, -4678.597],
                 cg_coeff=[0.026024508, 0.0364343]
                 ):
        # === surface specific data ===
        self.nSurface = nSurface
        # Reflection and transmission coefficients. See Fig. 1 in supplemental material
        self.t1 = 2. * nSurface / (nSurface + 1)
        self.t2 = 2. / (nSurface + 1)
        self.r1 = (nSurface - 1) / (nSurface + 1)
        self.r2 = (1 - nSurface) / (nSurface + 1)

        # === atom specific data ===
        self.kProbe = 2. * np.pi / wavelength  # wavector in vacuum
        self.gamma0 = gamma0
        self.dipole = 3.0 * np.sqrt(epsilon_0 * hbar * (2.0 * gamma0 * (10.0**6))
                                 * (wavelength**3) / (8.0 * np.pi))
        self.mass = atomMass * atomic_mass
        self.energyLevelsF = np.array(energyLevelsF)
        self.cg_coeff = np.array(cg_coeff)

    def getTransmission(self,
                        laserDetuning,
                        vapourTemperature,
                        broadening,
                        resonanceShift,
                        C3,
                        cellLength,
                        velocityStepLargeRelative=20,
                        velocityStepSmallAbsolute=2,
                        smallVelocity=10,
                        dz1=None):
        r"""
        Calculate thin-cell transmission spectra in presence.

        Args:
            laserDetuning: laser detuning (MHz)
            vapourTemperature: vapour temperature (:math:`^\circ\mathrm{C}`)
            broadening: additional Lorentzian broadening that
                accounts for collisions for example (MHz)
            resonanceShift: additional offset of transition resonance that
                accounts for collisions for example (MHz).
            C3: atom-surface van der Waals coefficient
                (:math:`\mathrm{kHz}.\mu\mathrm{m}^3`)
            cellLength: cell thickness (m)
            velocityStepLargeRelative: (Optional) defines velocity steps
                used for integration far from zero velocity defined **relative**
                to the mean 1D speed  of atoms. Default value is 20, resulting
                in steps of 1/20th of mean 1D speed of atoms.
            velocityStepSmallAbsolute: (Optional) defines velocity steps
                for small velocities around 0 in **absolute** unites (m/s).
            smallVelocity: (Optional) defines what is *small velocity* for
                velocityStepLargeRelative and velocityStepSmallAbsolute in
                units of m/s. Default is 10 m/s.
            dz1: (Optional) integration step in space in nm. Default value
                is None, and function will set this step correctly to
                integrate data from experiment (c.f. paper where results are
                published). Outside the range of the experiment considered,
                this parameter might need to be adjusted.

        Returns:
            normalised transmission
            :math:`T = |E(\textrm{with atoms})/E(\textrm{witout atoms})|^2`
        """

        # ===================== Put numbers in S.I units =====================

        detuning = (laserDetuning - resonanceShift) * 1e6 * 2 * np.pi

        # print detune# Delta in S.I unit and angular pulsation for calculation
        # Broadening in S.I unit and angular pulsation for calculation
        gammaTotal = 2 * np.pi * (self.gamma0 + broadening) * 1e6
        # Put C3 atom-surface coefficient in S.I starting from kHz.\mum^3 units
        C3 = 2 * np.pi * 1e-15 * C3

        temperature = vapourTemperature + 273.15  # Temperature in K

        # ========= Compute atomic density from vapour pressure curves =========

        if temperature < 301.65:
            vapourPressure = 10.0**(4.711 - 3999. / temperature)
        else:
            vapourPressure = 10.0**(8.232 - 4062. / temperature
                                    - 1.3359 * np.log10(temperature))
        N0 = 101325.0 * vapourPressure / (C_kb * temperature)


        # =============== Define vectors for future integration ===============

        # Define the 1D most probable velocity
        meanVelocity = np.sqrt(2 * C_kb * temperature / self.mass)

        velmax = meanVelocity
        dv1 = velocityStepSmallAbsolute
        dv2 = meanVelocity / velocityStepLargeRelative  # Integration step for velocities.
        velocityRange1 = np.arange(-2 * velmax, -smallVelocity, dv2)
        velocityRange2 = np.arange(-smallVelocity, smallVelocity, dv1)
        velocityRange3 = np.arange(smallVelocity, 2 * velmax, dv2)
        velocityRange = np.concatenate((velocityRange1,
                                        velocityRange2,
                                        velocityRange3))

        # set  integration step inside the cell
        if dz1 is None:
            if cellLength < 80e-9:
                dz1 = 5 * 1e-10
            elif cellLength > 80e-9 and cellLength < 110e-9:
                dz1 = 1 * 1e-9
            else:
                dz1 = 2 * 1e-9
        dz2 = dz1

        # Boundary vector for integral over dz', L=0, L=lcell
        # avoided due to vdW potential divergence
        zList = np.arange(1e-10, cellLength - 1e-10, dz1)



        '''                          ### Initialisation of iterable quantities ###
        '''

        sum_dz1 = 0  # Intitialise sum over dz1 to 0
        sum_dz2 = 0  # Intitialise sum over dz2 to 0
        # Initialise Resonant Field computed partially analytically to 0.
        E_AP = 0. + 0.j     # Initialise Resonant Field computed numerically to 0

        '''                          ###           Integrals Calculation      ###
        '''

        for f in range(self.energyLevelsF.size):  # sum over hyperfine transitions F -> F'
            delta = detuning - self.energyLevelsF[f] * 2 * np.pi * 1e6
            # prefactor
            pf = self.t1 * self.t2 / (1 - self.r2**2 * np.exp(2j * self.kProbe
                                                             * cellLength)) \
                * self.kProbe * N0 \
                * self.dipole**2 / (2 * epsilon_0 * hbar) * self.cg_coeff[f]

            for vIndex, v in enumerate(velocityRange[1:-1]):  # Loop over velocities
                dv = (velocityRange[vIndex + 1] - velocityRange[vIndex - 1]) / 2

                lambda0 = gammaTotal / 2 - 1j * (delta - self.kProbe * v)
                lambda0P = gammaTotal / 2 - 1j * (delta - self.kProbe * v)
                lambda0M = gammaTotal / 2 - 1j * (delta + self.kProbe * v)

                if v > 5:
                    EaP = 0   # field E_{A+}
                    r2EaM = 0   # field r2 * E_{A-}

                    for z1 in zList:
                        atomSurface1 = 1j * (C3 / (2 * z1**3)) - 1j * \
                              C3 / (2 * z1 * (cellLength - z1)**2)

                        # Correction of Lambda0 by atom-surface int in dz' integral
                        # Define Boundary for dz'' integration
                        zIntegrationList = np.arange(zList[0], z1, dz2)
                        sumEaP_rhoP = 0. + 0. * 1j  # Reinitialize interable to 0
                        sumEaP_rhoM = 0. + 0. * 1j
                        sumEaM_rhoP = 0. + 0. * 1j
                        sumEaM_rhoM = 0. + 0. * 1j

                        for z2 in zIntegrationList:
                            atomSurface2 = 1j * \
                                (C3 / (2 * z2**3)) - 1j * C3 / \
                                (2 * z2 * (cellLength - z2)**2)

                            # Correction of Lambda0 by atom-surface int in integral over z2
                            sumEaP_rhoP += 1. / v\
                                * np.exp(((z2 * (lambda0P + atomSurface2))
                                          - z1 * (lambda0P + atomSurface1)) / v) \
                                * dz2
                            sumEaP_rhoM += self.r2 \
                                * np.exp(2 * 1j * self.kProbe * cellLength)\
                                * np.exp(-2 * 1j * self.kProbe * z1) / v\
                                * np.exp(((z2 * (lambda0M + atomSurface2)) - \
                                        z1 * (lambda0M + atomSurface1)) / v)\
                                * dz2  # Sum over dz''
                            sumEaM_rhoP += 1. / v \
                                * np.exp(2 * 1j * self.kProbe * z1)\
                                * np.exp(((z2 * (lambda0P + atomSurface2)) - \
                                        z1 * (lambda0P + atomSurface1)) / v) \
                                * dz2
                            sumEaM_rhoM += self.r2 \
                                * np.exp(2 * 1j * self.kProbe * cellLength) / v \
                                * np.exp(((z2 * (lambda0M + atomSurface2)) - \
                                        z1 * (lambda0M + atomSurface1)) / v)\
                                * dz2

                        EaP += pf * (sumEaP_rhoP + sumEaP_rhoM) * dz1
                        r2EaM += self.r2 * pf * (sumEaM_rhoP + sumEaM_rhoM) * dz1

                    # Integrate over Maxwell Botlzman velocities
                    E_AP += -(EaP + r2EaM) \
                            * getBoltzmannDistribution(v, meanVelocity) * dv

                elif v < -5:
                    EaP = 0  # E_{A+}
                    r2EaM = 0 # r2* E_{A-}

                    for z1 in zList:
                        atomSurface1 = 1j * (C3 / (2 * z1**3)) - 1j * \
                              C3 / (2 * z1 * (cellLength - z1)**2)

                        # Correction of Lambda0 by atom-surface int in dz' integral
                        # Define Boundary for dz'' integration
                        zIntegrationList = np.arange(z1, zList[-1], dz2)

                        sumEaP_rhoP = 0. + 0. * 1j  # Reinitialize interable to 0
                        sumEaP_rhoM = 0. + 0. * 1j
                        sumEaM_rhoP = 0. + 0. * 1j
                        sumEaM_rhoM = 0. + 0. * 1j
                        for z2 in zIntegrationList:
                            atomSurface2 = 1j * (C3 / (2 * z2**3)) \
                                  - 1j * C3 / (2 * z2 *
                                               (cellLength - z2)**2)

                            sumEaP_rhoP -= 1. / v\
                                * np.exp(((z2 * (lambda0P + atomSurface2)) -
                                         z1 * (lambda0P + atomSurface1)) / v) \
                                * dz2
                            sumEaP_rhoM -= self.r2\
                                * np.exp(2 * 1j * self.kProbe * cellLength)\
                                * np.exp(-2 * 1j * self.kProbe * z1) / v\
                                * np.exp(((z2 * (lambda0M + atomSurface2)) -
                                          z1 * (lambda0M + atomSurface1)) / v)\
                                * dz2  # Sum over dz''
                            sumEaM_rhoP -= 1. / v\
                                * np.exp(2 * 1j * self.kProbe * z1)\
                                * np.exp(((z2 * (lambda0P + atomSurface2)) -
                                         z1 * (lambda0P + atomSurface1)) / v)\
                                * dz2
                            sumEaM_rhoM -= self.r2\
                                * np.exp(2 * 1j * self.kProbe * cellLength) / v \
                                * np.exp(((z2 * (lambda0M + atomSurface2)) -
                                         z1 * (lambda0M + atomSurface1)) / v)\
                                * dz2

                        EaP += pf * (sumEaP_rhoP + sumEaP_rhoM) * dz1
                        r2EaM += self.r2 * pf * (sumEaM_rhoP + sumEaM_rhoM) * dz1

                    # Integrate over Maxwell Botlzman velocities
                    E_AP += -(EaP + r2EaM) \
                            * getBoltzmannDistribution(v, meanVelocity) * dv

                else:  # Case where v=0
                    sumEaP = 0
                    sumEaM = 0

                    Lambda0a = gammaTotal / 2 - 1j * delta

                    for z1 in zList:
                        sumEaP += (1. / (gammaTotal / 2
                                           - 1j * (delta + C3 / z1**3
                                                   + C3 / (cellLength - z1)**3
                                                   - self.kProbe * v))\
                            + self.r2 \
                            * np.exp(2 * 1j * self.kProbe * (cellLength - z1)) / \
                                     (gammaTotal / 2
                                      - 1j * (delta + C3 / z1**3
                                              + C3 / (cellLength - z1)**3
                                              + self.kProbe * v)))\
                            * dz1
                        sumEaM += (np.exp(2 * 1j * self.kProbe * cellLength) / \
                            (gammaTotal / 2 - 1j * (delta + C3 / z1**3 +
                                           C3 / (cellLength - z1)**3 -
                                           self.kProbe * v))\
                            + self.r2\
                            * np.exp(2 * 1j * self.kProbe * (cellLength)) \
                            /(gammaTotal / 2 - 1j * (delta + C3 / z1**3
                                             + C3 / (cellLength - z1)**3
                                             + self.kProbe * v)))\
                            * dz1
                    E_AP += -1. * pf * (sumEaP + self.r2 * sumEaM) * \
                            getBoltzmannDistribution(v, meanVelocity) * dv

        E_0P = self.t1 * self.t2 / (1 - self.r2**2 * np.exp(2 * 1j
                                                            * self.kProbe
                                                            * cellLength))
        # Return the absolute value of the sum of incoherent field (1 here since it is basically )
        transmissionNormalised = np.abs((E_0P + E_AP) / (E_0P))**2
        # E_0 e^{ikz} that appears in both terms, normalized by the non resonant field E_0^{ikz}

        return transmissionNormalised


def getBoltzmannDistribution(velocity, meanVelocity):
    """
    Generate the Boltzmann distribution for the given mean velocity.

    Args:
        velocity: TO-DO
        meanVelocity: TO-DO
    """
    norm = 1. / (meanVelocity * np.sqrt(np.pi))
    return norm * (np.exp(-velocity**2 / meanVelocity**2))

# plt.savefig('C3_m10.pdf')
