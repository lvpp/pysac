#!/usr/bin/env python
'''Provides the basic implementation for the pairwise interaction
segment activity coefficient (SAC) equations.

This code favors readability, please check the accompayning article
for the equation numbers cited here.
'''

from abc import abstractmethod
import math
from copy import deepcopy
from numpy import zeros

RGAS = 0.001985875279
'''Gas constant in kcal/mol/K.'''

RGAS_SI = 8.314462618153
'''Gas constant in J/mol/K.'''

DEFAULT_Q_EFF = math.pi * 1.2**2
'''Default contact surface area in AA2'''

class SAC:
    '''
    This is a basic code for a system consisting in pairwise interacting surface segments.

    The user should derive this class and provide an implementation for the interaction energies,
    `calc_u` method. Optionally, a `calc_du_dT` should also be provided.
    '''

    def __init__(self, Q_eff=DEFAULT_Q_EFF, tol=1e-8, max_iter=400):
        self.tol = tol;
        self.max_iter = max_iter
        self.Q_eff = Q_eff;
    
    @abstractmethod
    def calc_u(self, T, i, j, m, n):
        '''
        This method should calculate the interaction energy for the pair `m-n`.
        The segment `m` is from molecule `i` and `n` from molecule `j`.

        Different model variants should implement this method. The returned value
        is expected to be in `kcal/mol`.
        '''
        pass

    def calc_du_dT(self, T, i, j, m, n):
        '''
        Override this method and calculate the derivative of the interaction
        energy with respect to the temperature (in case it depends on it).
        '''
        return 0

    def set_compounds(self, Q):
        '''
        Sets the compound list by actually sending the segment areas in AA2
        of each compound in the Q argument.

        An example on how to build the Q argument for a binary mixture with
        1 compound with 3 segments and another with just 1 segment follows:
        ```
        Q_1 = [20.167, 21.97,  80.23]
        Q_2 = [191.93]
        Q =   [Q_1, Q_2]
        ```
        '''
        self.ncomps = len(Q)
        self.Q = Q

        # Allocate the theta arrays, same size the given compound areas (Q)
        self.theta = deepcopy(Q)
        self.theta_pure = deepcopy(Q)
        # The segment gammas we initialize all as 1.0
        self.seg_gamma = [[1.0]*len(Q[i]) for i in range(self.ncomps)]
        self.seg_gamma_pure = [[1.0]*len(Q[i]) for i in range(self.ncomps)]
        # Psi and u have 4 dimensions [i][j][m][n], the first 2 we allocate here
        self.psi = [[0 for x in range(self.ncomps)] for y in range(self.ncomps)]
        self.u = [[0 for x in range(self.ncomps)] for y in range(self.ncomps)]
        # The compound total areas
        self.q = [0]*self.ncomps

        for i in range(self.ncomps):
            # Calculating the compound total areas
            for m in range(len(Q[i])):
                self.q[i] += Q[i][m]
            # Calculating the area fractions for pure compounds (theta_pure, Eq. 28 of Ref. 1)
            for m in range(len(self.theta_pure[i])):
                self.theta_pure[i][m] = Q[i][m]/self.q[i]

            # Allocate the psi array for [i][j][m][n]
            for j in range(self.ncomps):
                psi_ij = [[0] * len(Q[j]) for i in range(len(Q[i]))]
                self.psi[i][j] = psi_ij
                u_ij = [[0] * len(Q[j]) for i in range(len(Q[i]))]
                self.u[i][j] = u_ij
        
        # By default 298.15 K and equimolar (if not overriden later)
        self.set_temperature(298.15)
        self.set_composition([1.0/self.ncomps] * self.ncomps)

    def set_temperature(self, T):
        '''
        Sets the temperature for all future calculations.

        This will update all u and psi elements as well as solve all the
        pure compound gammas (as they depend only on temperature).
        '''
        self.T = T;
        self.inv_RT = 1/(RGAS * T)

        # Update the psi array
        for i in range(self.ncomps):
            for j in range(self.ncomps):
                psi_ij = self.psi[i][j]
                u_ij = self.u[i][j]
                for m in range(len(self.Q[i])):
                    for n in range(len(self.Q[j])):
                        u_ij[m][n] = self.calc_u(T, i, j, m, n)
                        # These are the Boltzmann factors, Eq. 2 of Ref. 1
                        psi_ij[m][n] = math.exp(-u_ij[m][n] * self.inv_RT)

        # Converge the pure compound segment gammas, they depend only on temperature
        self.__solve_gamma(True, self.theta_pure, self.seg_gamma_pure)

    def set_composition(self, z):
        '''
        Sets the molar fraction array z, the elements must sum 1.
        '''

        self.z = z
        # Reset the segment gammas we initialize all as 1.0
        self.seg_gamma = [[1.0]*len(self.theta[i]) for i in range(self.ncomps)]

        # First calculate the mixture average area
        q_avg = 0
        for i in range(self.ncomps):
            q_avg += z[i] * self.q[i]
        # Calculate theta for the mixture, Eq. 6 of Ref. 1
        for i in range(self.ncomps):
            for m in range(len(self.Q[i])):
                self.theta[i][m] = z[i] * self.Q[i][m] / q_avg
        self.q_avg = q_avg

    def calc_ln_gamma(self):
        '''
        Calculates the logarithm of the activity coefficients.
        
        Assumes the composition and temperature are already set.
        This method will also solve all the gamma's for the mixtre and has
        to be called before asking for helmholtz or internal energies.
        '''

        # First we converge the mixture gammas 
        self.__solve_gamma(False, self.theta, self.seg_gamma)

        # Initialize the ln gamma vector with all zeros and then add the
        # contribution of all segments for every compound.
        ln_gamma = [0] * self.ncomps;
        for i in range(self.ncomps):
            Q_i = self.Q[i]
            for m in range(len(Q_i)):
                ln_seg_gamma = math.log(self.seg_gamma[i][m])
                ln_seg_gamma_pure = math.log(self.seg_gamma_pure[i][m])
                nu_i_m = Q_i[m]/self.Q_eff

                # Compound activities, Eq. 27 of Ref. 1
                ln_gamma[i] += nu_i_m * (ln_seg_gamma - ln_seg_gamma_pure)

        return ln_gamma

    
    def __solve_gamma(self, for_pure, theta, seg_gamma):
        '''
        Method to solve the "self-consistency" equation (Eq. 10 of Ref. 1) by successive substitution.

        This method should not be called directly by the user.

        We use the same implementation for mixtures and pure compounds, given the theta and
        seg_gamma arrays. The only difference for pure compounds is that we skip the
        cross-compound interactions so we can have a single code base for this.
        '''
        niter = 0
        gamma_norm = 0
        while True:
            niter = niter+1
            gamma_norm_new = 0

            for i in range(self.ncomps):
                seg_gamma_i = seg_gamma[i]
                
                for m in range(len(seg_gamma_i)):
                    # The sum we need to invert in Eq 10 of Ref. 1
                    theta_gamma_psi_sum = 0.0

                    for j in range(self.ncomps):
                        if for_pure and i!=j:
                            continue # when converging the pure substance only i==j matters

                        theta_j = theta[j]
                        seg_gammaj = seg_gamma[j]
                        psi_ij = self.psi[i][j]

                        for n in range(len(theta_j)):
                            psi_mn = psi_ij[m][n]

                            if len(theta_j) == 1 and for_pure:
                                # Successive substitution can have problems to converge correctly
                                # for a pure component with single segment area, but for this case:
                                # seg_gamma = sqrt(1/psi_mn), so then the sum is as follows:
                                theta_gamma_psi_sum = 1.0/math.sqrt(1.0/psi_mn)
                            else:
                                theta_gamma_psi_sum += theta_j[n] * seg_gammaj[n] * psi_mn
						
                    # succesive substitution according to Eq. 10 of Ref. 1 damped with the latest value
                    seg_gamma_i[m] = (seg_gamma_i[m] + 1.0/theta_gamma_psi_sum)/2

                    # calculate a norm2 with our gamma's to check the convergence
                    gamma_norm_new += seg_gamma_i[m]**2

            # finish the norm calculation and check for convergence
            gamma_norm_new = math.sqrt(gamma_norm_new)
            if abs((gamma_norm_new - gamma_norm)/gamma_norm_new) <= self.tol:
                break
            gamma_norm = gamma_norm_new;

            if niter > self.max_iter:
                raise Exception(f"Maximum number of iterations when solving gamma: {self.max_iter}")
            
    def get_helmholtz(self):
        '''
        Returns the residual Helmholtz energy (divided by RT) for the mixture.
        '''
        a = 0
        for i in range(self.ncomps):
            seg_gamma_i = self.seg_gamma[i]
            Q_i = self.Q[i]
                    
            for m in range(len(seg_gamma_i)):
                nu_i_m = Q_i[m]/self.Q_eff
                # Eq. 23
                a += self.z[i] * nu_i_m * math.log(seg_gamma_i[m])

        return a
    
    def __get_partial_helmholtz(self, i, for_pure):
        '''
        Returns the residual partial Helmholtz energy (divided by RT) for the compound i.

        If the `for_pure` flag is `True`, then returns the pure compound version.
        '''
        ai = 0
        seg_gamma_i = self.seg_gamma_pure[i] if for_pure else self.seg_gamma[i]
        Q_i = self.Q[i]
                
        for m in range(len(seg_gamma_i)):
            nu_i_m = Q_i[m]/self.Q_eff
            # Eq. 25
            ai += nu_i_m * math.log(seg_gamma_i[m])

        return ai

    def get_helmholtz_partial(self, i):
        '''
        Returns the residual partial Helmholtz energy (divided by RT) for the compound i.
        '''
        return self.__get_partial_helmholtz(i, False);

    def get_helmholtz_pure(self, i):
        '''
        Returns the residual partial Helmholtz energy (divided by RT) for the compound i.
        '''
        return self.__get_partial_helmholtz(i, True);

    def get_energy_pure(self, i):
        '''
        Returns the residual internal energy (divided by RT) for the pure compound i.
        '''
        ui = 0
        seg_gamma_i = self.seg_gamma_pure[i]
        theta_i = self.theta_pure[i]
        nu_i = self.q[i]/self.Q_eff
        
        for j in range(self.ncomps):
            if i!=j:
                continue # it is pure, only i==j matters

            for m in range(len(seg_gamma_i)):
                theta_j = self.theta_pure[j]
                seg_gamma_j = self.seg_gamma_pure[j]

                psi_ij = self.psi[i][j]
                u_ij = self.u[i][j]

                for n in range(len(theta_j)):
                    # Eq. 7
                    theta_mn = theta_i[m]*seg_gamma_i[m] * theta_j[n]*seg_gamma_j[n] * psi_ij[m][n]
                    du_dT = self.calc_du_dT(self.T, i, j, m, n)

                    # Eq. 35 (see also Eq. 32)
                    ui += nu_i/2 * theta_mn * (u_ij[m][n] - self.T*du_dT)
        
        return ui * self.inv_RT

    def get_energy(self):
        '''
        Returns the residual internal energy (divided by RT) for the mixture.
        '''
        u = 0
        nu = self.q_avg/self.Q_eff
        for i in range(self.ncomps):
            seg_gamma_i = self.seg_gamma[i]
            theta_i = self.theta[i]
                    
            for j in range(self.ncomps):

                for m in range(len(seg_gamma_i)):
                    theta_j = self.theta[j]
                    seg_gamma_j = self.seg_gamma[j]

                    psi_ij = self.psi[i][j]
                    u_ij = self.u[i][j]

                    for n in range(len(theta_j)):
                        # Eq. 7 of Ref. 1
                        theta_mn = theta_i[m]*seg_gamma_i[m] * theta_j[n]*seg_gamma_j[n] * psi_ij[m][n]
                        du_dT = self.calc_du_dT(self.T, i, j, m, n)

                        # Eq. 35  of Ref. 1 (see also Eq. 32)
                        u += nu/2 * theta_mn * (u_ij[m][n] - self.T*du_dT)
        
        return u * self.inv_RT
    
    def get_nonrandom(self):
        '''
        Returns the nonrandom factors for every segment pair.
        '''
        alpha = deepcopy(self.psi)
        for i in range(self.ncomps):
            seg_gamma_i = self.seg_gamma[i]
            theta_i = self.theta[i]
            for j in range(self.ncomps):
                for m in range(len(seg_gamma_i)):
                    theta_j = self.theta[j]
                    seg_gamma_j = self.seg_gamma[j]

                    psi_ij = self.psi[i][j]

                    for n in range(len(theta_j)):
                        # Eq. 8
                        alpha[i][j][m][n] = seg_gamma_i[m]*seg_gamma_j[n]*psi_ij[m][n]
        
        return alpha
    
    def get_nonrandom_pure(self, i):
        '''
        Returns the nonrandom factors for every segment pair.
        '''
        alpha = deepcopy(self.psi)
        for i in range(self.ncomps):
            seg_gamma_i = self.seg_gamma_pure[i]
            theta_i = self.theta_pure[i]
            for j in range(self.ncomps):
                for m in range(len(seg_gamma_i)):
                    seg_gamma_j = self.seg_gamma_pure[j]
                    theta_j = self.theta_pure[j]

                    psi_ij = self.psi[i][j]

                    for n in range(len(theta_j)):
                        # Eq. 8 of Ref. 1
                        alpha[i][j][m][n] = seg_gamma_i[m]*seg_gamma_j[n]*psi_ij[m][n]
        
        return alpha
    
    def get_entropy_pure(self, i):
        '''
        Returns the residual entropy (divided by R) for the pure compound i.
        '''
        si = 0
        seg_gamma_i = self.seg_gamma_pure[i]
        theta_i = self.theta_pure[i]
        nu_i = self.q[i]/self.Q_eff
        
        for j in range(self.ncomps):
            if i!=j:
                continue # it is pure, only i==j matters

            for m in range(len(seg_gamma_i)):
                theta_j = self.theta_pure[j]
                seg_gamma_j = self.seg_gamma_pure[j]
                psi_ij = self.psi[i][j]

                for n in range(len(theta_j)):
                    du_dT = self.calc_du_dT(self.T, i, j, m, n)

                    # Eq. 7 of Ref. 1
                    theta_mn = theta_i[m]*seg_gamma_i[m] * theta_j[n]*seg_gamma_j[n] * psi_ij[m][n]

                    # Eq. 8 of Ref. 1
                    alpha_mn = seg_gamma_i[m] * seg_gamma_j[n] * psi_ij[m][n]

                    # Eq. 11 of Ref. 2
                    si -= nu_i/2 * theta_mn * (math.log(alpha_mn) + du_dT / RGAS )
        
        return si
    
    def get_entropy(self):
        '''
        Returns the entropy (divided by R) for the mixture.
        '''
        nu = self.q_avg/self.Q_eff
        s = 0

        for i in range(self.ncomps):
            theta_i = self.theta[i]
            seg_gamma_i = self.seg_gamma[i]

            for j in range(self.ncomps):
                theta_j = self.theta[j]
                seg_gamma_j = self.seg_gamma[j]
                psi_ij = self.psi[i][j]

                for m in range(len(theta_i)):
                    for n in range(len(theta_j)):
                        du_dT = self.calc_du_dT(self.T, i, j, m, n)

                        # Eq. 7 of Ref. 1
                        theta_mn = theta_i[m]*seg_gamma_i[m] * theta_j[n]*seg_gamma_j[n] * psi_ij[m][n]

                        # Eq. 8 of Ref. 1
                        alpha_mn = seg_gamma_i[m] * seg_gamma_j[n] * psi_ij[m][n]

                        # Eq. 11 of Ref. 2
                        s -= nu/2 * theta_mn * (math.log(alpha_mn) + du_dT / RGAS )
        
        return s

    def get_helmholtz2(self):
        '''
        Returns the Helmnotz (divided by RT) for the mixture, calculated with theta_mn.
        '''
        nu = self.q_avg/self.Q_eff
        a = 0

        for i in range(self.ncomps):
            theta_i = self.theta[i]
            seg_gamma_i = self.seg_gamma[i]

            for j in range(self.ncomps):
                theta_j = self.theta[j]
                seg_gamma_j = self.seg_gamma[j]
                psi_ij = self.psi[i][j]


                for m in range(len(theta_i)):
                    for n in range(len(theta_j)):
                        # Eq. 7 of Ref. 1
                        theta_mn = theta_i[m]*seg_gamma_i[m] * theta_j[n]*seg_gamma_j[n] * psi_ij[m][n]

                        # Eq. 8 of Ref. 2
                        a += nu/2 * theta_mn * math.log(seg_gamma_i[m] * seg_gamma_j[n])
        
        return a