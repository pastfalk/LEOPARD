LEOPARD - A dispersion solver for homogeneous plasmas with arbitrary gyrotropic velocity distributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a manual for the Fortran-90 code LEOPARD containing a short description of the program, an explanation of the input parameters, and some advice for the correct usage of the code. Note that this code version is adapted to the GNU Fortran Compiler and was tested with the compiler version gcc 6.1. For executing the code, first run 'make' and than execute './dsolve'.

General remarks
---------------
LEOPARD ('Linear Electromagnetic Oscillations in Plasmas with Arbitrary Rotationally-symmetric Distributions') is a dispersion relation solver which can derive the frequencies and growth rates of electromagnetic waves with arbitrary propagation angle in homogeneous plasmas with an arbitrary number of particle species. LEOPARD allows the computation of the (fully-kinetic) dielectric tensor components for both bi-Maxwellian or arbitrary gyrotropic velocity distributions. The implemented expressions for the dielectric tensor components are based on 'Kinetic Theory of Plasma Waves' by M. Brambilla (1998).

The solver is intended to enable a systematic study of parallel and oblique wave propagation in realistic plasmas. The velocity distribution can be provided to the code as a data set sampled on the 2d (vpara,vperp) velocity space. The distribution data may be derived from a parametric model distribution, from data obtained by kinetic simulations, or from spacecraft measurements (for more details on the required format of the velocity distribution, see section 'The velocity distribution').

You are welcome to use and distribute this code, and to adapt LEOPARD to your own purposes. If you publish some work which is based on results obtained from LEOPARD, please cite the corresponding paper Astfalk & Jenko, JGR 2017 and make sure that you mark any substantial changes in the code with respect to the original source code.


The velocity distribution
-------------------------
If you choose to include a particle species with a prescribed arbitrary gyrotropic distribution, switch the input parameter 'mode' to '1'. The code then reads the distribution from the file 'distribution/distributionX.dat' where 'X' is the index ('1', '2', '3', ...)  of the corresponding particle species where the numbering is according to all included species with arbitrary velocity distributions. The file is expected to sample the distribution F in parallel and perpendicular velocity space, v_para and v_perp,  with an equidistant grid in each direction. The data has to be arranged in three columns where column 1 = v_para, column 2 = v_perp, column 3 = F(vpara,vperp). The velocities are expected to be in ascending order and for each v_para, F is scanned for all v_perp, before proceeding with the next v_para. For reference, see the sample distribution file which is produced by 'distribution/print_bimax.py'.
The velocity components and the distribution values have to be normalized with respect to the Alfvén velocity of the first particle species. Furthermore, the velocity distribution has to be normalized with respect to the species density such that \int F(v_para,v_perp)*v_perp dv_perp dv_para = 1.0.


The program structure
---------------------
The core of the LEOPARD solver is an iterative root finding algorithm enclosed by a loop over the considered wavenumber interval. Before the loop is started, all required parameters are read from an input file 'input.dat' (see section 'The input file') and all necessary data structures are initialized.

The input provides the iterative root finding algorithm with the initial frequency guesses for the first three wavenumbers of the considered interval. For all subsequent wavenumbers, the initial guesses are determined by the routine 'polyfit()' which uses quadratic polynomials to interpolate previous solutions.

Starting from the supplied initial guess and the given wavenumber, the routine 'muller()' iterates a complex root of the dispersion relation using the Muller method. Thus, for every iteration the determinant of the dispersion tensor has to be evaluated which requires the computation of the dielectric tensor components. This is done by the routine 'disp_det()'. The dielectric tensor can be determined for a bi-Maxwellian distribution (mode='0') with prescribed beta parameters or for an arbitrary gyrotropic velocity distribution (mode='1') which is provided to the code as a file. The bi-Maxwellian scenario requires the evaluation of the plasma dispersion function which is done by the routine 'Z_func()'. For arbitrary distributions, the data is interpolated with cubic splines which allows a piecewise-analytic solution of the required velocity integrations (for more details, see Astfalk & Jenko, JGR 2017). 

After the loop successfully cycled through the wavenumber interval, all roots are printed to the output file 'omega.dat'.


The MPFUN package
-----------------
For the case of particles with arbitrary velocity distribution, the computation of the parallel velocity integrals and the determination of the generalized hypergeometric functions required for the computation of the perpendicular velocity integrals sometimes demands higher accuracy than double or quadruple precision could provide. Therefore, the MPFUN2015 package by David H. Bailey was included which allows thread-safe Fortran computations with arbitrary precision.
The package was found to sometimes cause a breakdown of the code. The problem seems to be related to issues with type conversion but the origin of this has not been located yet. If you find a fix for this issue, feel free to share your solution.



The input parameters
--------------------

&wavenumber

kstart - The lower border of the wavenumber interval the user intends to investigate.

kend   - The upper border of the wavenumber interval the user intends to investigate.

nk     - Here, the user can choose the required wavenumber resolution. nk determines at how many points the code will evaluate the dispersion relation within the chosen wavenumber interval.


Note:
All wavenumbers are given in units of the inertial length, d, of the first particle species.



&initial_guess

omega_r     - The initial guess for the real frequency from which the Muller method starts to find a root, omega(k), of the dispersion relation at the wavenumber k(1)=kstart.

omega_i     - The initial guess for the growth or damping rate from which the Muller method starts to find a root, omega(k), of the dispersion relation at the wavenumber k(1)=kstart.

increment_r - The frequency value by which the previously found root, omega(k), is incremented to provide the starting value for the next Muller iteration at the subsequent wavenumber, k(i+1)=k(i)+dk.

increment_i - The growth rate value by which the previously found root, omega(k), is incremented to provide the starting value for the next Muller iteration at the subsequent wavenumber, k(i+1)=k(i)+dk.


Note:
A proper initial guess is crucial for a successful root finding. If the guess lies too far away from the dispersion branch of interest, you may land on another branch. In general, LEOPARD always converges to a certain root. It has to be figured out by the user, whether this is the root he was searching for.

The increments are only necessary for the initial guesses for the root finding at k(2) and k(3). For subsequent wavenumbers a quadratic polynomial approximation determines all following initial guesses. If dk is not too high, this usually works well.

Both frequencies and growth rates are always given in units of the gyro frequency of the first particle species.



&setup

Nspecies - The number of particle species the user wants to include.

theta 	 - The propagation angle of the waves, i.e. the angle between the wave vector k and the background magnetic field (which is aligned with the z-axis in the chosen coordinate system).

delta 	 - Squared ratio of gyro frequency and plasma frequency of the first particle species.


Note:

The parallel and perpendicular wavenumbers are given as k_para=k*cos(theta) and k_perp=k*sin(theta).

Delta gives a measure for the magnetization of the plasma. Low delta corresponds to weak, high delta corresponds to strong magnetization.



&accuracy

rf_error    - The 'root finding error' gives the exit-condition for the Muller iteration. An error of 1.0d-2 or 1.0d-3 generally gives good results. But - of course - the choice depends on the accuracy requested by the user.

eps_error   - The 'epsilon error' gives the exit condition for the sum over the Bessel index n. Once the relative contribution of the computed dielectric tensor components for a given n gets smaller than the given eps_error, the code exits the loop.



Note:
If a solution seems fishy, play with these parameters and check whether the solution is numerically converged.

Choose the rf_error to be not too demanding, otherwise LEOPARD may run into convergence problems.



&species

mode_in      - Choose '0' for a bi-Maxwellian plasma or '1' for a plasma with arbitrary gyrotropic velocity distribution

q_in         - Charge of the particles in units of the charge of the first particle species.

mu_in 	     - Inverse mass of the particles in units of the inverse mass of the first particle species.

dens_in	     - Density of the particles in units of the density of the first particle species

drift_in     - This introduces a drift velocity to the bi-Maxwellian distribution (mode '0' only). The drift is normalized with respect to the Alfvén velocity.

beta_para_in - Beta parameter parallel to the background magnetic field (mode '0' only).

beta_perp_in - Beta parameter perpendicular to the background magnetic field (mode '0' only).


Note:
If you need more than the two default particle species, just add additional parameter blocks below the two default &species blocks. The choice, which particle species is declared in the first &species block, is of major importance since the normalization of all output data depends on this choice. E.g., if you choose protons to be the first particle species, then all frequencies and growth rates will be given in units of the proton gyrofrequency and the wavenumbers will be in units of the proton inertial length.

When including particle species with arbitrary velocity distribution, the size of the provided velocity grid will significantly affect the performance of LEOPARD. Good accuracy at sufficiently fast run times was found for distribution grids with 200 points in v_para and 100 points in v_perp - but in general the performance is highly dependent on the detailed velocity space structure of the distribution and the considered dispersion branch.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For more information on LEOPARD, see Astfalk & Jenko, JGR 2017. If you have further questions concerning the usage of the code or if you like to discuss some general issues, feel free to contact patrick.astfalk@gmail.com.
