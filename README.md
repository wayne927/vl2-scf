# vl2-scf
SCF code for Via Lactea II

This is the the code to quickly compute gravitational potential and accelerations inside the N-body data of a dark matter halo. Given the positions and masses of the N-body particles, the self-consistent field (SCF) method solves the Poisson equation by basis decomposition. The basis coefficients are first pre-computed, and then they can be used to evaluate the linear sum of analytical basis functions. The code package includes the coefficients of the main Via Lactea II halo and its subhalos, a driver code to demonstrate how to load these coefficients and setup an orbit, as well as the routines to generate the coefficients if you have your own N-body data.

Full details on the SCF method:
- Hernquist, L., & Ostriker, J. P. 1992, ApJ, 386, 375

Via Lactea II is a cosmological simulation, see
- http://www.ucolick.org/~diemand/vl/

The code here was used for the stream simulations in
- Ngan, W., Bozek, B., Carlberg, R. G., et al. 2015, ApJ, 803, 75
- Ngan, W., Carlberg, R. G., Bozek, B., et al. 2016, ApJ, 818, 194

