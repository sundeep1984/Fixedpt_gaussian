# Fixedpt_gaussian
Generate Gaussian binary vectors using the approach in R. Gutierrez, V. Torres, J.Valls, “Hardware Architecture of a Gaussian Noise Generator Based on Inversion Method, ” IEEE Transactions on Circuits &amp; Systems II, Vol.59, No.8, pp.501-505, 2012

In order to break the interval [0.5 1] into segments, you need to run segment.m, which implements the algorithm given in "Hierarchical Segmentation Schemes for Function Evaluation," D.U. Lee, W. Luk, J. Villasenor and P.Y.K. Cheung,  Proceedings. 2003 IEEE International Conference on Field-Programmable Technology (FPT). 

This script calls the function minimax and essentially does quadratic curve fitting across segments that it selects as per the minmax criterion described in the paper(minimizing max error of fit across all segments) and stores the coefficients for each segments in a file called coeffs.mat.

The main script for Gaussian number generation is Gauss_fixedpt.m, which takes as inputs the number of Gaussian samples to generate as well as wheter to write the test vectors. It also reads in coefficients from coeffs.txt to use for gaussian number generation
