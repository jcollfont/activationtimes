# Activation Times

This repository contains methods to compute the activation times from a series of potential recordings on the heart and gradient/hessian/laplacian estimators from a surface geometry.
This work has been developed during the research on ECGI of Jaume Coll-Font and Burak Erem and can be found in several publications:

* B. Erem, J. Coll-Font, R. Martinez Orellana, P. St’ovicek, and D. H. Brooks, “Using Transmural Regularization and Dynamic Modeling for Noninvasive Cardiac Potential Imaging of Endocardial Pacing With Imprecise Thoracic Geometry,” IEEE Trans. Med. Imaging, vol. 33, no. 3, pp. 726–738, Mar. 2014.
* M. J. M. Cluitmans, J. Coll-font, B. Erem, D. Brooks, P. Bonizzi, M. H. Karel, P. G. A. Volders, R. L. M. Peeters, and R. L. Westra, “Spatiotemporal Activation Time Estimation Improves Noninvasive Localization of Cardiac Electrical Activity,” in Computing in Cardiology, 2016, pp. 1185–1188.
* B. Erem, D. H. Brooks, P. M. van Dam, J. G. Stinstra, and R. S. MacLeod, “Spatiotemporal estimation of activation times of fractionated ECGs on complex heart surfaces,” in 2011 Annual International Conference of the IEEE Engineering in Medicine and Biology Society, 2011, pp. 5884–5887.

If you use this code in your research, we would appreciate that you reference them in your related publications. Please note that this is ongoing work, check out new developments in our webpages or in scholar Google!

The code has been developed in MATLAB and it has been mostly documented in the header files. However, this is research code and the documentation might, at times, be incomplete or poorly described.
If you have questions about the code or are interested in development, please contact Jaume Coll-Font at jcollfont ADD gmail.com.

## Contents:

### scripts
This folder contains scripts that run the functions in the other folders. We recommended to look at "RegMatrixEstimation_example.m" to have a reference on how to compute the estimators of gradient, Hessians and Laplacian functions on a surface geometry.

### computeGradientAndHessian
This folder contains all the functions that are necessary to compute the estimators of gradient, Hessians and Laplacian functions on a surface geometry. Please, look at the scripts folder for references on how to use them.
The most important functions are:
* **meshsurfdiffhessmatrix:** computes estimators using the neighbor nodes (defined by the surface graph).
* **meshVolDiffHessMatrix:** computes estimators that flexible definitions of neighborhood to compute the estimators.
* **LaplacianMatrixFromHessianMatrix:** takes the Hessian estimator and provides a Laplacian estimator.

### activationTimesFunctions
These functions are designed to compute the activation times. We recommend to start with the function "activationTimes_wrapper", which runs all the steps to do activation times.
The most important functions are:
* **activationTimes_wrapper:** runs the spatiotemporal activation times and the smoothing approach after.
* **spatiotemporalActtimes:** computes the spatiotemporal activation times method.
* **findMinDVDT:** computes activation times based on the classical definition of minimum dV/dT.
* **smoothactivationtimes:** applies and L-curve based method to smooth the activation times.
