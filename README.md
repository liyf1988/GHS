
Graphical Horseshoe (GHS) Version 1.0 07/06/2017

DESCRIPTION
-----------
Take Monte Carlo samples from the posterior distribution under the graphical horseshoe prior, to estimate the precision matrix for multivariate Gaussian data.

USAGE
-----
[GHS_omega_save,GHS_lambda_sq_save,GHS_tau_sq_save] = GHS(S,50,100,10000);
GHS_est = mean(GHS_omega_save,3)

ARGUMENTS
---------
S        Y'*Y : sample covariance matrix * n (symmetric)
n        Sample size
burnin   number of MCMC burnins
nmc      number of MCMC saved samples

VALUE
-----
omega_save       p by p by nmc matrices of saved posterior samples of the precision matrix
lambda_sq_save   p*(p-1)/2 by nmc vector of saved samples of lambda squared (local tuning parameter)
tau_sq_save      1 by nmc vector of saved samples of tau squared (global tuning parameter)

DETAILS
-------
Draw samples from the posterior distribution to estimate the precision matrix for multivariate Gaussian data, to generate the graphical horseshoe estimate by Li, Bhadra and Craig (2017).

The function uses matrix decomposition and variable change from the Bayesian graphical lasso by Wang (2012), and the variable augmentation for sampling under the horseshoe prior by Makalic and Schmidt (2016). 

REFERENCES
----------
Yunfan Li, Anindya Bhadra and Bruce A. Craig (2017). The graphical horseshoe estimator for inverse covariance matrices.
Wang, H. (2012). Bayesian graphical lasso models and efficient posterior computation. Bayesian Analysis, 7(4):867-886.
Makalic, E. and Schmidt, D. F. (2016). A simple sampler for the horseshoe estimator. IEE Signal Processing Letters, 23(1):179-182.

EXAMPLES
--------
%Generate eigenvalues for a positive definite matrix, p=50
a(1:20) = 1; a(21:35) = 0.75; a(36:50) = 0.25
%Generate a sparse random precision matrix with eigenvalues a
%with approximately density*30*30 nonzeros
A = sprandsym(50,0.05,a)

%Generate multivariate normal data, n=50
Mu = zeros(50,50)
Sigma = inv(A)
Y = mvnrnd(Mu,Sigma)

%Estimate the precision matrix by GHS
S = Y'*Y
%num_burntin=100, num_iteration=10000
[GHS_omega_save,GHS_lambda_sq_save,GHS_tau_sq_save] = GHS(S,50,100,5000);
GHS_est = mean(GHS_omega_save,3)
%plot the estimated precision matrix
figure
imagesc(GHS_est)


