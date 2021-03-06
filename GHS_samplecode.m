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
imagesc(GHS_est)