% Yunfan Li, March 2017
% Sampling for graphical horseshoe

function [omega_save,lambda_sq_save,tau_sq_save] = GHS(S,n,burnin,nmc)
% GHS MCMC sampler using data-augmented
% block (column-wise) Gibbs sampler
% Input:
%     S = Y'*Y : sample covariance matrix * n
%     n: sample size
%     burnin, nmc : number of MCMC burnins and saved samples

% Output:
%     omega_save: p by p by nmc matrices of saved posterior samples of
%     precision matrix
%     lambda_sq_save: p*(p-1)/2 by nmc vector of saved samples of lambda
%     squared (local tuning parameter)
%     tau_sq_save: 1 by nmc vector of saved samples of tau squared (global
%     tuning parameter)

[p] = size(S,1)
omega_save = zeros(p,p,nmc);
lambda_sq_save = zeros(p*(p-1)/2,nmc);
tau = zeros(nmc);

ind_all = zeros(p-1,p);
for i = 1:p
       if i==1  
       ind = [2:p]'; 
      elseif i==p
       ind = [1:p-1]'; 
      else
       ind = [1:i-1,i+1:p]';
       end
       
       ind_all(:,i) = ind;
end

% set initial values
Omega = eye(p); Sigma = eye(p);
Lambda_sq(1:p,1:p) = 1; Nu(1:p,1:p) = 1; tau_sq = 1; xi = 1;

for iter = 1: burnin+nmc    
          
%    if(mod(iter,1000)==0)
%     fprintf('iter = %d \n',iter);
%    end
       
%%% sample Sigma and Omega=inv(Sigma)
    for i = 1:p
      ind = ind_all(:,i);     
      Sigma_11 = Sigma(ind,ind); sigma_12 = Sigma(ind,i);
      sigma_22 = Sigma(i,i);
      s_21 = S(ind,i); s_22 = S(i,i);
      lambda_sq_12 = Lambda_sq(ind,i); nu_12 = Nu(ind,i);
      %% sample gamma and beta
      gamma = gamrnd((n/2+1),2/s_22);    % random gamma with shape=n/2+1, rate=s_22/2
      inv_Omega_11 = Sigma_11 - sigma_12*sigma_12'/sigma_22;
      inv_C = s_22*inv_Omega_11+diag(1./(lambda_sq_12*tau_sq));
      inv_C_chol = chol(inv_C);
      mu_i = -inv_C\s_21;
      beta = mu_i+ inv_C_chol\randn(p-1,1);
      omega_12 = beta; omega_22 = gamma + beta'*inv_Omega_11*beta;
      %% sample lambda_sq and nu
      rate = omega_12.^2/(2*tau_sq)+1./nu_12;
      lambda_sq_12 = 1./gamrnd(1,1./rate);    % random inv gamma with shape=1, rate=rate
      nu_12 = 1./gamrnd(1,1./(1+1./lambda_sq_12));    % random inv gamma with shape=1, rate=1+1/lambda_sq_12
      %% update Omega, Sigma, Lambda_sq, Nu
      Omega(i,ind) = omega_12; Omega(ind,i) = omega_12;
      Omega(i,i) = omega_22;
      temp = inv_Omega_11*beta;
      Sigma_11 = inv_Omega_11 + temp*temp'/gamma;
      sigma_12 = -temp/gamma; sigma_22 = 1/gamma;
      Sigma(ind,ind) = Sigma_11; Sigma(i,i) = sigma_22;
      Sigma(i,ind) = sigma_12; Sigma(ind,i) = sigma_12;
      Lambda_sq(i,ind) = lambda_sq_12; Lambda_sq(ind,i) = lambda_sq_12;
      Nu(i,ind) = nu_12; Nu(ind,i) = nu_12;
    end
    
%%% sample tau_sq and xi
    omega_vector = Omega(tril(true(size(Omega)),-1));
    lambda_sq_vector = Lambda_sq(tril(true(size(Lambda_sq)),-1));
    rate = 1/xi + sum(omega_vector.^2./(2*lambda_sq_vector));
    tau_sq = 1/gamrnd((p*(p-1)/2+1)/2, 1/rate);    % inv gamma w/ shape=(p*(p-1)/2+1)/2, rate=rate
    xi = 1/gamrnd(1,1/(1+1/tau_sq));    % inv gamma w/ shape=1, rate=1+1/tau_sq

%%% save Omega, lambda_sq, tau_sq
    if iter >burnin           
         omega_save(:,:,iter-burnin) = Omega;
         lambda_sq_save(:,iter-burnin) = lambda_sq_vector;
         tau_sq_save(iter-burnin) = tau_sq;
    end

end

end

