function [mu_peak, nu_peak, mu,nu] = Gaussian_randGen(n,num_peaks, sigma)
%function [mu_peak, nu_peak, mu,nu] = Gaussian_randGen(n,num_peaks, sigma)
%
% This function generates Guaussian mixtures 'mu' and 'nu' with 'num_peaks' number of
% peaks. 
%
% Input:
%                   n:  size of the background grid
%           num_peaks:  number of peaks in each of the two densities
%           sigma    :  variance in the covariance matrix 
%
% Output:
%     nu_peak, mu_peak:  peak points from in each of the two densities
%               nu, mu:  two densities to be generated
%
% generate peak point within the square [0.25,0.75]x[0.25,0.75]
%
% 
% This function is used in Section 3.1 Validation against the Hungarian Algorithm
% Chen, Z., Das, S., Chung, M.K. 2023, Sulcal Pattern Matching with the Wasserstein Distance, 
% International Symposium in Biomedcial Imaging (ISBI)
% https://github.com/laplcebeltrami/sulcaltree/blob/main/chen.2023.ISBI.pdf
%
% This function is downloaded from 
% https://github.com/laplcebeltrami/sulcaltree
%
%
%
% (C) 2023 Zijian Chen, Moo K. Chung
%     University of Wisconsin-Madison
%
%  History: Feb 13, 2023 created by Chen
%           March 14, 2023 checked and edied by Chung
%

mu_peak = 0.25 + 0.5 * rand([num_peaks,2]);  
nu_peak = 0.25 + 0.5 * rand([num_peaks,2]);


% generate mu and nu
Sigma = sigma*eye(2);

[x,y] = meshgrid(linspace(0,1,n));
mu=zeros(n^2,1); nu=zeros(n^2,1);

for i = 1:num_peaks
    mu = mu + mvnpdf([x(:) y(:)],mu_peak(i,:),Sigma);
    nu = nu + mvnpdf([x(:) y(:)],nu_peak(i,:),Sigma);
end


mu = reshape(mu,[n,n]);
mu = mu*n^2/sum(mu(:));

nu = reshape(nu,[n,n]);
nu = nu*n^2/sum(nu(:));

end