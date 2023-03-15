% This script is for Wassersten distancce based sulcal pattern registration 
% Section 3.1 Validation against the Hungarian Algorithm in
% 
% Chen, Z., Das, S., Chung, M.K. 2023, Sulcal Pattern Matching with the Wasserstein Distance, 
% International Symposium in Biomedcial Imaging (ISBI)
% https://github.com/laplcebeltrami/sulcaltree/blob/main/chen.2023.ISBI.pdf
%
% This script is downloaded from 
% https://github.com/laplcebeltrami/sulcaltree
%
%
%  Run the following scripts sequentially:
%   SCRIPT1_dataPreprocess 
%   SCRIPT2_registration 
%   SCRIPT3_validation (this script)
%
% (C) 2023 Zijian Chen, Moo K. Chung
%     University of Wisconsin-Madison
%
%  History: Feb 13, 2023 created by Chen
%           March 14, 2023 checked and edied by Chung
%
%
%

%% Generates Guaussian mixtures 'mu' and 'nu' with 'num_peaks' number of peaks. 
n = 100; %grid size
num_peaks=10; %number of peaks
sigma=0.001; %smoothness
[mu_peak, nu_peak, mu,nu] = Gaussian_randGen(n,num_peaks, sigma); %generate Gaussian mixtures
figure; subplot(1,2,1); imagesc(mu); subplot(1,2,2); imagesc(nu)

%% Wasserstin distance computation using gradient descent
[x1,y1] = meshgrid(0.5:n+0.5);
phi = (x1.^2 + y1.^2)/2; 
source = mat2gray(mu)+0.5*ones(n);
target = mat2gray(nu)+0.5*ones(n);
figure; subplot(1,2,1); imagesc(source); subplot(1,2,2); imagesc(target)
[~,~,phix,phiy]=wassersteinGD(target,source, phi,1, 200); %deformation field
DW_gd = wasserstein_distance(phix,phiy,mu) %Wasserstein distance. 

%% Hungarian algorithm
C = pdist2(mu_peak,nu_peak,"euclidean");
C = C.^2;
[assignment,~] = munkres(C);
tagetPt = [nu_peak(assignment,1),nu_peak(assignment,2)];
DW_hun = sqrt( 1/num_peaks* sum( vecnorm(tagetPt-mu_peak,2,2).^2  )  )

%% Releative error betwen the Hunagrain algorithm (baseline) and gradient descent
diff = (DW_hun - DW_gd)/DW_hun
%If we perform the above simulations 100 times and average the results, you
%obtain the results given in section 3.1. The relative error should be
%around 0.89 in this example. 


%% Generate Figure 4 in the paper
% If num_peaks=3 and sigma=0.01, you get similar picture

figure
subplot(1,2,1)
%display the displacement vector fields from the Huangarian algorithm
wassersteinH_display(mu_peak,nu_peak,tagetPt)

subplot(1,2,2)
%display the displacement vector fields from the gradient descent algorithm
wassersteinGD_display(mu,nu,phix,phiy)



