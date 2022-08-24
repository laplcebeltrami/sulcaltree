function resampled=SPHARMrepresent(surf,coeff,k,sigma)
% SPHARMsurf=SPHARMrepresent(surf,coeff,k,sigma)
%
% Reconstruct the k-th degree Weighted-SPHARM from given SPHARM coefficients. 
%  
% INPUT
% surf               : mesh for a sphere where the SPHARM is constructed.
%                
% coeff   : The estimated SPHARM coefficients (Fourier coefficients) given as a structured array
%                        coeff is the SPHARM coefficients given as (L+1)*(2L+1)*n_subject 
%                        coeff(3,:,10) is the 2nd degree SPHARM coefficients of all order for the 10th subject.  
%
% sigma            :  Smoothing bandwidth
%
% OUTPUT
% resampled data using coeff. 
%
%
% (C) 2007 Moo K. Chung 
%     Department of Biostatistics and Medical Informatics
%     Waisman Laboratory for Brain Imaging
%     University of Wisconsin-Maison
%  
% email://mkchung@wisc.edu
% http://www.stat.wisc.edu/softwares/weighted-SPHARM/weighted-SPHARM.html
%
% If you use this code, please reference the following paper. 
% You need to read the paper to understand the notations and the algorithm.
%
% Chung, M.K., Dalton, K.M., Shen, L., L., Evans, A.C., Davidson, R.J. 2007. 
% Weighted Fourier series representation and its application to quantifying 
% the amount of gray matter. IEEE Transactions on Medical Imaging, 26:566-581.
%
%
% July 6, 2007 created 
% Jan  9, 2008 modified
% August 24, 2022 modified further for scalar data
%----------------------------------------------------------------------------------------

n_subject = size(coeff,1);
n_vertex = size(surf.vertices,1);
[theta varphi]=EULERangles(surf);

% Initialization for coordinates
xcoord=zeros(n_vertex,1);
ycoord=zeros(n_vertex,1);
zcoord=zeros(n_vertex,1);

% 0-th degree construction

% real(Y_l) gives negative harmonics imag(Y_l) gives positive harmonics
% reading Y_lm and constructing design matrix. See the paper.
Y=Y_l(0,theta',varphi')';

betal=coeff(1,1);
xcoord=Y*betal;

% Iteratively add each degree up to the k-th degree.    

for l=1:k


    Y=Y_l(l,theta',varphi')';
    % real(Y_l) gives negative harmonics imag(Y_l) gives positive harmonics
    Y=[real(Y)   imag(Y(:,2:(l+1)))];
    
    betal=coeff(l+1,1:2*l+1);
    xsmooth=Y*betal';
    xcoord=xcoord+ exp(-l*(l+1)*sigma)*xsmooth;

end

%output the results in a proper shape

resampled=xcoord;


%---------------------------------------------------------------
function [theta,varphi]=EULERangles(surf);

n_vertex=size(surf.vertices,1);
c=mean(surf.vertices);  %mass center
surf.vertices=surf.vertices-kron(ones(n_vertex,1),c);  % translation

[theta,varphi,r] = cart2sph(surf.vertices(:,1),surf.vertices(:,2),surf.vertices(:,3));

% MATLAB coordinate systems are different from the convention used in the
% TMI paper.
temp = theta;
theta = pi/2 - varphi;
varphi = pi + temp;

%figure_wire(surf,'yellow')

%-----------------------------------------------------------------
function Y_l=Y_l(l,theta,varphi)
% computes spherical harmonics of degree l.
sz=length(theta);

m=0:l;
CLM=[];
exp_i_m=[];
sign_m=[];
SIGNM=[];
Pn=[];

for k = 0:(2*l)
    fact(k+1) = factorial(k);
end
clm = sqrt(((2*l+1)/(2*pi))*(fact(l-abs(m)+1)./fact(l+abs(m)+1)));
CLM=kron(ones(1,sz),clm');

for k = 0:l
    exp_i_m(k+1,:)= exp(i*k*varphi);
    sign_m(k+1) = (-1)^k;
end
exp_i_m(1,:)=exp_i_m(1,:)/sqrt(2);

SIGNM=kron(ones(1,sz),sign_m');
Pn=legendre(l,cos(theta));
Y_l=CLM.*SIGNM.*Pn.*exp_i_m;


