function [output, fourier]=SPHARMsmooth_signal(input,sphere, L,sigma)
%[output, fourier]=SPHARMsmooth_signal(input,sphere, L,sigma)
%
% Please read reference [1] to understand the input prameters, the notations 
% and the algorithm. If you use this code, please reference [1].
%
%
% input           : An input surface mesh format identical to the MATLAB format
% output          : An output surface mesh format identical to the MATLAB format
%                 
%  sphere      : Spherical mesh where the representation is constructed
%
%    L             : The maximal degree of the weighted-SPHARM representation.
%                      Read the paper below to find it optimally.
%
%  sigma        : bandwith of weighted-SPHARM representation
%                    It is the bandwidth of heat kernel on unit sphere.
%                    When sigma=0, it is the traditional SPHARM representation.
%                    range beween 0.0001 and 0.01 will be sufficient for cortical surfaces.
%
%
% surf_smooth     : The weighted-SPHARM result.
%
% fourier         : The estimated SPHARM coefficients (Fourier coefficients) given as a structured array
%                   containg coeff.x, coeff.y, coeff.z
%                   coeff.x is the SPHARM coefficients of x-cooridinate given as (L+1) by (2*L+1) matrix
%                   coeff.x(3,:) is the 2nd degree SPHARM coefficients of all order.
%
% NOTE: Modified from SPHARMsmooth2.m and handels multiple surfaces. This version is 
% different from SPHARMsmooth.m in http://www.stat.wisc.edu/~mchung/softwares/weighted-SPHARM/SPHARMsmooth.m
% SPHARMsmooth.m is based on a fixed connectivity across all surface.It
% also fixes the infinity error in Legendre function caused by new MATLAB
% version.
%
% Reference:
% [1] Chung, M.K., Dalton, K.M., Shen, L., L., Evans, A.C., Davidson, R.J. 2007.
% Weighted Fourier series representation and its application to quantifying
% the amount of gray matter. IEEE Transactions on Medical Imaging, 26:566-581.
%
%
%
% (C) 2006- Moo K. Chung, Shih-Gu Huang
%     Department of Biostatistics and Medical Informatics
%     Waisman Laboratory for Brain Imaging
%     University of Wisconsin-Maison
%
% email://mkchung@wisc.edu
% http://www.stat.wisc.edu/softwares/weighted-SPHARM/weighted-SPHARM.html
%
%
% Update history
%     Sept 19 2006; July 5, 2007; January 22, 2008; July 18, 2010
%     2019 May 28 Handles multiple surfaces with different mesh topology. Infinity error 
%                 in Legendre function caused by new MATLAB version was correcetd. By Shih-Gu Huang.
%     2019 Oct 07 Documents updated. 
%
%---------------------------------------------------------------------------

% coord=surf.vertices;
% n_vertex = size(coord,1);   % the number of vertices in a mesh.

[theta varphi]=EULERangles(sphere);

% x=coord(:,1);
% y=coord(:,2);
% z=coord(:,3);

x=input;
n_vertex = size(x,1); 
numx=size(x,2);

%INITIALIZATION

% xestiamte is the weighted-SPHARM of x-coordinate.
xestimate=zeros(n_vertex,numx);
% yestimate=zeros(n_vertex,1);
% zestimate=zeros(n_vertex,1);

% betax is the Fourier coefficients of x-coordinate.
betax=zeros(L+1,2*L+1,numx);
% betay=zeros(L+1,2*L+1);
% betaz=zeros(L+1,2*L+1);


%0-TH DEGREE. 
%Step 2 in the iterative resiual fitting (IRF) algorithm. See reference [1].

Y=Y_l(0,theta',varphi')';
Ycommon=inv(Y'*Y)*Y';

betal=Ycommon*x;
betax(1,1,:)=betal;
xsmooth=Y*betal;
xestimate=xsmooth;

% betal=Ycommon*y;
% betay(1,1)=betal';
% ysmooth=Y*betal;
% yestimate=ysmooth;
% 
% betal=Ycommon*z;
% betaz(1,1)=betal';
% zsmooth=Y*betal;
% zestimate=zsmooth;


%l-TH DEGREE ITERATION


for l=1:L
 
    % Step 4: residual. See the paper for detail
    x_j = x-xestimate;
%     y_j = y-yestimate;
%     z_j = z-zestimate;
 
    Y=Y_l(l,theta',varphi')';
    % real(Y_l) gives negative harmonics imag(Y_l) gives positive harmonics
   
    Y=[real(Y) imag(Y(:,2:(l+1)))];
    Y=real(Y);
    Ycommon=inv(Y'*Y)*Y';  %changed from Ycommon=inv(Y'*Y)*Y';

   %Y(1,:)
   
    % Step 5: refitting the residual. See the paper for detail
    betal=Ycommon*x_j;
    betax(l+1,1:2*l+1,:)=betal;
    xsmooth=Y*betal;
    xestimate=xestimate + exp(-l*(l+1)*sigma)*xsmooth;

%     betal=Ycommon*y_j;
%     betay(l+1,1:2*l+1)=betal';
%     ysmooth=Y*betal;
%     yestimate=yestimate + exp(-l*(l+1)*sigma)*ysmooth;
% 
%     betal=Ycommon*z_j;
%     betaz(l+1,1:2*l+1)=betal';
%     zsmooth=Y*betal;
%     zestimate=zestimate + exp(-l*(l+1)*sigma)*zsmooth;
    
% animation. remove next 4 lines if you don't want animation
%     temp=[xestimate; yestimate; zestimate];
%     surf_smooth.vertices=squeeze(reshape(temp,n_vertex,3));
%     surf_smooth.faces=surf.faces;
%     M(l) = getframe; 
%     figure_wire(surf_smooth,'yellow','white');

    
end;


%output the results in a proper shape
% temp=[xestimate; yestimate; zestimate];
% surf_smooth.vertices=squeeze(reshape(temp,n_vertex,3));
% surf_smooth.faces=surf.faces;

% fourier.x=betax;
% fourier.y=betay;
% fourier.z=betaz;
fourier=betax;
output=xestimate;


%---------------------------------------------------------------
function [theta,varphi]=EULERangles(surf);

n_vertex=size(surf.vertices,1);
c=mean(surf.vertices);  %mass center
surf.vertices=surf.vertices-kron(ones(n_vertex,1),c);  % translation

[theta,varphi,r] = cart2sph(surf.vertices(:,1),surf.vertices(:,2),surf.vertices(:,3));

% MATLAB coordinate systems are different from the convention used in the TMI paper.
temp = theta;
theta = pi/2 - varphi;
varphi = pi + temp;

%figure_wire(surf,'yellow')

%----------------------------------
function Y_l=Y_l(l,theta,varphi)
%Infinity error in legendre function is corrected and modified by Huang

% computes spherical harmonics of degree l.
% sz=length(theta);

% m=0:l;
% CLM=[];
exp_i_m=[];
% sign_m=[];
% SIGNM=[];
Pn=[];

% for k = 0:(2*l)
%     fact(k+1) = factorial(k);
% end
% clm = sqrt(((2*l+1)/(2*pi))*(fact(l-abs(m)+1)./fact(l+abs(m)+1))); % this clm will be close to 0 when l is large.
% clm=(-1).^m*sqrt(1/pi); % clm needs to be modified because Pn is normalized 
% CLM=kron(ones(1,sz),clm');    % CLM.*SIGNM=sqrt(1/pi)

for k = 0:l
    exp_i_m(k+1,:)= exp(i*k*varphi);
%     sign_m(k+1) = (-1)^k;
end

% c=sqrt(2);    % not used
exp_i_m(1,:)=exp_i_m(1,:)/sqrt(2);

% SIGNM=kron(ones(1,sz),sign_m'); % CLM.*SIGNM=sqrt(1/pi)
% Pn=legendre(l,cos(theta));     %this Pn goes to Inf when l is large
Pn=legendre(l,cos(theta),'norm'); % Pn is normalized to avoid Inf, and thus clm needs to be modified
% Y_l=CLM.*SIGNM.*Pn.*exp_i_m;  % CLM.*SIGNM=sqrt(1/pi)
Y_l=sqrt(1/pi)*Pn.*exp_i_m;
