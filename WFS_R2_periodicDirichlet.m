function [Phi,beta,est] = WFS_R2_periodicDirichlet(signal,vertices,deg,sigma)
%function [Phi,beta,est] = WFS_R2_periodicDirichlet(signal,vertices,deg,sigma)
% 
% The funciton performs heat kernel smoothing explained in
%
% Chen, Z., Das, S., Chung, M.K. 2023, Sulcal Pattern Matching with the Wasserstein Distance, 
% International Symposium in Biomedcial Imaging (ISBI)
% https://github.com/laplcebeltrami/sulcaltree/blob/main/chen.2023.ISBI.pdf
%
%
% The code is downloaded from 
% https://github.com/laplcebeltrami/sulcaltree
% If you are using the code, please reference the above paper
%
% (C) 2020- Zijian Chen, Moo K. Chung
% mkchung@wisc.edu
% Department of Biostatistics and Medical Informatics
% University of Wisconsin-Madison
%
%
% Update history:
%                 2023 Mar 15 Chung commented



L1=max(vertices(:,1)); L2=max(vertices(:,2));
Phi = []; lambda=[];

for j= 0:deg
    for k=1:deg
        psi_jk=eigenfunction_2dbox_perdir(vertices,j,k,L1,L2);
        Phi = [Phi,psi_jk];
        lambda = [lambda,(j*pi/L1)^2 + (k*pi/L2)^2,(j*pi/L1)^2 + (k*pi/L2)^2];
    end
end

lambda = lambda(:,any(Phi));
Phi=Phi(:,any(Phi));

beta = inv(Phi'*Phi)*Phi'*signal;

coef = exp(-sigma.*lambda);

est = coef.*Phi*beta;

end

function phi_jk=eigenfunction_2dbox_perdir(points,j,k,L1,L2)
    
    x=points(:,1);
    y=points(:,2);

    phi_jk1=cos(pi*j*x/L1).*sin(pi*k*y/L2);
    phi_jk2=sin(pi*j*x/L1).*sin(pi*k*y/L2);
    
    if j~=0 && k~=0
        phi_jk1 = 2/(L1*L2)*phi_jk1;
        phi_jk2 = 2/(L1*L2)*phi_jk2;
    elseif j==0 && k~=0
        phi_jk1 = 2/(L1*L2*2)*phi_jk1;
        phi_jk2 = 2/(L1*L2)*phi_jk2;
    end

    %lambda= (j*pi/L1)^2 + (k*pi/L2)^2; %eigenvalue

    %phi_jk=exp(-sigma*lambda)*[phi_jk1, phi_jk2];
    phi_jk=[phi_jk1, phi_jk2];

end

