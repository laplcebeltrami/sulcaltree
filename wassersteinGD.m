function [warpedImage,phi,phix,phiy,detD2] = wassersteinGD(I0,I1,phi,alpha,maxIter)
% function [warpedImage,phi,phi_dx_center,phi_dy_center,DW] = wassersteinGD(I0,I1,phi,alpha,maxIter)
% 
% The function solves the Kantorovich's dual problem of optimal
% transportation via gradient descent given in
% 
% Chartrand, R., Wohlberg, B., Vixie, K., Bollt, E.: A gradient descent 
%   solution to the monge- kantorovich problem. Applied Mathematical 
%   Sciences 3, 1071–1080 (2009)
% The codes for implementing the algorithm were originally provided by 
%   Pascal Getreuer (https://getreuer.info) in python.
%
%
% The function is written for  
% Chen, Z., Das, S., Chung, M.K. 2023, Sulcal Pattern Matching with the Wasserstein Distance, 
% International Symposium in Biomedcial Imaging (ISBI)
% https://github.com/laplcebeltrami/sulcaltree/blob/main/chen.2023.ISBI.pdf
%
%
% INPUT
%           I0:   target intensity matrix
%           I1:   original intensity matrix (to be deformed to I1)
%          phi:   initial map. One can take:
%                       [n1,m1] = meshgrid(0:m,0:m); 
%                       phi = ((m1-0.5).^2 + (n1-0.5).^2)/2;
%                 where m is the size of I0
%        alpha:   learning rate
%      maxIter:   maximum iteration times
%
% OUTPUT
%  warpedImage:   the deformed intensity matrix (of I1 towards I0)
%          phi:   the updated phi
%         phix:   x-derivative of phi
%         phiy:   y-derivative of phi
%
%
% This function is downloaded from 
% https://github.com/laplcebeltrami/sulcaltree
%   
%
%
%
%
% (C) 2022 Zijian Chen, Moo K. Chung
%     University of Wisconsin-Madison
%
%  History: Feb 03, 2023 created by Chen
%           Feb 07, 2023 checked by Chung



DW = zeros(maxIter,1);

[N,M]=size(I1); 


for i = 1 : maxIter


    phi_dx_vertices = phi(:,2:M+1)-phi(:,1:M);
    phi_dx_center = (phi_dx_vertices(1:end-1,:)+phi_dx_vertices(2:end,:))/2;

    phi_dy_vertices = phi(2:N+1,:)-phi(1:N,:);
    phi_dy_center = (phi_dy_vertices(:,1:end-1)+phi_dy_vertices(:,2:end))/2;

    phi_dx_center_r = [phi_dx_center(:,2:end),phi_dx_center(:,end)];
    phi_dx_center_l = [phi_dx_center(:,1),phi_dx_center(:,1:end-1)];
    phi_dx_center_d = [phi_dx_center(2:end,:);phi_dx_center(end,:)];
    phi_dx_center_u = [phi_dx_center(1,:);phi_dx_center(1:end-1,:)];

    phi_dy_center_d = [phi_dy_center(2:end,:);phi_dy_center(end,:)];
    phi_dy_center_u = [phi_dy_center(1,:);phi_dy_center(1:end-1,:)];

    phi_dxdx_center = (phi_dx_center_r-phi_dx_center_l)/2;
    phi_dxdy_center = (phi_dx_center_d-phi_dx_center_u)/2;
    phi_dydy_center = (phi_dy_center_d-phi_dy_center_u)/2;

    detD2 = phi_dxdx_center.*phi_dydy_center-phi_dxdy_center.*phi_dxdy_center;


    I1_of_nablaPhi = bilinear_interp(I1,phi_dx_center,phi_dy_center);

    warpedImage = detD2.*I1_of_nablaPhi;

    % update phi

    diff_I = I0 - warpedImage; 

    diff_I_u = [diff_I(1,:);diff_I];
    diff_I_ur = [diff_I_u,diff_I_u(:,end)];diff_I_ul = [diff_I_u(:,1),diff_I_u];

    diff_I_d = [diff_I;diff_I(end,:)];
    diff_I_dr = [diff_I_d,diff_I_d(:,end)];diff_I_dl = [diff_I_d(:,1),diff_I_d];

    diff_I_ver = (diff_I_ur+diff_I_ul+diff_I_dr+diff_I_dl)/4;

    phi = phi - alpha * diff_I_ver;

end

%%%%% for output
phi_dx_vertices = phi(:,2:M+1)-phi(:,1:M);
phix = (phi_dx_vertices(1:end-1,:)+phi_dx_vertices(2:end,:))/2;

phi_dy_vertices = phi(2:N+1,:)-phi(1:N,:);
phiy = (phi_dy_vertices(:,1:end-1)+phi_dy_vertices(:,2:end))/2;

%%%%%%

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function f_interp = bilinear_interp(f,x,y)

% This function is for bilinear interpretation of f at 4 surrounding points 
% whose x coordinates are given in 'x' and  y coordinates are given in 'y'

[N,M]=size(f);

x1 = floor(x);
x1(x1<=0)=0; x1(x1>=M-2)=M-2;

y1 = floor(y);
y1(y1<=0)=0; y1(y1>=N-2)=N-2;

wx = x - x1; wy = (y-y1);

% images in MATLAB have y as the first dimension and x as the second,
% so are accessed as I(y, x). For detailed description, see:
% https://www.mathworks.com/help/images/image-coordinate-systems.html

g11 = matInd(f,y1+1,x1+1);
g21 = matInd(f,y1+1,x1+1+1);
g12 = matInd(f,y1+1+1,x1+1);
g22 = matInd(f,y1+1+1,x1+1+1);

f_interp = (1-wy).*((1-wx).*g11+wx.*g21) + ...
    wy .*((1-wx).*g12+wx.*g22);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = matInd(f,mat1,mat2)

% please make sure mat1 and mat2 have the same size
[n,m]=size(mat1);

g = zeros(n,m);

for i = 1:n
    for j = 1:m
        idx = mat1(i,j);
        idy = mat2(i,j);
        g(i,j) = f(idx,idy);
    end
end

end



