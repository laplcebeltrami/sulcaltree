function DW = wasserstein_distance(phix, phiy, source)
%function DW = wasserstein_distance(phix, phiy, source)
%
% The function computes the Wasserstein distance using the deformation with
% the source image
%
%
% INPUT
%         phix  :   x-derivative of phi from wassersteinGD.m 
%         phiy  :   y-derivative of phi from wassersteinGD.m
%         source: source image
%
% OUTPUT
%         DW  :   Wasserstein distance
%
%
% The function is written for  
% Chen, Z., Das, S., Chung, M.K. 2023, Sulcal Pattern Matching with the Wasserstein Distance, 
% International Symposium in Biomedcial Imaging (ISBI)
% https://github.com/laplcebeltrami/sulcaltree/blob/main/chen.2023.ISBI.pdf
%
%
% This function is downloaded from 
% https://github.com/laplcebeltrami/sulcaltree
%   
%
%
% (C) 2022 Zijian Chen, Moo K. Chung
%     University of Wisconsin-Madison
%
%  History: Feb 03, 2023 created by Chen
%           Feb 07, 2023 checked by Chung


[m,n] = size(phix);
[x,y] = meshgrid(1:n,1:m);

DW = 0;

for i = 1:m
    for j = 1:n
        DW = DW + ( (x(i,j)/n-phix(i,j)/n)^2 + (y(i,j)/n-phiy(i,j)/n)^2 ) * source(i,j)  * 1/n * 1/m;
    end
end
DW = sqrt(DW);


end
