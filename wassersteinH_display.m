function wassersteinH_display(mu_peak,nu_peak,tagetPt)
% function wassersteinH_display(mu_peak,nu_peak,tagetPt)
% 
% The function displays the displacement vector from the Huangarian
% algoirthm. The function is written for  
% Chen, Z., Das, S., Chung, M.K. 2023, Sulcal Pattern Matching with the Wasserstein Distance, 
% International Symposium in Biomedcial Imaging (ISBI)
% https://github.com/laplcebeltrami/sulcaltree/blob/main/chen.2023.ISBI.pdf
%
%
% INPUT
%           mu_peak:   source scatter points
%           nu_peak:   target scatter points
%           tagetPt:  assignment from the Hungarian algorithm
%
% OUTPUT
%  warpedImage:   the deformed intensity matrix (of I1 towards I0)
%          phi:   the updated phi
%
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

hold on
scatter(mu_peak(:,1),mu_peak(:,2),'.r',SizeData = 300);
scatter(nu_peak(:,1),nu_peak(:,2),'.b',SizeData = 300);
set(gca,'YDir','normal')
quiver(mu_peak(:,1),mu_peak(:,2), ...
    tagetPt(:,1)-mu_peak(:,1),tagetPt(:,2)-mu_peak(:,2), ...
    0,MaxHeadSize=1,LineWidth=2);
hold off
xlim([0,1]);ylim([0,1]);
title('Hungarian algoirthm')