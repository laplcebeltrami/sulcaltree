function wassersteinGD_display(mu,nu,phix,phiy)
% function wassersteinGD_display(mu,nu,phix,phiy)
% 
% The function displays the displacement vector from the 
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
%           mu:   source image
%           nu:   target image
%         phix:   x-derivative of phi
%         phiy:   y-derivative of phi
%
% OUTPUT
%  warpedImage:   the deformed intensity matrix (of I1 towards I0)
%          phi:   the updated phi

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

n=size(mu,1);
[x,y] = meshgrid(linspace(0, 1, n));

s=4;
hold on;
imagesc(-nu, 'x',[0,1], 'y',[0,1],'AlphaData',0.9); %target
imagesc(mu, 'x',[0,1], 'y',[0,1],'AlphaData',0.7); %source
caxis([-max(nu(:)),max(nu(:))])
%colormap(rywb)
%colorbar
set(gca,'YDir','normal')
quiver(x(1:s:end,1:s:end), y(1:s:end,1:s:end), ...
       x(1:s:end,1:s:end)-phix(1:s:end,1:s:end)/n, y(1:s:end,1:s:end)-phiy(1:s:end,1:s:end)/n, ...
       0, 'Color','k',MaxHeadSize=0.05);
hold off
xlim([0,1]);ylim([0,1]);
title('Gradient descent')