% This script is for Wassersten distancce based sulcal pattern registration 
% Section 2.3 and 2.4 of the ISBI paper `Sulcal Pattern Matching with 
% the Wasserstein Distance`
%
%  Run the following scripts sequentially:
%   SCRIPT1_dataPreprocess 
%   SCRIPT2_registration (this script)
%   SCRIPT3_validation

% (C) 2022 Zijian Chen, Moo K. Chung
%     University of Wisconsin-Madison
%
%  History: Feb 03, 2023 created by Chen
%           Feb 07, 2023 checked by Chung



%% We validate between the proposed dual formulation against the Hungarian algorithm

% We transform them back into matrix form
%Generate [-pi,pi],[0,pi]

xnum = 800; ynum = 400;

subj1 = reshape(diffusion_map_gridPoint_vec_1,[ynum xnum]);
subj2 = reshape(diffusion_map_gridPoint_vec_2,[ynum xnum]);

% displaying the two images
figure; subplot(1,2,1)
imagesc([-pi,pi],[0,pi],subj1);caxis([0,0.1])
set(gca,'YDir','normal')
title("Source")

subplot(1,2,2)
imagesc([-pi,pi],[0,pi],subj2);caxis([0,0.1])
set(gca,'YDir','normal')
title("Target")
 
% 'subj1' and 'subj2' will serve as the two input images for registration


%% registration


subj1 = mat2gray(subj1); subj2 = mat2gray(subj2);
[n1,m1] = meshgrid(0:xnum,0:ynum);
phi = ((m1).^2 + (n1).^2)/2;
[warpedImage,phi_final,phix,phiy]=wassersteinGD(subj2, subj1, phi,1, 200);


% 1) displaying the two images

figure; subplot(1,3,1)
imagesc([-pi,pi],[0,pi],subj2);         % target image
set(gca,'YDir','normal'); colorbar; 
title("Target")

subplot(1,3,2)
imagesc([-pi,pi],[0,pi],warpedImage);   % warped image
set(gca,'YDir','normal'); colorbar;  
title("Warped image to Target")

subplot(1,3,3)
imagesc([-pi,pi],[0,pi],subj2-warpedImage);   % residual (error map)
set(gca,'YDir','normal'); colorbar
title("Residual Error (Target - Warped Image)")

% figure  % overlapping
% hold on
% imagesc([-pi,pi],[0,pi],subj2);         % target image
% imagesc([-pi,pi],[0,pi],warpedImage);   % warped image
% set(gca,'YDir','normal')
% hold off


% 2) deformation field 

[x,y] = meshgrid(0.5:799.5,0.5:399.5);

s=5; % display arrows for every 5 points.

figure   % on top of source image
subplot(1,2,1); 
hold on; imagesc(subj1);        
quiver(x(1:s:end,1:s:end), y(1:s:end,1:s:end), ...
       x(1:s:end,1:s:end)-phix(1:s:end,1:s:end), y(1:s:end,1:s:end)-phiy(1:s:end,1:s:end), ...
       1, 'Color','r');
set(gca,'YDir','normal')
title("Displacment field on top of Source")
xlim([0 xnum]); ylim([0 ynum])



subplot(1,2,2)   % on top of target image
hold on; imagesc(subj2);         
quiver(x(1:s:end,1:s:end), y(1:s:end,1:s:end), ...
       x(1:s:end,1:s:end)-phix(1:s:end,1:s:end), y(1:s:end,1:s:end)-phiy(1:s:end,1:s:end), ...
       1, 'Color','r');
set(gca,'YDir','normal')
title("Displacment field on top of Target")
xlim([0 xnum]); ylim([0 ynum])



% 3) Wasserstein distance
% Computed by the Monges's formulation. The deformation U is the gradient
% of phi ('phix' and 'phiy').
DW = 0;
intLength = 2*pi/xnum;
for i = 1:ynum
    for j = 1:xnum
        DW = DW + ((x(i,j)-phix(i,j))^2 + (y(i,j)-phiy(i,j))^2) * subj1(i,j) * (intLength)^2;
    end
end
DW = sqrt(DW);




