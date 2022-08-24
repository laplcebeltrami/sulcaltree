function sulcalcurve_display(mesh, scurve, c)
% function sulcalcurve_display(mesh, scurve, c)
% 
% Displays the edges of scurve on mesh with color c. 
%
% INPUT
%
% mesh   : surface mesh in MATLAB format
% scurve : Collection of curves with coordinates. scurve{i} is the i-th
%          curve data
% c      : color of curve
%
%
% The function is generated for study
%
% Huang, S.-G., Lyu, I., Qiu, A., Chung, M.K. 2020. Fast polynomial 
% approximation of heat kernel convolution on manifolds and its application 
% to brain sulcal and gyral graph pattern analysis, IEEE Transactions on 
% Medical Imaging 39:2201-2212 
% https://pages.stat.wisc.edu/~mchung/papers/huang.2020.TMI.pdf
%
% The code is downloaded from 
% https://github.com/laplcebeltrami/sulcaltree
% If you are using the code, please reference the above paper
%
% (C) 2020- Moo K. Chung, Ilwoo Lyu 
%
% mkchung@wisc.edu
% Department of Biostatistics and Medical Informatics
% University of Wisconsin-Madison
%
%
% Update history: 2020 August updated 
%                 2022 August 24 documented

nCurve = size(scurve,1);


for i = 1: nCurve
    index = scurve{i};
    curv = mesh.vertices(index, :);
    
    hold on; plot3(curv(:, 1), curv(:, 2), curv(:, 3),  c, 'LineWidth', 1);
end


% set the final visualization options
shading interp;
axis vis3d;
axis off;