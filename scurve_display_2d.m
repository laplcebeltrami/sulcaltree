function scurve_display_2d(mesh, scurve, c)
% function scurve_display_2d(mesh, scurve, c)
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
% (C) 2020- Zijian Chen, Ilwoo Lyu, Moo K. Chung
% mkchung@wisc.edu
% Department of Biostatistics and Medical Informatics
% University of Wisconsin-Madison
%
%
% Update history: 2020 August updated 
%                 2022 August 24 documented
%                 2023 Mar 15 Chung commented

nCurve = size(scurve,1);


for i = 1: nCurve
    index = scurve{i};
    curv = mesh.vertices(index, :);

    breakpoint = [];

    for j = 1:size(curv,1)-1
        if abs(curv(j,1)-curv(j+1,1))>=4
            breakpoint = [breakpoint;j];
        end
    end

    if isempty(breakpoint)
          hold on; plot(curv(:, 1), curv(:, 2),  c, 'LineWidth', 1);
    else
        breakpoint = [breakpoint;size(curv,1)];
        hold on; 
        for k = 1:size(breakpoint,1)-1
            plot(curv(1:breakpoint(k), 1), curv(1:breakpoint(k), 2),  c, 'LineWidth', 1);
            plot(curv(breakpoint(k)+1:breakpoint(k+1), 1), curv(breakpoint(k)+1:breakpoint(k+1), 2),  c, 'LineWidth', 1);
        end
    end
end

xlim([-pi pi]); ylim([0 pi]);