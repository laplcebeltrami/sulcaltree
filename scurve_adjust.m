function newcurve = scurve_adjust(mesh, scurve)
%function newcurve = scurve_adjust(mesh, scurve)
%
% adjust raw sulcal pattern to cut abnormal connection that cross the whole
% box (since we cannot display the periodicity)
% NOTE: This is for visualization only. Do not use in further processing.
%
%
% The function is written for  
% Chen, Z., Das, S., Chung, M.K. 2023, Sulcal Pattern Matching with the Wasserstein Distance, 
% International Symposium in Biomedcial Imaging (ISBI)
% https://github.com/laplcebeltrami/sulcaltree/blob/main/chen.2023.ISBI.pdf
%
%
% (C) 2022 Zijian Chen, Moo K. Chung
%     University of Wisconsin-Madison
%
%  History: Feb 03, 2023 created by Chen
%           Feb 07, 2023 checked by Chung


nCurve = size(scurve,1);
newcurve = {};

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
          newcurve{end+1} = scurve{i};
    else
        breakpoint = [breakpoint;size(curv,1)];
        for k = 1:size(breakpoint,1)-1
            newcurve{end+1} = index(1:breakpoint(k));
            newcurve{end+1} = index(breakpoint(k)+1:breakpoint(k+1));
        end
    end
end


end

