function scurvevalue = scurve_project(mesh, scurve, gcurve);
%scurvevalue = scurve_project(mesh, scurve, gcurve);
%
% 
% scurve: tree data that will be assigned value -1
% gcurve: tree data that will be assigned value 1
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
%
% (C) 2020- Moo K. Chung
%
% mkchung@wisc.edu
% Department of Biostatistics and Medical Informatics
% University of Wisconsin-Madison
%
% Update history: 2020 August updated 
%                 2022 August 24 documented

nVertices = size(mesh.vertices,1);
scurvevalue= zeros(nVertices,1);

nCurve = size(scurve,1);

for i = 1: nCurve
    index = scurve{i};
    scurvevalue(index)=-1;  
end



nCurve = size(gcurve,1);

for i = 1: nCurve
    index = gcurve{i};
    scurvevalue(index)=1;  
end

