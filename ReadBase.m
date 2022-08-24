function base = ReadBase(curve_name)
%function base = ReadBase(curve_name)
%
%The function is created for paper
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
% (C) 2020 Ilwoo Lyu
%
% mkchung@wisc.edu
% Department of Biostatistics and Medical Informatics
% University of Wisconsin-Madison
%
% Update history: 2020 August updated 
%                 2022 August 24 documented


line = textread(curve_name, '%s', 'delimiter', '\n');
n = size(line, 1);
base = cell(n, 1);
for i = 1: size(line, 1)
    base{i} = uint32(str2num(line{i}))' + 1;
end
