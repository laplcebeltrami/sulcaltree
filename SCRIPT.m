%SCRIPT generating sulcal tree patterns published in 
%
% Huang, S.-G., Lyu, I., Qiu, A., Chung, M.K. 2020. Fast polynomial 
% approximation of heat kernel convolution on manifolds and its application 
% to brain sulcal and gyral graph pattern analysis, IEEE Transactions on 
% Medical Imaging 39:2201-2212 
%
%subjects 130114 and 155938 are displayed
%------------
%subject 130114

%reading sphere for SPHARM modeling
sphereL1 = gifti('sphereL1.gii'); %spherical mesh corresponding to the left hemisphere of brain
figure; figure_wire(sphereL1); 

%Reading pial surface 
pialL1 = gifti('pialL1.gii'); 
figure; figure_wire(pialL1); 

%Reading white matter surface
whiteL1 = gifti('whiteL1.gii'); 
figure; figure_wire(whiteL1); 

%Reading sulcal trees defined on the pial surface
sulcalL1 = ReadBase('pialL1.scurve'); 

%Reading gyral trees defined on white matter surface
gyralL1 = ReadBase('whiteL1.gcurve'); 

%Put value -1 (heat sink) to sulcal curve and value 1 (heat source) to gyral curve
initialmap1 = scurve_project(whiteL1, sulcalL1,gyralL1); 


%Visualization of sulcal/gyral tree pattern on white matter surface as
%curves
figure; figure_trimesh(whiteL1,initialmap1,'rywb'); 
hold on; scurve_display(whiteL1, gyralL1,'r') %display gyral graph edges
hold on; scurve_display(whiteL1, sulcalL1, 'b') %display sulcal graphs
view([-90 0]);camlight
%print_pdf('sulcalcurve2')

%Visualization of sulcal/gyral tree pattern on white matter surface as
%graphs
figure; figure_trimesh(whiteL1,initialmap1,'rywb'); 
hold on; scurve_display(whiteL1, gyralL1,'r') %display gyral graph edges
hold on; scurve_displaynodes(whiteL1, gyralL1,'r') %display graph nodes
hold on; scurve_display(whiteL1, sulcalL1, 'b') %display sulcal graphs
hold on; scurve_displaynodes(whiteL1, sulcalL1,'b') %display graph nodes
view([-90 0]);camlight

%---------------
% Subject 155938

%reading sphere for SPHARM modeling
sphereL2 = gifti('sphereL2.gii'); %spherical mesh corresponding to the left hemisphere of brain
figure; figure_wire(sphereL2); 

%Reading pial surface 
pialL2 = gifti('pialL2.gii'); 
figure; figure_wire(pialL1); 

%Reading white matter surface 
whiteL2 = gifti('whiteL2.gii'); 
figure; figure_wire(whiteL1); 

%Reading sulcal trees defined on the pial surface
sulcalL2 = ReadBase('pialL2.scurve'); 

%Reading gyral trees defined on white matter surface
gyralL2 = ReadBase('whiteL2.gcurve'); 

%Put value -1 (heat sink) to sulcal curve and value 1 (heat source) to gyral curve
initialmap2 = scurve_project(whiteL2, sulcalL2,gyralL2); 

%Visualization of sulcal/gyral tree pattern on white matter surface as
%curves
figure; figure_trimesh(whiteL2,initialmap2,'rywb'); 
hold on; scurve_display(whiteL2, gyralL2,'r') %display gyral graph edges
hold on; scurve_display(whiteL2, sulcalL2, 'b') %display sulcal graphs
view([-90 0]);camlight

%Visualization of sulcal/gyral tree pattern on white matter surface as
%graphs
figure; figure_trimesh(whiteL2,initialmap2,'rywb'); 
hold on; scurve_display(whiteL2, gyralL2,'r') %display gyral graph edges
hold on; scurve_displaynodes(whiteL2, gyralL2,'r') %display graph nodes
hold on; scurve_display(whiteL2, sulcalL2, 'b') %display sulcal graphs
hold on; scurve_displaynodes(whiteL2, sulcalL2,'b') %display graph nodes
view([-90 0]);camlight


%---------------------------
%% Sulcal/gyral Trees on sphere

%%Visualization of sulcal/gyral tree pattern on sphere as curves

figure; subplot(1,2,1);
figure_trimesh(sphereL1,initialmap1,'rywb'); 
hold on; scurve_display(sphereL1, gyralL1,'r') %display gyral graph edges
hold on; scurve_display(sphereL1, sulcalL1, 'b') %display sulcal graphs
view([-90 0]);camlight; colorbar off

subplot(1,2,2); figure_trimesh(sphereL2,initialmap2,'rywb'); 
hold on; scurve_display(sphereL2, gyralL2,'r') %display gyral graph edges
hold on; scurve_display(sphereL2, sulcalL2, 'b') %display sulcal graphs
view([-90 0]);camlight; colorbar off

%----------------------
%heat kernel smoothing using SPHARM based on the study
%
%Chung, M.K., Dalton, K.M., Shen, L., L., Evans, A.C., Davidson, R.J. 2007. 
%Weighted Fourier series representation and its application to quantifying 
%the amount of gray matter. Special Issue of  IEEE Transactions on Medical Imaging, 
%on Computational Neuroanatomy. 26:566-581. 
%
%More codes
%https://pages.stat.wisc.edu/~mchung/softwares/weighted-SPHARM/weighted-SPHARM.html


[surfL1, fourierL1]=SPHARMsmooth_signal(initialmap1,sphereL1,40,0.0001); %heat kernel smoothing with bandwidth 0.001
[surfL2, fourierL2]=SPHARMsmooth_signal(initialmap2,sphereL2,40,0.0001); %heat kernel smoothing with bandwidth 0.001


figure; subplot(1,2,1);
figure_trimesh(sphereL1,surfL1,'rywb'); 
hold on; scurve_display(sphereL1, gyralL1,'r') %display gyral graph edges
hold on; scurve_display(sphereL1, sulcalL1, 'b') %display sulcal graphs
view([-90 0]);camlight; colorbar off; alpha(1)

subplot(1,2,2);
figure_trimesh(sphereL2,surfL2,'rywb'); 
hold on; scurve_display(sphereL2, gyralL2,'r') %display gyral graph edges
hold on; scurve_display(sphereL2, sulcalL2, 'b') %display sulcal graphs
view([-90 0]);camlight; colorbar off; alpha(1)

%Matching the sulcal/gyral trees might be computationally dtoo demanding. So
%we will match their heat kernel smoothing resampled at lower traingle
%resolutions.

%new spherical mesh with 
%    vertices: [2562×3 double]
%       faces: [5120×3 double]

sphere = sphere_tri('ico',4,1);
resampled=SPHARMrepresent(sphere, fourierL1,40, 0.0001); %resample heat kernel smoothing at 5120 triangle resolution

%Comparision of heat kernel smoothing in different resolutions
figure; subplot(1,2,1);
figure_trimesh(sphereL1,surfL1,'rywb'); caxis([-0.05 0.05])
title('mesh with 331078 triangles'); colorbar off;
subplot(1,2,2); figure_trimesh(sphere,resampled,'rywb'); caxis([-0.05 0.05])
title('mesh with 5120 triangles'); colorbar off;



