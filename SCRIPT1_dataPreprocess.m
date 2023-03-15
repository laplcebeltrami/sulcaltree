% This script is for preparing data for registration
% Section 2.1 and  2.2 of  
% Chen, Z., Das, S., Chung, M.K. 2023, Sulcal Pattern Matching with the Wasserstein Distance, 
% International Symposium in Biomedcial Imaging (ISBI)
% https://github.com/laplcebeltrami/sulcaltree/blob/main/chen.2023.ISBI.pdf
%
% This script is downloaded from 
% https://github.com/laplcebeltrami/sulcaltree
%
%
%  Run the following scripts sequentially:
%   SCRIPT1_dataPreprocess (this script)
%   SCRIPT2_registration 
%   SCRIPT3_validation
%
% (C) 2022 Zijian Chen, Moo K. Chung
%     University of Wisconsin-Madison
%
%    History: Feb 1, 2023 created by Chen
%             Feb 2, 2023 checked by Chung
%             Feb 3, 2023 changed faces in box mesh to to FREESURFER
%                         TRIANGULATION by Chen


%% Section 2.1.

%spherical mesh of two subjects
sphereL1 = gifti('sphereL1.gii');
sphereL2 = gifti('sphereL2.gii');

%sulcal curve of two subjects
sulcalL1 = ReadBase('pialL1.scurve');
sulcalL2 = ReadBase('pialL2.scurve');


%scurve_project will assign value -1 to sulcal patterns, we need it to be +1 => take -scurve_project(...)
initialmapL1 = -scurve_project(sphereL1, sulcalL1,[]); 
initialmapL2 = -scurve_project(sphereL2, sulcalL2,[]);    



% projection of the sulcal (Fig.1)

%----- subject 1 -----%

% create mesh in the rectagular box [-pi,pi]x[0,pi] using the coordinate
% just obtained
sphL_v_cart = sphereL1.vertices; % the vertices are given in Cartesian coordianate system (x,y,z)
sphL_v_sph_1=[];
[sphL_v_sph_1(:,1),sphL_v_sph_1(:,2),~] = cart2sph(sphL_v_cart(:,1),sphL_v_cart(:,2),sphL_v_cart(:,3)); % turn it into spherical coordinate (theta,varphi,1)
sphL_v_sph_1(:,2)=sphL_v_sph_1(:,2)+pi/2;

box_mesh_1.vertices = sphL_v_sph_1;
box_mesh_1.faces=sphereL1.faces;

% adjust raw sulcal pattern to cut abnormal connection that cross the whole box 
newSulcalL1 = scurve_adjust(box_mesh_1, sulcalL1)';
% display projection
figure; scurve_display_2d(box_mesh_1, newSulcalL1,'r');

%----- subject 2 -----%
sphL_v_cart = sphereL2.vertices; % the vertices are given in Cartesian coordianate system (x,y,z)
sphL_v_sph_2=[];
[sphL_v_sph_2(:,1),sphL_v_sph_2(:,2),~] = cart2sph(sphL_v_cart(:,1),sphL_v_cart(:,2),sphL_v_cart(:,3)); % turn it into spherical coordinate (theta,varphi,1)
sphL_v_sph_2(:,2)=sphL_v_sph_2(:,2)+pi/2;

box_mesh_2.vertices = sphL_v_sph_2;
box_mesh_2.faces=sphereL2.faces;

newSulcalL2 = curve_adjust(box_mesh_2, sulcalL2)';
figure; scurve_display_2d(box_mesh_2, newSulcalL2,'r');


%% Section 2.2. 


% Smoothing and resampling (Fig.2)

%----- Diffusion smoothing -----%

deg = 50;
sigma = 0.001;

[~,betaL_WFS_1,diffusion_map_orig_vec_1]=WFS_R2_periodicDirichlet(initialmapL1,sphL_v_sph_1,deg,sigma);
[~,betaL_WFS_2,diffusion_map_orig_vec_2]=WFS_R2_periodicDirichlet(initialmapL2,sphL_v_sph_2,deg,sigma);
 

%-----  resampling -----%
% we need to resample the two diffusion maps to the same regular rectnagular grid
% create regular grids over [-pi,pi]x[0,pi]

xnum = 800; ynum = xnum/2;

x = linspace(-pi,pi,xnum)';
y = linspace(0,pi,ynum)';
regGrid=[];
regGrid(:,1) = repelem(x,ynum);
regGrid(:,2) = repmat(y,xnum,1);

% resampling (approximately 150s runtime for each subject)
diffusion_map_gridPoint_vec_1 = WFS_R2_resampled(regGrid,betaL_WFS_1,deg,sigma);
diffusion_map_gridPoint_vec_2 = WFS_R2_resampled(regGrid,betaL_WFS_2,deg,sigma);

% Displaying the resampled smoothed sulcal pattern for registration.

% subject 1
figure; subplot(1,2,1); scurve_display_2d(box_mesh_1, newSulcalL1,'r');
hold on
patch('Faces',delaunay(regGrid),'Vertices',regGrid,'FaceVertexCData', ...
        diffusion_map_gridPoint_vec_1,'FaceColor','interp','EdgeColor','none')
scurve_display_2d(box_mesh_1, newSulcalL1,'r');
caxis([0,0.1])
hold off


% subject 2
subplot(1,2,2); scurve_display_2d(box_mesh_2, newSulcalL2,'r');
hold on
patch('Faces',delaunay(regGrid),'Vertices',regGrid,'FaceVertexCData', ...
        diffusion_map_gridPoint_vec_2,'FaceColor','interp','EdgeColor','none')
scurve_display_2d(box_mesh_2, newSulcalL2,'r');
caxis([0,0.1])
hold off




