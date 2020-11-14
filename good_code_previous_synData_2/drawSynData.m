% Example usage:3
clear all;
close all;
clc;
% mex CC=gcc-4.3 LD=g++-4.3 -O mexTopofix3d_synData.cpp
mex CC=gcc-4.3 LD=g++-4.3 -I/usr/local/Cellar/boost/1.66.0/include/ -O mexTopofix3d_good_chord.cpp
setDatapath;

NdataSet = 3;
selectDataSet = 3;  % 3->shrink.dz=58:good, 2:?, 1:?.
%% define the cut positions for different dataset.
% cutPosition = zeros(NdataSet,3,2);
% % Xudong: make sure the range(x0~x1, y0~y1, z0~z1) cover the mV and pM, 
% % or the program will crash.
% cutPosition(1,:,:) = [181,286;260,378;120,257]; %dz=138
% cutPosition(2,:,:) = [188,316;314,453;168,273];   %dz=106
% % cutPosition(2,:,:) = [188,200;314,350;168,200];   %shrink for checking.dz=32
% % cutPosition(3,:,:) = [157,300;280,411;170,247];   %shrink.dz=78
% cutPosition(3,:,:) = [180,286;280,400;170,227];   %shrink.dz=58, just find one chord.
% % cutPosition(3,:,:) = [147,316;270,411;163,258];   %dz=96
cutPosition(1,:,:) = [0,512;0,512;0,320];
cutPosition(2,:,:) = [0,512;0,512;0,320];
cutPosition(3,:,:) = [0,512;0,512;0,320];
%% obtian the coordinate range for different dataset.
Delta_X = cutPosition(selectDataSet,1,2) - cutPosition(selectDataSet,1,1);
Delta_Y = cutPosition(selectDataSet,2,2) - cutPosition(selectDataSet,2,1);
Delta_Z = cutPosition(selectDataSet,3,2) -  cutPosition(selectDataSet,3,1);

%% read 3D img, header, and mV&pM.
DataSet = cell(3,NdataSet);
DataSet(1,1:NdataSet) = {'CT4DStudy_1_phase_3.vol','CT4DStudy_1.hdr','Study01.20170830.07.6.14.1601.mat'};
DataSet(2,1:NdataSet) = {'CT4DStudy_2_phase_3.vol','CT4DStudy_2.hdr','Study02.20170911.05.3.16.641.mat'};
DataSet(3,1:NdataSet) = {'CT4DStudy_3_phase_3.vol','CT4DStudy_3.hdr','Study03.20170906.02.3.6.641.mat'};

% vol_name = fullfile(cd, DataSet{selectDataSet,1});
% hdr_name = fullfile(cd, DataSet{selectDataSet,2});
% mat_name = fullfile(cd, DataSet{selectDataSet,3});
vol_name = DataSet{selectDataSet,1};
hdr_name = DataSet{selectDataSet,2};
mat_name = DataSet{selectDataSet,3};
vol = load(vol_name, '-mat');
header = load(hdr_name, '-mat'); 
mVpM = load(mat_name, '-mat'); 

%%In original coord, draw mV and tips of pM and output VTK files.
% post = patch('vertices',mVpM.FC{1,1}.vertices,'faces',mVpM.FC{1,1}.faces);
% % post = patch('faces',mVpM.FC{1,T}.faces,'vertices',mVpM.FC{1,T}.vertices);
% set(post, 'FaceColor', 'red', 'EdgeColor', 'b','FaceAlpha',.3); 
% 
% ant = patch('vertices',mVpM.FC{2,1}.vertices,'faces',mVpM.FC{2,1}.faces);
% set(ant, 'FaceColor', 'g', 'EdgeColor', 'b','FaceAlpha',.3);
% hold on
% %show 9 critical points on the 3D plot.
% plot3(mVpM.PT(:,1,1),mVpM.PT(:,1,2),mVpM.PT(:,1,3),'yo','markerfacecolor','b')
% 
% %show corresponding names of 9 cirtical points.
% for q=1:length(mVpM.PTS)
%     text(mVpM.PT(q,1,1)+.5,mVpM.PT(q,1,2)+.5,mVpM.PT(q,1,3)+.5,mVpM.PTS{q}(12:end))
% end

% exportTriangulation2VTK('sampleExampleTri_1',mVpM.FC{1,1}.vertices,mVpM.FC{1,1}.faces)
% exportTriangulation2VTK('sampleExampleTri_2',mVpM.FC{2,1}.vertices,mVpM.FC{2,1}.faces)
%% modeify the intensity of original 3D intensity matrix.
% intensity3d = vol.CTVolume3D;   % size(intensity3d) = 512X512X320(int16) 3D matrix.
% new_intensity = double(intensity3d(cutPosition(selectDataSet,1,1):cutPosition(selectDataSet,1,2),...
%     cutPosition(selectDataSet,2,1):cutPosition(selectDataSet,2,2),...
%     cutPosition(selectDataSet,3,1):cutPosition(selectDataSet,3,2)));    %Specifying the cuting position.
% % figure(1),imshow3D(new_intensity,[-500,1000]);
% 
% %  new_intensity(new_intensity>550) = new_intensity(new_intensity>550)+100;
% % new_intensity(450>new_intensity>-50) = new_intensity(450>new_intensity>-50)-50;
% %  new_intensity(new_intensity<450) = new_intensity(new_intensity<450)-50;
% 
% %Normalized the 3d intensity matrix.
% % new_intensity3d = (new_intensity-min(new_intensity(:)))/(max(new_intensity(:))-min(new_intensity(:)))*40-20;
% % new_intensity3d = new_intensity - 500;
% new_intensity3d = new_intensity - repmat(500,size(new_intensity));   %||-->>
% % new_intensity3d = new_intensity3d - 400;

%%%%%%%%%% create new_intensity3d manually as synData %%%%%%%%%%
new_intensity3d = zeros(512,512,320);
size(new_intensity3d);
newOrigin = [cutPosition(selectDataSet,1,1),cutPosition(selectDataSet,2,1),cutPosition(selectDataSet,3,1)];



% FC{1,1}.vertices contain the coordinates of mitral valve_post.
% FC{2,1}.vertices contain the coordinates of mitral valve_ant.
% PT contains the 3d coordinates of 7 critical and 2 parpillary muscle tips
% points.
%%
dcmfield = header.dcmSeries;
resolutionVector = (reshape(dcmfield.Resolution,[1,3]));
mV = mVpM.FC;    %contains mesh of ant and post valve. 
pM = mVpM.PT;   % size(pM)    %contains 9 critical points.

%% dealing with pM.
% do not need to reverse z coordinates, and swap x and y.
pM_size = size(pM,1);

% pM(:,:,3) = repmat(320*resolutionVector(3),size(pM(:,:,3))) - pM(:,:,3);    % z is from inverse direction.||-->>>>
% key_points_1 = round(reshape(pM(8,1,:),[1,3])./resolutionVector);  %coord of one tip.||-->>>>
% key_points_1 = round(rdivide(reshape(pM(8,1,:),[1,3]),resolutionVector));
% key_points_2 = round(reshape(pM(9,1,:),[1,3])./resolutionVector);  %coord of another tip.||-->>>>
% key_points_2 = round(rdivide(reshape(pM(9,1,:),[1,3]),resolutionVector));
% temp = [key_points_1; key_points_2];
temppm = round(rdivide(reshape(pM,[pM_size,3]),resolutionVector));
% temp_pM_xy = temp(:,1);  % exchange x-y coordinate.
% temp(:,1) = temp(:,2);
% temp(:,2) = temp_pM_xy;

%%%%%%%%%%%%%%%%%%% Redefine the coord after cutting %%%%%%%%%%%%%%%%%%%
%temp(:,1),temp(:,2),temp(:,3) is the X,Y,Z in 3D img Volume.
% tips_all = temppm - newOrigin;  % size(halfVal_1)
tips_all = swapXY_320Z(temppm, 320);
size(tips_all);
% tips_all(tips_all<0) = 0;
% pM_X = tips_all(:,1);
% pM_Y = tips_all(:,2);
% pM_Z = tips_all(:,3);
% pM_X(pM_X>Delta_X) = Delta_X;
% pM_Y(pM_Y>Delta_Y) = Delta_Y;
% pM_Z(pM_Z>Delta_Z) = Delta_Z;
% tips_all(:,1) = pM_X;
% tips_all(:,2) = pM_Y;
% tips_all(:,3) = pM_Z;
% tips = tips_all(8:9,:); %only take the last two coordinates.
tips = tips_all;
size(tips);
% tip1:251   371   188 and tip2:238   345   232


%% dealing with mitral valve
% do not need to reverse z coordinates, and swap x and y.
% resolutionVector = [1,1,1];
size_halfVal_1 = size(mV{1,1}.vertices,1);
size_halfVal_2 = size(mV{2,1}.vertices,1);   % size(halfVal_2)
tempmv_ori = rdivide([mV{1,1}.vertices;mV{2,1}.vertices],resolutionVector);
tempmv = round(rdivide([mV{1,1}.vertices;mV{2,1}.vertices],resolutionVector)); 

% tempmv_xy = tempmv(:,1);  % exchange x-y coordinate.
% tempmv(:,1) = tempmv(:,2);
% tempmv(:,2) = tempmv_xy;
% tempmv(:,3) = repmat(320, size(tempmv(:,3))) - tempmv(:,3);  %invert z coordinate.||-->>

%%%%%%%%%%%%%%%%%%% Redefine the coord after cutting %%%%%%%%%%%%%%%%%%%
%tempmv(:,1),tempmv(:,2),tempmv(:,3) is the X,Y,Z in 3D img Volume.
% mitralValve = tempmv -newOrigin; %xudong: the results are good.
mitralValve_ori = swapXY_320Z(tempmv_ori, 320);
mitralValve = swapXY_320Z(tempmv, 320);
% mitralValve(mitralValve<0) = 0;

% mV_X = mitralValve(:,1);
% mV_Y = mitralValve(:,2);
% mV_Z = mitralValve(:,3);
% mV_X(mV_X>Delta_X) = Delta_X;
% mV_Y(mV_Y>Delta_Y) = Delta_Y;
% mV_Z(mV_Z>Delta_Z) = Delta_Z;
% mitralValve(:,1) = mV_X;
% mitralValve(:,2) = mV_Y;
% mitralValve(:,3) = mV_Z;
% size(mitralValve);


%% Redefine the coord after cutting.
% 
% disp('----Range of x from pM segmentation before chipping----');
% min(tips(:,1))
% max(tips(:,1))
% disp('----Range of y from pM segmentation before chipping----');
% min(tips(:,2))
% max(tips(:,2))
% disp('----Range of z from pM segmentation before chipping----');
% min(tips(:,3))
% max(tips(:,3))

% disp('----Range of x from mV segmentation before chipping----');
% min(mitralValve(:,1)) 
% max(mitralValve(:,1))
% disp('----Range of y from mV segmentation before chipping----');
% min(mitralValve(:,2)) 
% max(mitralValve(:,2))
% disp('----Range of z from mV segmentation before chipping----');
% min(mitralValve(:,3)) 
% max(mitralValve(:,3))
x_min = min(min(tips(:,1)),min(mitralValve(:,1)));
x_max = max(max(tips(:,1)),max(mitralValve(:,1)));
y_min = min(min(tips(:,2)),min(mitralValve(:,2)));
y_max = max(max(tips(:,2)),max(mitralValve(:,2)));
z_min = min(min(tips(:,3)),min(mitralValve(:,3)));
z_max = max(max(tips(:,3)),max(mitralValve(:,3)));
new_ori = [x_min, y_min, z_min];
translation = 10;
% [x_max, y_max, z_max];
tips_new = tips - new_ori + translation;
mitralValve_new = mitralValve - new_ori + translation;
new_intensity3d = zeros(max(cat(1,tips_new,mitralValve_new))+translation*2);


%% mexTopofix3d_good_chord() can find the chord %%
    chord_ub = 200; % used 800 before.
    chord_lb = 0;
    n_pix = 30;    % =1 -> 3 pix on the chord crosssection.
    Npix_mV = 3;
    Npix_pM = 4;    %=7
    changeToIntensity = -800;
    pM_tol = 15;
    continue_reduction = 1;
    
    changeToIntensity_mV = -800;    %800 and 1500-> same results.
    changeToIntensity_pM = -900;    %800 and 1500-> same results.
    ncomp_ub = 1;
    nhole_ub = -1;

%     figure(1),imshow3D(new_intensity3d,[-800,1500]);  %-->
    for death_1d=-500    %-200 find two chords

        [perturbM_0d,~,~,img3d_black_mV_pM,critM3d, dijM3d,...
        tempimg3d]=mexTopofix3d_good_chord(new_intensity3d,...
        ncomp_ub,chord_ub,n_pix,mitralValve_new,Npix_mV,tips_new,Npix_pM,...
        changeToIntensity, chord_lb,pM_tol,continue_reduction,...
        death_1d,changeToIntensity_mV, changeToIntensity_pM);
    end



%% drawing the results
% % fname = 'result.mat'
% % save(fname, 'new_intensity3d','perturbM_0d','img3d_black_mV_pM','critM3d');
 
name_1 = '1.img3d_black_mV_pM.vtk';
Mat2VTK(name_1,img3d_black_mV_pM,'ascii');

name_2 = '2.only2tips.vtk';
% % Mat2VTK(name_2,img3d_black_mV_pM - critM3d,'ascii');    %paraview: -500 show two tips and mV.
Mat2VTK(name_2,tempimg3d,'ascii');
name_3 = '3.markers.vtk';
Mat2VTK(name_3,critM3d,'ascii');

name_4 = '4.drawchords_w2w3_aveOn2Points.vtk';
Mat2VTK(name_4,dijM3d,'ascii');
figure(2), imshow3D(dijM3d,[-1500,1000]);
%% in real coordinate show the mitral valve and tips.
halfVal_1 = mitralValve_new(1:size_halfVal_1,:);
halfVal_2 = mitralValve_new(size_halfVal_1+1:size(mitralValve),:);
halfVal_ori_1 = mitralValve_ori(1:size_halfVal_1,:)- new_ori + translation;;
halfVal_ori_2 = mitralValve_ori(size_halfVal_1+1:size(mitralValve),:)- new_ori + translation;;
if(size(halfVal_1,1) ~= size_halfVal_1)
    disp("halfVal_1 ~= size_halfVal_1");
    size(halfVal_1,1)
    size_halfVal_1
end
if(size(halfVal_2,1) ~= size_halfVal_2)
    disp("halfVal_2 ~= size_halfVal_2");
    size(halfVal_2,1)
    size_halfVal_2
end

post = patch('vertices',halfVal_1,'faces',mVpM.FC{1,1}.faces);
% post = patch('vertices',mVpM.FC{1,T}.vertices,'faces',mVpM.FC{1,T}.faces);
set(post, 'FaceColor', 'red', 'EdgeColor', 'b','FaceAlpha',.3); 
ant = patch('vertices',halfVal_2,'faces',mVpM.FC{2,1}.faces);
set(ant, 'FaceColor', 'g', 'EdgeColor', 'b','FaceAlpha',.3);
hold on
tips_all = tips_all - new_ori + translation;
%show 9 critical points on the 3D plot.
plot3(tips_all(:,1),tips_all(:,2),tips_all(:,3),'yo','markerfacecolor','b')

%show corresponding names of 9 cirtical points.
for q=1:length(mVpM.PTS)
    text(tips_all(q,1)+.5,tips_all(q,2)+.5,tips_all(q,3)+.5,mVpM.PTS{q}(12:end))
end

% % when draw the mV, I should use the ori_sampleExampleTri_1 and _2.
% % exportTriangulation2VTK('ori_sampleExampleTri_1',(mV{1,1}.vertices)/resolutionVector,mVpM.FC{1,1}.faces)
% % exportTriangulation2VTK('ori_sampleExampleTri_2',mV{2,1}.vertices,mVpM.FC{2,1}.faces)
% exportTriangulation2VTK('ori_mV_1',halfVal_ori_1,mVpM.FC{1,1}.faces);
% exportTriangulation2VTK('ori_mV_2',halfVal_ori_2,mVpM.FC{2,1}.faces);
% % when calculate the chordae, I should use the the coord, halfVal_1 and _2.
exportTriangulation2VTK('mV_1',halfVal_1,mVpM.FC{1,1}.faces);
exportTriangulation2VTK('mV_2',halfVal_2,mVpM.FC{2,1}.faces);


%%
%%%%%%%%%%%%%%% show 3D img %%%%%%%%%%%%%%%%%%%%%
%%%%% C++(z, x, y) ---> matlab_Fig(z, y, x) %%%%%%
%%%%%%%%%%%%%%%% show 3D img %%%%%%%%%%%%%%%%%%%%%
% figure(13),imshow3D(new_intensity3d + perturbM_0d,[-1500,1000]);
% figure(14),imshow3D(img3d_black_mV_pM + perturbM_0d,[-1500,1000]);  %-->
% figure(15),imshow3D(img3d_black_mV_pM + critM3d,[-1500,1000]);  %-->
% [x,y,z] = ndgrid(-1:.05:1);
% figure(1), imshow3D(img3d_black_mV_pM,[-1500,1000]);
% figure(2), imshow3D(perturbM_0d,[-40,40]);
% figure(3), imshow3D(critM3d,[-40,40]);
% figure(4),isosurface(perturbM_0d,1/2);
% isosurface(new_intensity3d,1/2);
% daspect([1 1 1]); 
% axis vis3d; 

function result_vector_3d = swapXY_320Z(vector3d, zdim)
    result_vector_3d(:,1) = vector3d(:,2);
    result_vector_3d(:,2) = vector3d(:,1);
    result_vector_3d(:,3) = repmat(zdim, size(vector3d(:,3))) - vector3d(:,3);  %invert z coordinate.||-->>

end






