clear all;
clc;
mex CC=gcc-4.3 LD=g++-4.3 -O mexTopo3d_synData.cpp
setmypath;



% % figure(1),imshow(I4,[])
% % title('Original Image')
% 
% % figure(1),imshow3D(I4,[-20,20]);
% % figure(1),imshow3D(b,[0,255]);
% 
% A(:,:,1) = [0.1656    -0.2630    0.689;  
%             0.6020    0.6541    -0.7482];
%         
% A(:,:,2) = [0.4505    0.2290    0.1524;  
%             -0.0838    -0.9133    -0.8258];
%         
% A(:,:,3) = [0.5383    0.0782    0.1067;  
%             0.9961    -0.4427    -0.9619];
% 
% A(:,:,4) = [-0.0046    0.8173    0.0844;  
%             0.7749    0.8687    -0.3998];
%         
% % figure(1),imshow3D(A,[0,1]);  
% disp('size of A: ');size(A)
% 
% % B = cat()
% B = A(:,:,1);
% C = repmat(B,1,1,10);
% size(C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1 = imread('1.png');
b1 = double(a1);
c1 = b1(:,:,1);
% I1 = c1/255.0*40-20;
I1 = cat(3,c1,c1,c1,c1,c1,c1,c1,c1,c1,c1);


a2 = imread('2.png');
b2 = double(a2);
c2 = b2(:,:,1);
% I2 = c2/255.0*40-20;
I2 = cat(3,c2,c2,c2,c2,c2,c2,c2,c2,c2,c2);

a3 = imread('3.png');
b3 = double(a3);
c3 = b3(:,:,1);
% I3 = c3/255.0*40-20;
I3 = cat(3,c3,c3,c3,c3,c3,c3,c3,c3,c3,c3);

a4 = imread('4.png');
b4 = double(a4);
c4 = b4(:,:,1);
% I4 = c4/255.0*40-20;
I4 = cat(3,c4,c4,c4,c4,c4,c4,c4,c4,c4,c4);

a45 = imread('4-5.png');
b45 = double(a45);
c45 = b45(:,:,1);
% I5 = c5/255.0*40-20;
I45 = cat(3,c45,c45);

a451 = imread('4-5-1.png');
b451 = double(a451);
c451 = b451(:,:,1);
% I5 = c5/255.0*40-20;
I451 = cat(3,c451,c451);


a5 = imread('5.png');
b5 = double(a5);
c5 = b5(:,:,1);
% I5 = c5/255.0*40-20;
% I5 = cat(3,c5,c5,c5,c5,c5,c5,c5,c5,c5,c5);
I5 = cat(3,c5);

a6 = imread('6.png');
b6 = double(a6);
c6 = b6(:,:,1);
% I6 = c6/255.0*40-20;
I6 = cat(3,c6,c6);


% ze2dmatrix = ones(50,50)*20.0;
% ze3dmatrix = repmat(ze2dmatrix,1,1,24);
% imgGray_1 = repmat(I4,1,1,2);
% imgGray_2 = cat(3,ze3dmatrix, imgGray_1);
% imgGray_3 = cat(3,imgGray_2,ze3dmatrix);
imgGray_3 = cat(3,I1,I2,I3,I4,I451,I45,I5,I6,I5,I45, I451,I4,I3,I2,I1);
mitValve_Coord = I5;
papillaryMuscle_Coord = I5;
% imgGray_3 = round(imgGray_3);

figure(1),imshow3D(imgGray_3,[-40,40]);


name = 'synpic.vtk';
Mat2VTK(name,imgGray_3,'ascii');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% a = imread('15.png');
% b = double(a);
% c = b(:,:,1);
% I4 = c/255.0*40-20;
% ze2dmatrix = ones(50,50)*20.0;
% ze3dmatrix = repmat(ze2dmatrix,1,1,24);
% imgGray_1 = repmat(I4,1,1,2);
% imgGray_2 = cat(3,ze3dmatrix, imgGray_1);
% imgGray_3 = cat(3,imgGray_2,ze3dmatrix);
% mitValve_Coord = imgGray_2;disp('size of imgGray_2: ');size(imgGray_2)
% papillaryMuscle_Coord = imgGray_1;disp('size of imgGray_1: ');size(imgGray_1)
% % size(imgGray_3);
% % figure(1),imshow3D(imgGray,[-40,40]);

img_3d = double(imgGray_3);disp('size of img_3d: ');size(img_3d)
ncomp_ub = 1;
nhole_ub = -1;
Npix_increase_on_onePixchords =0;
Npix_increase_on_mV =2;
Npix_increase_on_pM =3;



[perturbM_0d,num_comp,num_hole]=mexTopo3d_synData(img_3d,ncomp_ub, nhole_ub, Npix_increase_on_onePixchords,mitValve_Coord, Npix_increase_on_mV,papillaryMuscle_Coord,Npix_increase_on_pM);


% figure(1),imshow3D(img_3d,[-40,40]);
figure(2), imshow3D(perturbM_0d,[-40,40]);
% % title(['perturbM_0d',int2str(i)]);
% figure(3), imshow3D(img_3d-perturbM_0d, [-40,40]); 
% % title(['img_3d-perturbM_0d-',int2str(i)]);
name = 'original_img.vtk';
Mat2VTK(name,img_3d,'ascii');
newname = 'perturbM_0d.vtk';
Mat2VTK(newname,perturbM_0d,'ascii');


fname = 'result.mat'
save(fname, 'img_3d','perturbM_0d')

[x,y,z] = ndgrid(-1:.05:1);
% v = perturbM_0d;
figure(4),isosurface(perturbM_0d,1/2);
isosurface(img_3d,1/2);
% daspect([1 1 1]); 
axis vis3d; 
camlight; 
camlight(100,30);
lighting gouraud;



