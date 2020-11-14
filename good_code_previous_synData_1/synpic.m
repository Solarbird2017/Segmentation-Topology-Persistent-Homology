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
imgGray_3 = cat(3,I1,I2,I3,I4,I45,I5,I6,I5,I45,I4,I3,I2,I1);
mitValve_Coord = I5;
papillaryMuscle_Coord = I5;
% imgGray_3 = round(imgGray_3);

figure(1),imshow3D(imgGray_3,[-40,40]);


name = 'synpic.vtk';
Mat2VTK(name,imgGray_3,'ascii');