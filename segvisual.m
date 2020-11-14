clear all;
% a = load('fcc.mat');
load Study01.20170830.07.6.14.1601.mat
% load Study02.20170911.05.3.16.641.mat
% load Study03.20170906.02.3.6.641.mat

%%=================== show mitral valve =======================
% %2X1 cell FC stroes the coordinates of the mesh.
% FC{1,T}.vertices contains the x,y,z coordinates
% FC{1,T}.faces has multiple rows, each row assigns 
% the which three vertices creates a triangel
% patch('Face',F,'Vertices',V) or patch('Vertices',V,'Face',F) creates multiple triangels
% set(??????1????1????2????2??)
% alpha(0.3)% set all patches transparency to 0.3

% T = 1;
post = patch('vertices',FC{1,1}.vertices,'faces',FC{1,1}.faces);
% post = patch('faces',FC{1,T}.faces,'vertices',FC{1,T}.vertices);
set(post, 'FaceColor', 'red', 'EdgeColor', 'b','FaceAlpha',.3); 

ant = patch('vertices',FC{2,1}.vertices,'faces',FC{2,1}.faces);
set(ant, 'FaceColor', 'g', 'EdgeColor', 'b','FaceAlpha',.3);

exportTriangulation2VTK('sampleExampleTri_1',FC{1,1}.vertices,FC{1,1}.faces)
exportTriangulation2VTK('sampleExampleTri_2',FC{2,1}.vertices,FC{2,1}.faces)
% daspect([1 1 1]); 
% axis vis3d; 
% camlight; 
% camlight(210,30)
% lighting gouraud

hold on
% %%======================= show tips ===========================
%show 9 critical points on the 3D plot.
plot3(PT(:,1,1),PT(:,1,2),PT(:,1,3),'yo','markerfacecolor','b')

%show corresponding names of 9 cirtical points.
for q=1:length(PTS)
    text(PT(q,1,1)+.5,PT(q,1,2)+.5,PT(q,1,3)+.5,PTS{q}(12:end))
end

%%==================== Additional info ========================
% FC{1,1}.vertices contain the coordinates of mitral valve_post.
% FC{2,1}.vertices contain the coordinates of mitral valve_ant.
% PT contains the 3d coordinates of 7 critical and 2 parpillary muscle tips
% points.


