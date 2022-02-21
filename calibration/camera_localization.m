close all
clear all
clc

%% Localization
% hardcode K of s20fe
K = [1.027599538696298e+03,0,3.113765846725814e+02;0,1.011012752325719e+03,5.483024790016393e+02;0,0,1];
% for s20fe use 1 service rectangle; the centre O is the centre of the
% court
O = [0, 0];
P1 = [0, 4.12];
P3 = [6.40, 0];
P2 = [6.40, 4.12];

img = imread('s20fe_single_photo.jpg');
imshow(img);
figure(1), imshow(img);
title('original image');
hold on;
disp('Select origin, P1 (right), P2 (right down), P3 (left)');
[x, y] = getpts();
imgO = [x(1); y(1); 1];
imgP1 = [x(2); y(2); 1];
imgP2 = [x(3); y(3); 1];
imgP3 = [x(4); y(4); 1];
hold off;
%% 
% fit the homography from scene to image 
H_tform = fitgeotrans([O; P1; P2; P3], [imgO(1:2).'; imgP1(1:2).'; imgP2(1:2).'; imgP3(1:2).'], 'projective');
H = H_tform.T.'; % get matrix 
% normalization parameter
lambda = 1 / norm(K \ H(:,1));
% compute the versors of the plane pi with regards to the camera (world)
% reference frame
i_pi = (K \ H(:,1)) * lambda;
j_pi = (K \ H(:,2)) * lambda;
k_pi = cross(i_pi,j_pi);
% find the translation vector of the plane reference frame with regards to
% the camera frame
o_pi = (K \ (lambda * H(:,3)));
% find the rotation matrix
R = [i_pi, j_pi, k_pi];
% exploit svd to cope with noise and retrieve an orthogonal matrix
[U, ~, V] = svd(R);
R = U * V';
% camera pose is the inverse
cameraOrientation = R.';
cameraPosition = -R.'*o_pi;
% plot the results
figure(2);
plotCamera('Location', cameraPosition, 'Orientation', cameraOrientation.', 'Size', 1, 'Color', [0,0,0]);
hold on
pcshow([[P1; P2; P3], zeros(size([P1; P2; P3],1), 1)],'green','VerticalAxisDir', 'up', 'MarkerSize', 150);
pcshow([[O], zeros(size([O],1),1)],'red','VerticalAxisDir', 'up', 'MarkerSize', 500);
set(gcf,'color','w');
set(gca,'color','w');
xlabel('X');
ylabel('Y');
zlabel('Z');