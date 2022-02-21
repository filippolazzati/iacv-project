close all
clear
clc
%% detect checkerboard points
n_images = 2;
imds = imageDatastore('calibration data\s20fe\*.jpg');
[imagePoints,boardSize,imagesUsed] = detectCheckerboardPoints(imds.Files(1:n_images),'HighDistortion',false);

for i = 1:n_images
  % Read image
  I = readimage(imds, i);
    figure(i);
  % Insert markers at detected point locations
  I = insertMarker(I, imagePoints(:,:,i), 'o', 'Color', 'red', 'Size', 10);

  imshow(I);
end
%% generate world points
squareSizeInMM = 19;
worldPoints = generateCheckerboardPoints(boardSize,squareSizeInMM);
[cameraParams,imagesUsed,estimationErrors] = estimateCameraParameters(imagePoints,worldPoints);
%% calibrate the camera
params = estimateCameraParameters(imagePoints,worldPoints);
%% remove distortion from ph1.jpg
I = imread('calibration data\s20fe\ph1.jpg');
J = undistortImage(I,cameraParams);
figure(1); imshow(I); title('original');
figure(2); imshow(J); title('corrected');
%% define camera calibration
% focal length
fx = params.FocalLength(1);
fy = params.FocalLength(2);
% principal point
x0 = params.PrincipalPoint(1);
y0 = params.PrincipalPoint(2);
% skew
s = params.Skew;
% calibration matrix
K = [fx, s, x0; ...
     0, fy, y0; ...
     0, 0, 1];
% aspect ratio
a = fx / fy;
% image of the absolute conic
iac = [a*a, 0, -x0*a*a; ...
      0, 1, -y0; ...
      -x0*a*a, -y0, fy*fy + a*a*x0*x0 + y0*y0];
iac2 = inv(K * K.'); 
%% Reconstruction of a vertical facade
img = imread('calibration data\s20fe\ph1.jpg');
imshow(img);
% select the number of directions and lines for each direction
directions = 2;
lines = 2;
% create a cell array in which store the parallel lines
parallelLines = cell(directions,1);
fprintf(['Draw ', num2str(directions) , ' families of parallel segments\n']);
col = 'rgbm';
% for each direction, draw all the parallel lines along it
for i = 1:directions
count = 1;
parallelLines{i} = nan(lines,3);
while(count <=lines)
figure(gcf);
title(['Draw ', num2str(lines),' segments: step ',num2str(count) ]);
segment1 = drawline('Color',col(i));
% segToLine takes two points of a line and returns the projective
% representation of the line
parallelLines{i}(count, :) = segToLine(segment1.Position);
count = count +1;
end
fprintf('Press enter to continue\n');
pause
end
V = nan(2,directions);
% for each direction, find the vanishing point by applying least squares
for i = 1:directions
A = parallelLines{i}(:,1:2);
b = - parallelLines{i}(:,3);
V(:,i) = A\b;
end
imLinfty = fitline(V);
imLinfty = imLinfty./(imLinfty(3));
%%
syms x1;
syms x2;
x=[x1,x2,1];
eq1 = imLinfty.'*x.'==0;
eq2= x*iac*x.'==0;
eqns = [eq1,eq2];
S = solve(eqns, [x1,x2],'ReturnConditions', true, 'Real', false);
sx1 = double(S.x1);
sx2 = double(S.x2);
% compute the circular points
I = [sx1(1), sx2(1), 1].';
J = [sx1(2), sx2(2), 1].';
% find the conic dual to the circular points
imDCCP = I*(J.') + J * (I.');
% normalize imDCCP
imDCCP = imDCCP ./ norm(imDCCP);
% find the rectification matrix
[U, S, ~] = svd(imDCCP);
Hrect = diag([sqrt(1/S(1, 1)), sqrt(1/S(2, 2)), 1])*(U.');
% rotate the image of 135°
H_rot1 = [0,-1,0;1,0,0;0,0,1]; % 90°
H_rot2 = [0.707,-0.707,0;0.707,0.707,0;0,0,1]; % 45°
H_rotation = H_rot1 * H_rot2;
% perform the rectification
tform = projective2d((H_rotation*Hrect).');
metric_rect_img = imwarp(img,tform);
figure;
imshow(metric_rect_img);
title('metric rectified image');


















