close all
clear
clc
%% detect checkerboard points
n_images = 4;
imds = imageDatastore('calibration/data/note8/*.jpg');
[imagePoints,boardSize,imagesUsed] = detectCheckerboardPoints(imds.Files(1:n_images), 'MinCornerMetric', 0.3, 'PartialDetections', false);

for i = 1:n_images
  % Read image
  I = readimage(imds, i);
  figure(i);
  % Insert markers at detected point locations
  I = insertMarker(I, imagePoints(:,:,i), 'o', 'Color', 'red', 'Size', 10);

  imshow(I);
end
%% generate world points
squareSizeInMM = 0.0185;
worldPoints = generateCheckerboardPoints(boardSize,squareSizeInMM);
[cameraParams,imagesUsed,estimationErrors] = estimateCameraParameters(imagePoints,worldPoints);

%%
figure;
u = undistortImage(imread('calibration/data/note8/1.jpg'), cameraParams);
imshow(u);