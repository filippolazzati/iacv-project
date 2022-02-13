close all;
clear;
clc;
%% 1) read video and compute background and filter by colour
% open the video
v1 = VideoReader('../mydata/s20fe.mp4'); % 405 frames, 720x1280
% read all the frames
frames = read(v1, [1 Inf]); % 4D array
background = median(frames, 4);
figure(1);
imshow(background);
%% read image to detect (frame number n)
% number of frame to detect
n = 100;
img = frames(:,:,:,n);
figure(2); imshow(img); title('img');
%%
J = imsubtract(img,background);
figure(3); imshow(J);
%%
se = strel('rectangle', [5 1]);
filtered = imopen(J,se);
figure(4); imshow(filtered);
%% mean
avg = mean(img, 3);
img_hsv = rgb2hsv(img);
figure(5); imshow(img_hsv);
V = img_hsv(:,:,3);
imshow(V);
%%
t = graythresh(V(:));
% compute the mask
mask = (V < t);
figure, imshow(mask);
%%
% show the image filtered by background subtraction
new_img = bsxfun(@times, filtered, cast(mask,class(filtered)));
figure(4);
imshow(new_img);
%%
% 1) one player is always on the top and the other below (filter)
% 2) draw bounding box around them (rectangle())













