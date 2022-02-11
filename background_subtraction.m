close all;
clear;
clc;
%% read video and compute background
% open the video
v1 = VideoReader('data/tennis-video.mp4'); % 405 frames, 720x1280
% read all the frames
frames = read(v1,[1 Inf]); % 4D array
% compute the background (median of all the frames)
background = median(frames, 4);
% check background is 720x1280x3
whos background
% show the background
figure(1); imshow(background); title('background');
%% read image to detect (frame number n)
% number of frame to detect
n = 100;
v2 = VideoReader('data/background-video.mp4');
% read all the frames
frames = read(v1,[1 Inf]); % 4D array
% take one image (the n-th frame) and show
img = frames(:,:,:,n);
figure(2); imshow(img); title('img');
%% select color
% take 100-th frame to get color of the ball
init_frame = frames(:,:,:,100);
% take points
figure(3); imshow(init_frame); title('select as many points in the ball as you can');
hold on;
fprintf('Zoom on the ball then press enter to continue\n');
pause
[x, y] = getpts;
% compute avg colour (remember that reference frame in images is inverted)
color = [0,0,0];
for i = 1 : length(x)
    new_color = init_frame(round(y(i)),round(x(i)),:);
    new_color = reshape(new_color,1,3);
    color = color + double(new_color);
end
color = color ./ length(x);
% rescale to 0-1
color = color ./ 255;
% close graph
hold off; close(3);
%% show an image full of the selected color
x = [0 1 1 0] ; y = [0 0 1 1] ;
figure(3);
fill(x,y,color);
color = color .* 255;
%% filter the image by the selected color (see pages 446-448 book digital image processing)
% use bounding box for improve efficiency; compute red standard deviation:
red_std = std(double(init_frame(:,:,1)),0,'all');
% compute green standard deviation:
green_std = std(double(init_frame(:,:,2)),0,'all');
% compute blue standard deviation:
blue_std = std(double(init_frame(:,:,3)),0,'all');
% decide percentage of std to use
k = 1.25;
% compute bounding box boundaries
red_min = max(0, color(1) - k*red_std);
red_max = min(255, color(1) + k*red_std);
green_min = max(0, color(1) - k*green_std);
green_max = min(255, color(1) + k*green_std);
blue_min = max(0, color(1) - k*blue_std);
blue_max = min(255, color(1) + k*blue_std);
% compute the mask
color_mask  = ( (img(:,:,1) >= red_min) & (img(:,:,1) <= red_max) ) & ...
    (img(:,:,2) >= green_min ) & (img(:,:,2) <= green_max) & ...
    (img(:,:,3) >= blue_min ) & (img(:,:,3) <= blue_max);
% show the image filtered by color (apply mask to rgb image)
figure(3);
% apply function @times (.*) to img and color_mask (casted to same type as
% img)
imshow(bsxfun(@times, img, cast(color_mask,class(img))));
title('Image filtered by color');
%% add morphological operations to color_mask
se = strel('square',3); % take a 3x3 morphological element
%mask_morph=imopen(color_mask,se);
color_mask=imclose(color_mask,se); % dilation + erosion
% show result
figure(3);
imshow(bsxfun(@times, img, cast(color_mask,class(img))));
title('Image filtered by color morph transformed');
%% filter by background subtraction (gray images - not sure it works)
% subtract the background
diff = abs(im2double(rgb2gray(img)) - im2double(rgb2gray(background)));
% show the result
figure; imshow(diff); title('diff');
% compute the threshold t through Otsu's method
t = graythresh(diff(:));
% compute the mask
diff_mask = (diff > t);
% show the image filtered by background subtraction
figure(4);
imshow(bsxfun(@times, img, cast(diff_mask,class(img))));
title('Image filtered by background subtraction');
%% combine the two masks to obtain an image filtered by both color and diff
% compute mask
mask = color_mask .* diff_mask;
% show the image filtered
figure(5);
imshow(bsxfun(@times, img, cast(mask,class(img))));
title('Image filtered');