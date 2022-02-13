close all;
clear;
clc;
%% read video and compute background
% open the video
v1 = VideoReader('mydata/s20fe.mp4'); % 405 frames, 720x1280
% read all the frames
frames = read(v1, [1 Inf]); % 4D array
background = median(frames, 4);

%%
img = background;
figure;
imshow(img);
%% color select
[x, y] = getpts;
color = rgb2hsv(img(round(y(1)), round(x(1)), :));
%%
addpath(genpath('functions'));

%%
S = size(img);
n = 50;
cx = ceil(S(2) / 2);
cy = ceil(S(1) / 2);
im2 = imcrop(img, [cx - n, cy - n, n, n]);
im2=rgb2hsv(im2);
%%
imgh = rgb2hsv(img);
color = median(im2, [1 2]);
%%
color = rgb2hsv([1 1 1]);
%%
% filter the image by the selected color (see pages 446-448 book digital image processing)
% use bounding box for improve efficiency; compute red standard deviation:
h_std = std(double(imgh(:,:,1)),0,'all');
s_std = std(double(imgh(:,:,2)),0,'all');
v_std = std(double(imgh(:,:,3)),0,'all');
% decide percentage of std to use
kh = 0.1;
ks = 1;
kv = 3;
% compute bounding box boundaries
h_min = max(0, color(1) - kh*h_std); h_max = min(255, color(1) + kh*h_std);
s_min = max(0, color(2) - ks*s_std); s_max = min(255, color(2) + ks*s_std);
v_min = max(0, color(3) - kv*v_std); v_max = min(255, color(3) + kv*v_std);
% compute the mask
color_mask  = ( (imgh(:,:,1) >= h_min) & (imgh(:,:,1) <= h_max) ) & ...
    (imgh(:,:,2) >= s_min ) & (imgh(:,:,2) <= s_max) & ...
    (imgh(:,:,3) >= v_min ) & (imgh(:,:,3) <= v_max);
% show the image filtered by color (apply mask to rgb image)
figure(3);
% apply function @times (.*) to img and color_mask (casted to same type as
% img)
res = bsxfun(@times, imgh, cast(color_mask,class(imgh)));
%res = rgb2gray(hsv2rgb(res));
%res = bwareaopen(res, 10);
imshow(res);
title('Image filtered by color');

resg = rgb2gray(hsv2rgb(res));
imeq = adapthisteq(resg);
imeq2 = histeq(resg);
figure;imshow(imeq);
figure;imshow(imeq2);
color_mask=imclose(res, se);
figure(3);
imshow(color_mask);
%imshow(bsxfun(@times, img, cast(color_mask,class(res))));
title('Image filtered by color morph transformed');

%t = graythresh(imeq(:));
% compute the mask
se = strel('square',2);
gmask = (imeq > .6);
bbb = bsxfun(@times, imeq, cast(gmask,class(imeq)));
s = imclose(bbb, se);
bbb = bsxfun(@times, bbb, cast(s,class(bbb)));
bbb = bwareaopen(bbb, 2);
figure; imshow(bbb);

%%
edges = edge(bbb, 'canny', .2);
figure;
imshow(edges);


%%
figure; imshow(bbb);
[H,o_pi,R] = hough(bbb);
P = houghpeaks(H, 4, 'threshold', 0.3*max(H(:)));
hlines = houghlines(edges, o_pi, R, P, 'FillGap', 3, 'MinLength', 40);
hold all;
lines = nan(length(hlines), 3);
for k = 1:length(hlines)
    seg = [hlines(k).point1; hlines(k).point2];
    lines(k,:) = segToLine(seg);
    plot(seg(:,1),seg(:,2),'Color', 'yellow', 'LineWidth', 2);
end
for i = 1:length(lines)
    for j = i:length(lines)
        p = cross(lines(i,:), lines(j, :));
        if p(3) == 0
            continue
        end
        p = p ./ p(3);
        if p(1) < 0 || p(1) > S(2) || p(2) < 0 || p(2) > S(1)
            continue
        end
        plot(p(1), p(2), 'ro', 'MarkerSize', 10);
    end
end