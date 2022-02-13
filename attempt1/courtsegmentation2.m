close all;
clear;
clc;
%% 1) read video and compute background and filter by colour
% open the video
v1 = VideoReader('../mydata/s20fe.mp4'); % 405 frames, 720x1280
% read all the frames
frames = read(v1, [1 Inf]); % 4D array
img = median(frames, 4);
figure(1);
imshow(img);
%% filter by size
S = size(img);
y_crop = [410 700];
y_mask = zeros(S(1), S(2));
y_mask(y_crop(1):y_crop(2),:) = 1;
y_mask = double(y_mask);
background = im2double(img) .* y_mask;
figure(2);
imshow(background);
%% filter by colour
n = 50;
cx = ceil(S(2) / 2);
cy = ceil(S(1) / 2);
im2 = imcrop(background, [cx - n, cy - n, n, n]);
% show median colour
color = median(im2, [1 2]);
x = [0 1 1 0] ; y = [0 0 1 1] ;
figure(2);
fill(x,y,color);
% move to hsv
im2=rgb2hsv(im2);
imgh = rgb2hsv(background);
color = rgb2hsv(color);
% filter by colour
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
%% 2) HSV
background_hsv = res;
figure(2);
imshow(background_hsv);
title('hsv');
%% 3) V
V = background_hsv(:,:,3);
figure(3);
imshow(V);
title('v');
%% 4) apply top-hat transform to V
se = strel('disk',7);
E = imtophat(V,se);
figure(4);
imshow(E);
title('E');
%% 5) threshold
%t = graythresh(E); % 0.0510
%t = 0.15;
t = adaptthresh(E);
B = imbinarize(E,t);
figure(5);
imshow(B);
title('B');
%% 6) extract horizontal lines
%se = strel('square',20);
se = strel('line', 20, 0);
marker = imerode(B,se);
figure(6);
imshow(marker);
title('marker hor');
%% reconstruct
hor_lines1 = imreconstruct(marker,B);
figure(7);
imshow(hor_lines1);
title('hor lines1');
%% opening hor
se = strel('line', 50, 0);
hor_lines2 = imopen(hor_lines1,se); % hor_lines1 !
figure(8);
imshow(hor_lines2);
title('hor lines2');
%% opening vert
se = strel('line', 1, 90);
hor_lines3 = imopen(hor_lines2,se); % hor_lines2 !
figure(9);
imshow(hor_lines3);
title('hor lines3');
%% 7) extract vertical lines
%se = strel('square',20);
se = strel('square', 3);
marker = imerode(B,se);
figure(10);
imshow(marker);
title('marker vert');
%% reconstruct
vert_lines1 = imreconstruct(marker,B);
figure(11);
imshow(vert_lines1);
title('vert lines1');
%% opening vert
se = strel('line', 5, 90);
vert_lines2 = imopen(vert_lines1,se); % vert_lines1 !
figure(12);
imshow(vert_lines2);
title('vert lines2');
%% opening hor
se = strel('line', 3, 0);
vert_lines3 = imopen(vert_lines2,se); % vert_lines2 !
figure(13);
imshow(vert_lines3);
title('vert lines3');

%% 8)
L = hor_lines3 + vert_lines3;
%% 9) reconstruction by dilation







% show
%montage(cat(3,background,background_hsv,V));




















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