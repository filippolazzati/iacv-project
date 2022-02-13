close all;
clear;
clc;
addpath(genpath('functions'));
%% read video and compute background
% open the video
v1 = VideoReader('mydata/s20fe.mp4'); % 405 frames, 720x1280
% read all the frames
frames = read(v1, [1 Inf]); % 4D array
background = median(frames, 4);

%%
S = size(background);
y_crop = [410 700];
y_mask = zeros(S(1), S(2));
y_mask(y_crop(1):y_crop(2),:) = 1;
y_mask = double(y_mask);
%%
img = im2double(background) .* y_mask;

figure;
imshow(img);
%%
n = 50;
cx = ceil(S(2) / 2);
cy = ceil(S(1) / 2);
im2 = imcrop(img, [cx - n, cy - n, n, n]);
im2=rgb2hsv(im2);
imgh=rgb2hsv(img);
color = median(im2, [1 2]);

kh=0.3;ks=1.5;kv=10;

h_std = std(double(imgh(:,:,1)),0,'all');
s_std = std(double(imgh(:,:,2)),0,'all');
v_std = std(double(imgh(:,:,3)),0,'all');

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
%%
se = strel('square',3);

x = imlocalbrighten(imtophat(res(:,:,3), se), .35);
figure;imshow(x);
x = im2double(imlocalbrighten(localcontrast(im2uint8(x), 0.3, -1)));

%figure;imshow(x);
%%
gmask = (x > .2);
masked = bsxfun(@times, x, cast(gmask,class(x)));
se2 = strel('square',3);
s = imclose(masked, se2);
masked = bsxfun(@times, masked, cast(s,class(masked)));
masked = bwareaopen(masked, 5);
%figure; imshow(masked);
%%
edges = edge(masked, 'canny');
%figure;
%imshow(edges);

figure; imshow(img);
[H,o_pi,R] = hough(edges);
P = houghpeaks(H, 30, 'threshold', 0.2*max(H(:)));
hlines = houghlines(edges, o_pi, R, P, 'FillGap', 10, 'MinLength', 100);
hold all;
lines = nan(length(hlines), 3);
for k = 1:length(hlines)
    seg = [hlines(k).point1; hlines(k).point2];
    l = segToLine(seg);
    lines(k,:) = l;
    %plot(seg(:,1),seg(:,2),'Color', 'yellow', 'LineWidth', 2);
    angle = rad2deg(atan2(-l(1), l(2)));
    if -5 <= angle && angle <= 5
        lc = 'b';
    else
        lc = 'r';
    end
    plotHCLine(l, lc);
end
%%
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


function plotHCLine(ll, color)
    l = ll ./ ll(3);
    disp(l);
    xlims = get(gca, 'XLim');
    t = linspace(-500, ceil(xlims(2))+500, 10);
    line(t, (-l(1) * t - 1)/l(2), 'Color', color, 'LineWidth', 2);
end