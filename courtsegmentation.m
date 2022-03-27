close all;
clear;
clc;
addpath(genpath('calibration'));
%% read video and compute background
v1 = VideoReader('videos/s20plus.mp4');
% read all the frames
frames = read(v1, [200 400]);
background = median(frames, 4);
%%
edges = edge(rgb2gray(background), 'canny', .2);
figure;
imshow(edges);
%%
figure; imshow(background);
[H,o_pi,R] = hough(edges);
P = houghpeaks(H, 8, 'threshold', .3*max(H(:)));
hlines = houghlines(edges, o_pi, R, P, 'FillGap', 20, 'MinLength', 100);
hold all;
lines = nan(length(hlines), 3);
for k = 1:length(hlines)
    seg = [hlines(k).point1; hlines(k).point2];
    l = segToLine(seg);
    lines(k,:) = l;
    angle = rad2deg(atan2(-l(1), l(2)));
    plotHCLine(l, 'r');
end

%%
function plotHCLine(ll, color)
    l = ll ./ ll(3);
    disp(l);
    xlims = get(gca, 'XLim');
    t = linspace(-500, ceil(xlims(2))+500, 10);
    line(t, (-l(1) * t - 1)/l(2), 'Color', color, 'LineWidth', 2);
end