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
background_bw = rgb2gray(background);

%%
Se = strel('disk', 5);
S1 = strel('square', 1);
hollow = ones(20);
hollow(2:19, 2:19) = 0;
S2 = strel('arbitrary', hollow);
init_frame = 9;
mask_thresh = 25;
player_thresh = 15;

F = [
    1 0 1 0;
    0 1 0 1;
    0 0 1 0;
    0 0 0 1
];

F_player = [
    1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1
];

n_particles = 2000;
frame_width = size(frames, 2);
frame_height = size(frames, 1);

Xstd_pos = 20;
Xstd_vel = 3;

Xstd_measure = 20;

window_semiwidth = ceil(frame_width / 2);
window_semiheight = ceil(frame_height / 2);

X = [
    randi(frame_width, 1, n_particles);  % position X
    randi(frame_height, 1, n_particles); % position Y
    zeros(2, n_particles)                % velocity
];

w = ones(1, n_particles) ./ n_particles;

nframe = init_frame;
pts = nan(size(frames, 4)-nframe, 2);
while nframe < size(frames, 4)
    j = nframe - init_frame + 1;
    frame = frames(:,:,:,nframe);
    frame_bw = rgb2gray(frame);

    bw_bg_sub = imclose(imabsdiff(frame_bw, background_bw), Se);
    mask_bg_sub = bw_bg_sub > mask_thresh;
    mask_bg_sub = mask_bg_sub & ~bwareaopen(mask_bg_sub, 50);

    mask_and = mask_bg_sub;

    for df = [-8, -6, -4, +4, +6, +8]
        prev_frame_bw = rgb2gray(frames(:,:,:,nframe-df));
        bw_frame_diff = imclose(imabsdiff(frame_bw, prev_frame_bw), Se);
        mask_frame_diff = bw_frame_diff > mask_thresh;
        mask_frame_diff = mask_frame_diff & ~bwareaopen(mask_frame_diff, 50);
        mask_and = mask_and & mask_frame_diff;
    end
   
    %mask_hit_miss = bwhitmiss(mask_and, S1, S2);
    mask_final = imdilate(mask_and, strel('disk', 2));

    figure(1);
    %imshow(frame);
    imshow([frame, gray2rgb(bw_bg_sub)]);
    hold all;

    %mask_players = bw_bg_sub > player_thresh;
    %mask_players = bwareaopen(imclose(mask_players, strel('disk', 10)), 200);
    %player_props = regionprops(mask_players, 'BoundingBox');
    %player_bbx = vertcat(player_props.BoundingBox);

    props = regionprops(mask_final, 'BoundingBox', 'Centroid');
    bbx = vertcat(props.BoundingBox);
    %ecc = vertcat(props.Eccentricity);
    cc = vertcat(props.Centroid);

    if (size(bbx, 1) > 0)
        for i = 1:size(bbx, 1)
            if bbx(i,3) > 0 && bbx(i,4) > 0
                rectangle('Position',[bbx(i,1),bbx(i,2),bbx(i,3),bbx(i,4)], 'EdgeColor','r','LineWidth',2);
                pts(j, :) = cc(i,:);
                %text(bbx(i,1), bbx(i,2) - 20, num2str(ecc(i)));
            end
        end
    end
    %{
    for p = 1:size(pts, 1)
        pt = pts(p,:);
        if not(isnan(pt(1)))
            plot(pt(1), pt(2), '.', 'MarkerSize', 10, 'Color', 'green');
        end
    end
    %}

    %%% PARTICLE FILTER %%%
    X = F * X;
    X(1:2,:) = X(1:2,:) + Xstd_pos * randn(2, n_particles);
    X(3:4,:) = X(3:4,:) + Xstd_vel * randn(2, n_particles);

    avg_pos_x = mean(X(1,:));
    avg_pos_y = mean(X(2,:));

    rectangle('Position',[ ...
        avg_pos_x - window_semiwidth, ...
        avg_pos_y - window_semiheight, ...
        2 * window_semiwidth, ...
        2 * window_semiheight ...
    ], 'EdgeColor','r','LineWidth',2);

    selected_cc = nan(1, 2);
    if size(cc, 1) > 0
        for c = 1:size(cc, 1)
            if (cc(c, 1) > avg_pos_x - window_semiwidth && cc(c, 1) < avg_pos_x + window_semiwidth) && (cc(c, 2) > avg_pos_y - window_semiheight && cc(c, 2) < avg_pos_y + window_semiheight)
                selected_cc = cc(c, :);
                window_semiwidth = ceil(frame_width / 10);
                window_semiheight = ceil(frame_height / 10);
                break;
            end
        end
    end

    if ~isnan(selected_cc(1))
        for p = 1:n_particles
            w(p) = max(0, abs(selected_cc(1) - X(1, p)) + abs(selected_cc(2) - X(2, p)) + randn * Xstd_measure);
        end
        w = tanh(1 ./ w);
        w = w ./ sum(w);
        pdf = cumsum(w);
        [~, ~, I] = histcounts(rand(1, n_particles), [0, pdf]);
        X = X(:, I); % +1 to avoid zero index
    end

    predicted = mean(X, 2);

    plot(X(1,:), X(2,:), '.', 'MarkerSize', 2, 'Color', 'green');

    plot(predicted(1), predicted(2), '.', 'MarkerSize', 15, 'Color', 'yellow');
    plot([predicted(1), predicted(1) + predicted(3) * 10], [predicted(2), predicted(2) + predicted(4) * 10], 'Color', 'yellow', 'LineWidth', 4);
    

    nframe = nframe + 1;

    pause(0.2);
end
%%
S = size(background);
y_crop = [420 700];
y_mask = zeros(S(1), S(2));
y_mask(y_crop(1):y_crop(2),:) = 1;
y_mask = double(y_mask);
%%
img = im2double(background) .* y_mask;

figure;
imshow(img);
%%
%{
n = 50;
cx = ceil(S(2) / 2);
cy = ceil(S(1) / 2);
im2 = imcrop(img, [cx - n, cy - n, n, n]);
im2=rgb2hsv(im2);
%}
imc = imcrop(img, [0, y_crop(1), S(2), y_crop(2)-y_crop(1)]);
imc = rgb2hsv(imc);
color = mean(imc, [1 2]);

imgh=rgb2hsv(img);

kh=.7;ks=1.5;kv=10;
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
%{
imgh_sm = imgh;
imgh_sm(:,:,3) = adapthisteq(imgh_sm(:,:,3));
figure;
imshow(hsv2rgb(imgh));
figure;
img_sm = hsv2rgb(imgh_sm);
imshow(img_sm);
%}
%%
se = strel('square', 3);

x = imbinarize(imtophat(res(:,:,3), se));
x = im2double(localcontrast(im2uint8(x), 0.3, -1));

figure;imshow(x);

%%
%gmask = (x > .1);
%masked = bsxfun(@times, x, cast(gmask,class(x)));
se2 = strel('square',3);
masked = x;
s = imclose(masked, se2);
masked = bsxfun(@times, masked, cast(s,class(masked)));
masked = bwareaopen(masked, 5);
figure; imshow(masked);
%%
edges = edge(masked, 'canny');
%figure;
%imshow(edges);

figure; imshow(img);
[H,o_pi,R] = hough(edges);
P = houghpeaks(H, 30, 'threshold', 0.3*max(H(:)));
hlines = houghlines(edges, o_pi, R, P, 'FillGap',25, 'MinLength', 100);
hold all;
lines = nan(length(hlines), 3);
for k = 1:length(hlines)
    seg = [hlines(k).point1; hlines(k).point2];
    l = segToLine(seg);
    lines(k,:) = l;
    plot(seg(:,1),seg(:,2),'Color', 'yellow', 'LineWidth', 5);
    angle = rad2deg(atan2(-l(1), l(2)));
    if -5 <= angle && angle <= 5
        lc = 'b';
    else
        lc = 'r';
    end
    plotHCLine(l, lc);
end
%%

% for i = 1:length(lines)
%     for j = i:length(lines)
%         p = cross(lines(i,:), lines(j, :));
%         if p(3) == 0
%             continue
%         end
%         p = p ./ p(3);
%         if p(1) < 0 || p(1) > S(2) || p(2) < 0 || p(2) > S(1)
%             continue
%         end
%         plot(p(1), p(2), 'ro', 'MarkerSize', 10);
%     end
% end


function plotHCLine(ll, color)
    l = ll ./ ll(3);
    disp(l);
    xlims = get(gca, 'XLim');
    t = linspace(-500, ceil(xlims(2))+500, 10);
    line(t, (-l(1) * t - 1)/l(2), 'Color', color, 'LineWidth', 2);
end

function res = imbin2rgb(img)
    res = 255 * repmat(uint8(img), 1, 1, 3);
end

function res = gray2rgb(img)
    res = cat(3, img, img, img);
end

function res = immask(img, mask)
    res = bsxfun(@times, img, cast(mask,class(img)));
end

function r = in_rect(pos, rect)
    r = pos(1) >= rect(1) && pos(1) <= rect(1) + rect(3) && ...
        pos(2) >= rect(2) && pos(2) <= rect(2) + rect(4);
end
