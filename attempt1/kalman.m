close all;
clear;
clc;
addpath(genpath('functions'));
%%

params_s20fe = load("calibration/params/s20fe.mat");
params_note8 = load("calibration/params/note8.mat");
params_s20plus = load("calibration/params/s20plus.mat");

%v1 = VideoReader('mydata/s20fe.mp4');
%v2 = VideoReader('mydata/note8.mp4');
v1 = VideoReader('videos/s20fe.mp4');
v2 = VideoReader('videos/note8.mp4');
v3 = VideoReader('videos/s20plus.mp4');

videos = {v1; v2; v3};

frames1 = read(v1, [200 400]);
frames2 = read(v2, [200 400]);
frames3 = read(v3, [200 400]);

background1_bw = rgb2gray(median(frames1, 4));
background2_bw = rgb2gray(median(frames2, 4));
background3_bw = rgb2gray(median(frames3, 4));
%%
field_fig = figure('Position', [10 10 1000 1000]);
draw_field(field_fig, params_s20fe, params_note8, params_s20plus);

cameraPoses = table;
cameraPoses.ViewId = uint32([1; 2; 3]);
cameraPoses.Orientation = {params_s20fe.orientation; params_note8.orientation; params_s20plus.orientation};
cameraPoses.Location = {params_s20fe.location; params_note8.location; params_s20plus.location};
intrinsics = [params_s20fe.cameraParams.Intrinsics, params_note8.cameraParams.Intrinsics, params_s20plus.cameraParams.Intrinsics];

frame_delta_note8 = +95;
frame_delta_s20plus = +43;

base_frame = 100;
init_frame_s20fe = base_frame;
init_frame_note8 = base_frame+frame_delta_note8;
init_frame_s20plus = base_frame+frame_delta_s20plus;

mask_thresh = 25;

nframe = init_frame_s20fe;

frame_buffer_size = 8;
frame_buffer = cell(frame_buffer_size, size(videos, 1));

for f = 1:frame_buffer_size
    for v = 1:size(videos, 1)
        frame_buffer{f, v} = rgb2gray(read(videos{v}, base_frame - frame_buffer_size + f - 1));
    end
end
%%
ransac_frames = 6;
points_buffer = cell(ransac_frames, 1);

% KALMAN
g = -9.81;
dt = 1 / v1.FrameRate;
mu_k0 = [5.485 11.89 1 0 0 0].';
sigma_k0 = eye(6) * 100;
Q = eye(6) * 0.1;
R = eye(3);
A_k = [
    1 0 0 dt 0 0;
    0 1 0 0 dt 0;
    0 0 1 0 0 dt;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
];
Bk_uk = [0, 0, 0.5 * g * dt^2, 0, 0, g * dt].';
C_k = [eye(3), zeros(3)];

max_lost_frames = 7;
lost_frames = 0;

mu_k = mu_k0; sigma_k = sigma_k0;

while nframe <= v1.NumFrames
    index0 = nframe - init_frame_s20fe;

    frame_s20fe = rgb2gray(read(v1, nframe));
    frame_note8 = rgb2gray(read(v2, nframe+frame_delta_note8));
    frame_s20plus = rgb2gray(read(v3, nframe+frame_delta_s20plus));

    [bbx1, centroids1] = ball_candidates(frame_s20fe, frame_buffer, index0, 1, background1_bw, mask_thresh);
    [bbx2, centroids2] = ball_candidates(frame_note8, frame_buffer, index0, 2, background2_bw, mask_thresh);
    [bbx3, centroids3] = ball_candidates(frame_s20plus, frame_buffer, index0, 3, background3_bw, mask_thresh);

    
    figure(100); hold all;
    imshow(frame_s20fe);
    draw_bbx(bbx1);
%{
    figure(200); hold all;
    imshow(frame_note8);
    draw_bbx(bbx2);

    figure(300); hold all;
    imshow(frame_s20plus);
    draw_bbx(bbx3);
%}
    buffer_idx = mod(index0, ransac_frames) + 1;
    triangulated_pts_count = 0;

    curr_pts = [];

    figure(field_fig); hold on;
    title(strcat('Frame: ', num2str(nframe)));
    for c1 = 1:size(centroids1, 1)
        for c2 = 1:size(centroids2, 1)
            for c3 = 1:size(centroids3, 1)
                p1 = centroids1(c1, :);
                p2 = centroids2(c2, :);
                p3 = centroids3(c3, :);
    
                track = pointTrack([1; 2; 3], [p1; p2; p3]);
                [wp, re] = triangulateMultiview(track, cameraPoses, intrinsics);
        
                if re < 50 && (-5 <= wp(1) && wp(1) <= 15.97) && (0 <= wp(2) && wp(2) <= 23.58)
                    triangulated_pts_count = triangulated_pts_count + 1;
                    points_buffer{buffer_idx}(triangulated_pts_count, :) = wp;
                end
            end
        end
    end

    if triangulated_pts_count == 0
        points_buffer{buffer_idx} = [];
    end

    all_pts = [];
    for i = 1:ransac_frames
         all_pts = [
             all_pts;
             [points_buffer{i}, ones(size(points_buffer{i}, 1), 1) * i]
         ];
    end

    [mu_k1, sigma_k1] = kalman_predict(A_k, mu_k, sigma_k, Bk_uk, Q);
    mu_k = mu_k1; sigma_k = sigma_k1;

    %{
    if (size(all_pts, 1) > 1)
        [plane, inlierIdx] = pcfitplane(pointCloud(all_pts(:, 1:3)), 0.05, [1 0 0], 40);
        inliers = all_pts(inlierIdx, :);

        min_dist_to_pred = Inf;
        selected_mu = mu_k;
        selected_sigma = sigma_k;
        for p = 1:size(inliers, 1)
            if inliers(p, 4) == buffer_idx
                [mu_k1, sigma_k1] = kalman_update(inliers(p, 1:3), mu_k, sigma_k, C_k, R);
                dist_to_pred = norm(mu_k(1:2) - mu_k1(1:2), 'fro');
                if dist_to_pred < min_dist_to_pred
                    selected_mu = mu_k1;
                    selected_sigma = sigma_k1;
                    min_dist_to_pred = dist_to_pred;
                end
            end
        end

        %if min_dist_to_prev < 10
            mu_k = selected_mu;
            sigma_k = selected_sigma;
        %else
            %fprintf('using prediction only: min_dist_to_prev=%d\n', min_dist_to_prev);
        %end


        %plot3(inliers(:, 1),inliers(:, 2),inliers(:, 3), '.', 'MarkerSize', 5, 'Color', 'green');
    end
    %}

    if (size(all_pts, 1) > 1)
        [plane, inlierIdx] = pcfitplane(pointCloud(all_pts(:, 1:3)), 0.05, [1 0 0], 40);
        inliers = all_pts(inlierIdx, :);
    
        min_dist_to_pred = Inf;
        selected_point = nan(3, 1);

        for p = 1:size(inliers, 1)
            if inliers(p, 4) == buffer_idx % if point from current frame
                dist_to_pred = norm(inliers(p, 1:3) - mu_k(1:3).');
                if dist_to_pred < min_dist_to_pred
                    selected_point = inliers(p, 1:3);
                    min_dist_to_pred = dist_to_pred;
                end
            end
        end

        if ~isnan(selected_point)
            window_k = 15;
    
            sigma_pos = sigma_k(1:3, 1:3);
            window = window_k * diag(sigma_pos);
    
            if (in_rect_centered(selected_point, [mu_k(1:3); window]))
                fprintf('Using measurement\n');
                lost_frames = 0;
                [mu_k1, sigma_k1] = kalman_update(selected_point, mu_k, sigma_k, C_k, R);
                mu_k = mu_k1; sigma_k = sigma_k1;
            else
                fprintf('Using prediction only\n');
                lost_frames = lost_frames + 1;
            end
        else
            fprintf('Using prediction only\n');
            lost_frames = lost_frames + 1;
        end

        %{
        for p = 1:size(pts, 1)
            [mu_k1, sigma_k1] = kalman_update(pts(p, :), mu_k, sigma_k, C_k, R);
            dist_to_pred = norm(mu_k(1:3) - mu_k1(1:3));
            if dist_to_pred < min_dist_to_pred
                selected_mu = mu_k1;
                selected_sigma = sigma_k1;
                min_dist_to_pred = dist_to_pred;
            end
        end
        %}

        %plot3(inliers(:, 1),inliers(:, 2),inliers(:, 3), '.', 'MarkerSize', 5, 'Color', 'green');
    end


    if sigma_k(1,1) < 20
        plot3(mu_k(1),mu_k(2),mu_k(3), '.', 'MarkerSize', 10, 'Color', 'green');
    
        %v_kmh = norm(mu_k(4:6))*3.6;
        %if v_kmh > 15 && v_kmh < 190
        %    c = 'green';
        %else
        %    c = 'blue';
        %end
        %plot3([mu_k(1) mu_k(1)+mu_k(4)*dt],[mu_k(2) mu_k(2)+mu_k(5)*dt],[mu_k(3) mu_k(3)+mu_k(6)*dt], 'Color', c);
        %text(mu_k(1),mu_k(2),mu_k(3)+0.1, num2str(norm(mu_k(4:6))*3.6), 'Color', 'yellow', 'FontSize', 8);

        %[sigma_x, sigma_y, sigma_z] = ellipsoid(mu_k(1),mu_k(2),mu_k(3),sigma_k(1,1),sigma_k(2,2),sigma_k(3,3));
        %surf(sigma_x, sigma_y, sigma_z, 'FaceAlpha',0.1);

        %[sigma_x, sigma_y, sigma_z] = ellipsoid(mu_k(1),mu_k(2),mu_k(3),sigma_k(4,4),sigma_k(5,5),sigma_k(6,6));
        %surf(sigma_x, sigma_y, sigma_z, 'FaceAlpha',0.1);
    end

    if lost_frames >= max_lost_frames
        mu_k = mu_k0;
        sigma_k = sigma_k0;
    end

    % frame overlay
    proj1 = worldToImage(params_s20fe.cameraParams, params_s20fe.rotationMatrix, params_s20fe.translationVector, mu_k(1:3).');
    figure(100); hold all; plot(proj1(1), proj1(2), '.', 'MarkerSize', 20, 'Color', 'yellow');

    frame_buffer{index0+1, 1} = frame_s20fe;
    frame_buffer{index0+1, 2} = frame_note8;
    frame_buffer{index0+1, 3} = frame_s20plus;
    nframe = nframe + 1;
    %pause(0.2);
end

%%

function [mu_k1, sigma_k1] = kalman_predict(A_k, mu_k, sigma_k, Bk_uk, Q)
    mu_k1 = A_k * mu_k + Bk_uk;
    if mu_k1(3) < 0
        mu_k1(3) = -mu_k1(3);
        mu_k1(6) = -mu_k1(6);
    end
    sigma_k1 = A_k * sigma_k * A_k.' + Q;
end

function [mu_k1, sigma_k1] = kalman_update(y_k, mu_k, sigma_k, C_k, R)
    K = sigma_k * C_k.' * inv(C_k * sigma_k * C_k.' + R);
    mu_k1 = mu_k + K * (y_k.' - C_k * mu_k);
    sigma_k1 = sigma_k - K * C_k * sigma_k;
end


function [bbx, cc] = ball_candidates(frame_bw, frame_buffer, index0, video_idx, background_bw, mask_thresh)
    bw_bg_sub = imclose(imabsdiff(frame_bw, background_bw), strel('disk', 5));
    %bw_bg_sub = imabsdiff(frame_bw, background_bw);
    mask_bg_sub = bw_bg_sub > mask_thresh;
    mask_bg_sub = mask_bg_sub & ~bwareaopen(mask_bg_sub, 200);

    mask_frame_diff = ones(size(frame_bw, 1), size(frame_bw, 2));

    for df = [8, 6, 4]
        frame_idx = mod(index0 - df, size(frame_buffer, 1)) + 1;
        other_frame_bw = frame_buffer{frame_idx, video_idx};
        bw_frame_diff = imclose(imabsdiff(frame_bw, other_frame_bw), strel('disk', 3));
        mask_other_frame_diff = bw_frame_diff > mask_thresh;
        mask_other_frame_diff = mask_other_frame_diff & ~bwareaopen(mask_other_frame_diff, 200);
        mask_frame_diff = mask_frame_diff & mask_other_frame_diff;
    end

    mask_and = mask_bg_sub & mask_frame_diff;
   
    %mask_hit_miss = bwhitmiss(mask_and, S1, S2);
    mask_final = imdilate(mask_and, strel('disk', 2));

    props = regionprops(mask_final, 'BoundingBox', 'Centroid');
    bbx = vertcat(props.BoundingBox);
    cc = vertcat(props.Centroid);
end

function valid = validate_ransac_poly(P, varargin)
disp(P);
    valid = size(P, 1) > 0 && size(P, 2) > 0 && ...
        P(1) < 0;
end

function draw_bbx(bbx)
    if (size(bbx, 1) > 0)
        for i = 1:size(bbx, 1)
            if bbx(i,3) > 0 && bbx(i,4) > 0
                rectangle('Position',[bbx(i,1),bbx(i,2),bbx(i,3),bbx(i,4)], 'EdgeColor','r','LineWidth',2);
            end
        end
    end
end

function draw_field(field_fig, params_s20fe, params_note8, params_s20plus)
    figure(field_fig);
    hold on;
    pcshow([[0 0], zeros(size([0 0],1),1)], 'red','VerticalAxisDir', 'up', 'MarkerSize', 5);
    
    plotCamera('Location', params_s20fe.location, 'Orientation', params_s20fe.orientation, 'Size', 0.2, 'Color', [1,0,0]);
    plotCamera('Location', params_note8.location, 'Orientation', params_note8.orientation, 'Size', 0.2, 'Color', [0,1,0]);
    plotCamera('Location', params_s20plus.location, 'Orientation', params_s20plus.orientation, 'Size', 0.2, 'Color', [0,0,1]);
    
    fieldSegs = [
    %   x1     x2     y1     y2
        0      10.97  0      0;
        10.97  10.97  0      23.78;
        10.97  0      23.78  23.78;
        0      0      23.78  0;
    
        1.37   1.37   0      23.78;
        9.6    9.6    0      23.78;
    
        0      10.97  11.89  11.89;
    
        1.37   9.6    5.49   5.49;
        1.37   9.6    18.29  18.29;
    
        5.485  5.485  5.49   18.29;
    
        5.485  5.485  0      0.5;
        5.485  5.485  23.28  23.78;
    ];
    
    plot(fieldSegs(:, 1:2).', fieldSegs(:, 3:4).', 'Color', 'b')
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

function r = in_rect_centered(pos, rect)
    % rect = [center_x center_y center_z halflength_x halflength_y halflength_z]

    r = pos(1) >= rect(1) - rect(4) && pos(1) <= rect(1) + rect(4) && ...
        pos(2) >= rect(2) - rect(5) && pos(2) <= rect(2) + rect(5) && ...
        pos(3) >= rect(3) - rect(6) && pos(3) <= rect(3) + rect(6);
end

