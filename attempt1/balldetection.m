close all;
clear;
clc;
addpath(genpath('functions'));
%%

params_s20fe = load("calibration/calibration data/s20fe/calibrationparams.mat");
params_note8 = load("calibration/calibration data/note8/calibrationparams.mat");

v1 = VideoReader('mydata/s20fe.mp4');
v2 = VideoReader('mydata/note8.mp4');

frames1 = read(v1, [1 Inf]);
frames2 = read(v2, [1 Inf]);

background1_bw = rgb2gray(median(frames1, 4));
background2_bw = rgb2gray(median(frames2, 4));

%%
draw_field(params_s20fe, params_note8);

camMatrix1 = cameraMatrix(params_s20fe.cameraParams, params_s20fe.rotationMatrix, params_s20fe.translationVector);
camMatrix2 = cameraMatrix(params_note8.cameraParams, params_note8.rotationMatrix, params_note8.translationVector);

frame_delta = 24;

init_frame_1 = 9 + frame_delta;
init_frame_2 = 9;

figure(1);
imshow(frames1(:,:,:,init_frame_1));
figure(2);
imshow(frames2(:,:,:,init_frame_2));

mask_thresh = 30;

nframe = init_frame_1;

while nframe < size(frames1, 4)
    [bbx1, centroids1] = ball_candidates(frames1, nframe, background1_bw, mask_thresh);
    [bbx2, centroids2] = ball_candidates(frames2, nframe - frame_delta, background2_bw, mask_thresh);

    figure(1); hold all;
    imshow(frames1(:,:,:,nframe));
    draw_bbx(bbx1);

    figure(2); hold all;
    imshow(frames2(:,:,:,nframe-frame_delta));
    draw_bbx(bbx2);

    if size(centroids1, 1) == 1 && size(centroids2, 1) == 1
        p1 = centroids1(1, :);
        p2 = centroids2(1, :);

        [wp] = triangulate(p1, p2, camMatrix1, camMatrix2);

        figure(3); hold on;
        plot3(wp(1), wp(2), wp(3), '.', 'MarkerSize', 20, 'Color', 'yellow');
    end

    nframe = nframe + 1;
end

%%
function [bbx, cc] = ball_candidates(frames, nframe, background_bw, mask_thresh)
    frame = frames(:,:,:,nframe);
    frame_bw = rgb2gray(frame);

    bw_bg_sub = imclose(imabsdiff(frame_bw, background_bw), strel('disk', 2));
    mask_bg_sub = bw_bg_sub > mask_thresh;
    mask_bg_sub = mask_bg_sub & ~bwareaopen(mask_bg_sub, 50);

    mask_and = mask_bg_sub;

    for df = [-8, -6, -4, +4, +6, +8]
        frame_idx = nframe + df;
        if frame_idx < 0 || frame_idx >= size(frames, 4)
            continue
        end
        prev_frame_bw = rgb2gray(frames(:,:,:,frame_idx));
        bw_frame_diff = imclose(imabsdiff(frame_bw, prev_frame_bw), strel('disk', 2));
        mask_frame_diff = bw_frame_diff > mask_thresh;
        mask_frame_diff = mask_frame_diff & ~bwareaopen(mask_frame_diff, 50);
        mask_and = mask_and & mask_frame_diff;
    end
   
    %mask_hit_miss = bwhitmiss(mask_and, S1, S2);
    mask_final = imdilate(mask_and, strel('disk', 2));

    props = regionprops(mask_final, 'BoundingBox', 'Centroid');
    bbx = vertcat(props.BoundingBox);
    cc = vertcat(props.Centroid);
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

function draw_field(params_s20fe, params_note8)
    figure(3);
    hold on;
    pcshow([[0 0], zeros(size([0 0],1),1)], 'red','VerticalAxisDir', 'up', 'MarkerSize', 5);
    
    plotCamera('Location', params_s20fe.location, 'Orientation', params_s20fe.orientation, 'Size', 0.2, 'Color', [1,0,0]);
    plotCamera('Location', params_note8.location, 'Orientation', params_note8.orientation, 'Size', 0.2, 'Color', [0,1,0]);
    
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

function r = in_rect(pos, rect)
    r = pos(1) >= rect(1) && pos(1) <= rect(1) + rect(3) && ...
        pos(2) >= rect(2) && pos(2) <= rect(2) + rect(4);
end