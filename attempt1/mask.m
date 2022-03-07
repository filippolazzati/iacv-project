close all;
clear;
clc;
addpath(genpath('functions'));
%%

v1 = VideoReader('mydata/s20fe.mp4');

frames1 = read(v1, [1 Inf]);

background1_bw = rgb2gray(median(frames1, 4));

%%
init_frame_1 = 9;

mask_thresh = 30;

nframe = init_frame_1;

while nframe < size(frames1, 4)
    [bbx1, centroids1, mask_bg_sub, mask_frame_diff, mask_final] = ball_candidates(frames1, nframe, background1_bw, mask_thresh);

    figure(1); hold all;
    %imshow([frames1(:,:,:,nframe), imbin2rgb(mask_bg_sub); imbin2rgb(mask_frame_diff), imbin2rgb(mask_final)]);
    imshow(frames1(:,:,:,nframe));
    draw_bbx(bbx1);

    nframe = nframe + 1;
end

%%
function [bbx, cc] = ball_candidates(frames, nframe, background_bw, mask_thresh)
    frame = frames(:,:,:,nframe);
    frame_bw = rgb2gray(frame);

    bw_bg_sub = imclose(imabsdiff(frame_bw, background_bw), strel('disk', 3));
    %bw_bg_sub = imabsdiff(frame_bw, background_bw);
    mask_bg_sub = bw_bg_sub > mask_thresh;
    mask_bg_sub = mask_bg_sub & ~bwareaopen(mask_bg_sub, 200);

    mask_frame_diff = ones(size(frame, 1), size(frame, 2));

    for df = [-8, -6, -4, +4, +6, +8]
        frame_idx = nframe + df;
        if frame_idx < 0 || frame_idx >= size(frames, 4)
            continue
        end
        other_frame_bw = rgb2gray(frames(:,:,:,frame_idx));
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

function draw_bbx(bbx)
    if (size(bbx, 1) > 0)
        for i = 1:size(bbx, 1)
            if bbx(i,3) > 0 && bbx(i,4) > 0
                rectangle('Position',[bbx(i,1),bbx(i,2),bbx(i,3),bbx(i,4)], 'EdgeColor','r','LineWidth',2);
            end
        end
    end
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
