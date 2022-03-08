close all
clear all
clc

%%
params_s20fe = load("calibration/params/s20fe.mat");
params_note8 = load("calibration/params/note8.mat");
params_s20plus = load("calibration/params/s20plus.mat");

v_s20fe = VideoReader('videos/s20fe.mp4');
v_note8 = VideoReader('videos/note8.mp4');
v_s20plus = VideoReader('videos/s20plus.mp4');
%%
base_frame = 303;
s20fe_to_s20plus = +43;
s20fe_to_note8 = +95;

img_s20fe = read(v_s20fe, base_frame);
figure(1); imshow(img_s20fe); hold on; [x, y] = getpts();
pts_s20fe = [x y];
img_note8 = read(v_note8, base_frame+s20fe_to_note8);
figure(2); imshow(img_note8); hold on; [x, y] = getpts();
pts_note8 = [x y];
img_s20plus = read(v_s20plus, base_frame+s20fe_to_s20plus);
figure(3); imshow(img_s20plus); hold on; [x, y] = getpts();
pts_s20plus = [x y];

%%
track = pointTrack([1; 2; 3], [pts_s20fe; pts_note8; pts_s20plus]);
poses = table;
poses.ViewId = uint32([1; 2; 3]);
poses.Orientation = {params_s20fe.orientation; params_note8.orientation; params_s20plus.orientation};
poses.Location = {params_s20fe.location; params_note8.location; params_s20plus.location};
intrinsics = [params_s20fe.cameraParams.Intrinsics, params_note8.cameraParams.Intrinsics, params_s20plus.cameraParams.Intrinsics];
[wp, re] = triangulateMultiview(track, poses, intrinsics);

%{
[wp, re] = triangulate(pts_note8, pts_s20plus, ...
    cameraMatrix(params_note8.cameraParams, params_note8.rotationMatrix, params_note8.translationVector), ...
    cameraMatrix(params_s20plus.cameraParams, params_s20plus.rotationMatrix, params_s20plus.translationVector) ...
);
%}
%%
figure(4);
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

plot3(wp(1), wp(2), wp(3), '.', 'MarkerSize', 30, 'Color', 'yellow');

set(gcf,'color','w');
set(gca,'color','w');
xlabel('X');
ylabel('Y');
zlabel('Z');

