close all
clear all
clc

%%
params_s20fe = load("calibration data/s20fe/calibrationparams.mat");
params_note8 = load("calibration data/note8/calibrationparams.mat");

%%
img_s20fe = imread('calibration data/s20fe_Moment.jpg');
figure(1); imshow(img_s20fe); hold on; [x, y] = getpts();
pts_s20fe = [x y];
img_note8 = imread('calibration data/note8_Moment.jpg');
figure(2); imshow(img_note8); hold on; [x, y] = getpts();
pts_note8 = [x y];

[wp, re] = triangulate(pts_s20fe, pts_note8, ...
    cameraMatrix(params_s20fe.cameraParams, params_s20fe.rotationMatrix, params_s20fe.translationVector), ...
    cameraMatrix(params_note8.cameraParams, params_note8.rotationMatrix, params_note8.translationVector) ...
);

%%
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

plot3(wp(1), wp(2), wp(3), '.', 'MarkerSize', 30, 'Color', 'yellow');

set(gcf,'color','w');
set(gca,'color','w');
xlabel('X');
ylabel('Y');
zlabel('Z');

