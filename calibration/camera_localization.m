%% Localization
close all

% Note8
P1 = [5.485 5.49];
P2 = [9.6 5.49];
P3 = [9.6 11.89];
P4 = [5.485 11.89];

% S20fe
%P1 = [5.485 5.49];
%P2 = [1.37 5.49];
%P3 = [1.37 11.89];
%P4 = [5.485 11.89];

% S20Plus
%P1 = [10.97 2*11.89];
%P2 = P1 - [1.37 0];
%P3 = P2 - [0 11.89];
%P4 = P1 - [0 11.89];

v = VideoReader('videos/note8.mp4');
frame = read(v, 1);

imshow(frame);

%%
hold on;
[x, y] = getpts();
hold off;

[rotationMatrix,translationVector] = extrinsics([x y], [P1; P2; P3; P4], cameraParams);

orientation = rotationMatrix';
location = -translationVector * orientation;