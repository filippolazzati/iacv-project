%% Localization

%O = [5.485 5.49];
%P1 = O + [4.115 0];
%P2 = O + [4.115 6.40];
%P3 = O + [0 6.40];

O = [10.97 11.89*2];
P1 = O + [-1.37 0];
P2 = O + [-1.37 -11.89];
P3 = O + [0 -11.89];

img = imread('calibration data/note8_single_photo.jpg');
imshow(img);
figure(1), imshow(img);
title('original image');
hold on;
[x, y] = getpts();
hold off;

[rotationMatrix,translationVector] = extrinsics([x y], [O; P1; P2; P3], cameraParams);

orientation = rotationMatrix';
location = -translationVector * orientation;