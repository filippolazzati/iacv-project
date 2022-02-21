close all
clear
clc

%% 1) Affine rectification ************************************************
img = imread('../mydata/calibration data/s20fe.jpg');
whos img
figure(1);
imshow(img);
% select the number of directions and lines for each direction
directions = 2;
lines = 2;
% create a cell array in which store the parallel lines
parallelLines = cell(directions,1);
fprintf(['Draw ', num2str(directions) , ' families of parallel segments\n']);
col = 'rgbm';
% for each direction, draw all the parallel lines along it
for i = 1:directions
count = 1;
parallelLines{i} = nan(lines,3);
while(count <=lines)
figure(gcf);
title(['Draw ', num2str(lines),' segments: step ',num2str(count) ]);
segment1 = drawline('Color',col(i));
% segToLine takes two points of a line and returns the projective
% representation of the line
parallelLines{i}(count, :) = segToLine(segment1.Position);
count = count +1;
end
fprintf('Press enter to continue\n');
pause
end
V = nan(2,directions);
% for each direction, find the vanishing point by applying least squares
for i = 1:directions
A = parallelLines{i}(:,1:2);
b = - parallelLines{i}(:,3);
V(:,i) = A\b;
end
%
imLinfty = fitline(V);
imLinfty = imLinfty./(imLinfty(3));
H_affine = [eye(2),zeros(2,1); imLinfty(:)'];
%use imwarp
tform = projective2d(H_affine');
aff_img_imwarp = imwarp(im2double(img),tform,'OutputView',imref2d(size(img)));
figure(1);
imshow(aff_img_imwarp);
% 2) Metric rectification
figure(1);
hold all;
fprintf('Draw pairs of orthogonal segments\n');
numConstraints = 2;
count = 1;
A = zeros(numConstraints,3);
% select pairs of orthogonal segments
while (count <=numConstraints)
    figure(gcf);
    title('Select pairs of orthogonal segments')
    col = 'rgbcmykwrgbcmykw';
    segment1 = drawline('Color',col(count));
    segment2 = drawline('Color',col(count));

    l = segToLine(segment1.Position);
    m = segToLine(segment2.Position);

    % each pair of orthogonal lines gives rise to a constraint on s
    % [l(1)*m(1),l(1)*m(2)+l(2)*m(1), l(2)*m(2)]*s = 0
    % store the constraints in a matrix A
     A(count,:) = [l(1)*m(1),l(1)*m(2)+l(2)*m(1), l(2)*m(2)];

    count = count+1;
end
%S = [x(1) x(2); x(2) 1];
[~,~,v] = svd(A);
s = v(:,end); %[s11,s12,s22];
S = [s(1),s(2); s(2),s(3)];
imDCCP = [S,zeros(2,1); zeros(1,3)]; % the image of the circular points
[U,D,V] = svd(S);
A = U*sqrt(D)*V';
H = eye(3);
H(1,1) = A(1,1);
H(1,2) = A(1,2);
H(2,1) = A(2,1);
H(2,2) = A(2,2);

Hrect = inv(H);
Cinfty = [eye(2),zeros(2,1);zeros(1,3)];

tform = projective2d(Hrect');

J = imwarp(aff_img_imwarp,tform);

figure;
imshow(J);
%%
close all
clear
clc

%% Camera calibration
linfty = [5.566577052383796e-04;1.438701790091674e-04;1];
Hrect = [2.802166479330810,0.623371618707219,0;0.623371618707221,1.197333776314424,0;0,0,1];
H = inv(Hrect);
h1 = H(:,1);
h2 = H(:,2);
h3 = H(:,3);
% extract the vertical vanishing point
img = imread('../mydata/calibration data/s20fe.jpg');
figure;
imshow(img);
title(['Draw 2 vertical parallel lines']);
segment1 = drawline('Color','r');
segment2 = drawline('Color','r');
line1 = segToLine(segment1.Position);
line2 = segToLine(segment2.Position);
% compute the vanishing point
v = cross(line1,line2);
v = v ./ v(3);
% take two lines to use to take two vanishing points on the horizontal plane:
title(['Draw 2 horizontal non-parallel lines']);
segment3 = drawline('Color','b');
segment4 = drawline('Color','b');
line3 = segToLine(segment3.Position);
line4 = segToLine(segment4.Position);
line3 = line3 ./ line3(3);
line4 = line4 ./ line4(3);
point1 = cross(line3,linfty);
point2 = cross(line4,linfty);
% solve the system
%%
syms x1;
syms x2;
syms x3;
syms x4;
omega = [x1,0,x2;0,1,x3;x2,x3,x4];
eq1 = point1.'*omega * v==0;
eq2 = point2.'*omega * v==0;
eq3 = h1'*omega*h2 == 0;
eq4 = h1'*omega*h1 == h2'*omega*h2;
eqns = [eq1,eq2,eq3,eq4];
S = solve(eqns, [x1,x2,x3,x4]);
sx1 = double(S.x1);
sx2 = double(S.x2);
sx3 = double(S.x3);
sx4 = double(S.x4);
IAC = [sx1,0,sx2;0,1,sx3;sx2,sx3,sx4]
% retrieve the parameters of K
a = sqrt(IAC(1,1));
u0 = - IAC(1,3) / (a*a);
v0 = - IAC(2,3);
fy = sqrt(IAC(3,3)-(a*a*u0*u0)-v0*v0);
fx = fy/a;
K = [fx,0,u0;0,fy,v0;0,0,1]


