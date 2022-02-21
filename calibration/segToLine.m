function [l] = segToLine(pts)
% pts = a 2-by-2 numeric matrix of the form [x1 y1; x2 y2] that represent a
%       line, where (x1;y1) = initial point, (x2;y2) = final point

% we take the points in homogeneous coordinates -> we simply add a 1 as
% last coordinate:

% initial point
a = [pts(1,:)';1];
% final point
b = [pts(2,:)';1];

% line through them is the cross product
l = cross(a,b);

% normalize it
l = l./norm(3);
end