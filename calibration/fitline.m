% FITLINE - Least squares fit of a line to a set of points
%
% XY = matrix whose columns are points of R^2 or P^2 (so XY has 2 or 3
% rows, by default if it has 2 rows this function adds a lot of 1s on the
% last row).

% what fitline does is imply applying the DLT algorithm but instead of
% looking for an homography (3x3) we are looking for a single vector (3x1)
% because we are looking for a line (3 parameters in homogeneous
% coordinates).

function [C, dist] = fitline(XY)
  
% rows = number of rows of XY (so if the points have 2 or 3 components);
% npts = number of columns of XY (so how many points we have).
  [rows,npts] = size(XY);    

  %we need at least 2 points to fit a line
  if npts < 2
      error('Too few points to fit line');
  end  
  
  % if the points have 2 coordinates -> they are in R^2 -> add 1 to move
  % them to P^2.
  if rows ==2    % Add homogeneous scale coordinate of 1 
      XY = [XY; ones(1,npts)];
  end
  
  if npts == 2    % Pad XY with a third column of zeros
    XY = [XY zeros(3,1)]; 
  end
  
  % Normalise points so that centroid is at origin and mean distance from
  % origin is sqrt(2). (it is the same of calling 
  % precond() function.
  [XYn, T] = normalise2dpts(XY);
  
  % IN THIS CASE, OUR MATRIX A IS THE MATRIX XY, BECAUSE XY * x = 0,
  % THEREFORE TRYING TO MINIMIZE THE L2 NORM OF argmin||XY*x|| (because we
  % want XY*x=0) WE OBTAIN THE SAME SOLUTION OF DLT, THAT IS THE LAST
  % COLUMN OF MATRIX V !!!!
  
  % Take C = [c1, c2, c3] that defines the line. Since all the points
  % belong to line C, we can write: XYn'*C = 0, namely 
  % c(1)*X + c(2)*Y + c(3) = 0, which gives an equation (a row) for each
  % point.
  [u d v] = svd(XYn',0);   % Singular value decomposition. The 0 means that
  % we want the thin SVD.
  
  C = v(:,3);              % Solution is last column of v.

  % Denormalise the solution
  C = T'*C;

  % If requested, calculate the distances from the fitted line to
  % the supplied data points 
  if nargout==2   
      dist = abs(C(1)*XY(1,:) + C(2)*XY(2,:) + C(3));
  end

end

