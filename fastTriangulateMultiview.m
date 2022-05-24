% triangulateMultiview triangulate 3-D locations of 2-D points matched across multiple views
%   worldPoints = triangulateMultiview(pointTracks, cameraPoses, intrinsics)
%   returns 3-D world points corresponding to points matched across
%   multiple images taken with a calibrated camera.
%
%   Inputs:
%   -------
%   pointTracks      - an N-element array of pointTrack objects, where each
%                      element contains two or more points matched across
%                      multiple images.
%
%   cameraPoses      - a table containing two columns: 'ViewId' and
%                      'AbsolutePose', typically produced by the poses
%                      method of imageviewset. The view IDs in cameraPoses
%                      refer to the view IDs in pointTracks.
%
%   intrinsics       - a scalar or an M-element array of cameraIntrinsics
%                      objects, where M is the number of camera poses. Use
%                      a scalar when images are captured using the same camera
%                      and a vector when images are captured by different cameras.
%
%   Output:
%   -------
%   worldPoints      - an N-by-3 array of [X Y Z] coordinates of the 3-D world
%                      points corresponding to the pointTracks.
%
%   [worldPoints, reprojectionErrors] = triangulateMultiview(...)
%   additionally returns reprojectionErrors, an N-by-1 vector containing
%   the mean reprojection error for each 3-D world point.
%
%   [worldPoints, reprojectionErrors, validIndex] = triangulateMultiview(...)
%   additionally returns the indices of valid world points that are located
%   in front of all the cameras. validIndex is an N-by-1 logical array
%   denoting the validity of each world point.
%
%   Notes
%   -----
%   - The validity of a world point with respect to a camera is determined
%     by projecting the world point onto the image using the camera matrix
%     and homogeneous coordinates. The world point is valid if the resulting
%     scale factor is positive.
%
%   Example 1: Scene reconstruction from multiple views
%   ---------------------------------------------------
%   % Load images
%   imageDir = fullfile(toolboxdir('vision'), 'visiondata', ...
%     'structureFromMotion');
%   images = imageDatastore(imageDir);
%
%   % Load precomputed camera parameters
%   data = load(fullfile(imageDir, 'cameraParams.mat'));
%
%   % Get camera intrinsic parameters
%   intrinsics = data.cameraParams.Intrinsics;
%
%   % Compute features for the first image
%   I = rgb2gray(readimage(images, 1));
%   I = undistortImage(I, intrinsics);
%   pointsPrev = detectSURFFeatures(I);
%   [featuresPrev, pointsPrev] = extractFeatures(I, pointsPrev);
%
%   % Load camera poses
%   load(fullfile(imageDir, 'cameraPoses.mat'));
%
%   % Create an imageviewset object
%   vSet = imageviewset;
%   vSet = addView(vSet, 1, absPoses(1), 'Points', pointsPrev, 'Features', featuresPrev);
%
%   % Compute features and matches for the rest of the images
%   for i = 2:numel(images.Files)
%     I = rgb2gray(readimage(images, i));
%     I = undistortImage(I, intrinsics);
%     points = detectSURFFeatures(I);
%     [features, points] = extractFeatures(I, points);
%     vSet = addView(vSet, i, absPoses(i), 'Points', points, 'Features', features);
%     pairsIdx = matchFeatures(featuresPrev, features, 'MatchThreshold', 5);
%     vSet = addConnection(vSet, i-1, i, 'Matches', pairsIdx);
%     featuresPrev = features;
%   end
%
%   % Find point tracks
%   tracks = findTracks(vSet);
%
%   % Get camera poses
%   cameraPoses = poses(vSet);
%
%   % Find 3-D world points
%   [worldPoints, reprojectionErrors, validIndex] = triangulateMultiview(...
%       tracks, cameraPoses, intrinsics);
%   idx = reprojectionErrors < 5 & worldPoints(:,3)< 20 & validIndex;
%   pcshow(worldPoints(idx, :), 'VerticalAxis', 'y', 'VerticalAxisDir', 'down', ...
%     'MarkerSize', 30);
%   hold on
%   plotCamera(cameraPoses, 'Size', 0.2);
%   hold off
%
%   Example 2: Structure from motion from multiple views
%   ----------------------------------------------------
%   % This example shows you how to estimate the poses of a calibrated
%   % camera from a sequence of views, and reconstruct the 3-D structure of
%   % the scene up to an unknown scale factor.
%   % <a href="matlab:helpview(fullfile(docroot,'toolbox','vision','vision.map'),'StructureFromMotionFromMultipleViewsExample')">View example</a>
%
%   See also triangulate, imageviewset, pointTrack, matchFeatures, rigid3d,
%       bundleAdjustment, bundleAdjustmentStructure, bundleAdjustmentMotion.

% Copyright 2015-2020 MathWorks, Inc.
%
% References:
%
% Hartley, Richard, and Andrew Zisserman. Multiple View Geometry in
% Computer Vision. Second Edition. Cambridge, 2000. p. 312

function [worldPoints, reprojectionErrors, validIndex] = fastTriangulateMultiview(pointTracks, ...
    camPoses, intrinsics)

%hasAbsolutePose = validateInputs(pointTracks, camPoses, intrinsics);
hasAbsolutePose = 0;

numCameras = size(camPoses, 1);
cameraMatrices = cell(numCameras, 1);
for i = 1:numCameras
    id = camPoses{i, 'ViewId'};
    if hasAbsolutePose
        R = camPoses{i, 'AbsolutePose'}.Rotation;
        t = camPoses{i, 'AbsolutePose'}.Translation;
    else
        R  = camPoses{i, 'Orientation'}{1};
        t  = camPoses{i, 'Location'}{1};
    end
    if numel(intrinsics) > 1 % cameraIntrinsics type
        cameraMatrices{id} = cameraMatrix(intrinsics(i), R', -t*R');
    else % cameraIntrinsics or cameraParameters type
        cameraMatrices{id} = cameraMatrix(intrinsics, R', -t*R');
    end
end

worldPoints = vision.internal.triangulateMultiViewPoints(pointTracks, cameraMatrices);
outputType = class(worldPoints);

if nargout > 1
    [reprojectionErrors, validIndex] = computeReprojectionErrors(worldPoints, cameraMatrices, pointTracks);
    reprojectionErrors = cast(reprojectionErrors, outputType);
end

%--------------------------------------------------------------------------
function hasAbsolutePose = validateInputs(pointTracks, camPoses, intrinsics)

validateattributes(pointTracks, {'pointTrack'}, {'nonempty'}, mfilename);

hasAbsolutePose = validateAbsolutePoses(camPoses);

validateCameraIntrinsics(intrinsics);

% Check the camera array
if ~isscalar(intrinsics) && numel(intrinsics) ~= height(camPoses)
    error(message('vision:sfm:unmatchedParamsPoses'));
end

%--------------------------------------------------------------------------
function hasAbsolutePose = validateAbsolutePoses(camPoses)

% Validate the table
validateattributes(camPoses, {'table'},{'nonempty'}, mfilename, 'cameraPoses');

% Check columns
hasAbsolutePose = ismember('AbsolutePose', camPoses.Properties.VariableNames);
if hasAbsolutePose
    vision.internal.inputValidation.validatePoseTableRigid3d(...
        camPoses, mfilename, 'cameraPoses');
else
    vision.internal.inputValidation.checkAbsolutePoses(...
        camPoses, mfilename, 'cameraPoses');
end

%--------------------------------------------------------------------------
function validateCameraIntrinsics(intrinsics)
% cameraParameters is supported but not recommended or documented and it
% can only be a scalar due to restriction in its constructor.
vision.internal.inputValidation.checkIntrinsicsAndParameters( ...
    intrinsics, false, mfilename);

%--------------------------------------------------------------------------
function [meanErrorsPerTrack, isInFrontOfAllCameras] = ...
    computeReprojectionErrors(points3d, cameraMatrices, tracks)

numPoints = size(points3d, 1);
points3dh = [points3d, ones(numPoints, 1)];
meanErrorsPerTrack = zeros(numPoints, 1);
isInFrontOfAllCameras = false(numPoints, 1);

for i = 1:numPoints
    p3d = points3dh(i, :);
    [reprojPoints2d, isInFrontOfAllCameras(i)] = reprojectPoint(...
        p3d, tracks(i).ViewIds, cameraMatrices);
    e = sqrt(sum((tracks(i).Points - reprojPoints2d).^2, 2));
    meanErrorsPerTrack(i) = mean(e);
end

%--------------------------------------------------------------------------
function [points2d, isInFrontOfAllCameras] = reprojectPoint(p3dh, viewIds, cameraMatrices)
num2DPoints = numel(viewIds);
points2d = zeros(num2DPoints, 2);
isInFrontOfCamera = false(num2DPoints, 1);
for i = 1:num2DPoints
    p2dh = p3dh * cameraMatrices{viewIds(i)};
    points2d(i, :) = p2dh(1:2) ./ p2dh(3);
    isInFrontOfCamera(i) = p2dh(3) > 0;
end
isInFrontOfAllCameras = all(isInFrontOfCamera);