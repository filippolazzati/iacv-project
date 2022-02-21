close all;
clear;
clc;
addpath(genpath('functions'));
%% read video and compute background
% open the video
v1 = VideoReader('../data/video2_Trim.mp4');
% read all the frames
frames = read(v1, [1 Inf]); % 4D array 720x1280x3x309 
% median of first n frames
n = 100;
mean = median(frames(:, :, :, 1:n),4);
%background_g = rgb2gray(background);
figure(1);
imshow(mean);
%
[height,width,N]= size(frames(:,:,:,1));
%%
% start stopwatch to measure time code needs to run
tic
% Variables
K= 3; % number of gaussian components (3-5)
M=3; %number of background components (B?)
alpha= 0.05; % learning rate % alpha =0.01
Thresh=0.25; % foreground threshold 0.25
var_init=900; % initial variance
rho= 0.03; % update coeff

w = nan(height,width,k);
var = nan(height,width,k);
% Initialize gaussian component means, standard deviation and weights
% randomly
    for i=1:height
        for j=1:width
            for k=1:K         
                w(i,j,k) = 1/K;                     % weights uniformly dist
                var(i,j,k) = var_init;                  % initialize variance
            end
        end
    end

for z=1:size(frames, 4)
    
   % take the current frame and convert to grayscale
    X_t = double(rgb2gray(frames(:,:,:,z)));
    
   % calculate the difference of pixel values from mean  
    for i=1:height
        for j=1:width
            for k=1:K
                diff(i,j,k)=abs(X_t(i,j) - mean(i,j,k));
            end
        end
    end
    % D = pdist2(X_t,X_t,'mahalanobis');

    for i=1:height
        for j=1:width
            match=0;
            for k=1:K
                if diff(i,j,k)<= 2.5*sqrt(var(i,j,k))
                    % update
                    match=1;
                    % I am increasing the probability (of belonging to BG)
                    w(i,j,k)=(1-alpha)*w(i,j,k)+alpha;
                    rho=alpha/w(i,j,k);
                    mean(i,j,k)=(1-rho)*mean(i,j,k)+rho*double(X_t(i,j));
                    var(i,j,k)=(1-rho)*var(i,j,k)+rho*((double(X_t(i,j)-mean(i,j,k)))^2);
                else
                    % decrease the probability (weight)
                    w(i,j,k)=(1-alpha)*w(i,j,k);
                end
            end

    % normalize weights to be interpreted as probabilities
    w(i,j,:)=w(i,j,:)./sum(w(i,j,:));
    % initialize this variable to the sum of all the means
    bg_bw(i,j)=0;
    for k=1:K
           bg_bw(i,j) = bg_bw(i,j)+ mean(i,j,k)*w(i,j,k);
    end

 % if no components match, create new component instead of the least
 % probable distribution, with current value as mean, initial variance and
 % low prior weight (low w)
    if (match == 0)
       % take the least probable distribution and lowest weight w (min
       % returns both the minimum of the pixel in the K gaussians and also
       % the minimum index
       [min_w, min_w_index] = min(w(i,j,:));  
       % set its mean to the current value
       mean(i,j,min_w_index) = double(X_t(i,j));
       % set variance to var_init
       var(i,j,min_w_index) = var_init;
    end

    % compute rank:=w/sigma
    rank = w(i,j,:)./sqrt(var(i,j,:));
    rank_ind = [1:1:K]; % array of numbers from 1 to K with step 1
            
    % sort rank values in decreasing order
    for k=2:K               
        for m=1:(k-1)
            
            if (rank(:,:,k) > rank(:,:,m))                     
                % swap max values
                rank_temp = rank(:,:,m);  
                rank(:,:,m) = rank(:,:,k);
                rank(:,:,k) = rank_temp;
                % swap max index values
                rank_ind_temp = rank_ind(m);  
                rank_ind(m) = rank_ind(k);
                rank_ind(k) = rank_ind_temp;    

            end
        end
    end

   % calculate foreground
    match = 0;
    k=1; % start from the maximum w/sigma
    % fg is the foreground
    fg(i,j) = 0;
    while ((match == 0)&&(k<=M))
    
        if (w(i,j,rank_ind(k)) >= Thresh)
            if (diff(i,j,rank_ind(k)) <= 2.5*sqrt(var(i,j,rank_ind(k))))
                fg(i,j) = 0;
                match = 1;
            else
                fg(i,j) = X_t(i,j);
            end
        end
        k = k+1;
    end

    end
    end

   %median filtering
   %fg3 = medfilt2(fg); % performs median filtering of the matrix A using the default 3-by-3 neighborhood,Padded with 0s. This is the default.

     %filtrage morphologique 
     se = strel('square',2); % take a 2x2 morphological element
     fg2=imopen(fg,se);
     fg3=imclose(fg2,se);
    
    %Otsu's threshold
    %fg3=fg2./255;
%     level = graythresh(fg3);
%     if level<0.05
%         level=0.10;
%     end
%     fg3 = im2bw(fg3,level);
 %fg3=fg;

  % figure(1);
   %subplot(2,1,1);
  % imshow(X_t);
   %subplot(2,1,2);
   %imshow(fg3);
   
% where you store results

%pathname = 'C:\Results all\Teste\';
disp('1')
   %imwrite(fg3,[pathname,['bin', num2str(z, '%.6d'),'.png']]);
   %hold on
   open(vidObj);
   % Add next frame to movie
   newFrameOut= getframe;
   writeVideo(vidObj,newFrameOut);
end
%%
close all;
clear;
clc;
addpath(genpath('functions'));
%% read video and compute background
% open the video
v1 = VideoReader('../mydata/s20fe.mp4'); % 405 frames, 720x1280
% read all the frames
frames = read(v1, [1 Inf]); % 4D array
background = median(frames, 4);
background_g = rgb2gray(background);
detector = vision.ForegroundDetector(...
       'NumTrainingFrames', 50, ...
       'NumGaussians', 5, ...
       'MinimumBackgroundRatio', 0.7, ...
       'InitialVariance', 400);
%%
Se = strel('disk', 10);
S1 = strel('square', 2);
hollow = ones(20);
hollow(2:19, 2:19) = 0;
S2 = strel('arbitrary', hollow);
init_frame = 9;
n_frame_diff = 2;
n_frame_diff2 = 1;
mask_thresh = 30;

nframe = init_frame;
pts = nan(size(frames, 4)-nframe, 2);
while nframe < size(frames, 4)
    j = nframe - init_frame + 1;
    f = frames(:,:,:,nframe);
    % background subtraction mask
    mask_bg_sub2 = detector(f);
    f_g = rgb2gray(f);
    f_prev_g = rgb2gray(frames(:,:,:,nframe-n_frame_diff));
    f_prev2_g = rgb2gray(frames(:,:,:,nframe-n_frame_diff2));

    bw_bg_sub = imclose(imabsdiff(f_g, background_g), Se);
    mask_bg_sub = bw_bg_sub > mask_thresh;
    
    mask_bg_sub = mask_bg_sub & ~bwareaopen(mask_bg_sub, 50);
    mask_and = mask_bg_sub & mask_bg_sub2;

    for df = [-8, -6, -4, +4, +6, +8]
        bw_frame_diff = imclose(imabsdiff(f_g, f_prev_g), Se);
        mask_frame_diff = bw_frame_diff > mask_thresh;
        mask_frame_diff = mask_frame_diff & ~bwareaopen(mask_frame_diff, 50);
        mask_and = mask_and & mask_frame_diff;
    end
    mask_hit_miss = bwhitmiss(mask_and, S1, S2);
    mask_final = imdilate(mask_hit_miss, strel('disk', 4));

    figure(1); imshow(f); hold all;
    props = regionprops(mask_final, 'BoundingBox', 'Centroid');
    bbx = vertcat(props.BoundingBox);
    %ecc = vertcat(props.Eccentricity);
    cc = vertcat(props.Centroid);

    if (size(bbx, 1) > 0)
        for i = 1:size(bbx, 1)
            if bbx(i,3) > 0 && bbx(i,4) > 0
                rectangle('Position',[bbx(i,1),bbx(i,2),bbx(i,3),bbx(i,4)], 'EdgeColor','r','LineWidth',2);
                pts(j, :) = cc(i,:);
                %text(bbx(i,1), bbx(i,2) - 20, num2str(ecc(i)));
            end
        end
    end

    for p = 1:size(pts, 1)
        pt = pts(p,:);
        if not(isnan(pt(1)))
            plot(pt(1), pt(2), '.', 'MarkerSize', 10, 'Color', 'green');
        end
    end
    
    nframe = nframe + 1;
end














