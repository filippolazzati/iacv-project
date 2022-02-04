clc 
clear all
close all
%%

% open the video
v1 = VideoReader('data/cetto.mp4'); % 405 frames, 720x1280
% read all the frames
frames = read(v1,[1 Inf]); % 4D array
% returns the width and height of frames in video and number of coor planes
[height,width,N]= size(frames(:,:,:,1));

% create avifile object to save results as video
 vidObj= VideoWriter('1deltafr with Th & cond.avi');

% start stopwatch to measure time code needs to run
tic

% Variables
    K= 3; % number of gaussian components (3-5)
    M=3; %number of background components
    alpha= 0.05; % learning rate % alpha =0.01
    Thresh=0.25; % foreground threshold 0.25
    var_init=900; % initial variance
    rho= 0.03; % update coeff

% Initialize gaussian component means, standard deviation and weights
% randomly
    for i=1:height
        for j=1:width
            for k=1:K

                mean(i,j,k) = rand*255;             % means random (0-255) 255=pixel range
                w(i,j,k) = 1/K;                     % weights uniformly dist
                var(i,j,k) = var_init;                % initialize variance

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
close(vidObj);
toc
%Framepersecond= (numFrames/toc)
%    figure(1),subplot(3,1,1),imshow(X_t)
%     subplot(3,1,2),imshow(uint8(bg_bw))
%     subplot(3,1,3),imshow(uint8(fg)) 
%     
%     Mov1(z)  = im2frame(uint8(fg),gray);           % put frames into movie
%     %Mov2(z)  = im2frame(uint8(bg_bw),gray);           % put frames into movie
%     
% end
%       
% movie2avi(Mov1,'GMM1','fps',25);           % save movie as avi 
% %movie2avi(Mov2,'MOG_background5','fps',25);           % save movie as avi 


% stop stopwatch
toc