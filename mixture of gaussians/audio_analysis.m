close all
clear
clc
%% read audio and video
% sample_rate = frequency with which samples are extracted [Hz]
% y = the audio data; it is a matrix mxn, where m is the number of samples read
% and n is the number of audio channels in the file. % 467966x2

% read audio
[y, sample_rate] = audioread('../mydata/note8.mp4');
info_audio = audioinfo('../mydata/note8.mp4');

% read video
v1 = VideoReader('../mydata/note8.mp4');
frames = read(v1,[1 Inf]); % 1280x720x3x292
info_video = info(vision.VideoFileReader('../mydata/note8.mp4'));

%% plot the channels by sampling instant
line_value = 0.2;
x = linspace(0, info_audio.Duration, size(y,1));
figure(1); plot(x, y(:, 1)); title('channel 1 (y = intensity, x = seconds)');
figure(2); plot(x, y(:, 2)); title('channel 2 (y = intensity, x = seconds)');
hold on;
plot(x, line_value * ones(1, size(y,1)), 'Color', 'red');
hold off;
%% Show frames corresponding to peaks in the audio signal
audio_intensity = 0.3;
ch1 = y(:, 1);
indexes = find(ch1 > audio_intensity); % find peaks' indexes
peaks_times = indexes / sample_rate; % peaks' instants
peak_frames = unique(floor(peaks_times * info_video.VideoFrameRate), 'stable');

imgs = frames(:, :, :, peak_frames(1));
for i = 2:size(peak_frames, 1)
    imgs = horzcat(imgs, frames(:, :, :, peak_frames(i)));
end

figure(3), imshow(imgs); title('images of impact');


%% read portion of audio (first 2 seconds)
samples = [1,2*sample_rate]; % read first 2 seconds since sample_rate = number of readings in 1s
[y,sample_rate] = audioread('../mydata/s20fe.mp4',samples);

%% play audio
sound(y, sample_rate);
