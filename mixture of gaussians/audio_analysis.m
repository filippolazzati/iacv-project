close all
clear
clc
%% read audio and video
% sample_rate = frequency with which samples are extracted [Hz]
% y = the audio data; it is a matrix mxn, where m is the number of samples read
% and n is the number of audio channels in the file. % 467966x2

% read audio
path = 'videos/note8.mp4';

[y, sample_rate] = audioread('videos/note8.mp4');
info_audio = audioinfo('videos/note8.mp4');

% read video
v1 = VideoReader('videos/note8.mp4');
frames = read(v1,[1 Inf]); % 1280x720x3x292
info_video = info(vision.VideoFileReader('videos/note8.mp4'));

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
peak_frames = unique(round(peaks_times * info_video.VideoFrameRate), 'stable');

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

%% example of signals synchronization through xcorr
n = 0:15;
x = 0.84.^n;
y = circshift(x,5);
[c,lags] = xcorr(x,y);
stem(lags,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% synchronize s20fe and note8
% s20fe
[y_s20fe, sample_rate_s20fe] = audioread('videos/s20fe.mp4');
info_audio_s20fe = audioinfo('videos/s20fe.mp4');
v_s20fe = VideoReader('videos/s20fe.mp4');
frames_s20fe = read(v_s20fe,[1 Inf]);
info_video_s20fe = info(vision.VideoFileReader('videos/s20fe.mp4'));

% note8
[y_note8, sample_rate_note8] = audioread('videos/note8.mp4');
info_audio_note8 = audioinfo('videos/note8.mp4');
v_note8 = VideoReader('videos/note8.mp4');
frames_note8 = read(v_note8,[1 Inf]); % 1280x720x3x292
info_video_note8 = info(vision.VideoFileReader('videos/note8.mp4'));

sample_rate = sample_rate_s20fe; % it is the same
frame_rate = info_video_s20fe.VideoFrameRate; % it is the same
%%
x1 = linspace(0, info_audio_s20fe.Duration, size(y_s20fe,1));
figure(1); plot(x1, y_s20fe(:, 2)); title('S20');
x2 = linspace(0, info_audio_note8.Duration, size(y_note8,1));
figure(2); plot(x2, y_note8(:, 2)); title('Note8');

%% take first channel of both and print xcorr
signal_s20fe = y_s20fe(:, 1);
signal_note8 = y_note8(:, 1);

[c,lags] = xcorr(signal_s20fe,signal_note8);
% NBBBBBB: I do the abs(c) (NOT SURE ABOUT THIS)
stem(lags,abs(c)) % show

%% take the maximum
[max_lag, index_max_lag] = max(abs(c)); % lag of maximum correlation

shift_between_signals = lags(index_max_lag); % the two signals are shifted by this term
% it is positive -> note8 è avanti rispetto a s20fe di
% shift_between_signals (sample_rate è lo stesso)

shift_in_seconds = shift_between_signals / sample_rate;

shift_in_frames = round(frame_rate * shift_in_seconds);
%% plot frames of the two videos synchronized
n_frame = 248;

figure(1), imshow(frames_s20fe(:,:,:,n_frame)), title('s20fe');
figure(2), imshow(frames_note8(:,:,:,n_frame - shift_in_frames-2)), title('note8');




%%

s20fe_to_s20plus = +43; % aggiungere a nframe dell's20fe per avere quello del s20plus

s20fe_to_note8 = +95;









