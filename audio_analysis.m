close all
clear
clc
%% Synchronize videos

% s20fe
[y_s20fe, sample_rate_s20fe] = audioread('videos/s20fe.mp4');
info_audio_s20fe = audioinfo('videos/s20fe.mp4');
v_s20fe = VideoReader('videos/s20fe.mp4');
frames_s20fe = read(v_s20fe,[1 400]);
info_video_s20fe = info(vision.VideoFileReader('videos/s20fe.mp4'));

% note8
[y_note8, sample_rate_note8] = audioread('videos/note8.mp4');
info_audio_note8 = audioinfo('videos/note8.mp4');
v_note8 = VideoReader('videos/note8.mp4');
frames_note8 = read(v_note8,[1 400]);
info_video_note8 = info(vision.VideoFileReader('videos/note8.mp4'));

sample_rate = sample_rate_s20fe; % it is the same
frame_rate = info_video_s20fe.VideoFrameRate; % it is the same
%%
figure(1);
x1 = linspace(0, size(frames_s20fe, 4) / v_s20fe.FrameRate, size(y_s20fe,1));
subplot(2,1,1);
plot(x1, y_s20fe(:, 2)); title('S20FE');
x2 = linspace(0, size(frames_note8, 4) / v_note8.FrameRate, size(y_note8,1));
subplot(2,1,2);
plot(x2, y_note8(:, 2)); title('Note8');

%% take first channel of both and print xcorr
signal_s20fe = y_s20fe(:, 1);
signal_note8 = y_note8(:, 1);

figure(3);
[c,lags] = xcorr(signal_s20fe,signal_note8);
stem(lags,abs(c))

%% take the maximum
[max_lag, index_max_lag] = max(abs(c)); % lag of maximum correlation

shift_between_signals = lags(index_max_lag); % the two signals are shifted by this term

shift_in_seconds = shift_between_signals / sample_rate;

shift_in_frames = round(frame_rate * shift_in_seconds);
%% plot frames of the two videos synchronized
n_frame = 248;

figure(1), imshow(frames_s20fe(:,:,:,n_frame)), title('s20fe');
figure(2), imshow(frames_note8(:,:,:,n_frame - shift_in_frames-2)), title('note8');

%% Results
s20fe_to_s20plus = +43; % nframe_s20plus = nframe_s20fe + 43
s20fe_to_note8 = +95;