%% All Parts of the term Project %%

clear all;
clc;

%%%%%%%%%%%%%%%%%%%%
%% %%% Part 1 %%% %%
%%%%%%%%%%%%%%%%%%%%

load('106m.mat'); 
% Loads the data from the '106m.mat' file,
% which contains raw ECG signal information for record '106' from the MIT-BIH Arrhythmia Database.

% Raw Signal Data
raw_signal = val(1, :);     % Extracts the raw signal from the loaded data.

% The following information are given in the info file along with the data
fs = 360;       % Given Sampling Frequency in Hz
gain = 200;     % Given Gain

% Calculating Amplitude to represent actual physiological signal in mV
% the amplitude values to represent the actual physiological signal in millivolts (mV)
% by dividing the raw signal by the gain.
amp_signal = raw_signal / gain;

% Calculating time for each data point
% Creates a time vector corresponding to each data point.
t = (0:(length(amp_signal)-1)) / fs; 

% data array containing time and corresponding amplitudes
data_array = [t', amp_signal'];

% Extracting the 2nd '15 minutes' of the dataset
start_time = 15 * 60;  % 15 minutes in seconds
end_time = 30 * 60;    % 30 minutes in seconds

% indices corresponding to the time range (second 15 minutes)
% conditions for the time interval of interest-30 minutes)

indices_15to30 = find((t >= start_time) & (t <= end_time));


% new data frame for the second 15 minutes (15-30 minutes)
df = data_array(indices_15to30, :);

%time = (0:(length(amplitude)-1)) / fs;
time = df(:, 1); % times (x)
amplitude = df(:, 2); % amplitudes
z_amp = amplitude - mean(amplitude);

%Plotting ECG signal for 15-30 minutes of data points

figure(1);
plot(time, amplitude);
xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('ECG Signal');

%%%%%%%%%%%%%%%%%%%%
%%%%%% Part 2 %%%%%%
%%%%%%%%%%%%%%%%%%%%

% Finding peaks to identify the number of heartbeats
% Number of Peaks = Number of Heartbeats (QRS Peaks)
% Using in-built function findpeaks

[peaks, peak_indices] = findpeaks(z_amp, time, 'MinPeakHeight', 0.2, 'MinPeakDistance', 0.4);
% 'MinPeakHeight', 0.2, 'MinPeakDistance', 0.5
% these thresholds can be changed
% I chose these, because mostly the QRS peaks have an amplitude of 0.2 mV and
% interval of at least 0.4 sec

heartbeats = length(peaks);
fprintf('Number of heartbeats (RAW ECG): %d\n', heartbeats);
avg_HBs_per_minute = heartbeats / (30 - 15);
fprintf('Average Number of heartbeats per minute (Raw ECG): %d\n', round(avg_HBs_per_minute));

%%%%%%%%%%%%%%%%%%%%
%%%%%% Part 3 %%%%%%
%%%%%%%%%%%%%%%%%%%%

% % Calculating the RR intervals
rr_intervals = diff(peak_indices); % difference for the adjacent QRS peaks = RR interval

% Displaying the RR intervals
%fprintf('RR Intervals (in seconds):\n');
%disp(rr_intervals);


%%%%%%%%%%%%%%%%%%%%
%%%%%% Part 4 %%%%%%
%%%%%%%%%%%%%%%%%%%%

% Plotting Histogram of RR intervals
figure(2);
nbins = 100;
histogram(rr_intervals, nbins);
xlabel('RR (seconds)');
ylabel('Number of Intervals');
title('Histogram of RR Intervals');
 
% %%%%%%%%%%%%%%%%%%%%
% %%%%%% Part 5 %%%%%%
% %%%%%%%%%%%%%%%%%%%%

% RR Interval as a function of time
% Calculating cumulative sum of RR intervals for the x-axis

new_time = cumsum(rr_intervals);
RR = [new_time,rr_intervals];

% Plotting RR-intervals over time
figure(3);
plot(new_time, rr_intervals);
title('RR-intervals Over Time');
xlabel('Time (s)');
ylabel('RR-interval (s)');

% %%%%%%%%%%%%%%%%%%%%
% %%%%%% Part 6 %%%%%%
% %%%%%%%%%%%%%%%%%%%%

% Interpolation
rr_interpolated_uniform = interp1(new_time, rr_intervals, 'cubic', 'pp');

% Uniformly sample the time vector
% number of samples = 900 (can be changed as needed)

uniformly_sampled_time = linspace(new_time(1), new_time(end), 900); 

% Evaluating the interpolated RR intervals at uniformly sampled time points

sampled_rr_evaluated = ppval(rr_interpolated_uniform, uniformly_sampled_time);

figure(4);
plot(uniformly_sampled_time, sampled_rr_evaluated);
title('Uniformly Interpolated RR (900 samples)');
xlabel('Time (s)');
ylabel('RR-interval (s)');

% %%%%%%%%%%%%%%%%%%%%
% %%%%%% Part 7 %%%%%%
% %%%%%%%%%%%%%%%%%%%%

% % Defining sampling rate and steps
fs2 = 128.0;
steps = 1 / fs2;

% Sampling from interpolation function
samples_128 = 1:steps:max(new_time);
rr_interpolated_128 = ppval(rr_interpolated_uniform, samples_128);

% Plotting original and interpolated RR intervals for whole 15 minutes
figure(5);
% Size
set(gcf, 'Position', [80, 80, 1000, 600]);
subplot(2, 1, 1);
plot(new_time, rr_intervals);
xlabel('Time (s)');
ylabel('RR-interval (s)');
title('Original RR Intervals');

subplot(2, 1, 2);
plot(samples_128, rr_interpolated_128);
xlabel('Time (s)');
ylabel('RR-interval (s)');
title('Interpolated RR Intervals (@ 128 samples/s');
sgtitle('Comparison of Original and Interpolated RR Intervals (Full 15 minutes)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The plot 15 minutes, it is difficult to distinguish the difference
% between the original and interpolated-sampled RR intervals plot
% Therefore, Plotting only for 100 seconds to compare the differences


% Plotting original and interpolated RR intervals for 100 seconds
%
figure(6);
% Size
set(gcf, 'Position', [80, 80, 1000, 600]);

subplot(2, 1, 1);
plot(new_time, rr_intervals);
xlabel('Time (s)');
ylabel('RR-interval (s)');
title('Original RR Intervals');
xlim([0, 100]);

subplot(2, 1, 2);
plot(samples_128, rr_interpolated_128);
xlabel('Time (s)');
ylabel('RR-interval (s)');
title('Interpolated RR Intervals (@ 128 samples/s');
xlim([0, 100]);

sgtitle('Comparison of Original and Interpolated RR Intervals (100 seconds only)');

% The interpolated-sampled rr-intervals plotted are smooth 
% when compared to original rr intervals


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
%%%%%% Part 8 %%%%%%
%%%%%%%%%%%%%%%%%%%%

% Calculating autocorrelation

% Shifting one of the signals by 100 samples
shifted_rr_intervals = circshift(rr_intervals, 100);

% Autocorelation using cross-corelation function
% Might Need to install Econometrics Toolbox
[Rxx, lags] = xcorr(rr_intervals,shifted_rr_intervals); 

% Plotting autocorrelation

figure(7);
plot(lags, Rxx);  
xlabel('Samples');
ylabel('Amplitude');
title('Autocorrelation of RR Intervals');
grid on;


% %%%%%%%%%%%%%%%%%%%%
% %%%%%% Part 9 %%%%%%
% %%%%%%%%%%%%%%%%%%%%

p = pspectrum(rr_interpolated_128, fs2);

% Plotting the power spectrum
figure(8);
plot(p);  % The pspectrum function already handles the plot
xlabel('Frequency (Hz)');
xlim([0, 50]);
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectrum of RR Intervals');
grid on;



% %%%%%%%%%%%%%%%%%%%%
% %%%%%% Part 10 %%%%%
% %%%%%%%%%%%%%%%%%%%%


% FIR filter parameters
cutoff_frequency = 50; % Hz
nyquist = fs / 2;
normalized_cutoff = cutoff_frequency / nyquist;
filter_order = 100; %  filter order

% Designing the FIR filter using the fir1 function
fir_filter = fir1(filter_order, normalized_cutoff, 'low');

% Applying the FIR filter to the raw ECG signal
filtered_signal = filter(fir_filter, 1, z_amp);

%Plotting raw and filtered ecg

figure(9);
subplot(2,1,1);
plot(time, z_amp);
xlabel('Time (s)');
ylabel('Amplitude');
title('Raw ECG Signal');

subplot(2,1,2);
plot(time, filtered_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered ECG Signal (Frequency < 50 Hz)');

%Frequency Response of the raw and filtered

N = length(z_amp);
frequencies = linspace(0, fs/2, N);

original_spectrum = abs(fft(z_amp)/N);
filtered_spectrum = abs(fft(filtered_signal)/N);

figure(10);
subplot(2,1,1);
plot(original_spectrum);
xlabel('Frequencies');
ylabel('Amplitude');
title('Frequency Spectrum Raw ECG Signal');

subplot(2,1,2);
plot(filtered_spectrum);
xlabel('Frequencies');
ylabel('Amplitude');
title('Frequency Specturm Filtered ECG Signal');

% Plotting the frequency response of the FIR filter
figure(11);
freqz(fir_filter, 1, 1024, fs); % fs is the sampling frequency of your signal
title('Frequency Response of FIR Filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[peaks2, peak_indices2] = findpeaks(filtered_signal, time, 'MinPeakHeight', 0.2, 'MinPeakDistance', 0.4);
heartbeats_filtered = length(peaks2);
fprintf('Number of heartbeats (filtered): %d\n', heartbeats_filtered);
avg_HBs_per_minute2 = heartbeats_filtered / (30 - 15);
fprintf('Average Number of heartbeats per minute (filtered ecg): %d\n', round(avg_HBs_per_minute2));

% Calculating the RR intervals
rr_intervals = diff(peak_indices); % difference for the adjacent QRS peaks = RR interval
% Displaying the RR intervals
fprintf('Number of RR Intervals (RAW ECG):%d\n',length(peak_indices));
rr_intervals_filtered = diff(peak_indices2);
fprintf('Number of RR Intervals (Filtered ECG):%d\n',length(peak_indices2));

%%%%%%%%%%%%%%%%%%%%
%% %%% Part 11 %%% %
%%%%%%%%%%%%%%%%%%%%

% Calculating crosscorrelation
% cross-corelation function
% Might Need to install Econometrics Toolbox

[Rxx2, lags2] = xcorr(filtered_signal,z_amp); 

% Plotting autocorrelation
figure(12);
plot(lags2, Rxx2);  
xlabel('Samples');
ylabel('Amplitude');
title('Crosscorrelation of Raw ECG & FIltered ECG');
grid on;

%%%%%%%%%%%%%%%%%%%%
%% %%% Part 12 %%% %
%%%%%%%%%%%%%%%%%%%%

% % Interpolation
new_time2 = cumsum(rr_intervals_filtered);
rr_filtered_interpolated_uniform = interp1(new_time2, rr_intervals_filtered, 'cubic', 'pp');

% Uniformly sample the time vector
% number of samples = 900 (can be changed as needed)
uniformly_filtered_sampled_time = linspace(new_time2(1), new_time2(end), 900); 

% Evaluating the interpolated RR intervals at uniformly sampled time points
sampled_rr_filtered_evaluated = ppval(rr_filtered_interpolated_uniform, uniformly_filtered_sampled_time);

figure(13);
plot(uniformly_filtered_sampled_time, sampled_rr_filtered_evaluated);
title('Uniformly Interpolated RR (Post Filtering)');
xlabel('Time (s)');
ylabel('RR-interval (s)');

% % Defining sampling rate and steps
fs2 = 128.0;
steps = 1 / fs2;

% Sampling from interpolation function
samples_128 = 1:steps:max(new_time2);
rr_filtered_interpolated_128 = ppval(rr_filtered_interpolated_uniform, samples_128);

% Plotting original and interpolated RR intervals for whole 15 minutes

figure(14);
set(gcf, 'Position', [80, 80, 1000, 600]);
subplot(2, 1, 1);
plot(new_time2, rr_intervals_filtered);
xlabel('Time (s)');
ylabel('RR-interval (s) (without interpolation)');
title('RR Intervals');

subplot(2, 1, 2);
plot(samples_128, rr_filtered_interpolated_128);
xlabel('Time (s)');
ylabel('RR-interval (s)');
title('Interpolated RR Intervals (@ 128 samples/s)');
sgtitle('Comparison of Original and Interpolated RR Intervals (Post Filtering)');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The plot 15 minutes, it is difficult to distinguish the difference
% between the original and interpolated-sampled RR intervals plot
% Therefore, Plotting only for 100 seconds to compare the differences

% Plotting original and interpolated RR intervals for 100 seconds

figure(15);

set(gcf, 'Position', [80, 80, 1000, 600]);
subplot(2, 1, 1);
plot(new_time2, rr_intervals_filtered);
xlabel('Time (s)');
ylabel('RR-interval (s)');
title('RR Intervals (without Interpolation)');
xlim([0, 100]);

subplot(2, 1, 2);
plot(samples_128, rr_filtered_interpolated_128);
xlabel('Time (s)');
ylabel('RR-interval (s)');
title('Interpolated RR Intervals (@ 128 samples/s');
xlim([0, 100]);

sgtitle('Comparison of Original and Interpolated RR Intervals (Post Filtering & 100 seconds only)');

% % The interpolated-sampled rr-intervals plotted are smooth 
% % when compared to original rr intervals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating autocorrelation

% Shifting one of the signals by 100 samples
shifted_rr_intervals_filtered = circshift(rr_intervals_filtered, 100);

% Autocorelation using cross-corelation function
% Might Need to install Econometrics Toolbox
[Rxx3, lags3] = xcorr(rr_intervals,shifted_rr_intervals_filtered); 

% Plotting autocorrelation
figure(16);
plot(lags, Rxx);  
xlabel('Samples');
ylabel('Amplitude');
title('Autocorrelation of RR Intervals (Post Filtering)');
grid on;

% Plotting the power spectrum
p2 = pspectrum(rr_filtered_interpolated_128, fs2);

figure(17);
plot(p2);  % The pspectrum function already handles the plot
xlabel('Frequency (Hz)');
xlim([0, 50]);
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectrum of RR Intervals (Post Filtering)');
grid on;

