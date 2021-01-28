function [output_signal, Fs, bpm] = get_BPMs(music, wn, showPlot)

[data, Fs] = audioread(music);              % read music and get data and sampling frequency
signal = data(:,2);                         % choose the data channel

%% resampling
ds = 6; 
signal = resample(signal, 1, ds);           % resample de 1/ds change from 48 kHz to 8 kHz
Fs = Fs/ds;                                 % new sampling frequency

%% number of points and time vector
N = length(signal);                         % number of points
t = (0 : N-1) / Fs;                         % time vector

%% filter design and signal filtering
if(wn ~= 0)
    order = 2;
    fnyquist = Fs/2;
    wn = wn / fnyquist;                    % normalized cutoff frequency (fnyquist)
    [b,a] = butter(order,wn);              % IIR butterworth filter
    output_signal = filtfilt(b,a,signal);  % applying the filter to original signal
else
    output_signal = signal;                % signal not filtered
end

%% autocorrelation
% compare the signal with it self with differents time shifts
 % this allows us to understand if the signal has a temporal structure
  % if it has we can determine the signal frequency
[r,lags] = xcorr(output_signal, 'coeff');

%% find spikes on the correlation graph
% 60%180 = 0.33 sec represents the minimum period of time to be analised
  % because a low period will result in higher bpm's
% 'MinPeakHeight' represents the minimum value of autocorrelation (4%) 
  % for which we want to count spikes
[spikes,loc] = findpeaks(r, 'MinPeakDistance', Fs*0.33, 'MinPeakHeight', 0.04);

%% obtain BPM
bpm = zeros(length(loc)-1,1);               % inicialize array of zeros
for i=1:length(loc)-1
    period = loc(i+1)/Fs - loc(i)/Fs;       % get the distance between spikes over time = period
    bpm(i) = (1/period) * 60;               % convert time in bpm's = frequency
end

%% find the most relevant interval of bpm 
[value,edges] = histcounts(bpm, 60:20:180);     % return edges (1st value of the interval range) values represents how much bmp are in each interval
[MaxBPM, iMaxBPM] = sort(value, 'descend');     % order descend by number of 'value' indexing the position of MaxBPM on array [value,edges]

% select the most relevant interval of bpms (the biggest value)
bestIntervalBPM = bpm >= edges(iMaxBPM(1)) & bpm <= edges(iMaxBPM(1)+1); %foreach bpm value we compare with the minimum and maximum value
                                                                            %of the interval (edges) and (edges+1) that has the biggest count of 'value' of bpm
                                                                            %if the condition is true the value of the position of 'bestIntervalBPM' is 1 else is 0
                                                                            %this array will be used to calculate the average bpm on the interval with the higher bpm count

imr = (MaxBPM(1)/sum(value))*100;              % get the accuracy of the most relevant interval
imr2 = (MaxBPM(2)/sum(value))*100;             % get the accuracy of the second most relevant interval

%% show results
fprintf('Interval: [%d %d]: %.1f%% Accuracy \n', edges(iMaxBPM(1)), edges(iMaxBPM(1)+1), imr);
fprintf('Interval: [%d %d]: %.1f%% Accuracy \n', edges(iMaxBPM(2)), edges(iMaxBPM(2)+1), imr2);
fprintf('Average BPM: %.1f \n\n', mean(bpm(bestIntervalBPM)));

%% show multiple plots
if(showPlot == 0)
    figure('Name',strcat('BPM Analysis_',music));
    %% filter analysis
    if(wn ~= 0)
        [H,W] = freqz(b, a, Fs, Fs);
        subplot(4,3,1);
        plot(W, abs(H)); 
        title('Filter response graph');
        xlabel('Frequency (Hz)'); 
        ylabel('Gain (dB)');
    end
    %% CALCÚLO E PLOT DE STFT
    nfft = Fs/2;                            % window size
    noverlap = round(nfft/2);               % number of overlaping points in adjacent segements
    
    [PS_original,f_original,t_original] = stft(signal, nfft, noverlap, Fs); 
    [PS_filtered,f_filtered,t_filtered] = stft(output_signal, nfft, noverlap, Fs); 
    
    res_espectral = Fs / nfft;              
    
    y = ceil( 1000 / res_espectral);        % calculation of the frequency range that we must use
    subplot(4,3,2); 
    mesh(f_original(2:y), t_original, PS_original(:, 2:y));     % 3D plot - STFT orignal signal
    title('STFT of original signal');
    xlabel('X - Freq. (Hz)'); 
    ylabel('Y - Tempo (s)'); 
    zlabel('Z - Potência');
    
    subplot(4,3,3); 
    mesh(f_filtered(2:y), t_filtered, PS_filtered(:, 2:y));     % 3D plot - STFT filtered signal
    title('STFT of filtered signal');
    xlabel('X - Freq. (Hz)'); 
    ylabel('Y - Tempo (s)'); 
    zlabel('Z - Potência');
    
    %% spectogram
    subplot(4,3,4);
    spectrogram(signal,1024,Fs/20,Fs/10,Fs,'yaxis');
    axis xy; axis tight; 
    title('spectogram of the original signal');
    xlabel('Time (sec)'); 
    ylabel('Frequency (Hz)');
    
    subplot(4,3,7); 
    spectrogram(output_signal,1024,Fs/20,Fs/10,Fs,'yaxis');
    axis xy; axis tight;
    title('spectogram of filtered signal');
    xlabel('Time (sec)'); 
    ylabel('Frequency (Hz)');
    
    %% plot signal in time (seconds)
    subplot(4,3,5); 
    plot(t, signal); 
    title('Original signal');
    xlabel('Time (sec)'); 
    ylabel('Amplitude');

    subplot(4,3,8); 
    plot(t, output_signal); 
    title('Filtered signal');
    xlabel('Time (sec)'); 
    ylabel('Amplitude');
    
    %% fft - Power Spectrum
    %shows the power
    f1 = (0:N-1) * (Fs/(N-1));                  % frequencies vector
    PS_before = abs(fft(signal)).^2;            % Power spectrum of the original signal
    PS_after = abs(fft(output_signal)).^2;      % Power spectrum of the signal after being filtered
    
    subplot(4,3,6); 
    plot(f1(1:end/2), PS_before(1:end/2)); 
    title('FFT of the original signal');
    xlabel('Frequency (Hz)'); 
    ylabel('abs(PS_before)');
    
    subplot(4,3,9); 
    plot(f1(1:end/2), PS_after(1:end/2)); 
    title('FFT of the filtered signal');
    xlabel('Frequency (Hz)'); 
    ylabel('abs(PS_after)');
    
    %% plot the correlation of the filtered signal
    subplot(4,3,11)
    plot(lags/Fs, abs(r)); 
    title('Correlation with filter');
    xlabel('Time(sec)'); 
    ylabel('Correlation')
    axis([-floor(t(end)) floor(t(end)) 0 1])
    hold on;   
    plot(lags(loc)/Fs, spikes, 'o');        % represent spikes on graph
    
end

%% show only the autocorrelation graph
if(showPlot == 1)
    figure('Name', strcat('Autocorrelation graph of the filtered signal_',music));
    plot(lags/Fs, abs(r)); 
    title('Correlation with filter');
    xlabel('Time(sec)'); 
    ylabel('Correlation')
    axis([-floor(t(end)) floor(t(end)) 0 1])
    hold on;   
    plot(lags(loc)/Fs, spikes, 'o');      % represent spikes on graph  
end