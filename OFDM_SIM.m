% Senjor Project:   OFDM Simulation using MATLAB
% Student:          Paul Lin
% Professor:        Dr. Cheng Sun
% Date:             June, 2010
% *************** MAIN PROGRAM FILE *************** %
% This is the OFDM simulation program's main file.
% It requires a 256-grayscale bitmap file (*.bmp image file) as data source
% and the following 5 script and function m-files to work:
%   ofdm_parameters.m, ofdm_base_convert.m, ofdm_modulate.m,
%   ofdm_frame_detect.m, ofdm_demod.m


% ####################################################### %
% ************* OFDM SYSTEM INITIALIZATION: ************* %
% **** setting up parameters & obtaining source data **** %
% ####################################################### %

% Turn off exact-match warning to allow case-insensitive input files
warning('off','MATLAB:dispatcher:InexactMatch'); 
close all;      % close all previously opened figures and graphs
clc             % clear the screen

fprintf('\n\n##########################################\n')
fprintf('#*********** OFDM Simulation ************#\n')
fprintf('##########################################\n\n')

% invoking ofdm_parameters.m script to set OFDM system parameters
ofdm_parameters;
load parameters.mat;
% save parameters for receiver 
% read data from input file
x = imread(file_in);

% arrange data read from image for OFDM processing 
h=size(x,1);
w=size(x,2);
t=size(x,3);
x = reshape(x, 1,[]);
baseband_tx = double(x);

% convert original data word size (bits/word) to symbol size (bits/symbol)
% symbol size (bits/symbol) is determined by choice of modulation method
baseband_tx = ofdm_base_convert(baseband_tx, word_size, symb_size);

% save original baseband data for error calculation later
save('err_calc.mat', 'baseband_tx');
% ####################################################### %
% ******************* OFDM TRANSMITTER ****************** %
% ####################################################### %

tic;    % start stopwatch

% generate header and trailer (an exact copy of the header)
f = 0.25;
header = sin(0:f*2*pi:f*2*pi*(head_len-1));
f=f/(pi*2/3);
header = header+sin(0:f*2*pi:f*2*pi*(head_len-1));

% arrange data into frames and transmit
frame_guard = zeros(1, symb_period);
time_wave_tx = [];
symb_per_carrier = ceil(length(baseband_tx)/carrier_count);
fig = 1;
if (symb_per_carrier > symb_per_frame)  % === multiple frames === %
    power = 0;
    while ~isempty(baseband_tx)
        % number of symbols per frame
        frame_len = min(symb_per_frame*carrier_count,length(baseband_tx));
        frame_data = baseband_tx(1:frame_len);
        % update the yet-to-modulate data
        baseband_tx(1:frame_len) = [];%baseband_tx((frame_len+1):(length(baseband_tx)));
        % OFDM modulation
        time_signal_tx = ofdm_modulate(frame_data,ifft_size,carriers,conj_carriers, carrier_count, symb_size, guard_time, fig);
        fig = 0; %indicate that ofdm_modulate() has already generated plots        
        % add a frame guard to each frame of modulated signal
        time_wave_tx = [time_wave_tx frame_guard time_signal_tx];
        frame_power = var(time_signal_tx);
    end
    % scale the header to match signal level
    power = power + frame_power;
    % The OFDM modulated signal for transmission
    time_wave_tx = [power*header time_wave_tx frame_guard power*header];
else                                    % === single frame === %
    % OFDM modulation
    time_signal_tx = ofdm_modulate(baseband_tx,ifft_size,carriers,conj_carriers, carrier_count, symb_size, guard_time, fig);
    
    % calculate the signal power to scale the header
    power = var(time_signal_tx);
    % The OFDM modulated signal for transmission
    time_wave_tx = [power*header frame_guard time_signal_tx frame_guard power*header];
end

% show summary of the OFDM transmission modeling
peak = max(abs(time_wave_tx(head_len+1:length(time_wave_tx)-head_len)));
sig_rms = std(time_wave_tx(head_len+1:length(time_wave_tx)-head_len));
peak_rms_ratio = (20*log10(peak/sig_rms));
fprintf('\nSummary of the OFDM transmission and channel modeling:\n')
fprintf('Peak to RMS power ratio at entrance of channel is: %f dB\n', peak_rms_ratio)
% ####################################################### %
% **************** COMMUNICATION CHANNEL **************** %
% ####################################################### %

% ===== signal clipping ===== %
clipped_peak = (10^(0-(clipping/20)))*max(abs(time_wave_tx));
time_wave_tx(abs(time_wave_tx)>=clipped_peak)...
    = clipped_peak.*time_wave_tx(abs(time_wave_tx)>=clipped_peak)...
    ./abs(time_wave_tx(abs(time_wave_tx)>=clipped_peak));

% ===== channel noise ===== %
power = var(time_wave_tx);  % Gaussian (AWGN)
SNR_linear = 10^(SNR_dB/10);
noise_factor = sqrt(power/SNR_linear);
noise = randn(1,length(time_wave_tx)) * noise_factor;
time_wave_rx = time_wave_tx + noise;

% show summary of the OFDM channel modeling
peak = max(abs(time_wave_rx(head_len+1:length(time_wave_rx)-head_len)));
sig_rms = std(time_wave_rx(head_len+1:length(time_wave_rx)-head_len));
peak_rms_ratio = (20*log10(peak/sig_rms));
fprintf('Peak to RMS power ratio at exit of channel is: %f dB\n', ...
    peak_rms_ratio)

% Save the signal to be received
save('received.mat', 'time_wave_rx', 'h', 'w','t');
fprintf('#******** OFDM data transmitted in %f seconds ********#\n\n', toc)


% ####################################################### %
% ********************* OFDM RECEIVER ******************* %
% ####################################################### %

disp('Press any key to let OFDM RECEIVER proceed...')
pause;
clear all;      % flush all data stored in memory previously
tic;    % start stopwatch

% invoking ofdm_parameters.m script to set OFDM system parameters
load parameters;

% receive data
load received.mat;
time_wave_rx = time_wave_rx.';

end_x = length(time_wave_rx);
start_x = 1;
data = [];
phase = [];
last_frame = 0;
unpad = 0;
if rem(w*h, carrier_count)~=0
    unpad = carrier_count - rem(w*h, carrier_count);
end
num_frame=ceil((h*w)*(word_size/symb_size)/(symb_per_frame*carrier_count));
fig = 0;
for k = 1:num_frame
    if k==1 || k==num_frame || rem(k,max(floor(num_frame/10),1))==0
        fprintf('Demodulating Frame #%d\n',k)
    end
    % pick appropriate trunks of time signal to detect data frame
    if k==1
        time_wave = time_wave_rx(start_x:min(end_x, ...
            (head_len+symb_period*((symb_per_frame+1)/2+1))));
    else
        time_wave = time_wave_rx(start_x:min(end_x, ...
            ((start_x-1) + (symb_period*((symb_per_frame+1)/2+1)))));
    end
    % detect the data frame that only contains the useful information
    frame_start = ...
        ofdm_frame_detect(time_wave, symb_period, envelope, start_x);
    if k==num_frame
        last_frame = 1;
        frame_end = min(end_x, (frame_start-1) + symb_period*...
            (1+ceil(rem(w*h,carrier_count*symb_per_frame)/carrier_count)));
    else
        frame_end=min(frame_start-1+(symb_per_frame+1)*symb_period, end_x);
    end
    % take the time signal abstracted from this frame to demodulate
    time_wave = time_wave_rx(frame_start:frame_end);
    % update the label for leftover signal
    start_x = frame_end - symb_period;
    
    if k==ceil(num_frame/2)
        fig = 1;
    end
    % demodulate the received time signal
    [data_rx, phase_rx] = ofdm_demod(time_wave, ifft_size, carriers, conj_carriers,guard_time, symb_size, word_size, last_frame, unpad, fig);
    if fig==1
        fig = 0;   % indicate that ofdm_demod() has already generated plots
    end
    
    phase = [phase phase_rx];
    data = [data data_rx];
end
phase_rx = phase;  % decoded phase
data_rx = data;    % received data

% convert symbol size (bits/symbol) to file word size (bits/byte) as needed
data_out = ofdm_base_convert(data_rx, symb_size, word_size);

fprintf('#********** OFDM data received in %f seconds *********#\n\n', toc)


% ####################################################### %
% ********************** DATA OUTPUT ******************** %
% ####################################################### %

% patch or trim the data to fit a w-by-h image
if length(data_out)>(w*h)        % trim extra data
    data_out = data_out(1:(w*h));
elseif length(data_out)<(w*h)    % patch a partially missing row
    buff_h = h;
    h = ceil(length(data_out)/w);
    % if one or more rows of pixels are missing, show a message to indicate
    if h~=buff_h
        disp('WARNING: Output image smaller than original')
        disp('         due to data loss in transmission.')
    end
    % to make the patch nearly seamless,
    % make each patched pixel the same color as the one right above it
    if length(data_out)~=(w*h)
        for k=1:(w*h-length(data_out))
            mend(k)=data_out(length(data_out)-w+k);
        end
        data_out = [data_out mend];
    end
end

% format the demodulated data to reconstruct a bitmap image
data_out = reshape(data_out, w, h)';
data_out = uint8(data_out);
% save the output image to a bitmap (*.bmp) file
imwrite(data_out, file_out, 'bmp');


% ####################################################### %
% ****************** ERROR CALCULATIONS ***************** %
% ####################################################### %

% collect original data before modulation for error calculations
load('err_calc.mat');

fprintf('\n#**************** Summary of Errors ****************#\n')

% Let received and original data match size and calculate data loss rate
if length(data_rx)>length(baseband_tx)
    data_rx = data_rx(1:length(baseband_tx));
    phase_rx = phase_rx(1:length(baseband_tx));
elseif length(data_rx)<length(baseband_tx)
    fprintf('Data loss in this communication = %f%% (%d out of %d)\n', ...
        (length(baseband_tx)-length(data_rx))/length(baseband_tx)*100, ...
        length(baseband_tx)-length(data_rx), length(baseband_tx))
end

% find errors
errors = find(baseband_tx(1:length(data_rx))~=data_rx);
fprintf('Total number of errors = %d (out of %d)\n', ...
    length(errors), length(data_rx))
% Bit Error Rate
fprintf('Bit Error Rate (BER) = %f%%\n',length(errors)/length(data_rx)*100)

% find phase error in degrees and translate to -180 to +180 interval
phase_tx = baseband_tx*360/(2^symb_size);
phase_err = (phase_rx - phase_tx(1:length(phase_rx)));
phase_err(phase_err>=180) = phase_err(phase_err>=180)-360;
phase_err(phase_err<=-180) = phase_err(phase_err<=-180)+360;
fprintf('Average Phase Error = %f (degree)\n', mean(abs(phase_err)))

% Error pixels
x = ofdm_base_convert(baseband_tx, symb_size, word_size);
x = uint8(x);
x = x(1:(size(data_out,1)*size(data_out,2)));
y = reshape(data_out', 1, length(x));
err_pix = find(y~=x);
fprintf('Percent error of pixels of the received image = %f%%\n\n', ...
    length(err_pix)/length(x)*100)

fprintf('##########################################\n')
fprintf('#******** END of OFDM Simulation ********#\n')
fprintf('##########################################\n\n')
[img_rece,map] = imread('qq_OFDM.bmp');
[img_sent,map0] = imread('qq.bmp');
figure(8)
subplot(121)
imshow(img_sent,map0)
subplot(122)
imshow(img_rece,map)
