% Senjor Project:   OFDM Simulation using MATLAB
% Student:          Paul Lin
% Professor:        Dr. Cheng Sun
% Date:             June, 2010
% ************* PARAMETERS INITIALIZATION ************* %
% This file configures parameters for the OFDM system.

% input/output file names
function  ofdm_parameters
file_in = 'image.bmp';
while isempty(file_in)
    file_in = input('source data filename: ', 's');
    if exist([pwd '/' file_in],'file')~=2
        fprintf ('"%s" does not exist in current directory.\n', file_in);
        file_in = [];
    end
end
file_out = [file_in(1:length(file_in)-4) '_OFDM.bmp'];
disp(['Output file will be: ' file_out])

% size of Inverse Fast Fourier Transform (must be power of 2)
ifft_size = 2048;            % force into the while loop below
while (isempty(ifft_size) || (rem(log2(ifft_size),1) ~= 0 || ifft_size < 8))
    ifft_size = input('IFFT size: ');
    if (isempty(ifft_size) || (rem(log2(ifft_size),1) ~= 0 || ifft_size < 8))
        disp('IFFT size must be at least 8 and power of 2.')
    end
end

% number of carriers
carrier_count = 1000;  % force into the while loop below
while (isempty(carrier_count) || (carrier_count>(ifft_size/2-2)) || carrier_count<2)
    carrier_count = input('Number of carriers: ');
    if (isempty(carrier_count) || (carrier_count > (ifft_size/2-2)))
        disp('Must NOT be greater than ("IFFT size"/2-2)')
    end
end

% bits per symbol (1 = BPSK, 2=QPSK, 4=16PSK, 8=256PSK)
symb_size = 4;              % force into the while loop below
while (isempty(symb_size) || (symb_size~=1 && symb_size~=2 && symb_size~=4 && symb_size~=8))
    symb_size = input('Modulation(1=BPSK, 2=QPSK, 4=16PSK, 8=256PSK): ');
    if (isempty(symb_size) || (symb_size~=1&&symb_size~=2&&symb_size~=4&&symb_size~=8))
        disp('Only 1, 2, 4, or 8 can be choosen')
    end
end

% channel clipping in dB
clipping = 2;
while isempty(clipping)
    clipping = input('Amplitude clipping introduced by communication channel (in dB):');
end

% signal to noise ratio in dB
SNR_dB = 10;
while isempty(SNR_dB)aa
    SNR_dB = input('Signal-to-Noise Ratio (SNR) in dB: ');
end

word_size = 8;          % bits per word of source data (byte)

guard_time = ifft_size/4; % length of guard interval for each symbol period
% 25% of ifft_size
% number of symbols per carrier in each frame for transmission
symb_per_frame = ceil(2^13/carrier_count);

% === Derived Parameters === %
% frame_len: length of one symbol period including guard time
symb_period = ifft_size + guard_time;
% head_len: length of the header and trailer of the transmitted data
head_len = symb_period*8;
% envelope: symb_period/envelope is the size of envelope detector
envelope = ceil(symb_period/256)+1;

% === carriers assigned to IFFT bins === %
% spacing for carriers distributed in IFFT bins
spacing = 0;
while (carrier_count*spacing) <= (ifft_size/2 - 2)
    spacing = spacing + 1;
end
spacing = spacing - 1;

% spead out carriers into IFFT bins accordingly
midFreq = ifft_size/4;
first_carrier = midFreq - round((carrier_count-1)*spacing/2);
last_carrier = midFreq + floor((carrier_count-1)*spacing/2);
carriers = (first_carrier:spacing:last_carrier) + 1;
conj_carriers = ifft_size - carriers + 2;
save parameters