% Senjor Project:   OFDM Simulation using MATLAB 
% Student:          Paul Lin 
% Professor:        Dr. Cheng Sun 
% Date:             June, 2010 
% ************* FUNCTION: ofdm_demod() ************* % 
% This function performs OFDM demodulation after data reception. 
  
function [decoded_symb, decoded_phase] = ofdm_demod... 
    (symb_rx, ifft_size, carriers, conj_carriers, ... 
    guard_time, symb_size, ~, last, unpad, fig) 
  
symb_period = ifft_size + guard_time; 
  
% reshape the linear time waveform into fft segments 
symb_rx_matrix = reshape(symb_rx(1:... 
    (symb_period*floor(length(symb_rx)/symb_period))), ... 
    symb_period, floor(length(symb_rx)/symb_period)); 
  
% ------------------------------------------ % 
% ##### remove the periodic time guard ##### % 
% ------------------------------------------ % 
symb_rx_matrix = symb_rx_matrix(guard_time+1:symb_period,:); 
  
% ------------------------------------------------------------------ % 
% ### take FFT of the received time wave to obtain data spectrum ### % 
% ------------------------------------------------------------------ % 
rx_spectrum_matrix = fft(symb_rx_matrix)'; 
  
% plot magnitude and phase of the received frequency spectrum 
if fig==1 
    limt = 1.1*max(abs(reshape(rx_spectrum_matrix',1,... 
        size(rx_spectrum_matrix,1)*size(rx_spectrum_matrix,2)))); 
    figure(5) 
    stem(0:ifft_size-1, abs(rx_spectrum_matrix(ceil... 
        (size(rx_spectrum_matrix,1)/2),1:ifft_size)),'b*-') 
    grid on 
    axis ([0 ifft_size -limt limt]) 
    ylabel('Magnitude') 
    xlabel('FFT Bin') 
    title('Magnitude of Received OFDM Spectrum') 
    figure(6) 
    plot(0:ifft_size-1, (180/pi)*angle(rx_spectrum_matrix(ceil... 
        (size(rx_spectrum_matrix,1)/2),1:ifft_size)'), 'go') 
    hold on 
    stem(carriers-1, (180/pi)*angle(rx_spectrum_matrix(2,carriers)'),'b*-') 
    stem(conj_carriers-1, (180/pi)*angle(rx_spectrum_matrix(ceil... 
        (size(rx_spectrum_matrix,1)/2),conj_carriers)),'b*-') 
     axis ([0 ifft_size -200 +200]) 
    grid on 
    ylabel('Phase (degrees)') 
    xlabel('FFT Bin') 
    title('Phase of Receive OFDM Spectrum') 
end 
  
% ----------------------------------------------------------------- % 
% ### extract columns of data on IFFT bins of all carriers only ### % 
% ----------------------------------------------------------------- % 
rx_spectrum_matrix = rx_spectrum_matrix(:,carriers); 
  
% --------------------------------------------- % 
% ### PSK (Phase Shift Keying) demodulation ### % 
% --------------------------------------------- % 
% calculate the corresponding phases from the complex spectrum 
rx_phase = angle(rx_spectrum_matrix)*(180/pi); 
% make negative phases positive 
rx_phase = rem((rx_phase+360), 360); 
  
% polar plot for the received symbols 
if fig==1 
    figure(7) 
    rx_mag = abs(rx_spectrum_matrix(ceil(size(rx_spectrum_matrix,1)/2),:)); 
    polar(rx_phase(ceil(size(rx_spectrum_matrix,1)/2),:)*(pi/180), ... 
        rx_mag, 'bd') 
    title('Received Phases') 
end 
  
% --------------------------------- % 
% ##### Differential Decoding ##### % 
% --------------------------------- % 
% reverse the differential coding 
decoded_phase = diff(rx_phase); 
% make negative phases positive 
decoded_phase = rem((decoded_phase+360), 360); 
  
% parellel to serial conversion of phases 
decoded_phase = reshape(decoded_phase', ... 
    1, size(decoded_phase,1)*size(decoded_phase,2)); 
  
% phase-to-data classification 
base_phase = 360/(2^symb_size); 
% phase-to-data translation 
decoded_symb = ... 
    floor(rem((decoded_phase/base_phase+0.5),(2^symb_size))); 
  
% obtain decoded phases for error calculations 
decoded_phase = rem(decoded_phase/base_phase+0.5, ... 
    (2^symb_size))*base_phase - 0.5*base_phase; 
  
% remove padded zeros during modulation 
if last==1 
    decoded_symb  = decoded_symb(1:(length(decoded_symb)-unpad)); 
    decoded_phase = decoded_phase(1:(length(decoded_phase)-unpad)); 
end 
 