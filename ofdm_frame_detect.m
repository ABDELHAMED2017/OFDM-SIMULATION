% Senjor Project:   OFDM Simulation using MATLAB 
% Student:          Paul Lin 
% Professor:        Dr. Cheng Sun 
% Date:             June, 2010 
% ************* FUNCTION: ofdm_frame_detect() ************* % 
% This function is to synchronize the received signal before demodulation 
% by detecting the starting point of a frame of received signal. 
  
function start_symb = ofdm_frame_detect(signal, symb_period, env, label) 
%   Find the approximate starting location 
  
signal = abs(signal); 
  
% ===== narrow down to an approximate start of the frame ===== % 
idx = 1:env:length(signal); 
samp_signal = signal(idx);  % sampled version of signal 
mov_sum = filter(ones(1,round(symb_period/env)),1,samp_signal); 
mov_sum = mov_sum(round(symb_period/env):length(mov_sum)); 
apprx = min(find(mov_sum==min(mov_sum))*env+symb_period); 
% move back by approximately 110% of the symbol period to start searching 
idx_start = round(apprx-1.1*symb_period); 
  
% ===== look into the narrow-downed window ===== % 
mov_sum = filter(ones(1,symb_period),1,signal(idx_start:round(apprx+symb_period/3))); 
mov_sum = mov_sum(symb_period:length(mov_sum)); 
null_sig = find(mov_sum==min(mov_sum)); 
  
start_symb = min(idx_start + null_sig + symb_period) - 1; 
% convert to global index 
start_symb = start_symb + (label - 1); 
 