%% Parameters %%

%% Parameters for convolutional codes simulation in BI-AWGN%%
EbNo=3:.2:11;    % Eb/No [dB] 
Nbits=1e4;     % Number of info bits per simulated block 

%% Error estimation parameters 
TargetErr=1e3;                  % Target number of (bit) errors 