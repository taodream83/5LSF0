%% Parameters for convolutional codes simulation in BI-AWGN%%


%%% Convolutional code parameters %%%%%%%%%%%
Nbits=1e5;     % Number of info bits per simulated block 
K=4;           % Convolutional encoder constraint length 
poly=[13,15];   % [5,7];   % Set of polynomials defining convolutional encoder
trellis=poly2trellis(K,poly);                        % Trellis structure
r=size(poly,1)/size(poly,2);    % Convolutional code rate 
tb=1e5;                         % Decoding window
DecType='unquant';
DecMode='trunc';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EbNo=2:0.2:7;    % Eb/No [dB] 


%% Error estimation parameters 
FERcheck=1;                  % Simulation controlled on frame error numbers 
BERcheck=0;                  % Simulation controlled on bit error numbers 
TargetErr=5e2;               % Target error numbers 