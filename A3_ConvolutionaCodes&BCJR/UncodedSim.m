%% Simulation of convolutional codes in BI-AWGN %%
clear all;
close all;
UncodedParameters; 
BER=ones(1,length(EbNo));

for nn=1:length(EbNo)
bit_err=0;
Nfr=0;
    while bit_err<=TargetErr
        %% TX %%
        info_bits=randi([0 1],1,Nbits); % Generate random input information bits %%

        %% Channel %%
        tx_sym=-2*info_bits+1;                             % Bit-to-symbol mapping  (binary antipodal modulation)
        SNR=EbNo(nn)+3;                             % Channel SNR
        rx_sym=awgn(tx_sym,SNR);                     % BI-AWGN channel

        %% RX %%%%%%
        snr=10^(SNR/10);
        llr=2*rx_sym*snr;                                   % llr computation
        hard_rxbits=(llr<0);                                %  Hard rx bits

 
        %% Bit error count
        bit_err=bit_err+sum(hard_rxbits~=info_bits);
        %berrframe=sum(dec_bits~=info_bits);                          % Number of bit errors per frame
        %bit_err=bit_err+berrframe;                                   % Total accumulated bit errors
        TotBits=Nfr*Nbits;                                           % Total number of bit errors
        BER(nn)=bit_err/TotBits;                                         % Current bit error rate
        Nfr=Nfr+1;
        
        %disp(['Current bit errors=' num2str(bit_err)]);
        %disp(['Current frame errors=' num2str(frame_err)]);
       
    
    end
   disp(['Uncoded BER=' num2str(BER(nn)) ' @EbNo=' num2str(EbNo(nn)) 'dB']);
  
end


semilogy(EbNo,BER,'sb--','MarkerFaceColor','w'); grid on;
xlabel('Eb/No [dB]');
ylabel('BER');


%% Save results %% 
Save=1;
if Save 
% % spath='/data/icave-nas/gabriele/ResearchData/TC for CV-QKD/OFC2024/';
% % if ~exist(spath,'dir')
% % mkdir(spath);
% % end

filename='BI-AWGN_uncoded_BERvsEbNo.mat';
filepath=filename;
save(filepath,'EbNo','BER'); 
end 
