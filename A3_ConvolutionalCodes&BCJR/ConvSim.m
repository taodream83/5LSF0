%% Simulation of convolutional codes in BI-AWGN %%
clear all;
close all;
Parameters_Conv; 
FER=ones(1,length(EbNo));
BER=ones(1,length(EbNo));

%%% Coded simulation 
for nn=1:length(EbNo)
Nfr=0;
frame_err=0;
bit_err=0;
NumErr=0;

    while NumErr<TargetErr
        %% TX %%
        info_bits=randi([0 1],1,Nbits); % Generate random input information bits %%

        %%% Convolutional encoding
        coded_bits=convenc(info_bits,trellis);               % Convolutional encoding
        Nc_bits=length(coded_bits);

        %% Channel %%
        tx_sym=-2*coded_bits+1;                             % Bit-to-symbol mapping  (binary antipodal modulation)
        SNR=EbNo(nn)+10*log10(r)+3;                             % Channel SNR
        rx_sym=awgn(tx_sym,SNR);                     % BI-AWGN channel

        %%% RX %%%
        snr=10^(SNR/10);
        llr=2*rx_sym*snr;                                   % llr computation
        hard_rxbits=(llr<0);                                %  Hard rx bits
        

        %% Convolutional decoding
        if strcmp(DecType,'hard')
            dec_bits=vitdec(hard_rxbits,trellis,tb,DecMode,DecType);
        else  
            dec_bits=vitdec(llr,trellis,tb,DecMode,DecType);
        end
        
        ErrFlag=~isequal(dec_bits,info_bits);
               
        if ErrFlag
            frame_err=frame_err+1;
            %% Update bit error count
        ErrVec=(dec_bits~=info_bits);
        %stem(ErrVec);
        berrframe=sum(ErrVec);                          % Number of bit errors per frame
        bit_err=bit_err+berrframe;                                   % Total accumulated bit errors
        end

        %% Frame error count
        Nfr=Nfr+1;
        FER(nn)=frame_err/Nfr;                                           % Current frame error rate
        
        %% Bit error count
        TotBits=Nfr*Nbits;                                               % Total number of transmitted information bits
        BER(nn)=bit_err/TotBits;                                         % Current bit error rate
        
        if ErrFlag
        fprintf(['Bit errors=' num2str(bit_err) '\n']); 
        fprintf(['Frame errors=' num2str(frame_err) '\n']);
        end

        if FERcheck
            NumErr=frame_err;
        elseif BERcheck
            NumErr=bit_err;
        end

    end
   
   disp(['Post-FEC BER=' num2str(BER(nn)) ' @EbNo=' num2str(EbNo(nn)) 'dB']);
   disp(['Post-FEC FER=' num2str(FER(nn)) ' @EbNo=' num2str(EbNo(nn)) 'dB']);

end

