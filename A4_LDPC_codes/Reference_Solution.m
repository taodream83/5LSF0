%% Defining te parity check matrix
% This parity check matrix is the rate 3/4 DVB-S2 code with blocklength 648
load('H.mat');

%% Answers to the question
%Degree distribution is 0.1L^8 + 0.1L^12 + 0.8L^1  and 0.5P^3 + 0.5P^4, the
%rate of the code is 1/5
% The girth of the code is 4, as a 4 cycle can be found in the parity check
% matrix
% The minimum distance of the code is 9, as the codeword with the lowest
% weight that is a non-zero codeword has weight 9


%% Doing Simulation 1
ncodewords = 4000;
SNR = 1:0.2:2.6;
blocklength = 648;
R = 0.75;
infolength = R*blocklength;
maxiter = 20;
alpha = 1;
encodercfg = ldpcEncoderConfig(sparse(logical(H)));
decodercfg1 = ldpcDecoderConfig(encodercfg,"norm-min-sum");
decodercfg2 = ldpcDecoderConfig(encodercfg);
BER = zeros(length(SNR),1);
BER2 = zeros(length(SNR),1);
for i  =1:length(SNR)
    for n = 1:ncodewords
        s = randi(2,infolength,1)-1;
        c = ldpcEncode(s,encodercfg);
        modSignal = pskmod(c,2,InputType='bit');
        [rxsig, noisevar] = awgn(modSignal,SNR(i));
        l = pskdemod(rxsig,2, ...
            OutputType='approxllr', ...
            NoiseVariance=noisevar);
        [s_hat,iterations] = ldpcDecode(l,decodercfg1,maxiter,"MinSumScalingFactor",alpha);
        [s_hat2,iterations2] = ldpcDecode(l,decodercfg2,maxiter);
        BER(i) = BER(i) + nnz(s_hat ~= s);
        BER2(i) = BER2(i) + nnz(s_hat2 ~= s);
    end
    BER(i) = BER(i)/(infolength*ncodewords);
    BER2(i) = BER2(i)/(infolength*ncodewords);
end

%% Doing Simulation 2

ncodewords = 4000;
SNR = 1:0.2:2.6;
blocklength = 648;
R = 0.75;
infolength = R*blocklength;
maxiter = 20;
alpha = 1:-0.2:0.2;
encodercfg = ldpcEncoderConfig(sparse(logical(H)));
decodercfg1 = ldpcDecoderConfig(encodercfg,"norm-min-sum");
decodercfg2 = ldpcDecoderConfig(encodercfg);
BER = zeros(length(SNR),length(alpha));
for i  =1:length(SNR)
    for n = 1:ncodewords
        s = randi(2,infolength,1)-1;
        c = ldpcEncode(s,encodercfg);
        modSignal = pskmod(c,2,InputType='bit');
        [rxsig, noisevar] = awgn(modSignal,SNR(i));
        l = pskdemod(rxsig,2, ...
            OutputType='approxllr', ...
            NoiseVariance=noisevar);
        for j = 1:length(alpha)
            [s_hat,iterations] = ldpcDecode(l,decodercfg1,maxiter,"MinSumScalingFactor",alpha(j));
            BER(i,j) = BER(i,j) + nnz(s_hat ~= s);
        end
    end
    BER(i,:) = BER(i,:)/(infolength*ncodewords);
end
