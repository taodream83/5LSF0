%% Script for coding assignment in Comms Theory %%
% Author: Gabriele Liga Dec 2023
clear all
close all

%% Code design 
% Extended Hamming code (8,4)
n=7;  k0=4;                                                                
H=de2bi(1:2^(n-k0)-1,n-k0).';                                              % Parity-check matrix of Hamming (7,4) 
Gh=[1,0,0,0,0,1,1; 0,1,0,0,1,0,1; 0,0,1,0,1,1,0; 0,0,0,1,1,1,1];           % Generator matrix for Hamming (7,4) in a systematic form
Heh=[H,zeros(3,1); ones(1,n+1)];                                           % Parity-check matrix of the extended Hamming (8,4)
Geh=[Gh,mod(Gh*ones(n,1),2)];                                              % Generator matrix of the extended Hamming (8,4) 
n=8;
r0=k0/n;

% (8,2) code with dmin=5 (t=2)  
G1=[1,0,1,1,1,1,0,0;0,1,1,1,0,0,1,1];                                      % Generator matrix of (8,2) code with dmin=5 (t=2)
k1=size(G1,1);
r1=k1/n;

% (8,1) code with dmin=8 (t=3)
G2=ones(1,n);                                                              %% Generator matrix off (8,1) code with dmin=8 (t=3)
k2=size(G2,1);
r2=k2/n;

% Generate codebooks  
C0=mod(de2bi(0:2^k0-1)*Geh,2);            % Hamming codebook
C1=mod(de2bi(0:2^k1-1)*G1,2);             % (8,2) codebook
C2=mod(de2bi(0:2^k2-1)*G2,2);             % (8,1) codebook
X0=1-2*C0;
X1=1-2*C1;
X2=1-2*C2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% BI-AWGN simulation  
Nfr=1e6;               % Number of transmitted codeword (frames)
SNR=-3:.5:2;
BER=zeros(1,length(SNR));
BER0=BER;
BER1=BER;
BER2=BER;

for ss=1:length(SNR)    
info_bits=randi([0 1],Nfr,n);    
info_bits0=randi([0 1],Nfr,k0);
info_bits1=randi([0 1],Nfr,k1);
info_bits2=randi([0 1],Nfr,k2);


%% Encoding 
c0=mod(info_bits0*Geh,2);
c1=mod(info_bits1*G1,2);
c2=mod(info_bits2*G2,2);

%% Modulation 
x=1-2*info_bits;  % Uncoded TX
x0=1-2*c0;
x1=1-2*c1;
x2=1-2*c2;

%% AWGN channel 
noise=sqrt(10^(-SNR(ss)/10))*randn(Nfr,n);
y=x+noise;
y0=x0+noise;                    
y1=x1+noise;                    
y2=x2+noise;                   

%% Hard-decision 
Ihat=(y<0);

%% ML decoding 
[~,mhat0]=max(squeeze(sum(repmat(y0,1,1,2^k0).*permute(repmat(X0,1,1,Nfr),[3,2,1]),2)),[],2);
[~,mhat1]=max(squeeze(sum(repmat(y1,1,1,2^k1).*permute(repmat(X1,1,1,Nfr),[3,2,1]),2)),[],2);
[~,mhat2]=max(squeeze(sum(repmat(y2,1,1,2^k2).*permute(repmat(X2,1,1,Nfr),[3,2,1]),2)),[],2);

Ihat0=C0(mhat0,1:k0);  
Ihat1=C1(mhat1,1:k1); 
Ihat2=C2(mhat2,k2); 

%% MC error counting 
BitErr=sum(info_bits~=Ihat,'all');
BER(ss)=BitErr/(Nfr*n);
BitErr=sum(info_bits0~=Ihat0,'all');
BER0(ss)=BitErr/(Nfr*k0);
BitErr=sum(info_bits1~=Ihat1,'all');
BER1(ss)=BitErr/(Nfr*k1);
BitErr=sum(info_bits2~=Ihat2,'all');
BER2(ss)=BitErr/(Nfr*k2);
end

figure(1) 
p=semilogy(SNR,BER,'--ks','MarkerFaceColor','w'); grid on; hold on;
p0=semilogy(SNR,BER0,'-o','MarkerFaceColor','w');
p1=semilogy(SNR,BER1,'-o','MarkerFaceColor','w'); 
p2=semilogy(SNR,BER2,'-o','MarkerFaceColor','w'); 
xlabel('SNR [dB]'); ylabel('BER'); 
legend([p,p0,p1,p2],'Uncoded','Hamming (8,2)','(8,2), d_{min}=5','(8,1), d_{min}=8');

figure(2) 
EbNo=SNR; EbNo0=SNR-10*log10(r0) ; EbNo1=SNR-10*log10(r1); EbNo2=SNR-10*log10(r2);
semilogy(EbNo,BER,'-o','MarkerFaceColor','w'); grid on; hold on;
semilogy(EbNo0,BER0,'-o','MarkerFaceColor','w'); 
semilogy(EbNo1,BER1,'-o','MarkerFaceColor','w');  
semilogy(EbNo2,BER2,'-o','MarkerFaceColor','w'); 
xlabel('Eb/No [dB]'); ylabel('BER'); 
legend([p,p0,p1,p2],'Uncoded','Hamming (8,2)','(8,2), d_{min}=5','(8,1), d_{min}=8');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%