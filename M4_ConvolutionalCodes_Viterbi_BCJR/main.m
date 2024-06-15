%% Main file to plot uncoded vs coded performance %%%%%%%%%%%%%%%%%%%%%%%%

ConvSim; 

% Plots results 
W=load("BI-AWGN_uncoded_BERvsEbNo.mat");  % loads uncoded results 
uncBER=W.BER;
uncEbNo=W.EbNo;
s1=semilogy(uncEbNo,uncBER,'or--','MarkerFaceColor','w'); hold on;            % uncoded performance 
s2=semilogy(EbNo,BER,'sb--','MarkerFaceColor','w'); grid on;                  % coded performance 
legend([s1 s2],'Uncoded BI-AWGN','Convolutional Code');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%