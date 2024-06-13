%% Reference code for the Coded Modulation module tutorial %%%%%%%%%%
% Computes mutual information of a given constellation in a Gaussian channel using Monte-Carlo sampling and plots it 
% as a function of the SNR %% 

%% Parameters 
X=-3:2:3;                % Constituent real constellation 
C=repmat(X,length(X),1)+1i*repmat(X.',1,length(X));
C=C(:);
C=C/sqrt(mean(abs(C).^2,1));      % Normalise to unitary energy

scatterplot(C);            % Plot constellation 

D=1e5;                      % Number of Monte-Carlo samples
SNR=-10:.5:25;                % Es/No 
m=log2(length(X)^2);

%% Computation
%Initialise quantities to compute
I=m*ones(1,length(SNR));
No=10.^(-SNR/10);  
M=size(C,1);               % Constellation size
m=log2(M);                 % Bits/symbol


for ss=1:length(SNR)        % SNR loop
    for i=1:M                  % Outer loop across constelletion symbols
         Z=sqrt(No(ss)/2)*(randn(1,D)+1i*randn(1,D));
         arg=0;
        for j=1:M              % Inner loop across constelletion symbols
            dij = C(i,:)-C(j,:);
            arg = exp(-(2*real(Z*dij)+sum(abs(dij)^2))/No(ss))+arg;
        end
        I(ss)= I(ss)-1/(D*M)*sum(log2(arg));
    end

end

figure;
plot(SNR,I,'LineWidth',2); grid on; xlabel('SNR [dB]'); ylabel('MI [bit/2D]');



