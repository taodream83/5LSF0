%% Script generating solution for assignment A1 - Linear block codes %%%

clear all, close all
ParMode=1;
MatlabWorkers=2;
Unc=0;
G= [1 1 1 0 1 0 0 0; 1 0 0 1 1 1 0 0 ; ...
    1 1 0 0 0 1 1 0; 0 1 1 0 0 0 1 1];       % Generator matrix (given in the assignment)
n=size(G,2);
k=size(G,1);
r=k/n;

%% Uncoded and coded transmission simulation
err_type='BER';
max_err_unc=1e3;
max_err=3e3;

EsNoUnc=1:.2:9;            % Es/No [dB]
EsNo=1:.2:5;              % Es/No [dB]
BERunc=zeros(1,length(EsNo));
FERunc=zeros(1,length(EsNo));
BER=zeros(1,length(EsNo));
FER=zeros(1,length(EsNo));

if Unc
    % Uncoded simulation
    for ee=1:length(EsNoUnc)
        err=0;
        b_err=0;
        esno=10^(EsNoUnc(ee)/10);
        tx_frames=0;
        tx_bits=0;

        while err<=max_err_unc

            u=randi([0  1], 1, n);      % Sequence of information bits
            x=1-2*u;                    % BPSK uncoded transmission
            sigma=sqrt(1/esno/2);             % Noise standard deviation (No/2)
            N=sigma*randn(1,n);
            y=x+N;                           % BI-AWGN channel
            % Hard decision for uncoded data
            u_hat=ones(1,n);
            u_hat(y>0)=0;
            berr_temp=sum(u~=u_hat);
            if berr_temp>0
                err=err+1;
                b_err=b_err+berr_temp;
            end

            tx_frames=tx_frames+1;
            tx_bits=tx_frames*n;
        end

        BERunc(ee)=b_err/tx_bits;
        FERunc(ee)=err/tx_frames;
        disp(['EsNo=' num2str(EsNoUnc(ee))  'dB_uncoded BER=' num2str(BERunc(ee))]);
    end

    semilogy(EsNoUnc,BERunc,'--or','MarkerFaceColor','w');  hold on; grid on;
end


%% Coded transmission
if ParMode
    Parpool=gcp('nocreate');
    if (~isempty(Parpool)  && Parpool.NumWorkers==MatlabWorkers)
    elseif ~isempty(gcp('nocreate'))
        delete(gcp);
        MatlabPool=parpool(MatlabWorkers);
    else
        MatlabPool=parpool(MatlabWorkers);
    end

    for ee=1:length(EsNo)
        dec_err=0;
        dec_berr=0;
        dec_frames=0;

        esno=10^(EsNo(ee)/10);
        while dec_err<=max_err
            parfor pp=1:MatlabWorkers       % Add parfor loop
                uc=randi([0  1], 1, k);      % Sequence of information bits
                c=mod(uc*G,2);                      % Encoding
                %c=[uc zeros(1, n-k)];                      % Encoding
                xc=1-2*c;
                sigma=sqrt(1/esno/2);             % Noise standard deviation (No/2)
                N=sigma*randn(1,n);
                yc=xc+N;                         % BI-AWGN channel

                %% Soft decision decoding
                x_test=1-2*mod(de2bi(0:2^k-1,k)*G,2);
                dm=(yc*x_test.');
                [dm_max,MLidx]=max(dm);
                u_hat_c=de2bi(MLidx-1,k);
                dec_berr_temp=sum(uc~=u_hat_c);
                if dec_berr_temp>0
                    dec_err=dec_err+1;
                    dec_berr=dec_berr+dec_berr_temp;
                end

                dec_frames=dec_frames+1;
                %disp(['Dec frames=' num2str(dec_frames)]);
                end
                %dec_err
                %dec_frames
        end

        dec_bits=dec_frames*k;
        BER(ee)=dec_berr/dec_bits;
        FER(ee)=dec_err/dec_frames;
        disp(['EsNo=' num2str(EsNo(ee))  'dB_coded BER=' num2str(BER(ee))]);
    end

else

    while dec_err<=max_err
        uc=randi([0  1], 1, k);      % Sequence of information bits
        c=mod(uc*G,2);                      % Encoding
        %c=[uc zeros(1, n-k)];                      % Encoding
        xc=1-2*c;
        sigma=sqrt(1/esno/2);             % Noise standard deviation (No/2)
        N=sigma*randn(1,n);
        yc=xc+N;                         % BI-AWGN channel

        %% Soft decision decoding
        x_test=1-2*mod(de2bi(0:2^k-1,k)*G,2);
        dm=(yc*x_test.');
        [dm_max,MLidx]=max(dm);
        u_hat_c=de2bi(MLidx-1,k);
        dec_berr_temp=sum(uc~=u_hat_c);
        if dec_berr_temp>0
            dec_err=dec_err+1;
            dec_berr=dec_berr+dec_berr_temp;
        end

        dec_frames=dec_frames+1;
        %disp(['Dec frames=' num2str(dec_frames)]);
        % end
        dec_err
        dec_frames
    end

    dec_bits=dec_frames*k;
    BER(ee)=dec_berr/dec_bits;
    FER(ee)=dec_err/dec_frames;
    disp(['EsNo=' num2str(EsNo(ee))  'dB_coded BER=' num2str(BER(ee))]);
end

EbNo=EsNo-10*log10(r);      % EbNo axis for simulated code
semilogy(EbNo,BER,'--ob','MarkerFaceColor','w');
delete(MatlabPool);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%