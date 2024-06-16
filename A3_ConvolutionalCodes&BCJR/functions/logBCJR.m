%% BCJR implementation (in log domain) 

%% Inputs: 
% yin - Soft channel output 
% FSM       - Finite state machine structure
%           .StateLaw   (SxM matrix where S is the number of states and M is the constellation cardinality. ij entry contains output state for input state i and input j)  
%           .OutputLaw  (SxS matrix ij entry contains FSM output for input state i and output state j)  
%           .State      (Sx1 string array containing state labels)    
%           .Input      (Mx1 Input vector)
%           .BranchMetrics (SxM matrix where ij entry contains branch metric corresponding to state i and input j)  


%% Outputs
%  yout     - Sequence of symbol-wise log-APP ratios (with respect to first input of the FSM) 


function yout=logBCJR(yin,FSM,channel)


M=length(FSM.Input);                         % Modulation cardinality
S=length(FSM.State);                         % Number of states in the trellis or trellis depth
T=length(yin);
N0=channel.N0;

%% Initialise alphas and betas (at time t=0 and t=T respectively)
Ceil=-1e9;
alpha=Ceil*ones(S,T+1); 
alpha(1,1)=0;
beta=Ceil*ones(S,T+1);  
beta(1,T+1)=0;
gamma=Ceil*ones(S,S,T+1);
gamma(:,1,T+1)=0;
yout=Ceil*ones(M,T);

%% Backward recursion (beta computation)
for tt=1:T
  
    for ss=1:S
        LinkedState=FSM.StateLaw(ss,:);           % States in section tt+1 linked to state ss in section tt
        for ll=1:length(LinkedState)
            branch=beta(LinkedState(ll)+1,T+2-tt)+gamma(ss,LinkedState(ll)+1,T+2-tt);  
            beta(ss,T+1-tt)=logsum(beta(ss,T+1-tt),branch);      % Add future beta value
        end
    end
    %% Compute branch metrics
        gamma(:,:,T+1-tt)=Gamma(yin(T+1-tt),FSM.OutputLaw,N0); 
end


%% Forward recursion (alpha computation)
for tt=1:T
    %% Trellis section iteration 
    for ss=1:S
        idx=(FSM.StateLaw==(ss-1));                          
        link=sum(idx,2);       
        LinkedState=find(link>0);             % States in section tt linked to state ss in section tt+1
    
        for ll=1:length(LinkedState)
            branch=alpha(LinkedState(ll),tt)+gamma(LinkedState(ll),ss,tt);   
            alpha(ss,tt+1)=logsum(alpha(ss,tt+1),branch);                % Add past alpha value  
        end
    end

    %% Compute log-APP
    for ii=1:M
        StateOut=FSM.StateLaw(:,ii);
        for jj=1:S
            branch=alpha(jj,tt)+gamma(jj,StateOut(jj)+1,tt)+beta(StateOut(jj)+1,tt);
            yout(ii,tt)=logsum(yout(ii,tt),branch);
        end
    end
    %yout(:,tt)=yout(:,tt)-yout(1,tt);      % log-APP ratios with respect to input 1
end


end



%% logsum function 
function ls=logsum(x,y)
ls=max(x,y)+log(1+exp(-abs(x-y)));
end

% Function computing gamma branch metrics under Gaussian assumption 
function gamma=Gamma(yin,u,N0)
gamma=-abs(yin-u).^2/(N0); 
end
