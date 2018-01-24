function [temporal_plv,Times]=moving_multitrial_pac2_sur(mode,signals,srate,Pf1,Pf2,Af1,Af2,window,step)
%estimation of phase-amplitude cross-frequency coupling
%
%[temporal_plv,Times]=moving_multitrial_pac2(mode,signals,srate,Pf1,Pf2,Af1,Af2,window,step)
% Current version: [20160322]
% 
% *** Inputs ***
% mode: use [0] for normal estimation and [1] for surrogate control.         
%       When [1], single-block resampling technique is applied for single-trial data and full-phase between-trials shuffling for multi-trial data. 
%       See J. Aru et al., "Untangling cross-frequency coupling in neuroscience", Current Opinion in Neurobiology 2015 
%
% signals = Single_trial- time series
% srate = frequency range
% Pf1,Pf2 = frequency range of the low-frequency 
% Af1,Af2 = frequency range of the high-frequency
%    window,step : width of temporal segments.... ; stepping window
%    Times -> the central-latency of each segment      
%          
% Outputs:
% temporal_plv = Phase-Amplitude Coupling array over the different temporal segments
%
%Cohen MX."Assessing transient cross-frequency coupling in EEG data"
% Journal of Neuroscience Methods 168 (2008) 494?499
%
%W.D. Penny et al., "Testing for nested oscillation"
%Journal of Neuroscience Methods 174 (2008) 50?61
%
% N.Laskaris 14/2/2015
%http://neuroinformatics.gr
%
% Rev.20151101, D.A. Adamos 01/Nov/2015, http://neuroinformatics.gr
%       -Code optimization
%       -Surrogate control feature added.



%get the phase of the low-frequency band
    if (nargin<9)  
    disp('Check number of arguments');
    return;
    end

    Ntrials=size(signals,1);
    [bb_p,aa_p]=butter(3,[Pf1/(srate/2),Pf2/(srate/2)]);
    low_filtered_signals=filtfilt(bb_p,aa_p,signals')';
    Phase_signals=angle(hilbert(low_filtered_signals'))'; % this is getting the phase time series

    if(mode)
        if(Ntrials==1)
            %P=randperm(size(Phase_signals,2),1);
            P=randsample([0.25*size(Phase_signals,2):0.75*size(Phase_signals,2)],1,true);
            surr_phase_signal=[Phase_signals(:,P+1:end) Phase_signals(:,1:P)];
            Phase_signals=surr_phase_signal;
        else
            P=randperm(Ntrials);
            surr_phase_signal=Phase_signals(P,:);
            Phase_signals=surr_phase_signal;
        end
            
    end
    %get the phase of high-frequency band

    %STEP 1: Filtering the original signal in the range of high-frequency range
                               % just filtering
    [bb,aa]=butter(3,[Af1/(srate/2),Af2/(srate/2)]);
    high_filtered_signals=filtfilt(bb,aa,signals')';

    %STEP 2:Get the analytic signal 
    Env_high_filtered_signals=abs(hilbert(high_filtered_signals'))'; % getting the amplitude envelope

    %STEP 3: Filtering the obrained envelope of the high-frequency range within the
    %frequency range of the low-frequency band
    lowfromhigh=filtfilt(bb_p,aa_p,Env_high_filtered_signals')'; 
    low_Env_high_filtered_signals=lowfromhigh;%-repmat(mean(lowfromhigh),Ntrials,1);

    %STEP 4:Get the phase
    Amp_phase_signals=angle(hilbert(low_Env_high_filtered_signals'))';

    Ntime=size(signals,2);

    temporal_plv=[]; 
    for ii=1:(Ntime-window)/step
        tt=[(ii-1)*step+1:(ii-1)*step+window];
        plv=abs(sum(sum(exp(1i*(Phase_signals(:,tt)-Amp_phase_signals(:,tt))))))/(Ntrials*(window));
        temporal_plv(ii)=plv;
        Times(ii)=tt(round(window/2));
    end
   

end