function [PLV_CP,IDX]=PLV_multichannel_signal(filtered_signals)

%[PLV_CP,IDX]=PLV_multichannel_signal(filtered_signals)
%                signals should be band-pass filtered 
%      PLV_CP: vectorized connectivity pattern , IDX: the indices to turn it to a Connectivity Matrix  
%
%  based on    LACHAUX "Measuring Phase Synchrony in Brain Signals" et al.1999

[Nsensors,Ntime]=size(filtered_signals);
A_fX=angle(hilbert(filtered_signals'))';  % no phase unwrap

IDX=[]; PLV_CP=[]; 
for i_sensor=1:Nsensors; for j_sensor=i_sensor+1:Nsensors; IDX=[IDX;i_sensor,j_sensor]; %i_sensor/Nsensors
phase1=A_fX(i_sensor,:);phase2=A_fX(j_sensor,:);
plv=abs(sum(exp(j*(phase1-phase2))))/Ntime ;
PLV_CP=[PLV_CP;plv];
    end ,end 




