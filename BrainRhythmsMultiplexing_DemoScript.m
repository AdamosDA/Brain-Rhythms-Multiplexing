clear all,close all
load single_subject.mat % load a single-subject data (average-rereferenced data)
                        % a music-listening and a resting-state dataset and
                        % along with the sensors' names and coordinates

FBANDS=[1 4; 4 8; 8 12; 12 20; 20 30; 30 45; 55 95]; % Define 7 Frequency bands
Fs=500; % Sampling Frequency


%%%STEP-1
%%%%%%%% Building Connectivity-Patterns; for each Brain Rhythm (frequency band) independently %%%%%%%%

% Rest
 REST_BAND_CPs=[];  % [(29x28/2)=406 pairs x 7 Bands]
 for i_band=1:size(FBANDS,1);
    [b,a]=butter(3,[FBANDS(i_band,1),FBANDS(i_band,2)]/(Fs/2));
    EEG=rest_eeg;
    filtered_EEG=filtfilt(b,a,EEG')';
    PLV_CP=[];[PLV_CP,IDX]=PLV_multichannel_signal(filtered_EEG); %1 column vector (vectorized WA-matrix)
     REST_BAND_CPs(:,i_band)=PLV_CP;    
 end    
 
% Music listening
 MUSIC_BAND_CPs=[];  % [(29x28/2)=406 pairs x 7 Bands]
 for i_band=1:size(FBANDS,1);
    [b,a]=butter(3,[FBANDS(i_band,1),FBANDS(i_band,2)]/(Fs/2));
    EEG=music_eeg;
    filtered_EEG=filtfilt(b,a,EEG')';
    PLV_CP=[];[PLV_CP,IDX]=PLV_multichannel_signal(filtered_EEG); %1 column vector (vectorized WA-matrix)
    MUSIC_BAND_CPs(:,i_band)=PLV_CP;    
 end    
 
 
 
%%%STEP-2
%%%%%%%% NetWork analysis - Defining Modular-organization based on Dominant Sets Algorithm %%%%%%%%

% The Dominant Sets algorithm ('iterative_dominant_set_extraction.m') follows an iterative, 
% subtractive-clustering approach that includes a random initialization. 
% Thus, the overall result may change from exacution to execution. 
% The ''official'' implementation is realized via 'recurrent_dominant_sets.m', 
% where after a number of repetition the ''best'' graph-partitioning result
% is finally outputted.

% Important Notice: For this demonstration, we will use a rng(seed) approach to minimize the execution time.
% For actual cases, edit the 'recurrent_dominant_sets.m' file by replacing the rng(seed) with the rng('shuffle'). 
% Then use a large ITERATION number below (e.g. 500) to guarantee optimal performance.  

ITERATIONS=1; 
MUSIC_DSgroups=[]; MUSIC_DSscores=[]; REST_DSgroups=[]; REST_DSscores=[];
MUSIC_DSscores=zeros(7,7);REST_DSscores=zeros(7,7);

tic
for i_band=1:size(FBANDS,1),
CP=REST_BAND_CPs(:,i_band);
WA=[];for i1=1:size(IDX,1); WA(IDX(i1,1),IDX(i1,2))=CP(i1);WA(IDX(i1,2),IDX(i1,1))=CP(i1); end
clear sorted_score sorted_cost_function;
[sorted_groups,sorted_cost_function]=recurrent_dominant_sets(WA,ITERATIONS);
REST_DSgroups(i_band,:)=sorted_groups;REST_DSscores(i_band,1:numel(sorted_cost_function))=sorted_cost_function;
end

for i_band=1:size(FBANDS,1),
CP=MUSIC_BAND_CPs(:,i_band);
WA=[];for i1=1:size(IDX,1); WA(IDX(i1,1),IDX(i1,2))=CP(i1);WA(IDX(i1,2),IDX(i1,1))=CP(i1); end
clear sorted_score sorted_cost_function;
[sorted_groups,sorted_cost_function]=recurrent_dominant_sets(WA,ITERATIONS);
MUSIC_DSgroups(i_band,:)=sorted_groups;MUSIC_DSscores(i_band,1:numel(sorted_cost_function))=sorted_cost_function;
end
toc

lim1=0; lim2=max(max([REST_DSscores MUSIC_DSscores]));

% Modular organization as a function of brain rhythm and recording condition. 
% Nodes belonging to the same functional community share color and size. Size indicates the level of within-group cohesiveness. 
% Presentation-format is similar to Figure 4 of the paper. 
% (i.e. upper line: Rest, lower line: Music-listening;  the brain-rhythms are presented in order of frequency,from left-to-right  ).

figure(1),clf,set(gcf,'color','w');clrs=colormap(lines(7));clrs=[0 0 0; clrs]
set(gcf,'units','points','position',[500,500,1600,400])
for j=1:7
        subplot(2,7,j),hold
        groups=REST_DSgroups(j,:);        
        spX=sensor_coordinates(:,1);spY=sensor_coordinates(:,2); 
        for i=1:29,if(groups(i)),size_index=max(1,round(((REST_DSscores(j,groups(i))-lim1)/(lim2-lim1))*ceil(lim2*10)));plot(spX(i),spY(i),'bo','markersize',10+1.5*size_index,'MarkerFaceColor',clrs(groups(i)+1,:)),axis([-200 200 -200 200]),end;end; %for ii=1:29; text(spX(ii)+5,spY(ii),sensors_names(ii,:)), end             
        for ii=1:29; text(spX(ii)-0.025*range(spX),spY(ii)+0.001*range(spY),num2str(ii),'fontsize',8), end, axis off;         

        subplot(2,7,7+j),hold
        groups=MUSIC_DSgroups(j,:);
        for i=1:29,if(groups(i)),size_index=max(1,round(((MUSIC_DSscores(j,groups(i))-lim1)/(lim2-lim1))*ceil(lim2*10)));plot(spX(i),spY(i),'bo','markersize',10+1.5*size_index,'MarkerFaceColor',clrs(groups(i)+1,:)),axis([-200 200 -200 200]),end;end; %for ii=1:29; text(spX(ii)+5,spY(ii),sensors_names(ii,:)), end             
        for ii=1:29; text(spX(ii)-0.025*range(spX),spY(ii)+0.001*range(spY),num2str(ii),'fontsize',8), end, axis off;         
        
end

%%%STEP-3
%%%%%%%%%  Exploring the change in modular organization of a brain rhythm by means of VI distance.

FRR_names=['delta '; 'theta ';'alpha ';'betaL ';'betaH ';'gammaL';'gammaH']; 
VI_DM=zeros(7,48);VI_DM=zeros(size(FBANDS,1),1);
for i_band=1:size(FBANDS,1);
            VI_DM(i_band)=partition_distance(REST_DSgroups(i_band,:),MUSIC_DSgroups(i_band,:));                    
end
figure(2),clf,bar(VI_DM,'k'); set(gca, 'XTick', 1:length(FRR_names),'XTickLabel',FRR_names);grid
title('VI distance of MusicListening modular pattern from Rest')



%%%STEP-4
%%%%%  Estimating Phase-Amplitude Coupling (PAC) for DELTA->BETA_high interaction

STEP=1; WINDOW=size(music_eeg,2)-1;% PAC Window = 8s
Pf1=1;Pf2=4; % Delta band
Af1=20;Af2=30; % Beta-high band


%% REST PAC
    for I=1:length(Pf1)       
          for ELECTRODE=1:size(rest_eeg,1)
                  clear REST_curve;REST_curve=squeeze(rest_eeg(ELECTRODE,:));
                  [REST_pac(I,ELECTRODE,:),EEG_REST_Times]=moving_multitrial_pac2_sur(0,REST_curve,Fs,Pf1,Pf2,Af1,Af2,WINDOW,STEP) ; 
              end
    end

%% MUSIC PAC
    for I=1:length(Pf1)       
          for ELECTRODE=1:size(music_eeg,1)
                  clear MUSIC_curve;MUSIC_curve=squeeze(music_eeg(ELECTRODE,:));
                  [MUSIC_pac(I,ELECTRODE,:),EEG_MUSIC_Times]=moving_multitrial_pac2_sur(0,MUSIC_curve,Fs,Pf1,Pf2,Af1,Af2,WINDOW,STEP) ; 
          end       
    end
    
figure(3),clf;set(gcf,'color','w');
bar(MUSIC_pac-REST_pac,'k');title('Increase of PAC (MusicListening - Rest)')
set(gca, 'XTick', 1:length(sensors_names),'XTickLabel',sensors_names);grid,


