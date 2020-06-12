clear all
close all
clc

%%
% Modified sLORETA, MNE, and MUSIC algorithms test code
% This code generates, tests, and saves randomly simulated regions of silence.
%
% See: README.txt and [1] for more info.

% [1] A. Chamanzar, M. Behrmann, and P. Grover,
%  "Neural silences can be localized using noninvasive scalp EEG",
%   To be submitted to Nature BME, 2020.


% L(nxp): leadfield/forward matrix
% S(pxt): source matrix
% eeg(nxt): sensor matrix measured
% E(nxt): measurement noise matrix

% eeg = LS+E: forward model

% Author: Alireza 	Date: 2020/05/18 12:00:48 	Revision: 0.1
% Copyright: This project is licensed under the MIT License - see the LICENSE.md file for details

%% Simulations of regions of silence on a real headmodel:

for repp = 1:70 % number of repetition of the simulations
    
    
    Delta_K_tot = [];
    Bias_tot = [];
    JI_tot = [];
    
    
    cd ..
    addpath('./MRI_prep_leadfield_ext')
    addpath('./SilenceMap')
    load OT_leadfield_symmetric_1662-128.mat % extrcated headmodel
    load OT_headmodel_symmetric_1662-128.mat % calculated leadfield matrix
    
    cortex = Cortex.Pial;
    cortex1 = cortex;
    L1 = L;
    sulc = Cortex.Sulc;
    
    
    source_loc = Cortex.Pial.vertices;
    
    t = 10000;
    
    k_original = 50;
    p = size(L1,2);
    n = size(L1,1);
    
    Gp = 5;% z^{gap} = 1cm
    
    ind_Long_fiss = find(source_loc(:,1)>-Gp & source_loc(:,1)<Gp);
    
    %initialization:
    GT_dist_to_sensor = 40;
    silence_indices = ind_Long_fiss;
    
    while(mean(GT_dist_to_sensor)>30 || ~isempty(intersect(silence_indices,ind_Long_fiss)))
        
        rng('shuffle')
        indx_center = randperm(size(source_loc,1),1)
        dist_vec = sum((source_loc- repmat(source_loc(indx_center,:),size(source_loc,1),1)).^2,2);
        [~,I_dist] = sort(dist_vec);
        
        
        silence_indices = (I_dist(1:k_original));
        
        X_act = zeros(size(L1,2),1);
        X_act(silence_indices) = 1;
        
        %Calculate the distance of simulated region of silence to the sensor space:
        GT_loc = cortex1.vertices(find(X_act),:);
        GT_dist_to_sensor = [];
        for j=1:size(GT_loc,1)
            GT_dist_to_sensor = cat(1,GT_dist_to_sensor,min(sqrt(sum((GT_loc(j,:)-sensor_locs).^2,2)),[],1));
        end
        
    end
    
    %Plot the simulated region of silence:
    plot_source_space_signal_vF(X_act,sulc,cortex1,cortex1);
    title('Red region = Actual region of silence','FontSize',20)
    
    gamma = 0.12;% Exp decay coeff of the source cov matrix (Cs)
    C_Exp = zeros(p);
    for i=1:p
        for j=i:p
            C_Exp(i,j) = gamma*(-norm(source_loc(i,:)-source_loc(j,:),2));
        end
    end
    C_Exp = C_Exp + permute(C_Exp,[2,1]);
    Sigma_S = exp(C_Exp);%Source cov matrix (Cs^{full}) without silence
    
    Noise_pow = 0.5e-7;% maximum noise power which controls the SNR of the simulated scalp signals
    Sigma_n = diag(Noise_pow*rand(n,1));%Noise cov matrix (Cz)
    mu = zeros(p,1);
    mu_n = zeros(n,1);
    
    S = mvnrnd(mu,Sigma_S,t)';
    E = mvnrnd(mu_n,Sigma_n,t)';
    S_silence = S;
    S_silence(silence_indices,:) = 0;%Simulation of the silent sources in S
    eeg = L1*S_silence;
    
    % Simulation of the eeg signals with frequency components upto 90Hz:
    Fs = 512;
    dd = designfilt('lowpassiir','FilterOrder',8, ...
        'PassbandFrequency',90,'PassbandRipple',0.2, ...
        'SampleRate',Fs);
    eeg = filter(dd,(eeg)')';
    eeg = eeg + E;%the simulated scalp eeg signals
    
    %Calculation of the SNR:
    eeg_wo_silence = L1*S;
    eeg_wo_silence = filter(dd,(eeg_wo_silence)')';
    SNR_M = mean(10*log10(var(eeg_wo_silence,[],2)./Noise_pow));%the avreage SNR of the simulated scalp signals
    
    %% %%%%%%%%%%%%%%%%%%%%%%% Modified source localization begins here %%%%%%%%%%%%%%%%%%
    %% Ref lookup table:
    
    i_ref_tot = [40;50;56;63;64;65;68;73;84;95];
    Ref_names = {{'CPz'}    {'Pz'}    {'POz'}    {'Oz'}    {'Iz'}    {'Fpz'}    {'AFz'}    {'Fz'}    {'FCz'}    {'Cz'}};
    % Ordered back to front: {'Iz'} {'Oz'} {'POz'} {'Pz'} {'CPz'} {'Cz'} {'FCz'} {'Fz'} {'AFz'} {'Fpz'}
    
    %%%%%%%%%%%%Modified source localization in the low-res source grid%%%%%%%%%%%%%%%%%%%%%%
    
    load OT_leadfield_symmetric_818-128.mat % extrcated headmodel
    load OT_headmodel_symmetric_818-128.mat % calculated leadfield matrix
    
    cortex = Cortex.Pial;
    cortex1 = cortex;
    L1 = L;%forward matrix
    sulc = Cortex.Sulc;%sulci depth to plot the cortical surface
    
    n = size(L1,1);%number of EEG electrodes
    p = size(L1,2);%number of sources in the brain
    
    %% Referencing and noise estimation:
    
    eeg = double(eeg);
    M = eye(n-1);
    
    i_ref = i_ref_tot(10);
    M = [M(:,1:i_ref-1),-ones(n-1,1),M(:,i_ref:end)];
    Y = M*eeg;
    L1 = M*L1;
    %%
    
    XMNE_initial = [];
    XMUSIC_initial = [];
    XLoreta_initial = [];
    
    LV = 1;
    
    [XMNE, XLORETA,  XMUSIC, I_MNE, I_sLoreta, I_MUSIC] = modified_src_loc(LV,L1,cortex.vertices,Y...
        ,XMNE_initial,XLoreta_initial,XMUSIC_initial);
    
    XMNE_initial = mean(XMNE,1);
    XMUSIC_initial = mean(XMUSIC,1);
    XLoreta_initial = mean(XLORETA,1);
    
    X_det_MUSIC = zeros(size(L1,2),1);
    X_det_MUSIC(I_MUSIC) = 1;
    
    X_det_MNE = zeros(size(L1,2),1);
    X_det_MNE(I_MNE) = 1;
    
    X_det_Loreta = zeros(size(L1,2),1);
    X_det_Loreta(I_sLoreta) = 1;
    
    %% Plots the output of the low-res localization:
    plot_source_space_signal_vF(X_det_MUSIC,sulc,cortex1,cortex1);
    plot_source_space_signal_vF(X_det_MNE,sulc,cortex1,cortex1);
    plot_source_space_signal_vF(X_det_Loreta,sulc,cortex1,cortex1);
    
    %% %%%%%%%%%%Modified source localization in the high-res source grid%%%%%%%%%%%%%%%%%%%%%%
    
    load OT_leadfield_symmetric_1662-128.mat % extrcated headmodel
    load OT_headmodel_symmetric_1662-128.mat % calculated leadfield matrix
    
    cortex1 = Cortex.Pial;
    L1 = L;
    sulc = Cortex.Sulc;
    
    n = size(L1,1);
    p = size(L1,2);
    
    %The effect of reference electrode:
    L1 = M*L1;
    %%
    LV = 0;
    ind = 1:size(L1,2);
    [XMNE, XLORETA,  XMUSIC, I_MNE, I_sLoreta, I_MUSIC] = modified_src_loc(LV,L1,cortex1.vertices,Y...
        ,XMNE_initial,XLoreta_initial,XMUSIC_initial);
    
    X_det_MUSIC = zeros(size(L1,2),1);
    X_det_MUSIC(I_MUSIC) = 1;
    
    X_det_MNE = zeros(size(L1,2),1);
    X_det_MNE(I_MNE) = 1;
    
    X_det_Loreta = zeros(size(L1,2),1);
    X_det_Loreta(I_sLoreta) = 1;
    
%     close all
    
    %% Plots the output of the high-res localization:
    plot_source_space_signal_vF(X_det_MUSIC,sulc,cortex1,cortex1);
    plot_source_space_signal_vF(X_det_MNE,sulc,cortex1,cortex1);
    plot_source_space_signal_vF(X_det_Loreta,sulc,cortex1,cortex1);
    
    %% Calculate the performance, i.e., \Delta COM & JI:
    
    Bias = [];
    Delta_K = [];
    JI = [];
    JI = cat(2,JI,sum(X_act & X_det_MUSIC)/sum(X_act | X_det_MUSIC));
    JI = cat(2,JI,sum(X_act & X_det_MNE)/sum(X_act | X_det_MNE));
    JI = cat(2,JI,sum(X_act & X_det_Loreta)/sum(X_act | X_det_Loreta));
    
    
    COM = mean(cortex1.vertices(find(X_det_MUSIC),:),1);
    COM_act = mean(cortex1.vertices(find(X_act),:),1);
    dcenters = sqrt(sum((COM_act-COM).^2));
    Bias = cat(1,Bias,dcenters);
    
    COM = mean(cortex1.vertices(find(X_det_MNE),:),1);
    COM_act = mean(cortex1.vertices(find(X_act),:),1);
    dcenters = sqrt(sum((COM_act-COM).^2));
    Bias = cat(1,Bias,dcenters);
    
    COM = mean(cortex1.vertices(find(X_det_Loreta),:),1);
    COM_act = mean(cortex1.vertices(find(X_act),:),1);
    dcenters = sqrt(sum((COM_act-COM).^2));
    Bias = cat(1,Bias,dcenters);
    
    Delta_K = cat(1,Delta_K,abs(sum(X_act,1)-sum(X_det_MUSIC,1))/sum(X_act,1));
    Delta_K = cat(1,Delta_K,abs(sum(X_act,1)-sum(X_det_MNE,1))/sum(X_act,1));
    Delta_K = cat(1,Delta_K,abs(sum(X_act,1)-sum(X_det_Loreta,1))/sum(X_act,1));
    
    Delta_K_tot = cat(2,Delta_K_tot,Delta_K);
    Bias_tot = cat(2,Bias_tot,Bias);
    JI_tot = cat(2,JI_tot,JI');
    
end

save MNE_MUSIC_LORETA JI_tot Bias_tot Delta_K_tot

%%

load MNE_MUSIC_LORETA

mean(JI_tot,2)
std(JI_tot,0,2)/sqrt(size(JI_tot,2))

mean(Bias_tot,2)
std(Bias_tot,0,2)/sqrt(size(JI_tot,2))

mean(Delta_K_tot,2)
std(Delta_K_tot,0,2)/sqrt(size(JI_tot,2))







