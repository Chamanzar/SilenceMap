clear all
close all
clc
%%
% SilenceMap with baseline algorithm test code
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

% Author: Alireza 	Date: 2020/05/03 17:00:48 	Revision: 0.1 
% Copyright: This project is licensed under the MIT License - see the LICENSE.md file for details

%% Simulations of regions of silence on a real headmodel:
cd ..
addpath('./MRI_prep_leadfield_ext')
addpath('./SilenceMap')
load OT_leadfield_symmetric_1662-128.mat % extrcated headmodel
load OT_headmodel_symmetric_1662-128.mat % calculated leadfield matrix

cortex = Cortex.Pial;
cortex1 = cortex;
L1 = L;
sulc = Cortex.Sulc  ;


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

%% %%%%%%%%%%%%%%%%%%%%%%% SilenceMap algorithm begins here %%%%%%%%%%%%%%%%%%
%% Ref lookup table:
i_ref_tot = [40;50;56;63;64;65;68;73;84;95];
Ref_names = {{'CPz'}    {'Pz'}    {'POz'}    {'Oz'}    {'Iz'}    {'Fpz'}    {'AFz'}    {'Fz'}    {'FCz'}    {'Cz'}};

% Ordered back to front: {'Iz'} {'Oz'} {'POz'} {'Pz'} {'CPz'} {'Cz'} {'FCz'} {'Fz'} {'AFz'} {'Fpz'}
%% Grid-search to choose the best reference electrode:

for ref_ind = 1:10
    
    %%%%%%%%%%%%SilenceMap in the low-res source grid%%%%%%%%%%%%%%%%%%%%%%
    
    load OT_leadfield_symmetric_818-128.mat % extrcated headmodel
    load OT_headmodel_symmetric_818-128.mat % calculated leadfield matrix
    
    cortex = Cortex.Pial;
    cortex1 = cortex;
    L1 = L;%forward matrix
    sulc = Cortex.Sulc  ;%sulci depth to plot the cortical surface
    
    n = size(L1,1);%number of EEG electrodes
    p = size(L1,2);%number of sources in the brain
    
    %% Referencing and noise estimation:
    
    eeg = double(eeg);
    M = eye(n-1);
    i_ref = i_ref_tot(ref_ind);
    M = [M(:,1:i_ref-1),-ones(n-1,1),M(:,i_ref:end)];
    Y = M*eeg;
    L1 = M*L1;
    
    %Noise estimation based on PSD:
    w_length = min(floor(0.5*size(Y,2)),256);
    w_over = floor(0.5*w_length);
    [pxx,f] = pwelch(Y',w_length,w_over,256,Fs);
    eta = mean(pxx(find(f==90):find(f==100),:),1)';%average over 90-100Hz freq band to obtain an est of flat PSD for noise
    sigma_z_sqrd = double(eta) * (100-0.1);
    Cz = diag(sigma_z_sqrd);%Estimated noise cov matrix
    %% Calculate variance reduction:
    M_Silence = Y;
    LL = L1'*L1;
    LL = (LL).^2;
    Var_norm_fact = sum(LL,2);
    Mu_tilda_Silence = (L1)'*M_Silence;
    %% Estimation of Var(\mu) = Var(Mu_tilda_Silence):
    [pxx,f] = pwelch(Mu_tilda_Silence',w_length,w_over,256,Fs);
    df = f(2)-f(1);
    sigma_mu_sqrd = sum(pxx,1)*df-mean(Mu_tilda_Silence,2)'.^2;%integration over psd to obtain sigma^2 = R_{mu,mu}(0)-mu^2
    sigma_mu_sqrd = double(sigma_mu_sqrd');
    
    %% Estimation of Var(eeg) = Var(M_Silence):
    [pxx,f] = pwelch(M_Silence',w_length,w_over,256,Fs);
    df = f(2)-f(1);
    sigma_eeg_sqrd = sum(pxx,1)*df-mean(M_Silence,2)'.^2;%integration over psd to obtain sigma^2 = R_{mu,mu}(0)-mu^2
    sigma_eeg_sqrd = double(sigma_eeg_sqrd');
    %%
    P_M = double((sigma_eeg_sqrd))-double((sigma_z_sqrd));%Estimated Scalp power wo noise
    sigma_mu_sqrd_wo_noise = sigma_mu_sqrd-diag(L1'*Cz*L1);
    Betta = sigma_mu_sqrd_wo_noise./Var_norm_fact;% The source contribution meausre (Beta_q)
    %% Finding silence center using low-res source grid:
    
    Gp = 5;% z^{gap} = 1cm
    source_loc = cortex1.vertices;
    % Calculate the source contribution meausre (\tilde{Beta}_q) for SilenceMap with baseline:
    Betta_new = hemispheric_base(Betta,source_loc,Gp);
    
    % figure;
    % h1 = plot_source_space_signal_vIV(Betta_new,sulc,cortex1,cortex1);
    % axis off
    % axis equal
    
    pow_flag = 0;% No power constraint in the low-res CSpeC
    k_search_grid = floor(linspace(2,200,20));
    Err = [];
    x_tot = [];
    
    P_normalized = P_M/max(P_M);
    for i=1:size(k_search_grid,2)
        k = k_search_grid(i);
        Sigma_s = ones(p,1);
        x = CSpeC(L1,cortex1,Betta_new,[],[],k,[],pow_flag,[]);% CSpeC function
        x_tot = cat(2,x_tot,x);
        [~,idx_x] = sort(x,'ascend');
        Sigma_s(idx_x(1:k)) = 0;
        P_hat = diag(L1*diag(Sigma_s)*L1');
        P_hat = P_hat/max(P_hat);
        Err = cat(1,Err,sum((P_normalized-P_hat).^2));%power mismatch (\Delta Pow)
    end
    
    [~,indx] = min(Err,[],1);
    k = k_search_grid(indx);
    
    x = CSpeC(L1,cortex1,Betta_new,[],[],k,[],pow_flag,[]);% CSpeC function
    X_det = zeros(p,1);
    X_det_L = X_det;
    
    [~,idx_x] = sort(x,'ascend');
    X_det(idx_x(1:k)) = 1;
    silence_center = mean(source_loc(X_det==1,:),1);
    % Plot function to show the localized region of silence:
    plot_source_space_signal_vF(X_det,sulc,cortex1,cortex1);
    
    %% %%%%%%%%%%SilenceMap in the high-res source grid%%%%%%%%%%%%%%%%%%%%%%
    
    P_M = diag(P_M);
    load OT_leadfield_symmetric_1662-128.mat % extrcated headmodel
    load OT_headmodel_symmetric_1662-128.mat % calculated leadfield matrix
    
    cortex1 = Cortex.Pial;
    L1 = L;
    sulc = Cortex.Sulc  ;
    
    n = size(L1,1);
    p = size(L1,2);
    
    %The effect of reference electrode:
    L1 = M*L1;
    
    
    % Calculate variance reduction:
    LL = L1'*L1;
    LL = (LL).^2;
    Mu_tilda_Silence = (L1'*M_Silence);
    
    
    %% Estimation of Var(\mu) = Var(Mu_tilda_Silence):
    
    [pxx,f] = pwelch(Mu_tilda_Silence',w_length,w_over,256,Fs);
    df = f(2)-f(1);
    sigma_mu_sqrd = sum(pxx,1)*df-mean(Mu_tilda_Silence,2)'.^2;%integration over psd to obtain sigma^2 = R_{mu,mu}(0)-mu^2
    sigma_mu_sqrd = double(sigma_mu_sqrd');
    sigma_mu_sqrd_wo_noise = (sigma_mu_sqrd-diag(L1'*Cz*L1));
    
    %% Beginning of iterations to estimate the source cov matrix Cs and region of silence alternatively:
    
    repp = 0;
    sigma_sq_hat = 1;
    gamma_p = 1;
    sigma_sq_hat_old = 0;
    silence_center_old = [0, 0, 0];
    silence_center_old_old = [0, 0, 0];
    Conv_time = [];
    
    Src_loc = cortex1.vertices;
    C_Exp = zeros(p);
    for i=1:p
        for j=i:p
            C_Exp(i,j) = (-norm(Src_loc(i,:)-Src_loc(j,:),2));
        end
    end
    C_Exp = C_Exp + permute(C_Exp,[2,1]);
    
    
    % Use the COM of the localized region in the low-res as an initial guess in the high-res grid:
    
    Cs = exp(C_Exp);
    dist_silence = sum((Src_loc-repmat(silence_center,p,1)).^2,2);
    [~,i_silence] = min(dist_silence);
    
    
    while((sum((silence_center_old-silence_center).^2,2)>100 ||...
            sum((silence_center_old_old-silence_center_old).^2,2)>100) && (repp<100))
        repp = repp + 1;
        P_e = diag(P_M);
        
        %% exponential decay in source correlation--identical dist:
        
        n = size(L1,1);
        p = size(L1,2);
        
        
        Cs_silent = Cs;
        Cs_silent(i_silence,:) = 0;
        Cs_silent(:,i_silence) = 0;
        elec_rank = diag(L1*Cs_silent*L1')./diag(L1*Cs*L1');
        [~,i_rank] = sort(elec_rank,'descend');
        
        % Finding the least square solutions for gamma and sigma_s:
        kk = 90;
        stp = 1;
        p = 1;
        lb = [zeros(1,p),0.1];
        ub = [];
        x0 = [1*ones(1,p),exp(1)];%initializations
        eqns = @(x)diag(L1(i_rank(1:stp:kk),:)*(diag(x(1:p))*(x(p+1).^C_Exp))*L1(i_rank(1:stp:kk),:)')...
            - diag(P_M(i_rank(1:stp:kk),i_rank(1:stp:kk)));
        options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective',...
            'MaxFunctionEvaluations', 1e6,'TolFun', 1e-30, 'TolX', 1e-15, 'MaxIterations', 5e2...
            ,'Display','iter');
        [xx,resnorm,res,eflag,output2] = lsqnonlin(eqns, x0, lb, ub, options);
        output2.funcCount
        sqrt(resnorm/length(1:stp:kk))
        
        sigma_sq_hat_old = sigma_sq_hat;
        sigma_sq_hat = diag(double(xx(1:p)));
        gamma_p = double(xx(p+1));
        Cs = sigma_sq_hat*(gamma_p.^C_Exp);%estimated source cov matrix (wo silence)
        
        
        P_base_hat = trace(L1*Cs*L1');
        P_ee = diag(L1*Cs*L1');
        
        
        %%
        p = size(L1,2);
        elec_index = i_rank(1:stp:kk);
        Betta = sigma_mu_sqrd_wo_noise./diag((L1'*L1)*Cs*(L1'*L1));
        Gp = 5;% z^{gap} = 1cm
        Betta_new = hemispheric_base(Betta,Src_loc,Gp);
        
        %     figure;
        %     h1 = plot_source_space_signal_vIV(Betta_new,sulc,cortex1,cortex1);
        %     axis off
        %     axis equal
        
        k_search_grid = floor(linspace(2,100,20));
        Err = [];
        x_tot = [];
        
        Epsil = (res.^2);
        pow_flag = 1;% power constraints in the high-res CSpeC
        
        for i=1:size(k_search_grid,2)
            k = k_search_grid(i);
            x = CSpeC(L1,cortex1,Betta_new,Cs,P_M,k,Epsil,pow_flag,elec_index);% CSpeC function
            x_tot = cat(2,x_tot,x);
            [~,idx_x] = sort(x,'ascend');
            Cs_silence = Cs;
            Cs_silence(idx_x(1:k),:)=0;
            Cs_silence(:,idx_x(1:k))=0;
            P_hat = diag(L1*Cs_silence*L1');
            P_hat = P_hat/max(P_hat);
            Err = cat(1,Err,sum((P_normalized-P_hat).^2));%power mismatch (\Delta Pow)
        end
        %
        [~,indx] = min(Err,[],1);
        k = k_search_grid(indx);

        x = CSpeC(L1,cortex1,Betta_new,Cs,P_M,k,Epsil,pow_flag,elec_index);% CSpeC function
        
        %% K-means in 1D:
        X_det = zeros(p,1);
        
        [~,idx_x] = sort(x,'ascend');
        X_det(idx_x(1:k)) = 1;
        % Plot function to show the localized region of silence:
        plot_source_space_signal_vF(X_det,sulc,cortex1,cortex1);
        i_silence = find(X_det);
        silence_center_old_old = silence_center_old;
        silence_center_old = silence_center;
        silence_center = mean(Src_loc(X_det==1,:),1);
        d_center = sum((silence_center_old-silence_center).^2,2);
        Result_file = sprintf('Noisy_OT_wHB_kk_%d_%d_indx_center_%d_i_ref_%d.mat', kk, repp, indx_center, i_ref)
        save(Result_file, 'X_det', 'X_act', 'd_center', 'sigma_sq_hat', 'gamma_p', 'X_det_L', 'SNR_M', 'Err', '-v7.3');
        
        
    end
end



