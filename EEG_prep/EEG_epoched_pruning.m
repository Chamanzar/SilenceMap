close all
clear all
clc
%% 
% This code does the preprocessing steps for the EEG data in Silence
% localization
%
% See: README.txt and [1] for more info.

% [1] A. Chamanzar, M. Behrmann, and P. Grover,
%  "Neural silences can be localized using noninvasive scalp EEG",
%   To be submitted to Nature BME, 2020.

% Author: Alireza 	Date: 2020/05/19 12:00:48 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details

%% Preprocessing steps for EEG data

current_path = pwd;
% eegpath = '/Path-to-EEGLab/';
eegpath = '/Users/Alireza/Downloads/eeglab2019_0';

cd(eegpath)
eeglab
cd(current_path)

%% Read in bdf file, average data to mastoids (electrodes 129 130), bandpass filter 1-100Hz, and assign channel names
% read in bdf
EEG = pop_biosig('SN_SSVEP.bdf', 'ref',[129 130] ,'refoptions',{'keepref' 'off'});
EEG=pop_chanedit(EEG, 'lookup',[eegpath '/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp']);

EEG.setname='SN';
EEG = eeg_checkset( EEG );
% filter
EEG = pop_eegfiltnew(EEG, 1,100,16896,0,[],1);
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename','SN_filtered.set','filepath',current_path);
close all
% STOP
%% visually check electrodes and interpolate if necessary. If all are fine, then run...
EEG = pop_loadset('filename','SN_filtered.set','filepath',current_path);        % attach electrode names
EEG = eeg_checkset( EEG );
pop_eegplot(EEG, 1, 1, 1);
CH_ind_vis = [128];% index of bad channel through visual inspection
EEG = eeg_interp(EEG, CH_ind_vis);
EEG = pop_saveset( EEG, 'filename','SN_filtered.set','filepath',current_path);
% STOP
%% eye electrodes are rereferenced as bipolar electrodes and relabelled. Then run ICA        
EEG = pop_loadset('filename','SN_filtered.set','filepath',current_path);      
% make eye electrodes bipolar electrodes and relabel
EEG = pop_eegchanoperator( EEG, {  'ch142 = ch129 - ch130 label hEOG', 'ch143 = ch132 - ch131 label vEOG'} , 'ErrorMsg', 'popup', 'Warning','on' );
EEG = pop_select( EEG,'nochannel',{'EXG3' 'EXG4' 'EXG5' 'EXG6' 'EXG7' 'GSR1' 'GSR2' 'Erg1' 'Erg2' 'Resp' 'Plet' 'Temp'});
pop_eegplot(EEG, 1, 1, 1);
EEG = eeg_checkset( EEG );
% run ICA (takes a while)
EEG = pop_runica(EEG, 'extended',1,'interupt','on');
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename','SN_SSVEP_postica.set','filepath',current_path);

%% Artifact removal based on x-correlation:
removed_ica = [];
EEG = pop_loadset('filename','SN_SSVEP_postica.set','filepath',current_path);
Y = EEG.data;   % set channel data to matrix Y
ica_weights = EEG.icaweights; % copy icaweights matrix
Y_ICA = ica_weights*Y;
Roh = corr(Y(129:end,:)',Y_ICA'); % xcorrelation of the artifact electrodes (lEOG, hEOG, ECG) with the ICA components
[~,Art_ind] = find(abs(Roh)>0.8*max(abs(Roh),[],2)); % finding the artifact ICA components--manually adjust the threshold
ind = EEG.icachansind;
ind(unique(Art_ind)) = [];
%         EEG_clean = icaproj(Y,ica_weights,ind);
%         EEG.data = EEG_clean;
EEG = pop_subcomp(EEG,Art_ind,1); % Removing the artifact ICA components
% removed_ica.ind = Art_ind; 
EEG = pop_saveset( EEG, 'filename','SN_SSVEP_postica_pruned.set','filepath',current_path);
EEG = eeg_checkset( EEG );

%% Epoch extraction:
EEG = pop_loadset('filename','SN_SSVEP_postica_pruned.set','filepath',current_path);
EEG = pop_epoch( EEG, { '1','2' }, [-1 2], 'newname', 'Continuous EEG Data epochs', 'epochinfo', 'yes'); % Epoch dataset using 'visual and rest' tasks
EEG = pop_rmbase( EEG, [-1000 0]); % remove baseline
EEG = pop_saveset( EEG, 'filename','SN_SSVEP_postica_pruned_epoched.set','filepath',current_path);
EEG = eeg_checkset( EEG );
eeglab redraw % redraw eeglab to show the new epoched dataset

%% Check electrodes and interpolate if necessary:
EEG = pop_loadset('filename','SN_SSVEP_postica_pruned_epoched.set','filepath',current_path);
EEG = eeg_checkset( EEG );
CH_stat = [];
for cnum=1:128
    [M,SD,sk,k,med,zlow,zhi,tM,tSD,tndx,ksh] = pop_signalstat( EEG, 1, cnum, 5);
    CH_stat = cat(1,CH_stat,[sk,k,SD]);
    close all
end
CH_ind = find(abs(CH_stat(:,1))>0.45 | CH_stat(:,2)>18 | CH_stat(:,3)>16); % manually adjust the thresholds
CH_ind = unique([CH_ind;CH_ind_vis]); % concatenate the indices of bad electrodes 
% detected by visual inspection and the elec statistics

% Removing the bad electrodes and their corresponding mirrored electrodes 
load SN_headmodel_symmetric_1780-128.mat
sensor_loc = headmodel.elec.chanpos;
sensor_loc_Mirrored = sensor_loc;
sensor_loc_Mirrored(:,1) = -sensor_loc_Mirrored(:,1);
CH_ind_ext = [];
for i=1:size(CH_ind,1)
    [~,ind_temp] = sort(sum((sensor_loc(CH_ind(i),:)-sensor_loc_Mirrored).^2,2),1);
    CH_ind_ext = cat(1,CH_ind_ext,ind_temp(1));
end

CH_ind = unique([CH_ind;CH_ind_ext]);  
EEG = eeg_interp(EEG, CH_ind); % interpolating the bad electrodes
EEG = pop_saveset( EEG, 'filename','SN_SSVEP_postica_pruned_epoched_Channel.set','filepath',current_path);
% Saving the indices of the removed/interp electrodes
save removed_channel_SN_SSVEP.mat CH_ind

%% Prune epochs:

EEG = pop_loadset('filename','SN_SSVEP_postica_pruned_epoched_Channel.set','filepath',current_path);
[EEG,~] = pop_eegthresh( EEG, 1, 1:128, -75, 75, -1000, 1998,0,1);% manually adjust the thresholds
EEG = pop_rejtrend( EEG, 1, 1:128, 1024, 50, 0.3, 0, 1,0);% manually adjust the thresholds
[EEG, locthresh, globthresh, nrej] = pop_jointprob( EEG, 1, 1:128, 4, 4, 0, 1, 0);% manually adjust the thresholds
[EEG, locthresh, globthresh, nrej] = pop_rejkurt( EEG, 1, 1:128, 5, 5, 0, 1, 0);% manually adjust the thresholds
[EEG,~] = pop_rejspec( EEG, 1, 'elecrange', 1:128,'threshold', [-50 50],...
    'freqlimits', [0 2], 'eegplotreject', 1);% manually adjust the thresholds
[EEG,~] = pop_rejspec( EEG, 1, 'elecrange', 1:128,'threshold', [-100 25],...
    'freqlimits', [20 40], 'eegplotreject', 1);% manually adjust the thresholds


EEG = pop_saveset( EEG, 'filename','SN_SSVEP_postica_pruned_epoched_rejected.set','filepath',current_path);






