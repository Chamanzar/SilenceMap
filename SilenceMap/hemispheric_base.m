function [Betta_new] = hemispheric_base(Betta,source_loc,Gp)
%%
% This code calculates the source contribution measure (Beta_new) for 
% silence localization with baseline.
%
% See: README.txt and [1] for more info.

% [1] A. Chamanzar, M. Behrmann, and P. Grover,
%  "Neural silences can be localized using noninvasive scalp EEG",
%   To be submitted to Nature BME, 2020.


% inputs:

% Betta: source contribution measure without baseline 
% source_loc: 3D locations of sources in the discretized brain model
% Gp: the width of the gap between two hemispheres (z^{gap} in mm)

% output:

% Betta_new: source contribution measure with baseline

% usage example: please see SilenceMap_with_baseline.m 

% Author: Alireza 	Date: 2020/05/21 11:30:25 	Revision: 0.1
% Copyright: This project is licensed under the MIT License - see the LICENSE.md file for details
%%
R_ind = find(source_loc(:,1)>Gp);
L_ind = find(source_loc(:,1)<-Gp);
voxel_loc_L = source_loc(L_ind,:);
voxel_loc_R = source_loc(R_ind,:);
voxel_loc_R(:,1) = -voxel_loc_R(:,1);
voxel_R_ind = [];
for i=1:size(voxel_loc_L,1)
    [~,ind_temp] = min(sum((voxel_loc_R-voxel_loc_L(i,:)).^2,2),[],1);
    voxel_R_ind = cat(1,voxel_R_ind,ind_temp);
end

sigma_mu_L = (Betta(L_ind));
sigma_mu_R = (Betta(R_ind));
delta_sigma_mu = (sigma_mu_L)./sigma_mu_R(voxel_R_ind);
Betta_L = delta_sigma_mu;
Betta_L(Betta_L>1) = 1;
Betta_R = delta_sigma_mu;
Betta_R(Betta_R<1) = 1;
Betta_R(Betta_R>1) = Betta_R(Betta_R>1).^(-1);

Betta_new = ones(size(Betta,1),1);
Betta_new(L_ind) = Betta_L;
Betta_new(R_ind(voxel_R_ind)) = Betta_R;
    
end

