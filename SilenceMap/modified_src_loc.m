function [XMNE, XLORETA,  XMUSIC, I_MNE, I_sLoreta, I_MUSIC] = modified_src_loc(LV,L,cortex,M_Silence_noisy...
    ,XMNE_initial, XLoreta_initial, XMUSIC_initial)

%%
% Modified sLORETA, MNE, and MUSIC algorithms
% This code implements the modified sLORETA, MNE, and MUSIC algorithms for 
% silence localization.
%
% See: README.txt and [1] for more info.

% [1] A. Chamanzar, M. Behrmann, and P. Grover,
%  "Neural silences can be localized using noninvasive scalp EEG",
%   To be submitted to Nature BME, 2020.


% inputs:

% LV: the level in the successive refinement approach. LV=1 for the low-res
% and LV=0 for the high-res grid.

% L: leadfield matrix

% cortex: discretized cortex vertices and faces

% M_Silence_noisy: scalp noisy signals

% XMNE_initial: COM of the localized region using modified MNE in the low-res. Leave it empty
% when using this function for the low-res grid.

% XLoreta_initial: COM of the localized region using modified sLORETA in the low-res. Leave it empty
% when using this function for the low-res grid.

% XMUSIC_initial: COM of the localized region using modified MUSIC in the low-res. Leave it empty
% when using this function for the low-res grid.

% output:

% XMNE: 3D coordinates of the sources in the localized region of silence
% using modified MNE

% XLORETA: 3D coordinates of the sources in the localized region of silence
% using modified sLORETA

% XMUSIC: 3D coordinates of the sources in the localized region of silence
% using modified MUSIC

% I_MNE: indices of the sources in the localized region of silence
% using modified MNE

% I_sLoreta: indices of the sources in the localized region of silence
% using modified sLORETA

% I_MUSIC: indices of the sources in the localized region of silence
% using modified MUSIC


% usage example: please see sLoreta_MNE_MUSIC_simulation_vF.m 

% Author: Alireza 	Date: 2020/05/19 11:48:00 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details


%% MNE & L-Curve:

lambdaa = logspace(-2,2,40);
Residual = [];
Solutn_norm = [];
b = M_Silence_noisy;
for lamb=1:size(lambdaa,2)
    S = (L'/(L*L'+lambdaa(lamb)*eye(size(L,1))))*b;
    Residual = cat(1,Residual,norm((L*S-b),'fro'));
    Solutn_norm = cat(1,Solutn_norm,norm(S,'fro'));
end

curvatur = Curvature_vI(Residual,Solutn_norm);
[~,I_lambda] = max(curvatur);

% figure;loglog(Residual(1:end),Solutn_norm(1:end));
% figure;plot((Residual(1:end)),(Solutn_norm(1:end)));
% hold on
% plot((Residual(I_lambda)),(Solutn_norm(I_lambda)),'r*');


S = (L'/(L*L'+lambdaa(I_lambda)*eye(size(L,1))))*b;
[~,Is] = sort(S.^2,1);

% S_var = var(S,0,2);

%Initial guess on k_MNE:
k_MNE = floor(size(L,2)/10);

Is = Is(1:k_MNE,:);
votes = [];
for i=1:size(L,2)
   votes = cat(2,votes,size(find(Is==i),1));
end
[votes_sorted,Iss] = sort(votes,'descend');

%Estimate k based on distance to origin in \beta_{q}-index curve:
if (LV==0)
    f_dist = (votes_sorted/max(votes_sorted)).^2+((1:size(L,2))/size(L,2)).^2;
    [~,k_hat] = min(f_dist);
    k_MNE = k_hat;
else
    k_MNE = 1;
end


I_MNE = Iss(1:min(2*k_MNE,size(L,2)));

% MNE exploiting the knowledge of contiguous region:
if (~isempty(XMNE_initial))
    Silence_center = XMNE_initial;
else
    Silence_center = mean(cortex(I_MNE,:),1);
end
dist_temp = sum((cortex-repmat(Silence_center,size(cortex,1),1)).^2,2);
[~,indx_new] = sort(dist_temp);
indx_new = indx_new(1:min(k_MNE,size(L,2)));
I_MNE = indx_new;

%% sLORETA & L-Curve:

lambdaa = logspace(-2,2,40);
Residual = [];
Solutn_norm = [];
b = M_Silence_noisy;
for lamb=1:size(lambdaa,2)
    S = (L'/(L*L'+lambdaa(lamb)*eye(size(L,1))))*b;
    Residual = cat(1,Residual,norm((L*S-b),'fro'));
    Solutn_norm = cat(1,Solutn_norm,norm(S,'fro'));
end

curvatur = Curvature_vI(Residual,Solutn_norm);
[~,I_lambda] = max(curvatur);

% figure;loglog(Residual(1:end),Solutn_norm(1:end));
% figure;plot((Residual(1:end)),(Solutn_norm(1:end)));
% hold on
% plot((Residual(I_lambda)),(Solutn_norm(I_lambda)),'r*');


S = (L'/(L*L'+lambdaa(I_lambda)*eye(size(L,1))))*b;
Var_S = diag((L'/(L*L'+lambdaa(I_lambda)*eye(size(L,1))))*L);

S = (S.^2)./Var_S;

[~,Is] = sort(S,1);

% S_var = var(S,0,2);

%Initial guess on k_MNE:
k_sLoreta = floor(size(L,2)/10);

Is = Is(1:k_sLoreta,:);
votes = [];
for i=1:size(L,2)
   votes = cat(2,votes,size(find(Is==i),1));
end
[votes_sorted,Iss] = sort(votes,'descend');

%Estimate k based on distance to origin in \beta_{q}-index curve:
if (LV==0)
    f_dist = (votes_sorted/max(votes_sorted)).^2+((1:size(L,2))/size(L,2)).^2;
    [~,k_hat] = min(f_dist);
    k_sLoreta = k_hat;
else
    k_sLoreta = 1;
end


I_sLoreta = Iss(1:min(2*k_sLoreta,size(L,2)));

% sLoreta exploiting the knowledge of contiguous region:
if (~isempty(XLoreta_initial))
    Silence_center = XLoreta_initial;
else
    Silence_center = mean(cortex(I_sLoreta,:),1);
end
dist_temp = sum((cortex-repmat(Silence_center,size(cortex,1),1)).^2,2);
[~,indx_new] = sort(dist_temp);
indx_new = indx_new(1:min(k_sLoreta,size(L,2)));
I_sLoreta = indx_new;

%% MUSIC:

% [U,S,V] = svd(M_Silence_noisy(:,1:1000));
[U,S,V] = svd(M_Silence_noisy);
S = diag(S);
S_tot = sum(S);
for i=2:size(S,1)
    S(i) = S(i) + S(i-1);
end
S = S/S_tot;
I_S = find(S>0.99);
U_s = U(:,1:I_S(1));


P_s = eye(size(U,1))-(U_s*U_s');

J = zeros(size(L,2),1);
for i=1:size(L,2)
    J(i) = norm(P_s*L(:,i))/norm(L(:,i));
end

[J_sorted,I_J] = sort(J,'descend');
% I_J = I_J(1:k);


%Estimate k based on distance to origin in \beta_{q}-index curve:
if (LV==0)
    f_dist = (J_sorted/max(J_sorted)).^2+((1:size(L,2))'/size(L,2)).^2;
    [~,k_hat] = min(f_dist);
    k_MUSIC = k_hat;
else
    k_MUSIC = 1;
end
I_MUSIC = I_J(1:min(2*k_MUSIC,size(L,2)));


% MUSIC exploiting the knowledge of contiguous region:
if (~isempty(XMUSIC_initial))
    Silence_center = XMUSIC_initial;
else
    Silence_center = mean(cortex(I_MUSIC,:),1);
end
dist_temp = sum((cortex-repmat(Silence_center,size(cortex,1),1)).^2,2);
[~,indx_new] = sort(dist_temp);
indx_new = indx_new(1:min(k_MUSIC,size(L,2)));
I_MUSIC = indx_new;

%% Output of detector as location of silence voxels:

%MNE:
XMNE = cortex(I_MNE,:);
%sLoreta:
XLORETA = cortex(I_sLoreta,:);
%MUSIC:
XMUSIC = cortex(I_MUSIC,:);

end