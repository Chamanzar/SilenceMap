clear all
close all
clc
%% 
% This code extracts the 3D head model from the prerpocessed MRI scans and 
% prepare the headmodel for the leadfield extraction in FieldTrip
%
% See: README.txt and [1] for more info.

% [1] A. Chamanzar, M. Behrmann, and P. Grover,
%  "Neural silences can be localized using noninvasive scalp EEG",
%   To be submitted to Nature BME, 2020.

% Author: Alireza 	Date: 2020/05/19 12:00:48 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details
%% Create a headmodel for preparing a leadfield matrix with fieldtrip

subjects_dir = './OT_defaced';% Please put the preprocessed MRI scan in the current folder
subject = 'OT_defaced';

% Load the cerebral cortex (gray matter-Pial) from the processed MRI scans:
fprintf('Loading the Pial surfaces\n');
[Pial_L.vertices, Pial_L.faces] = read_surf([subjects_dir , '/surf/lh.Pial']);
Pial_L.faces = Pial_L.faces + 1; % fixing the indices for Matlab
[Pial_R.vertices, Pial_R.faces] = read_surf([subjects_dir , '/surf/rh.Pial']);
Pial_R.faces = Pial_R.faces + 1;
Sulc_LH = read_curv([subjects_dir , '/surf/lh.sulc']);
Sulc_RH = read_curv([subjects_dir , '/surf/rh.sulc']);

% % Subsample the surfaces to reduce computational difficulty

% Constructing the whole brain model by concatenation of the LH and RH Pial surfaces
Pial.vertices = [Pial_L.vertices; Pial_R.vertices];
Pial_R.faces = Pial_R.faces + size(Pial_L.vertices,1);
Pial.faces = [Pial_L.faces; Pial_R.faces];
Sulc = [Sulc_LH; Sulc_RH];

cortex = Pial;
cortex1 = reducepatch(cortex,(1/5));
%Create an index of the preserved vertex.
[ind,loc] = ismember(cortex.vertices,cortex1.vertices,'rows'); 
ind = find(ind);
duplicated_faces = [];
%Remove repeated faces after subsampling:
faces_temp = sort(cortex1.faces,2);
for i=1:size(faces_temp,1)
    dface = faces_temp(i,:)-faces_temp;
    dface = sum(abs(dface),2);
    face_ind = find(dface==0);
    face_ind(face_ind<=i) = []; 
    duplicated_faces = [duplicated_faces;face_ind];
end
cortex1.faces(duplicated_faces,:) = [];
Sulc = Sulc(ind);
Pial.faces = cortex1.faces;
Pial.vertices = cortex1.vertices;


% Rotate Surfaces if necessary (visual check):

% phi_y = -12;
% phi = phi_y*pi/180;
% 
% Pial.vertices(:,1) = Pial.vertices(:,1);
% Pial.vertices(:,2) = cos(phi)*Pial.vertices(:,2)-sin(phi)*Pial.vertices(:,3);
% Pial.vertices(:,3) = sin(phi)*Pial.vertices(:,2)+cos(phi)*Pial.vertices(:,3);
% 
% 
% phi_x = -12;
% phi = phi_x*pi/180;
% 
% Pial.vertices(:,2) = Pial.vertices(:,2);
% Pial.vertices(:,1) = cos(phi)*Pial.vertices(:,1)-sin(phi)*Pial.vertices(:,3);
% Pial.vertices(:,3) = sin(phi)*Pial.vertices(:,1)+cos(phi)*Pial.vertices(:,3);
% 
% phi_z = +2;
% phi = phi_z*pi/180;
% 
% Pial.vertices(:,3) = Pial.vertices(:,3);
% Pial.vertices(:,1) = cos(phi)*Pial.vertices(:,1)-sin(phi)*Pial.vertices(:,2);
% Pial.vertices(:,2) = sin(phi)*Pial.vertices(:,1)+cos(phi)*Pial.vertices(:,2);


% Plot the Pial surface:
plot_source_space_signal_vF(ones(size(Pial.vertices,1),1),Sulc,Pial,Pial);
title('Red region = Actual region of silence','FontSize',20)


% Extraction of different layers of the brain:

% Scalp
[bnd(1).pos, bnd(1).tri] = read_surf([subjects_dir, '/bem/watershed/', [subject, '_outer_skin_surface']]);
bnd(1).tri = bnd(1).tri + 1;

% Skull
[bnd(2).pos, bnd(2).tri] = read_surf([subjects_dir, '/bem/watershed/', [subject, '_outer_skull_surface']]);
bnd(2).tri = bnd(2).tri + 1;

% CSF
[bnd(3).pos, bnd(3).tri] = read_surf([subjects_dir, '/bem/watershed/', [subject, '_inner_skull_surface']]);
bnd(3).tri = bnd(3).tri + 1;

% Brain
[bnd(4).pos, bnd(4).tri] = read_surf([subjects_dir, '/bem/watershed/', [subject, '_brain_surface']]);
bnd(4).tri = bnd(4).tri + 1;

% Plot the scalp, skull, CSF and brain surfaces
figure;
trisurf(bnd(1).tri, bnd(1).pos(:, 1), bnd(1).pos(:, 2), bnd(1).pos(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.6, 0.6, 0.6], 'FaceAlpha', 0.1);
hold on
trisurf(bnd(2).tri, bnd(2).pos(:, 1), bnd(2).pos(:, 2), bnd(2).pos(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 0.1);
hold on
trisurf(bnd(3).tri, bnd(3).pos(:, 1), bnd(3).pos(:, 2), bnd(3).pos(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.2, 0.2, 0.2], 'FaceAlpha', 0.1);
hold on
trisurf(bnd(4).tri, bnd(4).pos(:, 1), bnd(4).pos(:, 2), bnd(4).pos(:, 3), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.6, 0.4], 'FaceAlpha', 0.3);
axis equal;

% Construct the BEM model
vol.cfg.checksize = 100000;
vol.cfg.showcallinfo = 'yes';
vol.cfg.tissue = {'scalp', 'skull', 'csf', 'brain'};
vol.cfg.method = 'bem';
vol.cfg.trackconfig = 'off';
vol.cfg.checkconfig = 'loose';
vol.bnd = bnd;
vol.type = 'bemcp';
vol.unit = 'mm';
vol.mat = 1;
vol.cond = [1, 1/15, 5, 1];    % Conductivities


cfg = [];
cfg.channel = {'EEG'};                    
cfg.grid.pos = Pial.vertices;               % use the Gray matter as the source points
cfg.grid.inside = 1:size(Pial.vertices, 1); % all source points are inside of the brain
cfg.vol = vol;                            % BEM model

headmodel = cfg;
Cortex.Pial = Pial;
Cortex.Sulc = Sulc;
save([subject, '_headmodel.mat'], 'headmodel', 'Cortex', '-v7.3');
