% This code projects the electrode locations on the scalp surface of the headmodel 
% obtained from the MRI_headmodel.m and prepare the electrode info for the leadfield 
% extraction in FieldTrip
%
% See: README.txt and [1] for more info.

% [1] A. Chamanzar, M. Behrmann, and P. Grover,
%  "Neural silences can be localized using noninvasive scalp EEG",
%   To be submitted to Nature BME, 2020.

% Author: Alireza 	Date: 2020/05/20 12:08:11 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details
%%
% Set up the FieldTrip software (make sure to add the FieldTrip software to the matlab path)
ft_defaults;

% Load the headmodel
load('OT_defaced_headmodel.mat');
% Load the sensor locations based on 10-5 standard locations
load('Biosemi_128_ABC_standard.mat');

scalp_mesh = headmodel.vol.bnd(1);

% Stop here:
% You need to manually adjust the radius of the electrode grid to be located
% outside the scalp mesh:
max_radius = max(sqrt(sum(scalp_mesh.pos.^2, 2)));
radius_adjust = max_radius-1.5;
grid_elec = proj_elecs * (max_radius - radius_adjust);

% Rotate the electrode grid to match the scalp orientation:

phi_z = +90;
phi = phi_z*pi/180;

elec_grid_old = grid_elec;

grid_elec(:,3) = elec_grid_old(:,3);
grid_elec(:,1) = cos(phi)*elec_grid_old(:,1)-sin(phi)*elec_grid_old(:,2);
grid_elec(:,2) = sin(phi)*elec_grid_old(:,1)+cos(phi)*elec_grid_old(:,2);


phi_x = +10;
phi = phi_x*pi/180;

elec_grid_old = grid_elec;

grid_elec(:,1) = elec_grid_old(:,1);
grid_elec(:,2) = cos(phi)*elec_grid_old(:,2)-sin(phi)*elec_grid_old(:,3);
grid_elec(:,3) = sin(phi)*elec_grid_old(:,2)+cos(phi)*elec_grid_old(:,3);


% Change the position of the electrode grid if necessary:
pos_adjust = [0, -20, +5];
grid_elec = grid_elec + repmat(pos_adjust,size(grid_elec,1),1);


% plot the electrodes and the scalp mesh before projection:
figure;
hold on;
axis equal;

% Plot the scalp mesh
trisurf(scalp_mesh.tri, scalp_mesh.pos(:, 1), scalp_mesh.pos(:, 2), ...
        scalp_mesh.pos(:, 3), 'EdgeColor', [1, 1, 1], 'EdgeAlpha', 0.3, ...
		'FaceColor', [0.6, 0.6, 0.6], 'FaceAlpha', 0.5);

% Plot the electrode grid
scatter3(grid_elec(:, 1), grid_elec(:, 2), grid_elec(:, 3), 18, 'b', 'filled');

%% %%%%%%%%%%%%%Project the adjusted electrode grid on the scalp%%%%%%%%%%%%

n_elec = 128;

% Setup the elec structure as required by fieldtrip
elec.chanunit = cell(n_elec, 1);
elec.chanunit(:) = {'V'};    % Set all elements of cell array to V
elec.chanpos = grid_elec(1:n_elec,:);
elec.elecpos = grid_elec(1:n_elec,:);
elec.chantype = cell(n_elec, 1);
elec.chantype(:) = {'eeg'};  % Set all elements of cell array to EEG
elec.label = cell(n_elec, 1);
for i = 1:n_elec
	elec.label{i} = sprintf('Biosemi%d', i);
end
elec.type = sprintf('Biosemi_standard%d', n_elec);
elec.unit = 'mm';


% interactively coregister the electrodes to the BEM head model
% this is a visual check and refinement step
cfg = [];
cfg.method    = 'interactive';
cfg.elec      = elec;
cfg.headshape = scalp_mesh;
elec = ft_electroderealign(cfg);

% projection of the refined electrode grid onto the scalp surface
cfg = [];
cfg.method    = 'project';
cfg.elec      = elec;
cfg.headshape = scalp_mesh;
elec_proj = ft_electroderealign(cfg);


% Place the sensor configuration back into the headmodel
headmodel.elec = elec_proj;

% Plot to verify
plot_sensor_space_signal(randn(size(elec_proj.chanpos(:,:), 1), 1), ...
                         headmodel.vol.bnd, elec_proj.chanpos(:,:));

save('OT_defaced_headmodel.mat', 'headmodel', 'Cortex', ...
     '-v7.3');





