close all
clear all
clc
%% 
% This code extracts the leadfield matrix L using FieldTrip. We use the symmetric 
% brain model contructed in Symmetric_Leadfield_prep.m
%
% See: README.txt and [1] for more info.

% [1] A. Chamanzar, M. Behrmann, and P. Grover,
%  "Neural silences can be localized using noninvasive scalp EEG",
%   To be submitted to Nature BME, 2020.

% Author: Alireza 	Date: 2020/05/19 12:00:48 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details
%%
% Set up the FieldTrip software (make sure to add the FieldTrip software to the matlab path)
ft_defaults;

load('OT_headmodel_symmetric_818-128.mat');
    
%% %%%%%%Leadfield BEM calculation based on the source and sensor grid in the headmodel%%%%%

lf = ft_prepare_leadfield(headmodel); %Let FieldTrip take care of it!

n_elec = size(lf.leadfield{1}, 1);
n_dipoles = size(lf.leadfield, 2);

% Leadfield matrix
L = zeros(n_elec, n_dipoles);


% Compute normals to the cortical surface
pial_tri = triangulation(Cortex.Pial.faces(1:0.5*size(Cortex.Pial.faces,1),:), ...
                      Cortex.Pial.vertices(1:0.5*size(Cortex.Pial.vertices,1),:));
normals_L = vertexNormal(pial_tri);
normals_R = normals_L;
normals_R(:,1) = -normals_R(:,1);
vertices_norm = -[normals_L;normals_R];

%Calculating norm using fieldtrip function:
% [vertices_norm] = normals(Cortex.Pial.vertices, Cortex.Pial.faces, 'vertices');

% Ploting the calculated normals:
figure;
plot_source_space_signal_vIV(zeros(size(Cortex.Sulc)),Cortex.Sulc,Cortex.Pial,Cortex.Pial);
axis off
axis equal
hold on
quiver3(Cortex.Pial.vertices(:,1),Cortex.Pial.vertices(:,2),Cortex.Pial.vertices(:,3), ...
     vertices_norm(:,1),vertices_norm(:,2),vertices_norm(:,3),1,'Color','b');


% @ each source, project the calculated leadfield vector along the normal
% vector:
for j = 1:n_dipoles
	% Compute radial component of sensor values for the j'th unit dipole
	fwd_j = lf.leadfield{j};
	fwd_j = fwd_j * vertices_norm(j, :)';
	L(:, j) = fwd_j;
end

% Double check and zero any possible NaN values in L:
if any(isnan(L(:)))
	L(isnan(L)) = 0;
end
sensor_locs = headmodel.elec.chanpos;
% Save the normal leadfield matrix
leadfield_file = sprintf('OT_leadfield_symmetric_%d-%d.mat', n_dipoles, n_elec)
save(leadfield_file, 'L', 'sensor_locs', '-v7.3');