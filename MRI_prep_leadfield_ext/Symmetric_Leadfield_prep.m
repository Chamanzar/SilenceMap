close all
clear all
clc
%% 
% This code construct the symmetric brain models and prepare the headmodel
% for leadfield extraction using FieldTrip. We use the headmodel and the 
% sensor locations which are saved in Project_and_setup_sensors.m
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

load('OT_defaced_headmodel.mat');

%% %%%%%%%Mirroring the intact brain hemisphere to contruct the symmetric brain model%%%%%%%%

% Subsampling the Cortex surface:
cortex = Cortex.Pial;
cortex1 = reducepatch(cortex,(1/75));
%Create an index of the preserved vertex.
[ind,loc] = ismember(cortex.vertices,cortex1.vertices,'rows'); 
ind = find(ind);
Sulc = Cortex.Sulc;
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

Cortex.Sulc = Sulc(ind);
Cortex.Pial.faces = cortex1.faces;
Cortex.Pial.vertices = cortex1.vertices;

cortexR.vertices = Cortex.Pial.vertices(floor(size(Cortex.Pial.vertices,1)*0.5)-50:end,:);
sulc_R = Cortex.Sulc(floor(size(Cortex.Pial.vertices,1)*0.5)-50:end,:);

%Pruning the redundant faces:
ind_temp = ones(size(cortex1.faces,1),1);
for i=1:size(cortex1.faces,1)
    for j=1:size(cortex1.faces,2)
        if(cortex1.faces(i,j)<=floor(size(Cortex.Pial.vertices,1)*0.5)-51)
            ind_temp(i) = 0;
        end
    end
end
cortexR.faces = Cortex.Pial.faces((ind_temp==1),:);
cortexR.faces = cortexR.faces - floor(size(Cortex.Pial.vertices,1)*0.5)+51;

%Pruning the isolated vertices:
ind_temp = ones(size(cortexR.vertices,1),1);
vert_in_cortexR = unique(cortexR.faces);
for i=size(cortexR.vertices,1):-1:1
    if(size(find(i==vert_in_cortexR),1)==0)
       i
       ind_temp(i) = 0;
       cortexR.faces(cortexR.faces>i) = cortexR.faces(cortexR.faces>i)-1;
    end
end

cortexR.vertices = cortexR.vertices((ind_temp==1),:);
sulc_R = sulc_R((ind_temp==1));

cortexLR.faces = cortexR.faces;
cortexLR.vertices = cortexR.vertices;
sulc_LR = sulc_R;
cortexLR.vertices(:,1) = -cortexLR.vertices(:,1);

figure;
plot_source_space_signal_vIV(zeros(size(sulc_R)),sulc_R,cortexR,cortexR);
hold on
plot_source_space_signal_vIV(zeros(size(sulc_LR)),sulc_LR,cortexLR,cortexLR);
axis off
axis equal

%Brain location adjustments:

cortexLR.vertices(:,1) = cortexLR.vertices(:,1)+2;
cortexR.vertices(:,1) = cortexR.vertices(:,1)+2;


headmodel.grid.pos = cat(1,cortexLR.vertices,cortexR.vertices); % use the Gray matter as the source points
headmodel.grid.inside = 1:size(headmodel.grid.pos,1); % some source points are inside of the brain


% Construct the symmetric brain by concatination of the mirrored hemisphere and 
% the intact hemisphere:
Cortex.Pial.faces = cat(1,cortexLR.faces,cortexR.faces+size(cortexLR.vertices,1));
Cortex.Pial.vertices = cat(1,cortexLR.vertices,cortexR.vertices);
Cortex.Sulc = cat(1,sulc_LR,sulc_R);

cortex.faces = Cortex.Pial.faces;
cortex.vertices = Cortex.Pial.vertices;

figure;
plot_source_space_signal_vIV(zeros(size(Cortex.Sulc)),Cortex.Sulc,cortex,cortex);
axis off
axis equal


headmodel_file = sprintf('OT_headmodel_symmetric_%d-%d.mat', size(headmodel.grid.pos,1), size(headmodel.elec.chanpos,1))
save(headmodel_file,'headmodel','Cortex')


