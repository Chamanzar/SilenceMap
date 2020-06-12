function h1 = plot_source_space_signal_vF(x, sulc, cortex, cortex2)
%%
% This function plots the cortical surface as well as the region of silence
% on it.
%
% See: README.txt and [1] for more info.

% [1] A. Chamanzar, M. Behrmann, and P. Grover,
%  "Neural silences can be localized using noninvasive scalp EEG",
%   To be submitted to Nature BME, 2020.


% inputs:

% x: column vector with binary values, 0 indicating the active source and 1
% indicating the silent source

% sulc: Sulci depth in the cortex for plotting purposes 

% cortex2: discretized cortex vertices (can have high/low res)

% cortex: discretized cortex vertices for ploting x on the crortical
% surface (should have exactly the same number of vertices as the dimension
% of x

% output:

% h1: a handle to the graphical object in the figure

% usage example: please see sLoreta_MNE_MUSIC_simulation_vF.m,
% SilenceMap_with_baseline.m, and SilenceMap_without_baseline.m

% Author: Alireza 	Date: 2020/05/21 11:00:25 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details
%%

ind = find(x);
cortex_detected.vertices = cortex.vertices(ind,:);
Face_temp = [];
% Paint each triangle in the tessellated brain if at least on vertex in it
% belongs to the region of silence
for i=1:size(cortex.faces,1)
    if(((ismember(cortex.faces(i,1),ind)) ...
            || (ismember(cortex.faces(i,2),ind))...
            || (ismember(cortex.faces(i,3),ind))))
       Face_temp = cat(1,Face_temp,cortex.faces(i,:));      
    end
end
Ind = reshape(unique(Face_temp),[],1);
cortex_detected.vertices = cortex.vertices(Ind,:);

for i=1:size(Face_temp,1)
    for j=1:3
        Face_temp(i,j) = find(Face_temp(i,j)==Ind)';
    end
end

cortex_detected.faces = Face_temp;

y = zeros(size(Ind));
figure;
% h1 = plot_source_space_signal_vIV(zeros(size(Ind)),sulc,cortex2,cortex_detected);
hold on;

m = 64;  % 64-elements is each colormap
colormap([flip(gray(m)); flip(autumn(m))]);
h1 = trisurf(cortex2.faces, cortex2.vertices(:, 1), cortex2.vertices(:, 2), cortex2.vertices(:, 3), sulc, 'EdgeColor', 'none');
h2 = trisurf(cortex_detected.faces, cortex_detected.vertices(:, 1), cortex_detected.vertices(:, 2), cortex_detected.vertices(:, 3), 'EdgeColor', 'none', 'FaceColor', [1, 0, 0]);

c1 = min(m, round((m - 1) * (sulc - min(sulc)) / (max(sulc) - min(sulc))) + 1);
c2 = m + min(m, round((m - 1) * (y - min(y)) / (max(y) - min(y))) + 1);

set(h1, 'CData', c1);
set(h2, 'CData', c2);
caxis([min([c1; c2]), max([c1; c2])]);
h1.FaceAlpha = 0.5;

shading interp;
axis off
axis equal

end