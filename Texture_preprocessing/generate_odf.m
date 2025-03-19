% Generate synthetic ODF
%% load MTEX:
addpath '/Users/user/Documents/MATLAB/mtex-5.7.0'
startup_mtex
clear all

%% === load phase Crystal Symmetry: === %
% change to 'Ti64_alpha' or 'Ti64_beta'

% Ti64_alpha conventions
% DAMASK: X||a Y||b*
% Rolling: X||a* Y||b
% Chris uses X||b* Y||a Z||c*

% cs = crystalSymmetry('6/mmm', [3 3 4.7], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Ti-Hex', 'color', [0.53 0.81 0.98]);
% % odf.SS=specimenSymmetry('triclinic');
% hkil = [Miller(0,0,0,2,cs), Miller(1,0,-1,0,cs), Miller(1,1,-2,0,cs)];

% Ti64_beta conventions
cs = crystalSymmetry('m-3m', [3.2 3.2 3.2], 'mineral', 'Titanium cubic', 'color', [0.56 0.74 0.56]);
% odf.SS=specimenSymmetry('orthorombic');
hkil = [Miller(0,0,1,cs), Miller(1,1,0,cs), Miller(1,1,1,cs)];


% Al FCC conventions
% cs = crystalSymmetry('m-3m', [3.2 3.2 3.2], 'mineral', 'Titanium cubic', 'color', [0.56 0.74 0.56]);
% odf.SS=specimenSymmetry('orthorombic');
% hkil = [Miller(0,0,1,cs), Miller(1,1,0,cs), Miller(1,1,1,cs)];


% W BCC conventions
% cs = crystalSymmetry('m-3m', [3.2 3.2 3.2], 'mineral', 'Titanium cubic', 'color', [0.56 0.74 0.56]);
% hkil = [Miller(0,0,1,cs), Miller(1,1,0,cs), Miller(1,1,1,cs)];

%% define orientations...

eulers = [
    [0 0 0],
    ];

% quat = [
%     [1 0 0 0], % {100}
%     [0.888 0.0 0.325 -0.325] % {111}
%     ];

% convert orientations to matlab quaterions...
q = euler2quat(eulers(1,1)*degree,eulers(1,2)*degree,eulers(1,3)*degree,'Bunge'); % converted from bunge
disp(q)
% q = quaternion(transpose(quat)); % defined from quat

%% Estimate an ODF from the orientations...
% calcdensity method
ori = orientation(q, cs);
odf = calcDensity(ori, 'halfwidth', 0.2);
% for IPF colouring
IPFkey = ipfColorKey(odf.CS);
IPFkey.inversePoleFigureDirection = vector3d.(upper('Z')); % define direction for IPF colours
oriColors = IPFkey.orientation2color(ori);

%% plot 3D ODF...
newMtexFigure('figSize', 'large')
plot3d(odf)


%% plot pole figures...
newMtexFigure('figSize', 'large', 'layout', [1,3])
setMTEXpref('FontSize', 32);

% plotting convention for rolling pole figures
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','OutofPlane');

% plot
mrd_max = 3;
plotPDF(odf, hkil,'antipodal', 'contourf', 0:0.25:mrd_max, 'minmax') % plot with contouring
% plotPDF(ori, hkil, 'antipodal', 'property', oriColors, 'minmax'); % plot IPF color points

% Nicer annotations
text(vector3d.X,'X','VerticalAlignment','bottom'); % moving the vector3d axis labels outside of the hemisphere boundary
text(vector3d.Y,'Y','HorizontalAlignment','left');
f = gcm; % removes min label and moves max label
for i = 1:length(f.children)
    f.children(i).Title.Position=[1,1.25,1]; % use [-1,1.25,1] to make layout similar to Channel 5 and turn off minmax if you want to
    ax = getappdata(f.children(i),'sphericalPlot');
    ax.TL.Position=[0, 1.1, 1];           % max label posn
    ax.BL.Visible='off';                  % get rid of min
end

% mtexColorbar ('location', 'southoutside', 'title', 'mrd', 'FontSize', 30);

%% save odf...

save('ODF_0002ND1010RD.mat', 'odf');