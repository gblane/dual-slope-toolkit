%% Example: Dual-Slope Array Discovery
% This script demonstrates how to use DSdisc.m to identify valid 
% Single-Distance (SD), Single-Slope (SS), and Dual-Slope (DS) pairings
% for a given optode arrangement.
%
% Based on:
%   Blaney, G., Sassaroli, A., & Fantini, S. (2020). Design of a 
%   source-detector array for dual-slope diffuse optical imaging. 
%   Review of Scientific Instruments. https://doi.org/10.1063/5.0015512

clear; home;

%% Load Optode Coordinates
% Coors_HEXarray.mat contains AllSrcs and AllDets matrices [mm]
load('Coors_HEXarray.mat');

coors.AllSrcs = AllSrcs;
coors.AllDets = AllDets;
coors.name = 'Hexagonal_Array';

%% Define Discovery Parameters
% rRng: Range for SD distances [mm]
% lRng: Range for midpoints [mm]
% lTol: Tolerance for matching distances [mm]
nva.rRng = [0, 45]; 
nva.lRng = [5, 25];
nva.lTol = 1;
nva.plot = 'summary'; % Options: 'none', 'summary', 'full'

%% Run Array Discovery
% DSdisc identifying SD, SS, and DS pairings.
[armt, fh] = DSdisc(coors, nva);

%% Display Results
fprintf('Discovered:\n');
fprintf('  - %d SD pairs\n', size(armt.SDprs, 1));
fprintf('  - %d SS pairs\n', size(armt.SSprs, 3));
fprintf('  - %d DS sets\n', size(armt.DSprs, 4));

% The 'armt' struct can now be used with parseArrayData.m and arraySensMap.m
