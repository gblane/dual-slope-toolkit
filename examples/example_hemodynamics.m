%% Example: Dual-Slope Hemodynamics Calculation
% This script demonstrates how to calculate absorption changes and 
% hemoglobin concentrations using the standalone Dual-Slope (DS) functions.
%
% Based on:
%   Blaney, G, Sassaroli, A, Pham, T, Fernandez, C, Fantini, S. Phase
%   dual-slopes in frequency-domain near-infrared spectroscopy for enhanced
%   sensitivity to brain tissue: First applications to human subjects. J.
%   Biophotonics. 2019;e201960018. https://doi.org/10.1002/jbio.201960018

clear; home;

%% Load Example Data
% exampleData.mat contains raw intensity (I) and phase (phi) data
load('exampleData.mat');

%% Calculate Absorption Changes (dmua) with DS
% Define distances and baseline indices
opts.rho = [25, 35]; % [Short, Long] distance in mm
opts.blInds = 1:100; % Indices for baseline (rest) period

% Calculate dmua for Intensity and Phase Dual-Slopes
% (Assuming 2 wavelengths: 690nm and 830nm)
dmua_DSI = [];
dmua_DSphi = [];

lambda = [690, 830];

for lInd = 1:length(lambda)
    % Intensity DS (dmua in 1/mm)
    % Y_I should be [Short1, Long1, Short2, Long2]
    Y_I = [A_I(:, (lInd-1)*2+1), B_I(:, (lInd-1)*2+1), ...
           B_I(:, (lInd-1)*2+2), A_I(:, (lInd-1)*2+2)];
    dmua_DSI(:, lInd) = DSdmua(Y_I, 'intensity', opts);
    
    % Phase DS (dmua in 1/mm)
    % Y_phi should be in radians
    Y_phi = [A_phi(:, (lInd-1)*2+1), B_phi(:, (lInd-1)*2+1), ...
             B_phi(:, (lInd-1)*2+2), A_phi(:, (lInd-1)*2+2)] * pi/180;
    dmua_DSphi(:, lInd) = DSdmua(Y_phi, 'phase', opts);
end

%% Convert dmua to Hemodynamics (HbO, HbR)
% mua2OandD expects mua in 1/cm
[HbO_I, HbR_I] = mua2OandD(dmua_DSI * 10, lambda);
[HbO_P, HbR_P] = mua2OandD(dmua_DSphi * 10, lambda);

%% Plot Results
figure(1);
subplot(2,1,1);
plot(t, HbO_I, 'r', t, HbR_I, 'b');
title('Hemodynamics from Intensity Dual-Slope');
xlabel('Time (s)'); ylabel('\Delta Concentration (\muM)');
legend('HbO', 'HbR');

subplot(2,1,2);
plot(t, HbO_P, 'r', t, HbR_P, 'b');
title('Hemodynamics from Phase Dual-Slope');
xlabel('Time (s)'); ylabel('\Delta Concentration (\muM)');
legend('HbO', 'HbR');
