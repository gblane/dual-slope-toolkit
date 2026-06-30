function aux = calibrateAUX(data, varargin, nameValArgs)
% calibrateAUX Calibrate auxiliary channels from ISS/Imagent data.
%
% aux = calibrateAUX(data, Name=Value)
%
% Written by Cristianne Fernandez (Fall 2022; comments August 2023).
% Moved from DOIT-Toolbox into dual-slope-toolkit (2026).
%
% Inputs:
%   data - Struct containing AUX, timemat, and fs fields.
%
% Name-Value Inputs:
%   label   - String labels for AUX channels.
%   UnCal   - If true, apply calibration factors. If false, pass raw values.
%   *_cal   - Two-element calibration vectors [scale; offset].
%   Dthresh - Threshold for digital auxiliary channels.
%
% Outputs:
%   aux - Struct containing calibrated auxiliary channels plus t, fs, Labels.

arguments
    data struct
end
arguments (Repeating)
    varargin
end
arguments
    nameValArgs.label (:,1) string = [...
        "ABP"; "Cuff"; "Resp"; "HR"; ...
        "NotUsed"; "D3"; "D5"; "D13"]
    nameValArgs.UnCal (1,1) logical = true
    nameValArgs.ABP_cal (:,1) double = [-1.195e-1; -4.13239e1]
    nameValArgs.Cuff_cal (:,1) double = [6.1e-2; 7.743e-1]
    nameValArgs.Resp_cal (:,1) double = [1e1; 0]
    nameValArgs.HR_cal (:,1) double = [6.22e-2; 2.979e-1]
    nameValArgs.Dthresh (1,1) double = 3500
end

label = nameValArgs.label;
if numel(label) < size(data.AUX, 2)
    error('calibrateAUX:TooFewLabels', ...
        'Provide at least one label for each AUX channel.');
end

for ii = 1:size(data.AUX, 2)
    fieldName = label(ii);
    if nameValArgs.UnCal
        switch lower(fieldName)
            case "abp"
                aux.(fieldName)(:, 1) = ...
                    data.AUX(:, ii)*nameValArgs.ABP_cal(1) + ...
                    nameValArgs.ABP_cal(2);
            case "cuff"
                aux.(fieldName)(:, 1) = ...
                    data.AUX(:, ii)*nameValArgs.Cuff_cal(1) + ...
                    nameValArgs.Cuff_cal(2);
            case "resp"
                aux.(fieldName)(:, 1) = ...
                    data.AUX(:, ii)*nameValArgs.Resp_cal(1) + ...
                    nameValArgs.Resp_cal(2);
            case "hr"
                aux.(fieldName)(:, 1) = ...
                    data.AUX(:, ii)*nameValArgs.HR_cal(1) + ...
                    nameValArgs.HR_cal(2);
            case {"d3", "d5", "d13", "notused"}
                aux.(fieldName)(:, 1) = data.AUX(:, ii) >= nameValArgs.Dthresh;
            otherwise
                warning('calibrateAUX:UnknownLabel', ...
                    'Unknown AUX label %s; leaving channel uncalibrated.', ...
                    fieldName);
                aux.(fieldName)(:, 1) = data.AUX(:, ii);
        end
    else
        aux.(fieldName)(:, 1) = data.AUX(:, ii);
    end
end

aux.t = data.timemat;
aux.fs = data.fs;
aux.Labels = label;
end
