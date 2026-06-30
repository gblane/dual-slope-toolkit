function [datStruct] = rmBadChans(datStruct, NVA)
% rmBadChans Removes data channels/sets that exceed noise thresholds.
%
% [datStruct] = rmBadChans(datStruct, NVA)
%
% Written by Giles Blaney (Winter 2022; Ph.D. awarded May 2022)
% Modified by Cristianne Fernandez (August 2023)
%
% This function identifies and masks (sets to NaN) measurement channels 
% with total hemoglobin (HbT) noise levels exceeding a specified threshold.
%
% Inputs:
%   datStruct - Struct containing NIRS data (SD, SS, or DS format).
%
% Optional Name-Value-Arguments (NVA):
%   BLinds   - (t x 1) Logical vector defining baseline indices (Default: all true)
%   Tthresh  - Noise threshold for Total Hemoglobin [uM] (Default: 1)
%   window   - Window size for moving standard deviation [sec] (Default: 10)
%   fHP      - High-pass filter cut-off frequency [Hz] (Default: 100/60)
%   L_C, L_I, L_P - Differential Pathlength Factors for SD data types.
%   DSF_C, DSF_I, DSF_P - Differential Slope Factors for DS/SS data types.
%   SD, armt - Optional structs for advanced cross-validation between SD and DS.
%
% Outputs:
%   datStruct - Updated data struct with bad channels masked and 'useSet' flags.
%
% Shared-repo dependencies:
%   makeE is provided by ../dos-inverse-models.
%   circ_mean is provided by ../my-matlab.

    %% Parse Input
    arguments
        datStruct struct;

        NVA.BLinds (:,1) logical = true(size(datStruct.I, 1), 1);

        NVA.Tthresh double = 1; %uM
        NVA.window double = 10; %sec
        NVA.fHP double = 100/60; %Hz

        % From https://doi.org/10.1002/jbio.201960018 appendix
        NVA.L_C double = 31*7.50; %mm
        NVA.L_I double = 31*7.17; %mm
        NVA.L_P double = 31*1.22; %mm
        NVA.DSF_C double = 8.26;
        NVA.DSF_I double = 7.90;
        NVA.DSF_P double = 1.50;

        NVA.SD struct = [];
        NVA.armt struct = [];
    end
    
    if ~isempty(NVA.armt)
        SDinds=findDS_SDinds(NVA.armt.DSprs,...
            NVA.armt.SDprs);
    end

    lambda=datStruct.lambda;
    fs=datStruct.fs;
    
    datTyp=datStruct.typ;
    nNoise=NVA.window*fs;
    
    %% Calculate Noise Thresholds
    Text=sum(makeE('OD', lambda), 2); % ../dos-inverse-models; 1/(mm uM)
    muaNoiseThresh=NVA.Tthresh*Text; %1/mm
    
    fNy=fs/2;
    if NVA.fHP>=fNy
        error('fHP above Nyquist f');
    end
    noiseScl=fNy/(fNy-NVA.fHP);

    if strcmp(datTyp, 'SD') 
        CnoiseThresh=muaNoiseThresh*NVA.L_C;
        InoiseThresh=muaNoiseThresh*NVA.L_I;
        PnoiseThresh=muaNoiseThresh*NVA.L_P;
    elseif strcmp(datTyp, 'SS') || strcmp(datTyp, 'DS')
        CnoiseThresh=muaNoiseThresh*NVA.DSF_C;
        InoiseThresh=muaNoiseThresh*NVA.DSF_I;
        PnoiseThresh=muaNoiseThresh*NVA.DSF_P;
    end

    %% Loop Through Sets
    datStruct.useSet_C=true(size(datStruct.I, 2), length(lambda));
    datStruct.useSet_I=true(size(datStruct.I, 2), length(lambda));
    datStruct.useSet_P=true(size(datStruct.I, 2), length(lambda));
    for i=1:size(datStruct.I, 2)
        for Lind=1:length(lambda)
            for mTyp=['C', 'I', 'P']
                switch mTyp
                    case 'C'
                        noiseThresh=CnoiseThresh;
                    case 'I'
                        noiseThresh=InoiseThresh;
                    case 'P'
                        noiseThresh=PnoiseThresh;
                end

                SDoverRide=false;
                if ~isempty(NVA.SD) && strcmp(datTyp, 'DS')
                    if all(NVA.SD.(['useSet_' mTyp])(SDinds(i, :), Lind))
                        SDoverRide=true;
%                         noiseThresh(Lind)=noiseThresh(Lind)*2;
                    end
                end
                
                tmp=datStruct.(mTyp)(:, i, Lind);
                if NVA.fHP==0
                    if strcmp(datTyp, 'SD') && strcmp(mTyp, 'P')
                        dtmp=detrend(wrapToPi(tmp(NVA.BLinds)-...
                            circ_mean(tmp(NVA.BLinds))));
                    else
                        dtmp=detrend(tmp(NVA.BLinds));
                    end

                    datStruct.([mTyp 'noise'])(i, Lind)=median(movstd(dtmp,...
                        nNoise))*noiseScl;
                else
                    if strcmp(datTyp, 'SD') && strcmp(mTyp, 'P')
                        dtmp=detrend(wrapToPi(tmp-...
                            circ_mean(tmp(NVA.BLinds))));
                    else
                        dtmp=detrend(tmp);
                    end
                    try
                        tmpHP=highpass(dtmp, NVA.fHP, fs);
    
                        datStruct.([mTyp 'noise'])(i, Lind)=median(movstd(...
                            tmpHP(NVA.BLinds), nNoise))*noiseScl;
                    catch ME
                        warning('rmBadChans:HighpassFailed', ...
                            'High-pass filtering failed for %s set %d lambda %d: %s', ...
                            mTyp, i, Lind, ME.message);
                        datStruct.([mTyp 'noise'])(i, Lind)=Inf;
                    end
                end

                if datStruct.([mTyp 'noise'])(i, Lind)>noiseThresh(Lind) && ...
                        ~SDoverRide
%                 if datStruct.([mTyp 'noise'])(i, Lind)>noiseThresh(Lind)
                    tmp=NaN(size(tmp));
                    datStruct.(['useSet_' mTyp])(i, Lind)=false;
                end
                datStruct.(mTyp)(:, i, Lind)=tmp;
            end
        end
    end
end
