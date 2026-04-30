function [datStruct] =...
    rmBadChans(datStruct, NVA)
% Written by Giles Blaney Winter 2022
% Comments by Cristianne Fernandez August 2023
% This function removes sets that are above the set noise thresholds, determined 
% by channels with a Total Hemoglobin noise above the set threshold and
% frequency

%%%%%%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%
% datStruct - either SD, SS, or DS struct that has gone through
%           parseArrayData.m 
% OPTIONAL
% BLinds - Logical vector where true indicies will be considered baseline
% Tthresh - (microM) Threshold for the noise in Total Hemoglobin
% window - (sec) window size 
% fHP - (Hz) highpass filter frequency 
% L_C - path length for DC
% L_I - path length for AC 
% L_P - path length for phase
% DSF_C - for DC
% DSF_I - for AC
% DSF_P - for phase
%%%%%%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%
% datStruct - either SD, SS, or DS struct with missing 

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
    Text=sum(makeE('OD', lambda), 2); %1/(mm uM)
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
                eval(['noiseThresh=' mTyp 'noiseThresh;']);

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
                    catch
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