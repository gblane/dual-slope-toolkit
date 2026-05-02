function [SD, SS, DS] =...
    calAndAddDatatypes(SD, SS, DS, armt, absp, NVA)
% calAndAddDatatypes Calibrates raw data and calculates processed NIRS data types.
%
% [SD, SS, DS] = calAndAddDatatypes(SD, SS, DS, armt, absp, NVA)
%
% Written by Giles Blaney, Ph.D.
%
% This function performs intensity and phase calibration on raw SD data 
% using baseline absolute optical properties, then propagates these 
% calibrated values to calculate SS and DS data types (C, I, P).
%
% Inputs:
%   SD, SS, DS - Structs containing NIRS data/metadata for each category.
%   armt       - Struct defining the array arrangement and optode pairings.
%   absp       - Struct containing absolute optical property maps.
%
% Optional Name-Value-Arguments (NVA):
%   datTyps - String array of data types to calculate (Default: ["C", "I", "P"])
%   BLinds  - (nTime x 1) Logical vector defining baseline indices (Default: all true)
%   fmod    - Modulation frequency [Hz] (Default: 140.625e6)
%   nin     - Index of refraction (Default: 1.4)
%
% Outputs:
%   SD, SS, DS - Updated structs with calibrated data types.

    %% Parse Input
    arguments
        SD struct;
        SS struct;
        DS struct;
        armt struct;
        absp struct;

        NVA.datTyps = ["C", "I", "P"];

        NVA.BLinds (:,1) logical = true(size(SD.I, 1), 1);

        NVA.fmod double = 140.625e6; %Hz
        NVA.nin double = 1.4;
    end
    NVA.datTyps=unique([NVA.datTyps, "C", "I", "P"]);
    
    omega=2*pi*NVA.fmod;
    optProp.nin=NVA.nin;
    optProp.nout=1;
    optProp.musp=NaN; %1/mm
    optProp.mua=NaN; %1/mm
    
    %% Cal
    for i=1:size(SD.R, 2)
        rCent=sqrt(...
            (absp.XX-SD.loc(1, i)).^2+...
            (absp.YY-SD.loc(2, i)).^2);
        rCent(isnan(sum(absp.muaMap, 3)+sum(absp.muspMap, 3)))=Inf;
        [~, mapInd]=min(rCent, [], 'all');
        
        for Lind=1:length(SD.lambda)
            tmp=absp.muaMap(:, :, Lind);
            optProp.mua=tmp(mapInd);
            tmp=absp.muspMap(:, :, Lind);
            optProp.musp=tmp(mapInd);

            rs=[0, 0, 1/optProp.musp];
            rd=[SD.rho(i), 0, 0];
            
            Rtheo=complexReflectance(rs, rd, omega, optProp);
            RcwTheo=complexReflectance(rs, rd, 0, optProp);
            
            Rmeas=mean(SD.R(NVA.BLinds, i, Lind));
            RcwMeas=mean(SD.Rcw(NVA.BLinds, i, Lind));

            SD.Rcal(i, Lind)=Rmeas./Rtheo;
            SD.RcwCal(i, Lind)=RcwMeas./RcwTheo;

            SD.R(:, i, Lind)=SD.R(:, i, Lind)./...
                SD.Rcal(i, Lind);
            SD.Rcw(:, i, Lind)=SD.Rcw(:, i, Lind)./...
                SD.RcwCal(i, Lind);
        end
    end
    SD.isCal=true;
    
    %% Data-Types
    for mTyp=NVA.datTyps
        % SD
        for SDind=1:size(SD.R, 2)
            for Lind=1:length(SD.lambda)
                SD.(mTyp)(:, SDind, Lind)=...
                    calcData_datTyp(SD.rho(SDind),...
                    SD.R(:, SDind, Lind),...
                    SD.Rcw(:, SDind, Lind), mTyp);
            end
        end

        for Lind=1:length(SD.lambda)
            if ~any(strcmp(mTyp, ["C", "I", "P"]))
                tmpNm=join(['useSet_' mTyp], '');
                SD.(tmpNm)(:, Lind)=all([...
                    SD.useSet_C(:, Lind),...
                    SD.useSet_I(:, Lind),...
                    SD.useSet_P(:, Lind)], 2);
            end
        end

        % SS
        if ~isempty(SS)
            SS_SDinds=findSS_SDinds(armt.SSprs, armt.SDprs);
    
            for SSind=1:length(SS.drho)
                rhosTmp=SD.rho(SS_SDinds(SSind, :));
                drTmp=rhosTmp(2)-rhosTmp(1);

                for Lind=1:length(SS.lambda)
                    Rtemp=SD.R(:, SS_SDinds(SSind, :), Lind);
                    RcwTemp=SD.Rcw(:, SS_SDinds(SSind, :), Lind);
    
                    Ytmp=calcData_datTyp(rhosTmp, Rtemp, RcwTemp, mTyp);
    
                    if strcmp(mTyp, "P")
                        nom=wrapToPi(Ytmp(:, 2)-Ytmp(:, 1));
                    else
                        nom=Ytmp(:, 2)-Ytmp(:, 1);
                    end
                    SS.(mTyp)(:, SSind, Lind)=nom./drTmp;
                end
            end
            SS.isCal=true;

            for Lind=1:length(SS.lambda)
                if ~any(strcmp(mTyp, ["C", "I", "P"]))
                    tmpNm=join(['useSet_' mTyp], '');
                    SS.(tmpNm)(:, Lind)=all([...
                        SS.useSet_C(:, Lind),...
                        SS.useSet_I(:, Lind),...
                        SS.useSet_P(:, Lind)], 2);
                end
            end
        end

        % DS
        if ~isempty(DS)
            DS_SDinds=findDS_SDinds(armt.DSprs, armt.SDprs);
            
            for DSind=1:size(DS.drhos, 2)
                [rhosTmp, inds]=sort(SD.rho(DS_SDinds(DSind, :)));
                DS_SDinds(DSind, :)=DS_SDinds(DSind, inds);
                drTmp1=rhosTmp(3)-rhosTmp(1);
                drTmp2=rhosTmp(4)-rhosTmp(2);

                for Lind=1:length(DS.lambda)
                    Rtemp=SD.R(:, DS_SDinds(DSind, :), Lind);
                    RcwTemp=SD.Rcw(:, DS_SDinds(DSind, :), Lind);
    
                    Ytmp=calcData_datTyp(rhosTmp, Rtemp, RcwTemp, mTyp);

                    if strcmp(mTyp, "P")
                        nom1=wrapToPi(Ytmp(:, 3)-Ytmp(:, 1));
                        nom2=wrapToPi(Ytmp(:, 4)-Ytmp(:, 2));
                    else
                        nom1=Ytmp(:, 3)-Ytmp(:, 1);
                        nom2=Ytmp(:, 4)-Ytmp(:, 2);
                    end
                    DS.(mTyp)(:, DSind, Lind)=...
                        mean([nom1./drTmp1, nom2./drTmp2], 2);
                end
            end
            DS.isCal=true;

            for Lind=1:length(DS.lambda)
                if ~any(strcmp(mTyp, ["C", "I", "P"]))
                    tmpNm=join(['useSet_' mTyp], '');
                    DS.(tmpNm)(:, Lind)=all([...
                        DS.useSet_C(:, Lind),...
                        DS.useSet_I(:, Lind),...
                        DS.useSet_P(:, Lind)], 2);
                end
            end
        end
    end
end