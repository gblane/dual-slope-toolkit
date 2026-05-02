function [SD, SS, DS, aux] = parseArrayData(data, armt, ISSmap, NVA)
% parseArrayData Parses raw ISS NIRS data into SD, SS, and DS structures.
%
% [SD, SS, DS, aux] = parseArrayData(data, armt, ISSmap, NVA)
%
% Written by Giles Blaney, Ph.D. (Winter 2023)
%
% This function processes raw frequency-domain data from an ISS instrument
% and organizes it into Single-Distance (SD), Single-Slope (SS), and 
% Dual-Slope (DS) structures based on a specified array arrangement.
% 
% INPUTS:
%   data - Struct output from load_Imagent() or load_ISS(), containing:
%          - timemat: (t x 1) Vector of time [s]
%          - fs: Sampling rate [Hz]
%          - A, B, ...: Detector structs with AC, DC, and Ph matrices.
%          - AUX: (t x nAux) Matrix of auxiliary signals.
%   armt - Struct from DSdisc() defining the array geometry:
%          - rSrc, rDet: Optode coordinates [mm]
%          - SDprs, SSprs, DSprs: Optode pairing definitions.
%   ISSmap - Struct mapping ISS channels to the physical array:
%          - numLam: Number of wavelengths.
%          - lambda: Vector of wavelengths [nm].
%          - lets: Detector letter names.
%          - inds: Mapping array (lambda x source x detector).
%
% Optional Name-Value-Arguments (NVA):
%   BLinds - (t x 1) Logical vector defining baseline indices (Default: all true)
%   BLaux  - Name of the auxiliary channel to use for defining BLinds.
%   doAbs  - Logical flag to calculate absolute optical properties (Default: true)
%   muaBnd - [min, max] Bounds for absolute absorption recovery [1/mm].
%   muspBnd - [min, max] Bounds for absolute reduced scattering recovery [1/mm].
%   calibrateAUX_NVA - Cell array of arguments for calibrateAUX().
%
% OUTPUTS:
%   SD, SS, DS - Structs containing processed data types (C, I, P), 
%                coordinates, and metadata for each measurement category.
%   aux - Struct containing calibrated auxiliary signals.

    %% Parse Input
    arguments
        data struct;
        armt struct;
        ISSmap struct;
        
        NVA.BLinds (:,1) logical = true(length(data.timemat), 1);
        NVA.BLaux string = [];
        
        NVA.doAbs logical = true;
        
        NVA.muaBnd (1,2) double = [0, 0.1]; %1/mm
        NVA.muspBnd (1,2) double = [0, 10]; %1/mm

        NVA.calibrateAUX_NVA (1,:) cell = {};
    end

    tN=length(data.timemat);
    tInds=1:tN;
    
    tmp=struct2pairs(NVA);
    tmp=[tmp, NVA.calibrateAUX_NVA];
    aux=calibrateAUX(data, tmp{:});

    if ~isempty(NVA.BLaux)
        NVA.BLinds=aux.(NVA.BLaux);
    end

    %% SD
    if ~isempty(armt.SDprs)
        for SDind=1:size(armt.SDprs, 1)
            SrcInd=armt.SDprs(SDind, 1);
            DetInd=armt.SDprs(SDind, 2);
            
            DetISSnam=ISSmap.lets(DetInd);
            
            for Lind=1:ISSmap.numLam
                SrcISSind=ISSmap.inds(Lind, SrcInd, DetInd);
                
                if ~isnan(SrcISSind)
                    % [t, SDind, lambda]
                    SD.Rcw(:, SDind, Lind)=data.(DetISSnam).DC(:, SrcISSind);
                    ACtmp=data.(DetISSnam).AC(:, SrcISSind);
                    
                    Ptmp=wrapTo2Pi(data.(DetISSnam).Ph(:, SrcISSind)*pi/180);
                    if max(diff(Ptmp))>pi
                        Ptmp=wrapToPi(Ptmp);
                    end
                    indsTemp=find(diff(Ptmp)>pi/16);
                    PtrshInds=unique([indsTemp; indsTemp+1]);
                    PtrshInds(PtrshInds>tN)=[];
                    Ptmp(PtrshInds)=NaN;
                    Ptmp=interp1(tInds(~isnan(Ptmp)), Ptmp(~isnan(Ptmp)), tInds,...
                        'linear', 'extrap');
                    
                    SD.R(:, SDind, Lind)=ACtmp.*exp(1i*Ptmp.');
                else
                    SD.Rcw(:, SDind, Lind)=NaN(tN, 1);
                    SD.R(:, SDind, Lind)=NaN(tN, 1);
                end
            end
            
            SD.rho(1, SDind)=...
                vecnorm(armt.rSrc(SrcInd, :)-armt.rDet(DetInd, :));
            SD.loc(:, SDind)=...
                mean([armt.rSrc(SrcInd, :); armt.rDet(DetInd, :)]);
        end
        
        for mTyp=["C", "I", "P"]
            SD.(mTyp)=calcData_datTyp(SD.rho, SD.R, SD.Rcw, mTyp);
        end
        SD.typ='SD';
        SD.t=data.timemat;
        SD.fs=data.fs;
        SD.lambda=ISSmap.lambda;
    else
        SD=[];
    end

    %% SS
    if ~isempty(armt.SSprs)
        for SSind=1:size(armt.SSprs, 3)
            SrcInds=armt.SSprs(:, 1, SSind);
            DetInds=armt.SSprs(:, 2, SSind);
            
            for Lind=1:ISSmap.numLam
                rhosTmp=NaN(size(SrcInds.'));
                locsTmp=NaN(length(SrcInds), 3);
                RcwTemp=NaN(tN, length(SrcInds));
                Rtemp=NaN(tN, length(SrcInds));
                
                for SDind=1:length(SrcInds)
                    SrcInd=SrcInds(SDind);
                    DetInd=DetInds(SDind);
                    
                    DetISSnam=ISSmap.lets(DetInd);
                    
                    rhosTmp(SDind)=vecnorm(armt.rSrc(SrcInd, :)-...
                        armt.rDet(DetInd, :));
                    locsTmp(SDind, :)=mean([armt.rSrc(SrcInd, :);...
                        armt.rDet(DetInd, :)]);
                    
                    SrcISSind=ISSmap.inds(Lind, SrcInd, DetInd);

                    if ~isnan(SrcISSind)
                        % [t, SDind]
                        RcwTemp(:, SDind)=data.(DetISSnam).DC(:, SrcISSind);
                        ACtmp=data.(DetISSnam).AC(:, SrcISSind);
                        
                        Ptmp=wrapTo2Pi(data.(DetISSnam).Ph(:, SrcISSind)*pi/180);
                        if max(diff(Ptmp))>pi
                            Ptmp=wrapToPi(Ptmp);
                        end
                        indsTemp=find(diff(Ptmp)>pi/8);
                        PtrshInds=unique([indsTemp; indsTemp+1]);
                        PtrshInds(PtrshInds>tN)=[];
                        Ptmp(PtrshInds)=NaN;
                        Ptmp=interp1(...
                            tInds(~isnan(Ptmp)), Ptmp(~isnan(Ptmp)), tInds,...
                            'linear', 'extrap');

                        Rtemp(:, SDind)=ACtmp.*exp(1i*Ptmp.');
                    else
                        RcwTemp(:, SDind)=NaN(tN, 1);
                        Rtemp(:, SDind)=NaN(tN, 1);
                    end
                end
                
                drTmp=rhosTmp(2)-rhosTmp(1);
                for mTyp=["C", "I", "P"]
                    Ytmp=calcData_datTyp(rhosTmp, Rtemp, RcwTemp, mTyp);

                    if strcmp(mTyp, "P")
                        nom=wrapToPi(Ytmp(:, 2)-Ytmp(:, 1));
                    else
                        nom=Ytmp(:, 2)-Ytmp(:, 1);
                    end
                    SS.(mTyp)(:, SSind, Lind)=nom./drTmp;
                end
                
                SS.rhos(1:2, SSind)=rhosTmp;
                SS.drho(1, SSind)=rhosTmp(2)-rhosTmp(1);
                SS.loc(1:3, SSind)=mean(locsTmp);
            end
        end
        SS.typ='SS';
        SS.t=data.timemat;
        SS.fs=data.fs;
        SS.lambda=ISSmap.lambda;
    else
        SS=[];
    end
    
    %% DS
    if ~isempty(armt.DSprs)
        for DSind=1:size(armt.DSprs, 4)
            SSprsTmp=armt.DSprs(:, :, :, DSind);
            
            for Lind=1:ISSmap.numLam
                SSrhosTmp=NaN(2, size(SSprsTmp, 3));
                SSdrhoTmp=NaN(size(SSprsTmp, 3), 1);
                SSlocTmp=NaN(size(SSprsTmp, 3), 3);
                SSCtmp=NaN(tN, size(SSprsTmp, 3));
                SSItmp=NaN(tN, size(SSprsTmp, 3));
                SSPtmp=NaN(tN, size(SSprsTmp, 3));
                R=NaN(size(SSprsTmp, 3), size(SSprsTmp, 1));
                for SSind=1:size(SSprsTmp, 3)
                    SrcInds=SSprsTmp(:, 1, SSind);
                    DetInds=SSprsTmp(:, 2, SSind);
                    
                    rhosTmp=NaN(size(SrcInds.'));
                    locsTmp=NaN(length(SrcInds), 3);
                    RcwTemp=NaN(tN, length(SrcInds));
                    Rtemp=NaN(tN, length(SrcInds));
                    for SDind=1:length(SrcInds)
                        SrcInd=SrcInds(SDind);
                        DetInd=DetInds(SDind);
                        
                        DetISSnam=ISSmap.lets(DetInd);
                        
                        rhosTmp(SDind)=vecnorm(armt.rSrc(SrcInd, :)-...
                            armt.rDet(DetInd, :));
                        locsTmp(SDind, :)=mean([armt.rSrc(SrcInd, :);...
                            armt.rDet(DetInd, :)]);
        
                        SrcISSind=ISSmap.inds(Lind, SrcInd, DetInd);

                        if ~isnan(SrcISSind)
                             % [t, SDind]
                            RcwTemp(:, SDind)=data.(DetISSnam).DC(:, SrcISSind);
                            ACtmp=data.(DetISSnam).AC(:, SrcISSind);
            
                            Ptmp=wrapTo2Pi(data.(DetISSnam).Ph(:, SrcISSind)*pi/180);
                            if max(diff(Ptmp))>pi
                                Ptmp=wrapToPi(Ptmp);
                            end
                            indsTemp=find(diff(Ptmp)>pi/8);
                            PtrshInds=unique([indsTemp; indsTemp+1]);
                            PtrshInds(PtrshInds>tN)=[];
                            Ptmp(PtrshInds)=NaN;
                            Ptmp=interp1(...
                                tInds(~isnan(Ptmp)), Ptmp(~isnan(Ptmp)), tInds,...
                                'linear', 'extrap');

                            Rtemp(:, SDind)=ACtmp.*exp(1i*Ptmp.');
                        else
                            RcwTemp(:, SDind)=NaN(tN, 1);
                            Rtemp(:, SDind)=NaN(tN, 1);
                        end
                        R(SSind, SDind)=mean(Rtemp(NVA.BLinds, SDind));
                    end
                    if  SSind == 1
                        DS.R(:, 1:2, DSind, Lind) = Rtemp;

                    elseif SSind==2
                        DS.R(:, 3:4, DSind, Lind) = Rtemp;
                    end
                    drTmp=rhosTmp(2)-rhosTmp(1);
                    for mTyp=["C", "I", "P"]
                        Ytmp=calcData_datTyp(rhosTmp, Rtemp, RcwTemp, mTyp);
    
                        if strcmp(mTyp, "P")
                            nom=wrapToPi(Ytmp(:, 2)-Ytmp(:, 1));
                        else
                            nom=Ytmp(:, 2)-Ytmp(:, 1);
                        end
                        eval(['SS' convertStringsToChars(mTyp)...
                            'tmp(:, SSind)=nom./drTmp;']);
                    end
                    
                    SSrhosTmp(1:2, SSind)=rhosTmp;
                    SSdrhoTmp(SSind, 1)=rhosTmp(2)-rhosTmp(1);
                    SSlocTmp(SSind, 1:3)=mean(locsTmp);
                end
                
                DS.rhos(1:4, DSind)=[SSrhosTmp(1:2, 1); SSrhosTmp(1:2, 2)];
                DS.drhos(1:2, DSind)=SSdrhoTmp;
                DS.loc(1:3, DSind)=mean(SSlocTmp);
                
                % [t, DSind, lambda]
                DS.C(:, DSind, Lind)=mean(SSCtmp, 2);
                DS.I(:, DSind, Lind)=mean(SSItmp, 2);
                DS.P(:, DSind, Lind)=mean(SSPtmp, 2);
                
                if NVA.doAbs
                    RRtmp=[R(1, :), R(2, :)];
                    rhosTmp=[SSrhosTmp(1:2, 1); SSrhosTmp(1:2, 2)]';
                    [DS.mua(DSind, Lind), DS.musp(DSind, Lind), ~]=...
                        DSR2muamuspEB_iterRecov(rhosTmp, RRtmp);                    

                    if DS.musp(DSind, Lind)<NVA.muspBnd(1) ||...
                            DS.mua(DSind, Lind)<NVA.muaBnd(1) ||...
                            DS.musp(DSind, Lind)>NVA.muspBnd(2) ||...
                            DS.mua(DSind, Lind)>NVA.muaBnd(2)
                        DS.musp(DSind, Lind)=NaN;
                        DS.mua(DSind, Lind)=NaN;
                    end
                end
            end
        end
        DS.typ='DS';
        DS.t=data.timemat;
        DS.fs=data.fs;
        DS.lambda=ISSmap.lambda;
    else
        DS=[];
    end
end