function [SD, SS, DS, aux] = parseArrayData(data, armt, ISSmap, NVA)
% [SD, SS, DS, aux] = parseArrayData(data, armt, ISSmap)
% Giles Blaney Ph.D. Winter 2023
% This function takes in raw data loaded from an ISS instrument output file
% (data) along with the source-detector arrangement (armt) and map of 
% plugged in sources and detectors (ISSmap). The output is a struct for 
% each measurement type (SD, SS, and DS) along with auxiliary signals
% (aux).
% 
% INPUTS:   - data: Struct output from the load_Imagent() or load_ISS
%                   functions, containing:
%               - timemat: (t,1) Vector of time [s]
%               - fs: Sample rate [Hz]
%               - A, B,... struct with detector letter name containing:
%                   - AC: (t,ISSsrcIdx) Matrix of FD amplitude [arb.]
%                   - DC: (t,ISSsrcIdx) Matrix of CW signal [arb.]
%                   - Ph: (t,ISSsrcIdx) Metrix of FD phase [deg]
%               - AUX: (t,auxIdx) Matrix of auxiliary signals [various]
%           - armt: Struct output from the DSdisc() function, containing 
%                   (SDprs, SSprs, or DSprs can be empty):
%               - rSrc: (srcIdx,3) Matrix of source coordinates [mm]
%               - rDet: (detIdx,3) Matrix of detector coordinates [mm]
%               - SDprs: (SDidx,2) Single-Distance pairs in the format:
%                       [srcIdx1, detIdx1;
%                        srcIdx2, detIdx2;
%                          ...      ...   ]
%               - SSprs: (SDidx,2,SSidx) Single-Slope pairs in the format:
%                       SSprs(:,:,SSidx)=[srcIdx1, detIdx1;
%                                         srcIdx2, detIdx2]
%               - DSprs: (SDidx,2,SSidx,DSidx) Dual-Slope pairs in the
%                       format:
%                       DSprs(:,:,1,DSidx)=[srcIdx11, detIdx11;
%                                           srcIdx12, detIdx12]
%                       DSprs(:,:,2,DSidx)=[srcIdx21, detIdx21;
%                                           srcIdx22, detIdx22]
%           - ISSmap: Struct containing:
%               - numLam: Number of wavelengths
%               - lambda: (1,lamIdx) Vector of wavelengths [nm]
%               - lets: (1,detIdk) Char vector of detector letter names
%               - inds: (lamIdx,srcIdx,detIdx) Array containing the
%                       ISSsrcIndx (used in data struct) for each 
%                       combination of lamIdx, srcIdx, and detIdx
%           - Optional Name-Value-Arguments (NVAs):
%               - 'BLinds': (t,1) Logical vector which is true during
%                           baseline. Baseline data is used for calculation
%                           of absolute optical properties for each DS set
%               - 'BLaux': (1,:) String containing the name of the aux 
%                           channel to be used for BLinds. This overrides 
%                           the above 'BLinds' input
%               - 'doAbs': Logical controlling whether or not absolute
%                           optical properties are calculated
%               - 'muaBnd': (1,2) Vector of upper and lower bound for
%                           recovery of the absolute absorption coefficient
%                           (mua), values outside of bound return NaN
%               - 'muspBnd': (1,2) Vector of upper and lower bound for
%                           recovery of the absolute reduced scattering 
%                           coefficient (musp), values outside of bound 
%                           return NaN
%               - 'calibrateAUX_NVA': (1,:) Cell array containing 
%                           Name-Value-Argument pairs to be passed to the
%                           calibrateAUX function
% OUTPUTS:  - SD: Struct for Single-Distance (SD) containing:
%               - typ: String of value 'SD' to identify struct type
%               - t: (t,1) Vector of time [s]
%               - fs: Sampling rate [Hz]
%               - lambda: (1,lamIdx) Vector of wavelengths [nm]
%               - rho: (1,SDidx) Vector of source-detector distances [mm]
%               - loc: (3,SDidx) Matrix of SD sets optode centroid [mm]
%               - R: (t,SDidx,lamIdx) Complex Reflectance [arb./mm^2]
%               - Rcw: (t,SDidx,lamIdx) Reflectance [arb./mm^2]
%               - C: (t,SDidx,lamIdx) CW data-type {ln(rho^2Rcw)}
%               - I: (t,SDidx,lamIdx) FD amplitude data-type {ln(rho^2|R|)}
%               - P: (t,SDidx,lamIdx) FD phase {angle(R)} [rad]
%           - SS: Struct for Single-Slope (SS) containing:
%               - typ: String of value 'SS' to identify struct type
%               - t: (t,1) Vector of time [s]
%               - fs: Sampling rate [Hz]
%               - lambda: (1,lamIdx) Vector of wavelengths [nm]
%               - rhos: (SDidx,SSidx) Source-detector distances [mm]
%               - drho: (1,SSidx) Differences in distances [mm]
%               - loc: (3,SSidx) Matrix of SS sets optode centroid [mm]
%               - C: (t,SSidx,lamIdx) Slope of CW data-type
%                       {Delta ln(rho^2Rcw) / Delta rho} [1/mm]
%               - I: (t,SSidx,lamIdx) Slope of FD amplitude data-type
%                       {Delta ln(rho^2|R|) / Delta rho} [1/mm]
%               - P: (t,SSidx,lamIdx) Slope of FD phase data-type
%                       {Delta angle(R) / Delta rho} [rad/mm]
%           - DS: Struct for Dual-Slope (DS) containing:
%               - typ: String of value 'DS' to identify struct type
%               - t: (t,1) Vector of time [s]
%               - fs: Sampling rate [Hz]
%               - lambda: (1,lamIdx) Vector of wavelengths [nm]
%               - rhos: (SDidx,DSidx) Source-detector distances [mm]
%               - drho: (SSidx,DSidx) Differences in distances [mm]
%               - loc: (3,DSidx) Matrix of DS sets optode centroid [mm]
%               - mua: (DSidx,lamIdx) Matrix of recovered absorption
%                       coefficient [1/mm]
%               - musp: (DSidx,lamIdx) Matrix of recovered reduced
%                       scattering coeficent [1/mm]
%               - C: (t,DSidx,lamIdx) Slope of CW data-type
%                       {Delta ln(rho^2Rcw) / Delta rho} [1/mm]
%               - I: (t,DSidx,lamIdx) Slope of FD amplitude data-type
%                       {Delta ln(rho^2|R|) / Delta rho} [1/mm]
%               - P: (t,DSidx,lamIdx) Slope of FD phase data-type
%                       {Delta angle(R) / Delta rho} [rad/mm]
%           - aux: Struct of auxiliary data containing:
%               - t: (t,1) Vector of time [s]
%               - fs: Sampling rate [Hz]
%               - Labels: (:,1) Vector of strings containing the names of
%                       the auxilary channels
%               - XXX, YYY, ZZZ,... (t,1) Fields with names of auxiliary
%                       channel labels [various]

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