function [absp] = arrayAbsMap(DS, NVA)
% arrayAbsMap Generates a 2D map of absolute optical properties from DS data.
%
% [absp] = arrayAbsMap(DS, NVA)
%
% Written by Giles Blaney, Ph.D. (Fall 2022)
% Modified by Cristianne Fernandez (Aug 2023)
%
% Inputs:
%   DS - Struct containing absolute optical properties for each Dual-Slope set.
%        Fields:
%        - mua: (nDS x nLambda) Absolute absorption [1/mm]
%        - musp: (nDS x nLambda) Absolute reduced scattering [1/mm]
%        - loc: (3 x nDS) Coordinates of DS set centroids [mm]
%        - lambda: (1 x nLambda) Wavelengths [nm]
%
% Optional Name-Value-Arguments (NVA):
%   method - Smoothing method to be used ('GaussCent')
%   len    - Characteristic length of Gaussian smoothing [mm] (Default: 30)
%   dr     - Spatial resolution of the coordinate system [mm] (Default: 1)
%   rHalo  - Radius to apply absolute properties around each centroid [mm] (Default: 30)
%
% Outputs:
%   absp - Struct containing the generated maps of absolute optical properties.
%        Fields:
%        - muaMap, muspMap: 2D maps for each wavelength.
%        - x, y, XX, YY: Grid coordinates.

%% Parse Input
    arguments
        DS struct;

        NVA.method string ...
            {mustBeMember(NVA.method,{'GaussCent'})} = 'GaussCent';
        NVA.len double = 30; %mm
        
        NVA.dr double = 1; %mm
        NVA.rHalo double = 30; %mm
    end

    absp.mua=DS.mua;
    absp.musp=DS.musp;
    absp.loc=DS.loc;
    absp.lambda=DS.lambda;

    %% Make Coor System
    absp.x=floor(min(absp.loc(1, :))-NVA.rHalo):NVA.dr:...
        ceil(max(absp.loc(1, :))+NVA.rHalo);
    absp.y=floor(min(absp.loc(2, :))-NVA.rHalo):NVA.dr:...
        ceil(max(absp.loc(2, :))+NVA.rHalo);
    [absp.XX, absp.YY]=meshgrid(absp.x, absp.y);
    
    %% Run
    switch NVA.method
        case 'GaussCent'
            for Lind=1:size(DS.mua, 2)
                wStack=NaN([size(absp.XX), size(absp.mua, 1)]);
                muaTmp=absp.mua(:, Lind);
                muspTmp=absp.musp(:, Lind);

                for DSind=1:size(DS.mua, 1)
                    rCentTmp=sqrt((absp.XX-absp.loc(1, DSind)).^2+...
                        (absp.YY-absp.loc(2, DSind)).^2);
                    wStack(:, :, DSind)=...
                        exp(-rCentTmp.^2/(2*NVA.len^2));
                end
                
                if ~all(isnan(absp.mua(:, Lind)))
                    wStack_mua=wStack./sum(wStack(:, :,...
                        ~isnan(muaTmp)), 3);
                    absp.muaMap(:, :, Lind)=sum(wStack_mua.*...
                        reshape(muaTmp,...
                        [1, 1, size(absp.mua, 1)]), 3, 'omitnan');
                else
                    mua_other=absp.mua(:, setdiff(1:end, Lind));
                    lam_other=absp.lambda(setdiff(1:end, Lind));
                    if ~all(isnan(mua_other(:)))
                        mua_other_avg=mean(mua_other, 'omitnan');
                        lam_other(isnan(mua_other_avg))=[];
                        mua_other_avg(isnan(mua_other_avg))=[];
                        
                        Eother=makeE('T', lam_other);
                        T=Eother\mua_other_avg.';

                        E=makeE('T', absp.lambda(Lind));
                        mua_ass=E*T;
                        warning('T from other lambda assumed for mua');
                    else
                        [mua_ass, ~, ~, ~]=...
                            assumeOptProp(absp.lambda(Lind));
                        warning('Default mua assumed');
                    end
                    absp.muaMap(:, :, Lind)=mua_ass*...
                        ones(size(absp.XX));
                end

                if ~all(isnan(absp.musp(:, Lind)))
                    wStack_musp=wStack./sum(wStack(:, :,...
                        ~isnan(muspTmp)), 3);
                    absp.muspMap(:, :, Lind)=sum(wStack_musp.*...
                        reshape(muspTmp,...
                        [1, 1, size(absp.musp, 1)]), 3, 'omitnan');
                else
                    musp_other=absp.musp(:, setdiff(1:end, Lind));
                    lam_other=absp.lambda(setdiff(1:end, Lind));
                    if ~all(isnan(musp_other(:)))
                        musp_other_avg=mean(musp_other, 'omitnan');
                        lam_other(isnan(musp_other_avg))=[];
                        musp_other_avg(isnan(musp_other_avg))=[];
                        
                        if length(musp_other_avg)==1
                            [~, ~, ~, b_ass]=...
                                assumeOptProp(absp.lambda(Lind));
                            
                            musp_ass=musp_other_avg*...
                                (absp.lambda(Lind)/lam_other)^-b_ass;
                            
                            warning('b assumed for musp');
                        elseif length(musp_other_avg)==2
                            b_ass=-log(musp_other_avg(1)/musp_other_avg(2))/...
                                log(lam_other(1)/lam_other(2));

                            musp_ass=musp_other_avg(1)*...
                                (absp.lambda(Lind)/lam_other(1))^-b_ass;

                            warning('b from other lambda assumed for musp');
                        else
                            x=lam_other;
                            y=musp_other_avg;
                            f=fit(x, y, 'a*(x/700)^-b');
                            musp_ass=f.a*(absp.lambda(Lind)/700)^-f.b;

                            warning('b from other lambda assumed for musp');
                        end
                    else
                        [~, musp_ass, ~, ~]=...
                            assumeOptProp(absp.lambda(Lind));
                        warning('Default musp assumed');
                    end
                    absp.muspMap(:, :, Lind)=musp_ass*...
                        ones(size(absp.XX));
                end
            end
        otherwise
            error('Unknown method');
    end

    %% Mask Halo
    mStack=false(size(wStack));
    for DSind=1:size(DS.mua, 1)
        rCentTmp=sqrt((absp.XX-absp.loc(1, DSind)).^2+...
            (absp.YY-absp.loc(2, DSind)).^2);
        mStack(:, :, DSind)=rCentTmp<=NVA.rHalo;
    end
    mask=any(mStack, 3);
    for Lind=1:size(DS.mua, 2)
        temp=absp.muaMap(:, :, Lind);
        temp(~mask)=NaN;
        absp.muaMap(:, :, Lind)=temp;

        temp=absp.muspMap(:, :, Lind);
        temp(~mask)=NaN;
        absp.muspMap(:, :, Lind)=temp;
    end
end