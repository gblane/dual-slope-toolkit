function [nois] = arrayNoiseMap(datStruct, nameValArgs)
% arrayNoiseMap Generates a 2D map of noise distribution across a NIRS array.
%
% [nois] = arrayNoiseMap(datStruct, nameValArgs)
%
% Written by Giles Blaney, Ph.D.
%
% This function interpolates discrete channel noise measurements onto a 
% continuous 2D spatial map using Gaussian smoothing.
%
% Inputs:
%   datStruct - Struct containing NIRS data and noise metrics (Cnoise, Inoise, Pnoise).
%
% Optional Name-Value-Arguments:
%   method - Smoothing method to be used ('GaussCent')
%   len    - Characteristic length of Gaussian smoothing [mm] (Default: 30)
%   dr     - Spatial resolution of the output grid [mm] (Default: 1)
%   rHalo  - Radius around centroids to apply noise values [mm] (Default: 30)
%
% Outputs:
%   nois - Struct containing the generated noise maps and grid coordinates.

    %% Parse Input
    arguments
        datStruct struct;
        
        nameValArgs.method string ...
            {mustBeMember(nameValArgs.method,{'GaussCent'})} = 'GaussCent';
        nameValArgs.len double = 30; %mm
        
        nameValArgs.dr double = 1; %mm
        nameValArgs.rHalo double = 30; %mm
    end
    
    nois.typ=datStruct.typ;
    nois.loc=datStruct.loc;
    nois.lambda=datStruct.lambda;

    datStruct.Cnoise(isinf(datStruct.Cnoise))=NaN;
    datStruct.Inoise(isinf(datStruct.Inoise))=NaN;
    datStruct.Pnoise(isinf(datStruct.Pnoise))=NaN;

    %% Make Coor System
    nois.x=floor(min(nois.loc(1, :))-nameValArgs.rHalo):nameValArgs.dr:...
        ceil(max(nois.loc(1, :))+nameValArgs.rHalo);
    nois.y=floor(min(nois.loc(2, :))-nameValArgs.rHalo):nameValArgs.dr:...
        ceil(max(nois.loc(2, :))+nameValArgs.rHalo);
    [nois.XX, nois.YY]=meshgrid(nois.x, nois.y);
    
    %% Run
    switch nameValArgs.method
        case 'GaussCent'
            for Lind=1:length(datStruct.lambda)
                wStack=NaN([size(nois.XX), size(nois.loc, 2)]);

                for setInd=1:size(nois.loc, 2)
                    rCentTmp=sqrt((nois.XX-nois.loc(1, setInd)).^2+...
                        (nois.YY-nois.loc(2, setInd)).^2);
                    wStack(:, :, setInd)=...
                        exp(-rCentTmp.^2/(2*nameValArgs.len^2));
                end
                wStack=wStack./sum(wStack, 3);

                nois.CnoiseMap(:, :, Lind)=sum(wStack.*...
                    reshape(datStruct.Cnoise(:, Lind),...
                    [1, 1, size(nois.loc, 2)]), 3, 'omitnan');
                nois.InoiseMap(:, :, Lind)=sum(wStack.*...
                    reshape(datStruct.Inoise(:, Lind),...
                    [1, 1, size(nois.loc, 2)]), 3, 'omitnan');
                nois.PnoiseMap(:, :, Lind)=sum(wStack.*...
                    reshape(datStruct.Pnoise(:, Lind),...
                    [1, 1, size(nois.loc, 2)]), 3, 'omitnan');
            end
        otherwise
            error('Unknown method');
    end

    %% Mask Halo
    mStack=false(size(wStack));
    for setInd=1:size(nois.loc, 2)
        rCentTmp=sqrt((nois.XX-nois.loc(1, setInd)).^2+...
            (nois.YY-nois.loc(2, setInd)).^2);
        mStack(:, :, setInd)=rCentTmp<=nameValArgs.rHalo;
    end
    mask=any(mStack, 3);
    for Lind=1:length(nois.lambda)
        temp=nois.CnoiseMap(:, :, Lind);
        temp(~mask)=NaN;
        nois.CnoiseMap(:, :, Lind)=temp;

        temp=nois.InoiseMap(:, :, Lind);
        temp(~mask)=NaN;
        nois.InoiseMap(:, :, Lind)=temp;

        temp=nois.PnoiseMap(:, :, Lind);
        temp(~mask)=NaN;
        nois.PnoiseMap(:, :, Lind)=temp;
    end
end