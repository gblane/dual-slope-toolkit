function [datStruct] = arrayData2dmua(datStruct, absp, NVA)
% arrayData2dmua Calculates absorption changes and hemodynamics for NIRS data.
%
% [datStruct] = arrayData2dmua(datStruct, absp, NVA)
%
% Written by Giles Blaney, Ph.D. (Winter 2023)
% Modified by Cristianne Fernandez (August 2023)
%
% Inputs:
%   datStruct - Struct containing NIRS data (SD, SS, or DS format).
%               Must have gone through parseArrayData.m and rmBadChans.m.
%   absp      - Struct of absolute optical properties maps from arrayAbsMap.m.
%
% Optional Name-Value-Arguments (NVA):
%   BLinds   - (nTime x 1) Logical vector defining baseline indices (Default: all true)
%   datTyps  - String array of data types to process (e.g., ["C", "I", "P"])
%   fmod     - Modulation frequency [Hz] (Default: 140.625e6)
%   nin      - Index of refraction of the medium (Default: 1.4)
%   opt2L    - Struct for two-layer pathlength calculations:
%              - totL: '2L' (true two-layer) or 'homoEff' (effective homogeneous)
%              - thk: Layer thickness [mm]
%              - nin: [n1, n2] Index of refraction (Default: [1.4, 1.4])
%              - nout: Index of refraction outside (Default: 1)
%              - musp: [musp1, musp2] Reduced scattering [1/mm]
%              - mua: [mua1, mua2] Absorption [1/mm]
%              - fmod: Modulation frequency [Hz]
%              - h_end: Number of Bessel function zeros (Default: 2000)
%              - B: Radius of cylindrical boundary [mm] (Default: 150)
%
% Outputs:
%   datStruct - Updated input struct with dmua, dO (HbO), dD (HbR), and dT (HbT).

    %% Parse Input
    arguments
        datStruct struct;
        absp struct;
        
        NVA.BLinds (:,1) logical = true(size(datStruct.I, 1), 1);

        NVA.datTyps = ["C", "I", "P"];

        NVA.fmod double = 140.625e6; %Hz
        NVA.nin double = 1.4;

        % add feature to get a 2L pathlength to get DPF and DSF by Cristy 
        NVA.opt2L struct 
    end
    NVA.datTyps=unique([NVA.datTyps, "C", "I", "P"]);
    
    if ~isfield(NVA, 'opt2L')
        opt2L.totL = 'homoEff';
    else
        opt2L = NVA.opt2L; 
    end 

    datTyp=datStruct.typ;
    switch datTyp
        case 'SD'
            coefNm='DPF';
        case {'SS', 'DS'}
            coefNm='DSF';
        otherwise
            error('Unknown data type');
    end

    BLinds=NVA.BLinds;
    
    lambda=datStruct.lambda;
    E=makeE('OD', lambda);

    %% Set Loop
    for i=1:size(datStruct.I, 2)
        %% Find mua and musp
        rCent=sqrt(...
            (absp.XX-datStruct.loc(1, i)).^2+...
            (absp.YY-datStruct.loc(2, i)).^2);
        rCent(isnan(sum(absp.muaMap, 3)+sum(absp.muspMap, 3)))=Inf;
        [~, mapInd]=min(rCent, [], 'all');
        
        muaTmp=NaN(1, length(lambda));
        muspTmp=NaN(1, length(lambda));
        for Lind=1:length(lambda)
            tmp=absp.muaMap(:, :, Lind);
            muaTmp(Lind)=tmp(mapInd);
            tmp=absp.muspMap(:, :, Lind);
            muspTmp(Lind)=tmp(mapInd);
        end

        %% Meas Typ Loop
        for mTyp=NVA.datTyps
            nm=join([coefNm '_' mTyp], '');
            nm0=join([mTyp '0'], '');
            nmdmua=join(['dmua_' mTyp], '');

            nmdO=join(['dO_' mTyp], '');
            nmdD=join(['dD_' mTyp], '');
            nmdT=join(['dT_' mTyp], '');
            nmdX=join(['dX_' mTyp], '');

            %% Calc DPF or DSF
            switch datTyp
                case 'SD'
                    datStruct.(nm)(i, :)=DPF_DSF_calc(...
                        [0, 0, 0], [datStruct.rho(i), 0, 0],...
                        muaTmp, muspTmp, nm,...
                        NVA.fmod, NVA.nin, opt2L);
                case 'SS'
                    rdTmp=[datStruct.rhos(1:2, i), zeros(2, 2)];
                    datStruct.(nm)(i, :)=DPF_DSF_calc(...
                        [0, 0, 0; 0, 0, 0], rdTmp,...
                        muaTmp, muspTmp, nm,...
                        NVA.fmod, NVA.nin, opt2L);
                case 'DS'
                    rdTmp=[datStruct.rhos(1:2, i), zeros(2, 2)];
                    datStruct.(nm)(i, :, 1)=DPF_DSF_calc(...
                        [0, 0, 0; 0, 0, 0], rdTmp,...
                        muaTmp, muspTmp, nm,...
                        NVA.fmod, NVA.nin, opt2L);
                    rdTmp=[datStruct.rhos(3:4, i), zeros(2, 2)];
                    datStruct.(nm)(i, :, 2)=DPF_DSF_calc(...
                        [0, 0, 0; 0, 0, 0], rdTmp,...
                        muaTmp, muspTmp, nm,...
                        NVA.fmod, NVA.nin, opt2L);
                otherwise
                    error('Unknown data type');
            end
    
            %% dmua
            for lInd=1:length(lambda)
                coef=datStruct.(nm)(i, lInd, :);

                switch datTyp
                    case 'SD'
                        if strcmp(mTyp, 'P')
                            inds=~isnan(datStruct.(mTyp)(:, i, lInd));
                            datStruct.(nm0)(i, lInd)=...
                                wrapToPi(circ_mean(...
                                datStruct.(mTyp)(and(inds, BLinds), i, lInd)));

                            datStruct.(nmdmua)(:, i, lInd)=...
                                -wrapToPi(datStruct.(mTyp)(:, i, lInd)-...
                                datStruct.(nm0)(i, lInd))/...
                                (datStruct.rho(i)*coef);
                        else
                            datStruct.(nm0)(i, lInd)=...
                                squeeze(mean(...
                                datStruct.(mTyp)(BLinds, i, lInd),...
                                'omitnan'));

                            datStruct.(nmdmua)(:, i, lInd)=...
                                -(datStruct.(mTyp)(:, i, lInd)-...
                                datStruct.(nm0)(i, lInd))/...
                                (datStruct.rho(i)*coef);
                        end
                    case {'SS', 'DS'}
                        datStruct.(nm0)(i, lInd)=...
                            squeeze(mean(...
                            datStruct.(mTyp)(BLinds, i, lInd),...
                            'omitnan'));

                        datStruct.(nmdmua)(:, i, lInd)=...
                            -(datStruct.(mTyp)(:, i, lInd)-...
                            datStruct.(nm0)(i, lInd))/mean(coef, 3);
                    otherwise
                        error('Unknown data type');
                end
            end

            %% dBlood
            dmuaTmp=squeeze(datStruct.(nmdmua)(:, i, :)).';
            tmp=(E\dmuaTmp).';
            datStruct.(nmdO)(:, i)=tmp(:, 1);
            datStruct.(nmdD)(:, i)=tmp(:, 2);
            datStruct.(nmdT)(:, i)=tmp(:, 1)+tmp(:, 2); % O + D
            datStruct.(nmdX)(:, i)=tmp(:, 1)-tmp(:, 2); % O - D

        end
    end
end