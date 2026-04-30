function [datStruct] =...
    arrayData2dmua(datStruct, absp, NVA)
% Giles Blaney Ph.D. Winter 2023
% Comments written by Cristianne Fernandez August 2023
%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datStruct - either SD, SS, or DS struct that has gone through
%           parseArrayData.m and rmBadChans.m
% absp  - absolute optical properties for each DS set output from arrayAbsMap
%%% Optional Inputs 
% BLinds - Logical vector where true indicies will be considered baseline
% datTypes - String aray of data types to run 
% fmod - modulation frequency (MHz)
% nin - index of refraction inside (-)
% opts2L         - OPTIONAL; struct - when you want to do the 2-Layer you
%                   need to have these inputs
%                   totL   - String controlling the total pathlength to be used
%                        with values:
%                        '2L'      - Use true two-layer total pathlengths
%                        'homoEff' - Use homogenous total pathlengths from
%                        effective homogenous optical properties, as-if
%                        effective homogenous optical properties were used
%                        to find DPF and DSF in a real-world measurment
%                   thk  - Layer thickness. (mm)
%                   nin  - (default=[1.4, 1.4]) Index of refraction inside. (-)
%                   nout - (default=1) Index of refraction outside. (-)
%                   musp - Size lambda x layer (default=[1.20, 0.25] 1/mm) 
%                           Reduced scattering. (1/mm)
%                   mua  - Size lambda x layer (default=[0.008, 0.020] 1/mm) Absorption. (1/mm) 
%                   fmod   - (default=1.40625 Hz) Modulation frequency {Hz}
%                   h_end  - (default=2000) Number of Bessel function zeros
%                   B      - (default=150 mm) Radius of cylindrical boundary {mm}
%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataStruct -  either SD, SS, or DS struct that now has dmua, Oxy and
%               DeOxy calculated for the datTypes


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