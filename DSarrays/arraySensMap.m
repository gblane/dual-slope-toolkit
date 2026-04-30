function [sens] = arraySensMap(datStruct, armt, absp, NVA)
% Giles Blaney Ph.D. Winter 2023
    %% Parse Input
    arguments
        datStruct struct;
        armt struct;
        absp struct;

        NVA.Sthresh (1,1) double = 5e-5;

        NVA.z_min (1,1) double = 0; %mm
        NVA.z_max (1,1) double = 30; %mm
        NVA.z_lay (1,1) double = 5; %mm

        NVA.haloRad (1,1) double = 35; %mm

        NVA.dr (1,3) double = [1, 1, 1]; %mm

        NVA.nin (1,1) double = 1.4;
        NVA.nou (1,1) double = 1;
        
        NVA.fmod (1,1) double = 140.625e6; %Hz

        NVA.datTyps = ["C", "I", "P"];
    end
    NVA.datTyps=unique([NVA.datTyps, "C", "I", "P"]);

    argNms=fieldnames(NVA);
    for i=1:length(argNms)
        eval(sprintf('%s=NVA.%s;',...
            argNms{i}, argNms{i}));
    end
    clear NVA;
    
    sens.V=prod(dr); %mm^3
    sens.z_all=z_min:dr(3):z_max;
    sens.z_lay2=sens.z_all>=z_lay;

    sens.typ=datStruct.typ;

    %% Make Sens Map Initial r
    for mTyp=datTyps
        if any(strcmp(mTyp, ["C", "I", "P"]))
            eval(sprintf('useSrc=all(armt.%suseSrc_%s, 2);',...
                datStruct.typ, mTyp));
            eval(sprintf('useDet=all(armt.%suseDet_%s, 2);',...
                datStruct.typ, mTyp));
        else
            useSrc=all([...
                armt.([datStruct.typ 'useSrc_C']),...
                armt.([datStruct.typ 'useSrc_I']),...
                armt.([datStruct.typ 'useSrc_P'])], 2);
            useDet=all([...
                armt.([datStruct.typ 'useDet_C']),...
                armt.([datStruct.typ 'useDet_I']),...
                armt.([datStruct.typ 'useDet_P'])], 2);
        end


        AllOpt=[armt.rSrc(useSrc, :); armt.rDet(useDet, :)];
        if ~isempty(AllOpt)
            shp=alphaShape(AllOpt(:, 1), AllOpt(:, 2), inf);
    
            x_min=min([round(AllOpt(:, 1), -1);...
                round(AllOpt(:, 1), -1)-10])-haloRad;
            x_max=max([round(AllOpt(:, 1), -1);...
                round(AllOpt(:, 1), -1)+10])+haloRad;
            y_min=min([round(AllOpt(:, 2), -1);...
                round(AllOpt(:, 2), -1)-10])-haloRad;
            y_max=max([round(AllOpt(:, 2), -1);...
                round(AllOpt(:, 2), -1)+10])+haloRad;
    
            [XX, YY]=meshgrid(x_min:dr(1):x_max, y_min:dr(2):y_max);
    
            r_xy=[XX(:), YY(:)];
            inds=or(...
                min(sqrt((AllOpt(:, 1).'-r_xy(:, 1)).^2+...
                (AllOpt(:, 2).'-r_xy(:, 2)).^2), [], 2)<=haloRad,...
                inShape(shp, r_xy(:, 1), r_xy(:, 2)));
            
            eval(sprintf('sens.rxy_%s=r_xy(inds, :);',...
                mTyp));
        else
            eval(sprintf('sens.rxy_%s=[];',...
                mTyp));
        end
    end

    %% Find mua and musp
    mua=NaN(size(datStruct.loc, 2), length(datStruct.lambda));
    musp=NaN(size(datStruct.loc, 2), length(datStruct.lambda));
    for prsInd=1:size(datStruct.loc, 2)
        rCent=sqrt(...
            (absp.XX-datStruct.loc(1, prsInd)).^2+...
            (absp.YY-datStruct.loc(2, prsInd)).^2);
        rCent(isnan(sum(absp.muaMap, 3)+sum(absp.muspMap, 3)))=Inf;
        [~, mapInd]=min(rCent, [], 'all');
        
        for Lind=1:length(datStruct.lambda)
            tmp=absp.muaMap(:, :, Lind);
            mua(prsInd, Lind)=tmp(mapInd);
            tmp=absp.muspMap(:, :, Lind);
            musp(prsInd, Lind)=tmp(mapInd);
        end
    end
    
    %% Pathlengths
    optProp.nin=nin;
    optProp.nout=nou;
    omega=2*pi*fmod;

    switch datStruct.typ
        case 'SD'
            for Lind=1:length(datStruct.lambda)
                for SDind=1:length(datStruct.rho)
                    SrcInd=armt.SDprs(SDind, 1);
                    DetInd=armt.SDprs(SDind, 2);
            
                    optProp.mua=mua(SDind, Lind);
                    optProp.musp=musp(SDind, Lind);
                    
                    rsTMP=armt.rSrc(SrcInd, :)+[0, 0, 1/optProp.musp];
                    rdTMP=armt.rDet(DetInd, :);
                    
                    [Ltmp, Rtmp]=...
                        complexTotPathLen(rsTMP, rdTMP, omega, optProp);
                    [LcwTmp, RcwTmp]=...
                        complexTotPathLen(rsTMP, rdTMP, 0, optProp);
                    
                    for mTyp=datTyps
                        % [r_xy, layer, SD, lambda]
                        Lnm=join(['L_' mTyp], '');
                        sens.(Lnm)(1, 1, SDind, Lind)=...
                            calcPathLen_datTyp(Rtmp, Ltmp, RcwTmp, LcwTmp,...
                            mTyp);

                        eval(sprintf(...
                            'r_xy=sens.rxy_%s;',...
                            mTyp));
                        
                        if ~isempty(r_xy)
                            lTmp_1=zeros(size(r_xy, 1), 1);
                            lcwTmp_1=zeros(size(r_xy, 1), 1);
                            lTmp_2=zeros(size(r_xy, 1), 1);
                            lcwTmp_2=zeros(size(r_xy, 1), 1);
                            for zInd=1:length(sens.z_all)
                                
                                lTmp=complexPartPathLen(...
                                    rsTMP,...
                                    [r_xy, ones(size(r_xy, 1), 1)*sens.z_all(zInd)],...
                                    rdTMP, sens.V, omega, optProp);
                                lcwTmp=complexPartPathLen(...
                                    rsTMP,...
                                    [r_xy, ones(size(r_xy, 1), 1)*sens.z_all(zInd)],...
                                    rdTMP, sens.V, 0, optProp);
            
                                lTmp(isnan(lTmp))=0;
                                lcwTmp(isnan(lcwTmp))=0;
            
                                if ~sens.z_lay2(zInd)
                                    lTmp_1=lTmp_1+lTmp;
                                    lcwTmp_1=lcwTmp_1+lcwTmp;
                                else
                                    lTmp_2=lTmp_2+lTmp;
                                    lcwTmp_2=lcwTmp_2+lcwTmp;
                                end
                            end
                            
                            % [r_xy, layer, SD, lambda]
                            lnm=join(['l_' mTyp], '');
                            sens.(lnm)(:, 1, SDind, Lind)=...
                                calcPathLen_datTyp(...
                                Rtmp, lTmp_1, RcwTmp, lcwTmp_1,...
                                mTyp);
                            sens.(lnm)(:, 2, SDind, Lind)=...
                                calcPathLen_datTyp(...
                                Rtmp, lTmp_2, RcwTmp, lcwTmp_2,...
                                mTyp);
                        else
                            eval(sprintf(...
                                'sens.l_%s=[];',...
                                mTyp));
                        end
                    end
                end
            end

        case 'SS'
            for Lind=1:length(datStruct.lambda)
                for SSind=1:size(datStruct.rhos, 2)
                    SrcInds=armt.SSprs(:, 1, SSind);
                    DetInds=armt.SSprs(:, 2, SSind);
                    
                    optProp.mua=mua(SSind, Lind);
                    optProp.musp=musp(SSind, Lind);
                    
                    for SDind=1:length(SrcInds)
                        SrcInd=SrcInds(SDind);
                        DetInd=DetInds(SDind);
                        
                        rsTMP=armt.rSrc(SrcInd, :)+[0, 0, 1/optProp.musp];
                        rdTMP=armt.rDet(DetInd, :);
                        
                        [Ltmp, Rtmp]=...
                            complexTotPathLen(rsTMP, rdTMP, omega, optProp);
                        [LcwTmp, RcwTmp]=...
                            complexTotPathLen(rsTMP, rdTMP, 0, optProp);

                        for mTyp=datTyps
                            % [r_xy, layer, SD, SS, lambda]
                            Lnm=join(['L_' mTyp], '');
                            sens.(Lnm)(1, 1, SDind, SSind, Lind)=...
                                calcPathLen_datTyp(Rtmp, Ltmp, RcwTmp, LcwTmp,...
                                mTyp);
                        
                            eval(sprintf(...
                                'r_xy=sens.rxy_%s;',...
                                mTyp));
                            
                            if ~isempty(r_xy)
                                lTmp_1=zeros(size(r_xy, 1), 1);
                                lcwTmp_1=zeros(size(r_xy, 1), 1);
                                lTmp_2=zeros(size(r_xy, 1), 1);
                                lcwTmp_2=zeros(size(r_xy, 1), 1);
                                for zInd=1:length(sens.z_all)
                                    
                                    lTmp=complexPartPathLen(...
                                        rsTMP,...
                                        [r_xy, ones(size(r_xy, 1), 1)*sens.z_all(zInd)],...
                                        rdTMP, sens.V, omega, optProp);
                                    lcwTmp=complexPartPathLen(...
                                        rsTMP,...
                                        [r_xy, ones(size(r_xy, 1), 1)*sens.z_all(zInd)],...
                                        rdTMP, sens.V, 0, optProp);
                
                                    lTmp(isnan(lTmp))=0;
                                    lcwTmp(isnan(lcwTmp))=0;
                
                                    if ~sens.z_lay2(zInd)
                                        lTmp_1=lTmp_1+lTmp;
                                        lcwTmp_1=lcwTmp_1+lcwTmp;
                                    else
                                        lTmp_2=lTmp_2+lTmp;
                                        lcwTmp_2=lcwTmp_2+lcwTmp;
                                    end
                                end
                                
                                % [r_xy, layer, SD, SS, lambda]
                                lnm=join(['l_' mTyp], '');
                                sens.(lnm)(:, 1, SDind, SSind, Lind)=...
                                    calcPathLen_datTyp(...
                                    Rtmp, lTmp_1, RcwTmp, lcwTmp_1,...
                                    mTyp);
                                sens.(lnm)(:, 2, SDind, SSind, Lind)=...
                                    calcPathLen_datTyp(...
                                    Rtmp, lTmp_2, RcwTmp, lcwTmp_2,...
                                    mTyp);
                            else
                                eval(sprintf(...
                                    'sens.l_%s=[];',...
                                    mTyp));
                            end
                            
                        end
                    end
                end
            end

        case 'DS'
            for Lind=1:length(datStruct.lambda)
                for DSind=1:size(datStruct.rhos, 2)
                    SSprsTmp=armt.DSprs(:, :, :, DSind);
                    
                    optProp.mua=mua(DSind, Lind);
                    optProp.musp=musp(DSind, Lind);
                    
                    for SSind=1:size(SSprsTmp, 3)
                        SrcInds=SSprsTmp(:, 1, SSind);
                        DetInds=SSprsTmp(:, 2, SSind);
                    
                        for SDind=1:length(SrcInds)
                            SrcInd=SrcInds(SDind);
                            DetInd=DetInds(SDind);
            
                            rsTMP=armt.rSrc(SrcInd, :)+[0, 0, 1/optProp.musp];
                            rdTMP=armt.rDet(DetInd, :);

                            [Ltmp, Rtmp]=...
                                complexTotPathLen(rsTMP, rdTMP, omega, optProp);
                            [LcwTmp, RcwTmp]=...
                                complexTotPathLen(rsTMP, rdTMP, 0, optProp);
            
                            for mTyp=datTyps
                                % [r_xy, layer, SD, SS, DS, lambda]
                                Lnm=join(['L_' mTyp], '');
                                sens.(Lnm)(1, 1, SDind, SSind, DSind, Lind)=...
                                    calcPathLen_datTyp(Rtmp, Ltmp, RcwTmp, LcwTmp,...
                                    mTyp);

                                eval(sprintf(...
                                    'r_xy=sens.rxy_%s;',...
                                    mTyp));
                                
                                if ~isempty(r_xy)
                                    lTmp_1=zeros(size(r_xy, 1), 1);
                                    lcwTmp_1=zeros(size(r_xy, 1), 1);
                                    lTmp_2=zeros(size(r_xy, 1), 1);
                                    lcwTmp_2=zeros(size(r_xy, 1), 1);
                                    for zInd=1:length(sens.z_all)
                                        
                                        lTmp=complexPartPathLen(...
                                            rsTMP,...
                                            [r_xy, ones(size(r_xy, 1), 1)*sens.z_all(zInd)],...
                                            rdTMP, sens.V, omega, optProp);
                                        lcwTmp=complexPartPathLen(...
                                            rsTMP,...
                                            [r_xy, ones(size(r_xy, 1), 1)*sens.z_all(zInd)],...
                                            rdTMP, sens.V, 0, optProp);
                    
                                        lTmp(isnan(lTmp))=0;
                                        lcwTmp(isnan(lcwTmp))=0;
                    
                                        if ~sens.z_lay2(zInd)
                                            lTmp_1=lTmp_1+lTmp;
                                            lcwTmp_1=lcwTmp_1+lcwTmp;
                                        else
                                            lTmp_2=lTmp_2+lTmp;
                                            lcwTmp_2=lcwTmp_2+lcwTmp;
                                        end
                                    end
                                    
                                    % [r_xy, layer, SD, SS, DS, lambda]
                                    lnm=join(['l_' mTyp], '');
                                    sens.(lnm)(:, 1, SDind, SSind, DSind, Lind)=...
                                        calcPathLen_datTyp(...
                                        Rtmp, lTmp_1, RcwTmp, lcwTmp_1,...
                                        mTyp);
                                    sens.(lnm)(:, 2, SDind, SSind, DSind, Lind)=...
                                        calcPathLen_datTyp(...
                                        Rtmp, lTmp_2, RcwTmp, lcwTmp_2,...
                                        mTyp);
                                else
                                    eval(sprintf(...
                                        'sens.l_%s=[];',...
                                        mTyp));
                                end
                            end
                        end
                    end
                end
            end

        otherwise
            error('Unknown Datatype');
    end

    %% Calc Sen
    for mTyp=datTyps
        lnm=join(['l_' mTyp], '');
        if ~isempty(sens.(lnm))
            switch datStruct.typ
                case 'SD'
                    for SDind=1:size(armt.SDprs, 1)
                        for Lind=1:length(datStruct.lambda)
                            Snm=join(['S_' mTyp], '');
                            lnm=join(['l_' mTyp], '');
                            Lnm=join(['L_' mTyp], '');
                            sens.(Snm)(:, :, SDind, Lind)=...
                                sens.(lnm)(:, :, SDind, Lind)./...
                                sens.(Lnm)(:, :, SDind, Lind);
                        end
                    end

                case 'SS'
                    for SSind=1:size(armt.SSprs, 3)
                        for Lind=1:length(datStruct.lambda)
                            Snm=join(['S_' mTyp], '');
                            lnm=join(['l_' mTyp], '');
                            Lnm=join(['L_' mTyp], '');
                            sens.(Snm)(:, :, SSind, Lind)=...
                                (sens.(lnm)(:, :, 2, SSind, Lind)-...
                                sens.(lnm)(:, :, 1, SSind, Lind))./...
                                ...
                                (sens.(Lnm)(:, :, 2, SSind, Lind)-...
                                sens.(Lnm)(:, :, 1, SSind, Lind));
                        end
                    end

                case 'DS'
                    for DSind=1:size(armt.DSprs, 4)
                        for Lind=1:length(datStruct.lambda)
                            Snm=join(['S_' mTyp], '');
                            lnm=join(['l_' mTyp], '');
                            Lnm=join(['L_' mTyp], '');
                            sens.(Snm)(:, :, DSind, Lind)=...
                                ((sens.(lnm)(:, :, 2, 1, DSind, Lind)-...
                                sens.(lnm)(:, :, 1, 1, DSind, Lind))+...
                                ...
                                (sens.(lnm)(:, :, 2, 2, DSind, Lind)-...
                                sens.(lnm)(:, :, 1, 2, DSind, Lind)))./...
                                ...
                                ...
                                ((sens.(Lnm)(:, :, 2, 1, DSind, Lind)-...
                                sens.(Lnm)(:, :, 1, 1, DSind, Lind))+...
                                ...
                                (sens.(Lnm)(:, :, 2, 2, DSind, Lind)-...
                                sens.(Lnm)(:, :, 1, 2, DSind, Lind)));
                        end
                    end

                otherwise
                    error('Unknown Datatype');
            end
        else
            Snm=join(['S_' mTyp], '');
            sens.(Snm)=[];
        end
    end

    %% Remove Weak Sen
    for mTyp=datTyps
        Snm=join(['S_' mTyp], '');
        lnm=join(['l_' mTyp], '');
        rxynm=join(['rxy_' mTyp], '');

        if ~isempty(sens.(rxynm))
            senAvgTmp=mean(abs(mean(sens.(Snm)(:, 2, :, :), 4)), 3);
    
            wkInds=senAvgTmp<Sthresh;
            
            sens.(rxynm)(wkInds, :)=[];
            sens.(lnm)(wkInds, :, :, :, :, :)=[];
            sens.(Snm)(wkInds, :, :, :)=[];
        end
    end

end