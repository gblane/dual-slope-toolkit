%% Setup
clear; home;

%% Find File
filesTMP=dir('*.set');
if length(filesTMP)>1
    error(['More than one .set file found, '...
        'place only one dataset in same folder']);
end

filename=filesTMP.name(1:(end-4));

load([filename '_analOutputA.mat']);

filesTMP=dir('Armt*.mat');
if length(filesTMP)>1
    error(['More than one Armt file found, '...
        'place only one in same folder']);
end
load(filesTMP.name);

z_min=0; %mm
z_max=30; %mm
z_lay=5; %mm

%% Make sens map r
dataTyps={'SDI', 'SDP', 'SSI', 'SSP', 'DSI', 'DSP'};

haloRad=10; %mm
dr=[1, 1, 1]; %mm
sen.V=prod(dr); %mm^3

sen.z_all=z_min:dr(3):z_max;
sen.z_lay2=sen.z_all>=z_lay;

for i=1:length(dataTyps)
    datGeo=dataTyps{i}(1:2);
    datTyp=dataTyps{i}(3);
    
    eval(sprintf('useSrc=%s.useSrcFold_%s;',...
        datGeo, datTyp));
    eval(sprintf('useDet=%s.useDetFold_%s;',...
        datGeo, datTyp));
    
    AllOpt=[armt.rSrc(useSrc, :); armt.rDet(useDet, :)];
    if ~isempty(AllOpt)
        shp=alphaShape(AllOpt(:, 1), AllOpt(:, 2), inf);

        x_min=min([round(AllOpt(:, 1), -1); round(AllOpt(:, 1), -1)-10])-haloRad;
        x_max=max([round(AllOpt(:, 1), -1); round(AllOpt(:, 1), -1)+10])+haloRad;
        y_min=min([round(AllOpt(:, 2), -1); round(AllOpt(:, 2), -1)-10])-haloRad;
        y_max=max([round(AllOpt(:, 2), -1); round(AllOpt(:, 2), -1)+10])+haloRad;

        [XX, YY]=meshgrid(x_min:dr(1):x_max, y_min:dr(2):y_max);

        r_xy=[XX(:), YY(:)];
        inds=or(min(sqrt((AllOpt(:, 1).'-r_xy(:, 1)).^2+...
            (AllOpt(:, 2).'-r_xy(:, 2)).^2), [], 2)<=haloRad,...
            inShape(shp, r_xy(:, 1), r_xy(:, 2)));
        
        eval(sprintf('sen.%s.rxy_%s=r_xy(inds, :);',...
            datGeo, datTyp));
    else
        eval(sprintf('sen.%s.rxy_%s=[];',...
            datGeo, datTyp));
    end
end

%% Pathlengths Setup
datTyps={'I', 'P'};

optProp.nin=1.4;
optProp.nout=1;

fmod=140.625e6; %Hz
omega=2*pi*fmod;

%% SD Pathlengths
for Lind=1:length(lambda)
    for SDind=1:size(SD.rho, 2)
        SrcInd=armt.SDprs(SDind, 1);
        DetInd=armt.SDprs(SDind, 2);
        
        DSlocDist=vecnorm(DS.loc-SD.loc(:, SDind), 2, 1);
        [~, DSinds]=sort(DSlocDist);

        optProp.mua=nanmean(DS.mua(DSinds(1:DSn_closeOptPropAvg), Lind));
        optProp.musp=nanmean(DS.musp(DSinds(1:DSn_closeOptPropAvg), Lind));
        if any(isnan(optProp.mua)) || any(isnan(optProp.musp))
            optProp.mua=nanmean(DS.mua(:, Lind));
            optProp.musp=nanmean(DS.musp(:, Lind));
        end
        
        rsTMP=armt.rSrc(SrcInd, :)+[0, 0, 1/optProp.musp];
        rdTMP=armt.rDet(DetInd, :);

        % [r_xy, layer, SD, lambda]
        [sen.SD.L(1, 1, SDind, Lind),...
            sen.SD.R(1, 1, SDind, Lind)]=...
            complexTotPathLen(rsTMP, rdTMP, omega, optProp);
        [sen.SD.Lcw(1, 1, SDind, Lind),...
            sen.SD.C(1, 1, SDind, Lind)]=...
            complexTotPathLen(rsTMP, rdTMP, 0, optProp);
        
        for DTind=1:length(datTyps)
            eval(sprintf(...
                'r_xy=sen.%s.rxy_%s;',...
                'SD', datTyps{DTind}));
            
            if ~isempty(r_xy)
                lTmp_1=zeros(size(r_xy, 1), 1);
                lcwTmp_1=zeros(size(r_xy, 1), 1);
                lTmp_2=zeros(size(r_xy, 1), 1);
                lcwTmp_2=zeros(size(r_xy, 1), 1);
                for zInd=1:length(sen.z_all)

                    lTmp=complexPartPathLen(...
                        rsTMP,...
                        [r_xy, ones(size(r_xy, 1), 1)*sen.z_all(zInd)],...
                        rdTMP, sen.V, omega, optProp);
                    lcwTmp=complexPartPathLen(...
                        rsTMP,...
                        [r_xy, ones(size(r_xy, 1), 1)*sen.z_all(zInd)],...
                        rdTMP, sen.V, 0, optProp);

                    lTmp(isnan(lTmp))=0;
                    lcwTmp(isnan(lcwTmp))=0;

                    if ~sen.z_lay2(zInd)
                        lTmp_1=lTmp_1+lTmp;
                        lcwTmp_1=lcwTmp_1+lcwTmp;
                    else
                        lTmp_2=lTmp_2+lTmp;
                        lcwTmp_2=lcwTmp_2+lcwTmp;
                    end
                end
                
                % [r_xy, layer, SD, lambda]
                eval(sprintf(...
                    'sen.%s.l_%s(:, 1, SDind, Lind)=lTmp_1;',...
                    'SD', datTyps{DTind}));
                eval(sprintf(...
                    'sen.%s.lcw_%s(:, 1, SDind, Lind)=lcwTmp_1;',...
                    'SD', datTyps{DTind}));
                eval(sprintf(...
                    'sen.%s.l_%s(:, 2, SDind, Lind)=lTmp_2;',...
                    'SD', datTyps{DTind}));
                eval(sprintf(...
                    'sen.%s.lcw_%s(:, 2, SDind, Lind)=lcwTmp_2;',...
                    'SD', datTyps{DTind}));
            else
                % [r_xy, layer, SD, lambda]
                eval(sprintf(...
                    'sen.%s.l_%s=[];',...
                    'SD', datTyps{DTind}));
                eval(sprintf(...
                    'sen.%s.lcw_%s=[];',...
                    'SD', datTyps{DTind}));
                eval(sprintf(...
                    'sen.%s.l_%s=[];',...
                    'SD', datTyps{DTind}));
                eval(sprintf(...
                    'sen.%s.lcw_%s=[];',...
                    'SD', datTyps{DTind}));
            end
        end
    end
end

%% SS Pathlengths
for Lind=1:length(lambda)
    for SSind=1:size(SS.rhos, 2)
        SrcInds=armt.SSprs(:, 1, SSind);
        DetInds=armt.SSprs(:, 2, SSind);
        
        DSlocDist=vecnorm(DS.loc-SS.loc(:, SSind), 2, 1);
        [~, DSinds]=sort(DSlocDist);
        
        optProp.mua=nanmean(DS.mua(DSinds(1:DSn_closeOptPropAvg), Lind));
        optProp.musp=nanmean(DS.musp(DSinds(1:DSn_closeOptPropAvg), Lind));
        if any(isnan(optProp.mua)) || any(isnan(optProp.musp))
            optProp.mua=nanmean(DS.mua(:, Lind));
            optProp.musp=nanmean(DS.musp(:, Lind));
        end
        
        for SDind=1:length(SrcInds)
            SrcInd=SrcInds(SDind);
            DetInd=DetInds(SDind);
            
            rsTMP=armt.rSrc(SrcInd, :)+[0, 0, 1/optProp.musp];
            rdTMP=armt.rDet(DetInd, :);
            
            % [r_xy, layer, SD, SS, lambda]
            [sen.SS.L(1, 1, SDind, SSind, Lind),...
                sen.SS.R(1, 1, SDind, SSind, Lind)]=...
                complexTotPathLen(rsTMP, rdTMP, omega, optProp);
            
            for DTind=1:length(datTyps)
                eval(sprintf(...
                    'r_xy=sen.%s.rxy_%s;',...
                    'SS', datTyps{DTind}));
                
                if ~isempty(r_xy)
                    lTmp_1=zeros(size(r_xy, 1), 1);
                    lTmp_2=zeros(size(r_xy, 1), 1);
                    for zInd=1:length(sen.z_all)

                        lTmp=complexPartPathLen(...
                            rsTMP,...
                            [r_xy, ones(size(r_xy, 1), 1)*sen.z_all(zInd)],...
                            rdTMP, sen.V, omega, optProp);

                        lTmp(isnan(lTmp))=0;

                        if ~sen.z_lay2(zInd)
                            lTmp_1=lTmp_1+lTmp;
                        else
                            lTmp_2=lTmp_2+lTmp;
                        end
                    end
                    
                    % [r_xy, layer, SD, SS, lambda]
                    eval(sprintf(...
                        'sen.%s.l_%s(:, 1, SDind, SSind, Lind)=lTmp_1;',...
                        'SS', datTyps{DTind}));
                    eval(sprintf(...
                        'sen.%s.l_%s(:, 2, SDind, SSind, Lind)=lTmp_2;',...
                        'SS', datTyps{DTind}));
                else
                    % [r_xy, layer, SD, SS, lambda]
                    eval(sprintf(...
                        'sen.%s.l_%s=[];',...
                        'SS', datTyps{DTind}));
                    eval(sprintf(...
                        'sen.%s.l_%s=[];',...
                        'SS', datTyps{DTind}));
                end
            end
        end
    end
end

%% DS Pathlengths
for Lind=1:length(lambda)
    for DSind=1:size(DS.rhos, 2)
        SSprsTmp=armt.DSprs(:, :, :, DSind);
        
        DSlocDist=vecnorm(DS.loc-DS.loc(:, DSind), 2, 1);
        [~, DSinds]=sort(DSlocDist);

        optProp.mua=nanmean(DS.mua(DSinds(1:DSn_closeOptPropAvg), Lind));
        optProp.musp=nanmean(DS.musp(DSinds(1:DSn_closeOptPropAvg), Lind));
        if any(isnan(optProp.mua)) || any(isnan(optProp.musp))
            optProp.mua=nanmean(DS.mua(:, Lind));
            optProp.musp=nanmean(DS.musp(:, Lind));
        end
        
        for SSind=1:size(SSprsTmp, 3)
            SrcInds=SSprsTmp(:, 1, SSind);
            DetInds=SSprsTmp(:, 2, SSind);
        
            for SDind=1:length(SrcInds)
                SrcInd=SrcInds(SDind);
                DetInd=DetInds(SDind);

                rsTMP=armt.rSrc(SrcInd, :)+[0, 0, 1/optProp.musp];
                rdTMP=armt.rDet(DetInd, :);

                % [r_xy, layer, SD, SS, DS, lambda]
                [sen.DS.L(1, 1, SDind, SSind, DSind, Lind),...
                    sen.DS.R(1, 1, SDind, SSind, DSind, Lind)]=...
                    complexTotPathLen(rsTMP, rdTMP, omega, optProp);

                for DTind=1:length(datTyps)
                    eval(sprintf(...
                        'r_xy=sen.%s.rxy_%s;',...
                        'DS', datTyps{DTind}));
                    
                    if ~isempty(r_xy)
                        lTmp_1=zeros(size(r_xy, 1), 1);
                        lTmp_2=zeros(size(r_xy, 1), 1);
                        for zInd=1:length(sen.z_all)

                            lTmp=complexPartPathLen(...
                                rsTMP,...
                                [r_xy, ones(size(r_xy, 1), 1)*sen.z_all(zInd)],...
                                rdTMP, sen.V, omega, optProp);

                            lTmp(isnan(lTmp))=0;

                            if ~sen.z_lay2(zInd)
                                lTmp_1=lTmp_1+lTmp;
                            else
                                lTmp_2=lTmp_2+lTmp;
                            end
                        end
                        
                        % [r_xy, layer, SD, SS, DS, lambda]
                        eval(sprintf(...
                            'sen.%s.l_%s(:, 1, SDind, SSind, DSind, Lind)=lTmp_1;',...
                            'DS', datTyps{DTind}));
                        eval(sprintf(...
                            'sen.%s.l_%s(:, 2, SDind, SSind, DSind, Lind)=lTmp_2;',...
                            'DS', datTyps{DTind}));
                    else
                        % [r_xy, layer, SD, SS, DS, lambda]
                        eval(sprintf(...
                            'sen.%s.l_%s=[];',...
                            'DS', datTyps{DTind}));
                        eval(sprintf(...
                            'sen.%s.l_%s=[];',...
                            'DS', datTyps{DTind}));
                    end
                end
            end
        end
    end
end

%% Calc Sen
for DTind=1:length(datTyps)
    DTnm=sprintf('%s', datTyps{DTind});
    
    lnm=sprintf('l_%s', datTyps{DTind});
    lcwnm=sprintf('lcw_%s', datTyps{DTind});
    Lnm='L';
    Lcwnm='Lcw';
    
    % SD
    % [r_xy, layer, SD, lambda]
    Ltmp=calcPathLen_datTyp(...
        sen.SD.R, sen.SD.(Lnm), sen.SD.C, sen.SD.(Lcwnm), datTyps{DTind});
    lTmp=calcPathLen_datTyp(...
        sen.SD.R, sen.SD.(lnm), sen.SD.C, sen.SD.(lcwnm), datTyps{DTind});
    if ~isempty(lTmp)
        for SDind=1:size(armt.SDprs, 1)
            for Lind=1:size(lambda)

                SDsen=lTmp(:, :, SDind, Lind)./Ltmp(:, :, SDind, Lind);

                sen.SD.(DTnm)(:, :, SDind, Lind)=SDsen;
            end
        end
    else
        sen.SD.(DTnm)=[];
    end
    
    % SS
    % [r_xy, layer, SD, SS, lambda]
    Ltmp=calcPathLen_datTyp(...
        sen.SS.R, sen.SS.(Lnm), [], [], datTyps{DTind});
    lTmp=calcPathLen_datTyp(...
        sen.SS.R, sen.SS.(lnm), [], [], datTyps{DTind});
    if ~isempty(lTmp)
        for SSind=1:size(armt.SSprs, 3)
            for Lind=1:size(lambda)

                SSsen=...
                    (lTmp(:, :, 2, SSind, Lind)-lTmp(:, :, 1, SSind, Lind))./...
                    (Ltmp(:, :, 2, SSind, Lind)-Ltmp(:, :, 1, SSind, Lind));

                sen.SS.(DTnm)(:, :, SSind, Lind)=SSsen;
            end
        end
    else
        sen.SS.(DTnm)=[];
    end
    
    % DS
    % [r_xy, layer, SD, SS, DS, lambda]
    Ltmp=calcPathLen_datTyp(...
        sen.DS.R, sen.DS.(Lnm), [], [], datTyps{DTind});
    lTmp=calcPathLen_datTyp(...
        sen.DS.R, sen.DS.(lnm), [], [], datTyps{DTind});
    if ~isempty(lTmp)
        for DSind=1:size(armt.DSprs, 4)
            for Lind=1:size(lambda)

                SSsen1=...
                    (lTmp(:, :, 2, 1, DSind, Lind)-lTmp(:, :, 1, 1, DSind, Lind))./...
                    (Ltmp(:, :, 2, 1, DSind, Lind)-Ltmp(:, :, 1, 1, DSind, Lind));
                SSsen2=...
                    (lTmp(:, :, 2, 2, DSind, Lind)-lTmp(:, :, 1, 2, DSind, Lind))./...
                    (Ltmp(:, :, 2, 2, DSind, Lind)-Ltmp(:, :, 1, 2, DSind, Lind));

                DSsen=(SSsen1+SSsen2)/2;

                sen.DS.(DTnm)(:, :, DSind, Lind)=DSsen;
            end
        end
    else
        sen.DS.(DTnm)=[];
    end
end

%% Save
save([filename '_analOutputB.mat'], '-v7.3');

%% Functions
function LY = calcPathLen_datTyp(R, L, Rcw, Lcw, datTyp)
    switch datTyp
        case 'I' % ln(r^2 |R|)
            LY=real(L);
            
        case 'P' % angle(R)
            LY=imag(L);
            
        case 'Re' % ln(r^2 Re(R))
            LY=real(L)-(imag(R)./real(R)).*imag(L);
            
        case 'Im' % ln(r^2 Im(R))
            LY=real(L)+(real(R)./imag(R)).*imag(L);
            
        case 'ReN' % Re(R/Rcw)
            L_Re=calcPathLen_datTyp(R, L, Rcw, Lcw, 'Re');
            LY=real(R./Rcw).*(L_Re-Lcw);
            
        case 'ImN' % Im(R/Rcw)
            L_Im=calcPathLen_datTyp(R, L, Rcw, Lcw, 'Im');
            LY=imag(R./Rcw).*(L_Im-Lcw);
            
        case 'ReNpImN' % Re(R/DC)+Im(R/DC)
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            L_ImN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ImN');
            LY=L_ReN+L_ImN;
            
        case 'ImNmReN' % Im(R/DC)-Re(R/DC)
            L_ImN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ImN');
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            LY=L_ImN-L_ReN;
            
        case 'ReNpP' % Re(R/Rcw)+phi
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            LY=L_ReN+L_P;
            
        case 'ReNmP' % Re(R/Rcw)-phi
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            LY=L_ReN-L_P;
            
        case 'PmReN' % phi-Re(R/Rcw)
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            LY=L_P-L_ReN;
            
        case 'ImNpP' % Im(R/Rcw)+phi
            L_ImN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ImN');
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            LY=L_ImN+L_P;
            
        case 'ImNmP' % Im(R/Rcw)-phi
            L_ImN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ImN');
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            LY=L_ImN-L_P;
            
        case 'PmImN' % phi-Im(R/Rcw)
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            L_ImN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ImN');
            LY=L_P-L_ImN;
            
        case 'ReNpP2' % Re(R/Rcw)+phi^2/2
            L_ReN=calcPathLen_datTyp(R, L, Rcw, Lcw, 'ReN');
            L_P=calcPathLen_datTyp(R, L, Rcw, Lcw, 'P');
            LY=L_ReN+angle(R).*L_P;
            
        otherwise
            warning('Unknown data type');
            LY=[];
    end
end

function Y = calcData_datTyp(rho, R, Rcw, datTyp)
    switch datTyp
        case 'I' % ln(r^2 |R|)
            Y=log(rho.^2.*abs(R));
            
        case 'P' % angle(R)
            Y=angle(R);
            
        case 'Re' % ln(r^2 Re(R))
            Y=log(rho.^2.*real(R));
            
        case 'Im' % ln(r^2 Im(R))
            Y=log(rho.^2.*imag(R));
            
        case 'ReN' % Re(R/Rcw);
            Y=real(R./Rcw);
            
        case 'ImN' % Im(R/Rcw)
            Y=imag(R./Rcw);
            
        case 'ReNpImN' % Re(R/DC)+Im(R/DC)
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            Y=ReN+ImN;
            
        case 'ImNmReN' % Im(R/DC)-Re(R/DC)
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            Y=ImN-ReN;
            
        case 'ReNpP' % Re(R/Rcw)+phi
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=ReN+P;
            
        case 'ReNmP' % Re(R/Rcw)-phi
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=ReN-P;
            
        case 'PmReN' % phi-Re(R/Rcw)
            P=calcData_datTyp(rho, R, Rcw, 'P');
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            Y=P-ReN;
            
        case 'ImNpP' % Im(R/Rcw)+phi
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=ImN+P;
            
        case 'ImNmP' % Im(R/Rcw)-phi
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=ImN-P;
            
        case 'PmImN' % phi-Im(R/Rcw)
            P=calcData_datTyp(rho, R, Rcw, 'P');
            ImN=calcData_datTyp(rho, R, Rcw, 'ImN');
            Y=P-ImN;
            
        case 'ReNpP2' % Re(R/Rcw)+phi^2/2
            ReN=calcData_datTyp(rho, R, Rcw, 'ReN');
            P=calcData_datTyp(rho, R, Rcw, 'P');
            Y=ReN+P.^2/2;
            
        otherwise
            warning('Unknown data type');
            Y=[];
    end
end