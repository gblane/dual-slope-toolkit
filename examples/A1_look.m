%% Setup
clear; home;

svFig=true;

rowStart=302;
swNum=16;

dtNoise=10; %sec
fHPnoise=100/60; %Hz
% TnoiseThresh=inf; %uM turn off noise removal
TnoiseThresh=1*sqrt(11); %uM

%% Find File
filesTMP=dir('Armt*.mat');
if length(filesTMP)>1
    error(['More than one Armt file found, '...
        'place only one in same folder']);
end
load(filesTMP.name);

filesTMP=dir('*.set');
if length(filesTMP)>1
    error(['More than one .set file found, '...
        'place only one dataset in same folder']);
end

filename=filesTMP.name(1:(end-4));

%% Load
if ~exist([filename '_data.mat'], 'file')
    [data, setups]=load_Imagent([filename '.txt'], rowStart, swNum, true);
    save([filename '_data.mat'], 'data', 'setups');
else
    load([filename '_data.mat']);
end
fs=setups.Fs;
t=timeAxis(fs, length(data.timemat));
lambda=ISSmap.lambda;

%Fix AUX
act=(data.AUX(:, 2)>150);
BLinds=(data.AUX(:, 1)<-100);

nNoise=dtNoise*fs;

%% Calculate Noise Thresholds
% Assume mua0=0.01 1/mm, musp0=1 1/mm, rho=35 mm, lambda=760 nm
Text=sum(makeE('OD', lambda), 2); %1/(mm uM)
muaNoiseThresh=TnoiseThresh*Text; %1/mm

% From https://doi.org/10.1002/jbio.201960018 appendix
L_Cass=35*7.50; %mm
L_Iass=35*7.17; %mm
L_Pass=35*1.22; %mm
DSF_Iass=7.90;
DSF_Pass=1.50;

CnoiseThresh=muaNoiseThresh*L_Cass;
InoiseThresh=muaNoiseThresh*L_Iass;
PnoiseThresh=muaNoiseThresh*L_Pass;

SInoiseThresh=muaNoiseThresh*DSF_Iass;
SPnoiseThresh=muaNoiseThresh*DSF_Pass;

%% Pull SD Sets
SD.useSet_C=true(size(armt.SDprs, 1), length(lambda));
SD.useSet_I=true(size(armt.SDprs, 1), length(lambda));
SD.useSet_P=true(size(armt.SDprs, 1), length(lambda));
for SDind=1:size(armt.SDprs, 1)
    SrcInd=armt.SDprs(SDind, 1);
    DetInd=armt.SDprs(SDind, 2);
    
    DetISSnam=char(DetInd+64);
    
    for Lind=1:length(lambda)
        SrcISSind=ISSmap.inds(Lind, SrcInd);
        
        % [t, SDind, lambda]
        Ctmp=data.(DetISSnam).DC(:, SrcISSind);
        if median(movstd(...
                highpass(...
                (Ctmp(BLinds)-mean(Ctmp(BLinds)))/mean(Ctmp(BLinds)),...
                fHPnoise, fs),...
                nNoise))>CnoiseThresh(Lind)
            Ctmp=NaN(size(Ctmp));
            SD.useSet_C(SDind, Lind)=false;
        end
        SD.C(:, SDind, Lind)=Ctmp;
        
        Itmp=data.(DetISSnam).AC(:, SrcISSind);
        if median(movstd(...
                highpass(...
                (Itmp(BLinds)-mean(Itmp(BLinds)))/mean(Itmp(BLinds)),...
                fHPnoise, fs),...
                nNoise))>InoiseThresh(Lind)
            Itmp=NaN(size(Itmp));
            SD.useSet_I(SDind, Lind)=false;
        end
        SD.I(:, SDind, Lind)=Itmp;
        
        Ptmp=wrapTo2Pi(data.(DetISSnam).Ph(:, SrcISSind)*pi/180);
        if max(diff(Ptmp))>pi
            Ptmp=wrapToPi(Ptmp);
        end
        indsTemp=find(diff(Ptmp)>pi/16);
        PtrshInds=unique([indsTemp; indsTemp+1]);
        PtrshInds(PtrshInds>length(t))=[];
        Ptmp(PtrshInds)=NaN;
        Ptmp=interp1(t(~isnan(Ptmp)), Ptmp(~isnan(Ptmp)), t,...
            'linear', 'extrap');
        if median(movstd(...
                highpass(...
                wrapToPi(Ptmp(BLinds)-circ_mean(Ptmp(BLinds))),...
                fHPnoise, fs),...
                nNoise))>PnoiseThresh(Lind)
            Ptmp=NaN(size(Ptmp));
            SD.useSet_P(SDind, Lind)=false;
        end
        SD.P(:, SDind, Lind)=Ptmp;
    end
    
    SD.rho(1, SDind)=vecnorm(armt.rSrc(SrcInd, :)-armt.rDet(DetInd, :));
    SD.loc(:, SDind)=mean([armt.rSrc(SrcInd, :); armt.rDet(DetInd, :)]);
end

% Find good optodes
SD.useSrc_C=false(size(armt.rSrc, 1), length(lambda));
SD.useSrc_I=false(size(armt.rSrc, 1), length(lambda));
SD.useSrc_P=false(size(armt.rSrc, 1), length(lambda));
SD.useDet_C=false(size(armt.rDet, 1), length(lambda));
SD.useDet_I=false(size(armt.rDet, 1), length(lambda));
SD.useDet_P=false(size(armt.rDet, 1), length(lambda));
for SDind=1:size(armt.SDprs, 1)
    SrcInd=armt.SDprs(SDind, 1);
    DetInd=armt.SDprs(SDind, 2);
    
    for Lind=1:length(lambda)
        if SD.useSet_C(SDind, Lind)
            SD.useSrc_C(SrcInd, Lind)=true;
            SD.useDet_C(DetInd, Lind)=true;
        end
        if SD.useSet_I(SDind, Lind)
            SD.useSrc_I(SrcInd, Lind)=true;
            SD.useDet_I(DetInd, Lind)=true;
        end
        if SD.useSet_P(SDind, Lind)
            SD.useSrc_P(SrcInd, Lind)=true;
            SD.useDet_P(DetInd, Lind)=true;
        end
    end
end

%% Pull SS Sets
SS.useSet_I=true(size(armt.SSprs, 3), length(lambda));
SS.useSet_P=true(size(armt.SSprs, 3), length(lambda));
for SSind=1:size(armt.SSprs, 3)
    SrcInds=armt.SSprs(:, 1, SSind);
    DetInds=armt.SSprs(:, 2, SSind);
    
    for Lind=1:length(lambda)
        rhosTmp=NaN(size(SrcInds));
        locsTmp=NaN(length(SrcInds), 3);
        Itemp=NaN(length(t), length(SrcInds));
        Ptemp=NaN(length(t), length(SrcInds));
        for SDind=1:length(SrcInds)
            SrcInd=SrcInds(SDind);
            DetInd=DetInds(SDind);
        
            DetISSnam=char(DetInd+64);
        
            rhosTmp(SDind)=vecnorm(armt.rSrc(SrcInd, :)-...
                armt.rDet(DetInd, :));
            locsTmp(SDind, :)=mean([armt.rSrc(SrcInd, :);...
                armt.rDet(DetInd, :)]);
        
            SrcISSind=ISSmap.inds(Lind, SrcInd);
    
            % [t, SDind]
            Itemp(:, SDind)=data.(DetISSnam).AC(:, SrcISSind);

            Ptmp=wrapTo2Pi(data.(DetISSnam).Ph(:, SrcISSind)*pi/180);
            if max(diff(Ptmp))>pi
                Ptmp=wrapToPi(Ptmp);
            end
            indsTemp=find(diff(Ptmp)>pi/8);
            PtrshInds=unique([indsTemp; indsTemp+1]);
            PtrshInds(PtrshInds>length(t))=[];
            Ptmp(PtrshInds)=NaN;
            Ptmp=interp1(t(~isnan(Ptmp)), Ptmp(~isnan(Ptmp)), t,...
                'linear', 'extrap');
            Ptemp(:, SDind)=Ptmp;
        end
        
        SS.rhos(1:2, SSind)=rhosTmp;
        SS.drho(1, SSind)=rhosTmp(2)-rhosTmp(1);
        SS.loc(1:3, SSind)=mean(locsTmp);
        
        % [t, SSind, lambda]
        SS.I(:, SSind, Lind)=...
            (log(rhosTmp(2)^2*Itemp(:, 2))-log(rhosTmp(1)^2*Itemp(:, 1)))/...
            (rhosTmp(2)-rhosTmp(1));
        SS.P(:, SSind, Lind)=...
            wrapToPi(Ptemp(:, 2)-Ptemp(:, 1))/...
            (rhosTmp(2)-rhosTmp(1));
        
        Itmp=SS.I(:, SSind, Lind);
        if median(movstd(...
                highpass(...
                Itmp(BLinds),...
                fHPnoise, fs),...
                nNoise))>SInoiseThresh(Lind)
            Itmp=NaN(size(Itmp));
            SS.useSet_I(SSind, Lind)=false;
        end
        SS.I(:, SSind, Lind)=Itmp;
        
        Ptmp=SS.P(:, SSind, Lind);
        if median(movstd(...
                highpass(...
                Ptmp(BLinds),...
                fHPnoise, fs),...
                nNoise))>SPnoiseThresh(Lind)
            Ptmp=NaN(size(Ptmp));
            SS.useSet_P(SSind, Lind)=false;
        end
        SS.P(:, SSind, Lind)=Ptmp;
    end
end

% Find good optodes
SS.useSrc_I=false(size(armt.rSrc, 1), length(lambda));
SS.useSrc_P=false(size(armt.rSrc, 1), length(lambda));
SS.useDet_I=false(size(armt.rDet, 1), length(lambda));
SS.useDet_P=false(size(armt.rDet, 1), length(lambda));
for SSind=1:size(armt.SSprs, 3)
    SrcInds=armt.SSprs(:, 1, SSind);
    DetInds=armt.SSprs(:, 2, SSind);
    
    for Lind=1:length(lambda)
        if SS.useSet_I(SSind, Lind)
            SS.useSrc_I(SrcInds, Lind)=true;
            SS.useDet_I(DetInds, Lind)=true;
        end
        if SS.useSet_P(SSind, Lind)
            SS.useSrc_P(SrcInds, Lind)=true;
            SS.useDet_P(DetInds, Lind)=true;
        end
    end
end

%% Pull DS Sets
DS.useSet_I=true(size(armt.DSprs, 4), length(lambda));
DS.useSet_P=true(size(armt.DSprs, 4), length(lambda));
for DSind=1:size(armt.DSprs, 4)
    SSprsTmp=armt.DSprs(:, :, :, DSind);
    
    for Lind=1:length(lambda)
        SSrhosTmp=NaN(2, size(SSprsTmp, 3));
        SSdrhoTmp=NaN(size(SSprsTmp, 3), 1);
        SSlocTmp=NaN(size(SSprsTmp, 3), 3);
        SSItmp=NaN(length(t), size(SSprsTmp, 3));
        SSPtmp=NaN(length(t), size(SSprsTmp, 3));
        R=NaN(2, 2);
        for SSind=1:size(SSprsTmp, 3)
            SrcInds=SSprsTmp(:, 1, SSind);
            DetInds=SSprsTmp(:, 2, SSind);
            
            rhosTmp=NaN(size(SrcInds));
            locsTmp=NaN(length(SrcInds), 3);
            for SDind=1:length(SrcInds)
                SrcInd=SrcInds(SDind);
                DetInd=DetInds(SDind);

                DetISSnam=char(DetInd+64);

                rhosTmp(SDind)=vecnorm(armt.rSrc(SrcInd, :)-...
                    armt.rDet(DetInd, :));
                locsTmp(SDind, :)=mean([armt.rSrc(SrcInd, :);...
                    armt.rDet(DetInd, :)]);

                SrcISSind=ISSmap.inds(Lind, SrcInd);

                % [t, SDind]
                Itemp(:, SDind)=data.(DetISSnam).AC(:, SrcISSind);

                Ptmp=wrapTo2Pi(data.(DetISSnam).Ph(:, SrcISSind)*pi/180);
                if max(diff(Ptmp))>pi
                    Ptmp=wrapToPi(Ptmp);
                end
                indsTemp=find(diff(Ptmp)>pi/8);
                PtrshInds=unique([indsTemp; indsTemp+1]);
                PtrshInds(PtrshInds>length(t))=[];
                Ptmp(PtrshInds)=NaN;
                Ptmp=interp1(t(~isnan(Ptmp)), Ptmp(~isnan(Ptmp)), t,...
                    'linear', 'extrap');
                Ptemp(:, SDind)=Ptmp;
                
                R(SSind, SDind)=mean(Itemp(BLinds, SDind))*exp(1i*...
                    wrapToPi(circ_mean(Ptemp(BLinds, SDind))));
            end
            
            SSrhosTmp(1:2, SSind)=rhosTmp;
            SSdrhoTmp(SSind, 1)=rhosTmp(2)-rhosTmp(1);
            SSlocTmp(SSind, 1:3)=mean(locsTmp);

            % [t, SSind]
            SSItmp(:, SSind)=...
                (log(rhosTmp(2)^2*Itemp(:, 2))-log(rhosTmp(1)^2*Itemp(:, 1)))/...
                (rhosTmp(2)-rhosTmp(1));
            SSPtmp(:, SSind)=...
                wrapToPi(Ptemp(:, 2)-Ptemp(:, 1))/...
                (rhosTmp(2)-rhosTmp(1));
        end
        
        DS.rhos(1:4, DSind)=[SSrhosTmp(1:2, 1); SSrhosTmp(1:2, 2)];
        DS.drhos(1:2, DSind)=SSdrhoTmp;
        DS.loc(1:3, DSind)=mean(SSlocTmp);
        
        % [t, DSind, lambda]
        DS.I(:, DSind, Lind)=mean(SSItmp, 2);
        DS.P(:, DSind, Lind)=mean(SSPtmp, 2);
        
        Itmp=DS.I(:, DSind, Lind);
        if median(movstd(...
                highpass(...
                Itmp(BLinds),...
                fHPnoise, fs),...
                nNoise))>SInoiseThresh(Lind)
            Itmp=NaN(size(Itmp));
            DS.useSet_I(DSind, Lind)=false;
        end
        DS.I(:, DSind, Lind)=Itmp;
        
        Ptmp=DS.P(:, DSind, Lind);
        if median(movstd(...
                highpass(...
                Ptmp(BLinds),...
                fHPnoise, fs),...
                nNoise))>SPnoiseThresh(Lind)
            Ptmp=NaN(size(Ptmp));
            DS.useSet_P(DSind, Lind)=false;
        end
        DS.P(:, DSind, Lind)=Ptmp;
        
        if ~DS.useSet_I(DSind, Lind) && ~DS.useSet_P(DSind, Lind)
            DS.mua(DSind, Lind)=NaN;
            DS.musp(DSind, Lind)=NaN;
        else
            RRtmp=[R(1, :), R(2, :)];
            rhosTmp=[SSrhosTmp(1:2, 1); SSrhosTmp(1:2, 2)]';
            [DS.mua(DSind, Lind), DS.musp(DSind, Lind), iterTMP(DSind, Lind)]=...
                DSR2muamuspEB_iterRecov(rhosTmp, RRtmp);
        end
        
        if DS.musp(DSind, Lind)<0 || DS.mua(DSind, Lind)<0
            DS.musp(DSind, Lind)=NaN;
            DS.mua(DSind, Lind)=NaN;
        end
    end
end

% Find good optodes
DS.useSrc_I=false(size(armt.rSrc, 1), length(lambda));
DS.useSrc_P=false(size(armt.rSrc, 1), length(lambda));
DS.useDet_I=false(size(armt.rDet, 1), length(lambda));
DS.useDet_P=false(size(armt.rDet, 1), length(lambda));
for DSind=1:size(armt.DSprs, 4)
    SSprsTmp=armt.DSprs(:, :, :, DSind);
    
    for Lind=1:length(lambda)
        SrcIndsTmp=SSprsTmp(:, 1, :);
        DetIndsTmp=SSprsTmp(:, 2, :);
        if DS.useSet_I(DSind, Lind)
            DS.useSrc_I(SrcIndsTmp(:), Lind)=true;
            DS.useDet_I(DetIndsTmp(:), Lind)=true;
        end
        if DS.useSet_P(DSind, Lind)
            DS.useSrc_P(SrcIndsTmp(:), Lind)=true;
            DS.useDet_P(DetIndsTmp(:), Lind)=true;
        end
    end
end

%% Save
save([filename '_analOutputA.mat']);

%% Plot SD I
figure(1); clf;
cols=turbo(length(SD.rho)*2);
for i=1:length(SD.rho)
    sp=(SD.rho(i)>=30)+1;
    
    for j=1:length(lambda)
        ci=i+length(SD.rho)*(j-1);
        
        subplot(2, 1, sp);
        semilogy(t/(60*60), SD.C(:, i, j), '-', 'color', cols(ci, :));
        hold on;
    end
end
for i=1:2
    subplot(2, 1, i); hold off;
    xlim([t(1), t(end)]/(60*60));
    
    xlabel('t (hr)');
    ylabel('I_{CW} (counts)');
    
    tmp=sum(SD.C);
    
    switch i
        case 1
            title(sprintf('\\rho < 30 mm\n%d/%d NaNs',...
                sum(isnan(tmp(:))), numel(tmp)));
        case 2
            title(sprintf('\\rho > 30 mm\n%d/%d NaNs',...
                sum(isnan(tmp(:))), numel(tmp)));
        otherwise
    end
end

sgtitle('All SD Sets (Both wavelengths)');
saveFig;


figure(2); clf;
cols=turbo(length(SD.rho)*2);
for i=1:length(SD.rho)
    sp=(SD.rho(i)>=30)+1;
    
    for j=1:length(lambda)
        ci=i+length(SD.rho)*(j-1);
        
        subplot(2, 1, sp);
        semilogy(t/(60*60), SD.I(:, i, j), '-', 'color', cols(ci, :)); hold on;
    end
end
for i=1:2
    subplot(2, 1, i); hold off;
    xlim([t(1), t(end)]/(60*60));
    
    xlabel('t (hr)');
    ylabel('I (counts)');
    
    tmp=sum(SD.I);
    
    switch i
        case 1
            title(sprintf('\\rho < 30 mm\n%d/%d NaNs',...
                sum(isnan(tmp(:))), numel(tmp)));
        case 2
            title(sprintf('\\rho > 30 mm\n%d/%d NaNs',...
                sum(isnan(tmp(:))), numel(tmp)));
        otherwise
    end
end

sgtitle('All SD Sets (Both wavelengths)');
saveFig;

%% Plot SD P
figure(3); clf;
cols=turbo(length(SD.rho)*2);
for i=1:length(SD.rho)
    sp=(SD.rho(i)>=30)+1;
    
    for j=1:length(lambda)
        ci=i+length(SD.rho)*(j-1);
        
        subplot(2, 1, sp);
        plot(t/(60*60), SD.P(:, i, j), '-', 'color', cols(ci, :)); hold on;
    end
end
for i=1:2
    subplot(2, 1, i); hold off;
    xlim([t(1), t(end)]/(60*60));
    
    xlabel('t (hr)');
    ylabel('\phi (rad)');
    
    tmp=sum(SD.P);
    
    switch i
        case 1
            title(sprintf('\\rho < 30 mm\n%d/%d NaNs',...
                sum(isnan(tmp(:))), numel(tmp)));
        case 2
            title(sprintf('\\rho > 30 mm\n%d/%d NaNs',...
                sum(isnan(tmp(:))), numel(tmp)));
        otherwise
    end
end

sgtitle('All SD Sets (Both wavelengths)');
saveFig;

%% Plot SD dI
figure(4); clf;
cols=turbo(length(SD.rho)*2);
for i=1:length(SD.rho)
    sp=(SD.rho(i)>=30)+1;
    
    for j=1:length(lambda)
        ci=i+length(SD.rho)*(j-1);
        
        subplot(2, 1, sp);
        I0=mean(SD.I(1:100, i, j));
        plot(t/(60*60), (SD.I(:, i, j)-I0)/I0, '-', 'color', cols(ci, :)); hold on;
    end
end
for i=1:2
    subplot(2, 1, i); hold off;
    xlim([t(1), t(end)]/(60*60));
    
    xlabel('t (hr)');
    ylabel('\DeltaI/I_0');
    
    tmp=sum(SD.I);
    
    switch i
        case 1
            title(sprintf('\\rho < 30 mm\n%d/%d NaNs',...
                sum(isnan(tmp(:))), numel(tmp)));
        case 2
            title(sprintf('\\rho > 30 mm\n%d/%d NaNs',...
                sum(isnan(tmp(:))), numel(tmp)));
        otherwise
    end
end

sgtitle('All SD Sets (Both wavelengths)');
saveFig;

%% Plot SD dP
figure(5); clf;
cols=turbo(length(SD.rho)*2);
for i=1:length(SD.rho)
    sp=(SD.rho(i)>=30)+1;
    
    for j=1:length(lambda)
        ci=i+length(SD.rho)*(j-1);
        
        subplot(2, 1, sp);
        P0=circ_mean(SD.P(1:100, i, j));
        plot(t/(60*60), wrapToPi(SD.P(:, i, j)-P0), '-', 'color', cols(ci, :)); hold on;
    end
end
for i=1:2
    subplot(2, 1, i); hold off;
    xlim([t(1), t(end)]/(60*60));
    
    xlabel('t (hr)');
    ylabel('\Delta\phi (rad)');
    
    tmp=sum(SD.P);
    
    switch i
        case 1
            title(sprintf('\\rho < 30 mm\n%d/%d NaNs',...
                sum(isnan(tmp(:))), numel(tmp)));
        case 2
            title(sprintf('\\rho > 30 mm\n%d/%d NaNs',...
                sum(isnan(tmp(:))), numel(tmp)));
        otherwise
    end
end

sgtitle('All SD Sets (Both wavelengths)');
saveFig;

%% Plot DS
figure(6); clf;
cols=turbo(size(DS.rhos, 2)*2);
for i=1:size(DS.rhos, 2)
    for j=1:length(lambda)
        ci=i+size(DS.rhos, 2)*(j-1);
        
        subplot(2, 1, 1);
        plot(t/(60*60), DS.I(:, i, j), '-', 'color', cols(ci, :)); hold on;
    end
end; hold off;
xlim([t(1), t(end)]/(60*60));
xlabel('t (hr)');
ylabel('DSI (1/mm)');
tmp=sum(DS.I);
title(sprintf('%d/%d NaNs',...
    sum(isnan(tmp(:))), numel(tmp)));

for i=1:size(DS.rhos, 2)
    for j=1:length(lambda)
        ci=i+size(DS.rhos, 2)*(j-1);
        
        subplot(2, 1, 2);
        plot(t/(60*60), DS.P(:, i, j), '-', 'color', cols(ci, :)); hold on;
    end
end; hold off;
xlim([t(1), t(end)]/(60*60));
xlabel('t (hr)');
ylabel('DS\phi (1/mm)');
tmp=sum(DS.P);
title(sprintf('%d/%d NaNs',...
    sum(isnan(tmp(:))), numel(tmp)));

sgtitle('All DS Sets (Both wavelengths)');
saveFig;

%% Plot Opt Props
figure(7); clf;

subplot(2, 2, 1);
lamInd=1;
useSetTmp=~isnan(DS.mua(:, lamInd));
plotAbsOptProp(DS, 'mua', useSetTmp, lamInd, 25);
axis equal;
xlim([min(DS.loc(1, :))-25, max(DS.loc(1, :))+25]);
ylim([min(DS.loc(2, :))-25, max(DS.loc(2, :))+25]);
shading flat;
cb=colorbar;
view(0, 90.0);
ylabel(cb, '\mu_a (1/mm)');
xlabel('x (mm)');
ylabel('y (mm)');
title('830 nm');

subplot(2, 2, 3);
lamInd=1;
useSetTmp=~isnan(DS.musp(:, lamInd));
plotAbsOptProp(DS, 'musp', useSetTmp, lamInd, 25);
axis equal;
xlim([min(DS.loc(1, :)), max(DS.loc(1, :))]);
xlim([min(DS.loc(1, :))-25, max(DS.loc(1, :))+25]);
ylim([min(DS.loc(2, :))-25, max(DS.loc(2, :))+25]);
shading flat;
cb=colorbar;
view(0, 90.0);
ylabel(cb, '\mu_s^'' (1/mm)');
xlabel('x (mm)');
ylabel('y (mm)');
title('830 nm');

subplot(2, 2, 2);
lamInd=2;
useSetTmp=~isnan(DS.mua(:, lamInd));
plotAbsOptProp(DS, 'mua', useSetTmp, lamInd, 25);
axis equal;
xlim([min(DS.loc(1, :))-25, max(DS.loc(1, :))+25]);
ylim([min(DS.loc(2, :))-25, max(DS.loc(2, :))+25]);
shading flat;
cb=colorbar;
view(0, 90.0);
ylabel(cb, '\mu_a (1/mm)');
xlabel('x (mm)');
ylabel('y (mm)');
title('690 nm');

subplot(2, 2, 4);
lamInd=2;
useSetTmp=~isnan(DS.musp(:, lamInd));
plotAbsOptProp(DS, 'musp', useSetTmp, lamInd, 25);
axis equal;
xlim([min(DS.loc(1, :))-25, max(DS.loc(1, :))+25]);
ylim([min(DS.loc(2, :))-25, max(DS.loc(2, :))+25]);
shading flat;
cb=colorbar;
view(0, 90.0);
ylabel(cb, '\mu_s^'' (1/mm)');
xlabel('x (mm)');
ylabel('y (mm)');
title('690 nm');

saveFig;

%% Functions
function h=plotAbsOptProp(DS, opPropNm, useSet, lamInd, blobRad)

    mu=DS.(opPropNm)(useSet, lamInd);
    xCent=DS.loc(1, useSet);
    yCent=DS.loc(2, useSet);
    
    x=floor(min(xCent)-blobRad):ceil(max(xCent)+blobRad);
    y=floor(min(yCent)-blobRad):ceil(max(yCent)+blobRad);
    
    [XX, YY]=meshgrid(x, y);
    
    mu_stack=NaN([size(XX), length(mu)]);
    
    for i=1:length(mu)
        muLay=NaN(size(XX));
        inds=sqrt((XX(:)-xCent(i)).^2+(YY(:)-yCent(i)).^2)<=blobRad;
        
        muLay(inds)=mu(i);
        
        mu_stack(:, :, i)=muLay;
        
    end
    
    muaMap=nanmean(mu_stack, 3);
    
    h=pcolor(x, y, muaMap); shading flat; hold on;
    
    for i=1:length(mu)
        if mu(i)>=mean(mu)
            col='k';
        else
            col='w';
        end
        
%         text(xCent(i)-blobRad/2, yCent(i),...
%             sprintf('%.1e',  mu(i)), 'color', col); hold on;
    end
end
