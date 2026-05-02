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
% abp=data.AUX(:, 1);
resp=data.AUX(:, 3);
hr=data.AUX(:, 4);

dtAct=15; %sec
dtPostAct=30; %sec
t0Act=t(diff(act)==1);

fLP=fHPnoise; %Hz

optsFold.DTbools=[true, false];
optsFold.nDT=10;

%% Do Fold
% sigLP=lowpass(detrend(abp), fLP, 1/median(diff(t)),...
%     'ImpulseResponse', 'iir');
% [tFold, abpFold, ~, ~]=foldingAvgASync(t, sigLP,...
%     t0Act, 0, dtAct+dtPostAct, optsFold);

sigLP=lowpass(detrend(resp), fLP, 1/median(diff(t)),...
    'ImpulseResponse', 'iir');
[tFold, respFold, ~, ~]=foldingAvgASync(t, sigLP,...
    t0Act, 0, dtAct+dtPostAct, optsFold);

sigLP=lowpass(detrend(hr), fLP, 1/median(diff(t)),...
    'ImpulseResponse', 'iir');
[~, hrFold, ~, ~]=foldingAvgASync(t, sigLP,...
    t0Act, 0, dtAct+dtPostAct, optsFold);

%SDI
SD.useSetFold_I=all(SD.useSet_I, 2);
for i=1:size(SD.dO_I, 2)
    for Lind=1:length(lambda)
        if ~isnan(sum(SD.dmua_I(:, i, Lind)))
            sigLP=lowpass(SD.dmua_I(:, i, Lind), fLP, 1/median(diff(t)),...
                'ImpulseResponse', 'iir');

            [~, SD.dmuaFold_I(:, i, Lind),...
                ~, SD.dmuadmuaFold_I(:, :, i, Lind)]=...
                foldingAvgASync(t, sigLP,...
                t0Act, 0, dtAct+dtPostAct, optsFold);
        else
            SD.dmuaFold_I(:, i, Lind)=NaN(size(tFold));
            SD.dmuadmuaFold_I(:, :, i, Lind)=...
                NaN(length(tFold), length(t0Act));
        end
    end
    if ~isnan(sum(SD.dO_I(:, i)))
        OsigLP=lowpass(SD.dO_I(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');
        DsigLP=lowpass(SD.dD_I(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');

        [~, SD.Ofold_I(:, i), sigmaTmpO, SD.OOfold_I(:, :, i)]=...
            foldingAvgASync(t, OsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
        [~, SD.Dfold_I(:, i), sigmaTmpD, SD.DDfold_I(:, :, i)]=...
            foldingAvgASync(t, DsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
    else
        SD.Ofold_I(:, i)=NaN(size(tFold));
        SD.Dfold_I(:, i)=NaN(size(tFold));
        
        SD.OOfold_I(:, :, i)=NaN(length(tFold), length(t0Act));
        SD.DDfold_I(:, :, i)=NaN(length(tFold), length(t0Act));
    end
end

%SDP
SD.useSetFold_P=all(SD.useSet_P, 2);
for i=1:size(SD.dO_P, 2)
    for Lind=1:length(lambda)
        if ~isnan(sum(SD.dmua_P(:, i, Lind)))
            sigLP=lowpass(SD.dmua_P(:, i, Lind), fLP, 1/median(diff(t)),...
                'ImpulseResponse', 'iir');

            [~, SD.dmuaFold_P(:, i, Lind),...
                ~, SD.dmuadmuaFold_P(:, :, i, Lind)]=...
                foldingAvgASync(t, sigLP,...
                t0Act, 0, dtAct+dtPostAct, optsFold);
        else
            SD.dmuaFold_P(:, i, Lind)=NaN(size(tFold));
            SD.dmuadmuaFold_P(:, :, i, Lind)=...
                NaN(length(tFold), length(t0Act));
        end
    end
    if ~isnan(sum(SD.dO_P(:, i)))
        OsigLP=lowpass(SD.dO_P(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');
        DsigLP=lowpass(SD.dD_P(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');

        [~, SD.Ofold_P(:, i), sigmaTmpO, SD.OOfold_P(:, :, i)]=...
            foldingAvgASync(t, OsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
        [~, SD.Dfold_P(:, i), sigmaTmpD, SD.DDfold_P(:, :, i)]=...
            foldingAvgASync(t, DsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
    else
        SD.Ofold_P(:, i)=NaN(size(tFold));
        SD.Dfold_P(:, i)=NaN(size(tFold));
        
        SD.OOfold_P(:, :, i)=NaN(length(tFold), length(t0Act));
        SD.DDfold_P(:, :, i)=NaN(length(tFold), length(t0Act));
    end
end

% Find good optodes for SD
SD.useSrcFold_I=false(size(armt.rSrc, 1), 1);
SD.useSrcFold_P=false(size(armt.rSrc, 1), 1);
SD.useDetFold_I=false(size(armt.rDet, 1), 1);
SD.useDetFold_P=false(size(armt.rDet, 1), 1);
for SDind=1:size(armt.SDprs, 1)
    SrcInd=armt.SDprs(SDind, 1);
    DetInd=armt.SDprs(SDind, 2);
    
    if SD.useSetFold_I(SDind)
        SD.useSrcFold_I(SrcInd)=true;
        SD.useDetFold_I(DetInd)=true;
    end
    if SD.useSetFold_P(SDind)
        SD.useSrcFold_P(SrcInd)=true;
        SD.useDetFold_P(DetInd)=true;
    end
end

%SSI
SS.useSetFold_I=all(SS.useSet_I, 2);
for i=1:size(SS.dO_I, 2)
    for Lind=1:length(lambda)
        if ~isnan(sum(SS.dmua_I(:, i, Lind)))
            sigLP=lowpass(SS.dmua_I(:, i, Lind), fLP, 1/median(diff(t)),...
                'ImpulseResponse', 'iir');

            [~, SS.dmuaFold_I(:, i, Lind),...
                ~, SS.dmuadmuaFold_I(:, :, i, Lind)]=...
                foldingAvgASync(t, sigLP,...
                t0Act, 0, dtAct+dtPostAct, optsFold);
        else
            SS.dmuaFold_I(:, i, Lind)=NaN(size(tFold));
            SS.dmuadmuaFold_I(:, :, i, Lind)=...
                NaN(length(tFold), length(t0Act));
        end
    end
    if ~isnan(sum(SS.dO_I(:, i)))
        OsigLP=lowpass(SS.dO_I(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');
        DsigLP=lowpass(SS.dD_I(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');

        [~, SS.Ofold_I(:, i), sigmaTmpO, SS.OOfold_I(:, :, i)]=...
            foldingAvgASync(t, OsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
        [~, SS.Dfold_I(:, i), sigmaTmpD, SS.DDfold_I(:, :, i)]=...
            foldingAvgASync(t, DsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
    else
        SS.Ofold_I(:, i)=NaN(size(tFold));
        SS.Dfold_I(:, i)=NaN(size(tFold));
        
        SS.OOfold_I(:, :, i)=NaN(length(tFold), length(t0Act));
        SS.DDfold_I(:, :, i)=NaN(length(tFold), length(t0Act));
    end
end

%SSP
SS.useSetFold_P=all(SS.useSet_P, 2);
for i=1:size(SS.dO_P, 2)
    for Lind=1:length(lambda)
        if ~isnan(sum(SS.dmua_P(:, i, Lind)))
            sigLP=lowpass(SS.dmua_P(:, i, Lind), fLP, 1/median(diff(t)),...
                'ImpulseResponse', 'iir');

            [~, SS.dmuaFold_P(:, i, Lind),...
                ~, SS.dmuadmuaFold_P(:, :, i, Lind)]=...
                foldingAvgASync(t, sigLP,...
                t0Act, 0, dtAct+dtPostAct, optsFold);
        else
            SS.dmuaFold_P(:, i, Lind)=NaN(size(tFold));
            SS.dmuadmuaFold_P(:, :, i, Lind)=...
                NaN(length(tFold), length(t0Act));
        end
    end
    if ~isnan(sum(SS.dO_P(:, i)))
        OsigLP=lowpass(SS.dO_P(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');
        DsigLP=lowpass(SS.dD_P(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');

        [~, SS.Ofold_P(:, i), sigmaTmpO, SS.OOfold_P(:, :, i)]=...
            foldingAvgASync(t, OsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
        [~, SS.Dfold_P(:, i), sigmaTmpD, SS.DDfold_P(:, :, i)]=...
            foldingAvgASync(t, DsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
    else
        SS.Ofold_P(:, i)=NaN(size(tFold));
        SS.Dfold_P(:, i)=NaN(size(tFold));
        
        SS.OOfold_P(:, :, i)=NaN(length(tFold), length(t0Act));
        SS.DDfold_P(:, :, i)=NaN(length(tFold), length(t0Act));
    end
end

% Find good optodes for SS
SS.useSrcFold_I=false(size(armt.rSrc, 1), 1);
SS.useSrcFold_P=false(size(armt.rSrc, 1), 1);
SS.useDetFold_I=false(size(armt.rDet, 1), 1);
SS.useDetFold_P=false(size(armt.rDet, 1), 1);
for SSind=1:size(armt.SSprs, 3)
    SrcInds=armt.SSprs(:, 1, SSind);
    DetInds=armt.SSprs(:, 2, SSind);
    
    if SS.useSetFold_I(SSind)
        SS.useSrcFold_I(SrcInds)=true;
        SS.useDetFold_I(DetInds)=true;
    end
    if SS.useSetFold_P(SSind)
        SS.useSrcFold_P(SrcInds)=true;
        SS.useDetFold_P(DetInds)=true;
    end
end

%DSI
DS.useSetFold_I=all(DS.useSet_I, 2);
for i=1:size(DS.dO_I, 2)
    for Lind=1:length(lambda)
        if ~isnan(sum(DS.dmua_I(:, i, Lind)))
            sigLP=lowpass(DS.dmua_I(:, i, Lind), fLP, 1/median(diff(t)),...
                'ImpulseResponse', 'iir');

            [~, DS.dmuaFold_I(:, i, Lind),...
                ~, DS.dmuadmuaFold_I(:, :, i, Lind)]=...
                foldingAvgASync(t, sigLP,...
                t0Act, 0, dtAct+dtPostAct, optsFold);
        else
            DS.dmuaFold_I(:, i, Lind)=NaN(size(tFold));
            DS.dmuadmuaFold_I(:, :, i, Lind)=...
                NaN(length(tFold), length(t0Act));
        end
    end
    if ~isnan(sum(DS.dO_I(:, i)))
        OsigLP=lowpass(DS.dO_I(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');
        DsigLP=lowpass(DS.dD_I(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');

        [~, DS.Ofold_I(:, i), sigmaTmpO, DS.OOfold_I(:, :, i)]=...
            foldingAvgASync(t, OsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
        [~, DS.Dfold_I(:, i), sigmaTmpD, DS.DDfold_I(:, :, i)]=...
            foldingAvgASync(t, DsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
    else
        DS.Ofold_I(:, i)=NaN(size(tFold));
        DS.Dfold_I(:, i)=NaN(size(tFold));
        
        DS.OOfold_I(:, :, i)=NaN(length(tFold), length(t0Act));
        DS.DDfold_I(:, :, i)=NaN(length(tFold), length(t0Act));
    end
end

%DSP
DS.useSetFold_P=all(DS.useSet_P, 2);
for i=1:size(DS.dO_P, 2)
    for Lind=1:length(lambda)
        if ~isnan(sum(DS.dmua_P(:, i, Lind)))
            sigLP=lowpass(DS.dmua_P(:, i, Lind), fLP, 1/median(diff(t)),...
                'ImpulseResponse', 'iir');

            [~, DS.dmuaFold_P(:, i, Lind),...
                ~, DS.dmuadmuaFold_P(:, :, i, Lind)]=...
                foldingAvgASync(t, sigLP,...
                t0Act, 0, dtAct+dtPostAct, optsFold);
        else
            DS.dmuaFold_P(:, i, Lind)=NaN(size(tFold));
            DS.dmuadmuaFold_P(:, :, i, Lind)=...
                NaN(length(tFold), length(t0Act));
        end
    end
    if ~isnan(sum(DS.dO_P(:, i)))
        OsigLP=lowpass(DS.dO_P(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');
        DsigLP=lowpass(DS.dD_P(:, i), fLP, 1/median(diff(t)),...
            'ImpulseResponse', 'iir');

        [~, DS.Ofold_P(:, i), sigmaTmpO, DS.OOfold_P(:, :, i)]=...
            foldingAvgASync(t, OsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
        [~, DS.Dfold_P(:, i), sigmaTmpD, DS.DDfold_P(:, :, i)]=...
            foldingAvgASync(t, DsigLP,...
            t0Act, 0, dtAct+dtPostAct, optsFold);
    else
        DS.Ofold_P(:, i)=NaN(size(tFold));
        DS.Dfold_P(:, i)=NaN(size(tFold));
        
        DS.OOfold_P(:, :, i)=NaN(length(tFold), length(t0Act));
        DS.DDfold_P(:, :, i)=NaN(length(tFold), length(t0Act));
    end
end

% Find good optodes for DS
DS.useSrcFold_I=false(size(armt.rSrc, 1), 1);
DS.useSrcFold_P=false(size(armt.rSrc, 1), 1);
DS.useDetFold_I=false(size(armt.rDet, 1), 1);
DS.useDetFold_P=false(size(armt.rDet, 1), 1);
for DSind=1:size(armt.DSprs, 4)
    SSprsTmp=armt.DSprs(:, :, :, DSind);
    SrcIndsTmp=SSprsTmp(:, 1, :);
    DetIndsTmp=SSprsTmp(:, 2, :);
    
    if DS.useSetFold_I(DSind)
        DS.useSrcFold_I(SrcIndsTmp(:))=true;
        DS.useDetFold_I(DetIndsTmp(:))=true;
    end
    if DS.useSetFold_P(DSind)
        DS.useSrcFold_P(SrcIndsTmp(:))=true;
        DS.useDetFold_P(DetIndsTmp(:))=true;
    end
end

%% Save
save([filename '_analOutputA.mat']);

%% Plot Folded AUX
figure(100); clf;
% plot(tFold, (abpFold-mean(abpFold))/(max(abs(abpFold-mean(abpFold))))); 
hold on;
plot(tFold, (respFold-mean(respFold))/(max(abs(respFold-mean(respFold))))); 
plot(tFold, (hrFold-mean(hrFold))/(max(abs(hrFold-mean(hrFold))))); 
yl=ylim;
plot([1, 1]*dtAct, yl, ':k');
ylim(yl);
hold off;
ylabel('(arb.)');
xlabel('t (sec)');
% legend('ABP', 'RESP', 'HR', 'location', 'best');
legend('RESP', 'HR', 'location', 'best');
title('Stim ----------> Rest');

saveFig;

%% Plot Folded Blood
dataTyps={'SDI', 'SDP', 'DSI', 'DSP'};

for dTypInd=1:length(dataTyps)
    eval(sprintf(...
        'inds=find(%s.useSetFold_%s);',...
        dataTyps{dTypInd}(1:2), dataTyps{dTypInd}(3)));
    
    eval(sprintf(...
        'Odata=%s.Ofold_%s;',...
        dataTyps{dTypInd}(1:2), dataTyps{dTypInd}(3)));
    eval(sprintf(...
        'Ddata=%s.Dfold_%s;',...
        dataTyps{dTypInd}(1:2), dataTyps{dTypInd}(3)));
    
    eval(sprintf(...
        'OOdata=%s.OOfold_%s;',...
        dataTyps{dTypInd}(1:2), dataTyps{dTypInd}(3)));
    eval(sprintf(...
        'DDdata=%s.DDfold_%s;',...
        dataTyps{dTypInd}(1:2), dataTyps{dTypInd}(3)));
    
    sp=numSubplots(length(inds));
    figure(100+dTypInd); clf;
    yl=[quantile([Odata(:); Ddata(:)], 0.01),...
        quantile([Odata(:); Ddata(:)], 0.99)];
    for i=1:length(inds)
        subplot(sp(1), sp(2), i); hold on;
        for j=1:size(OOdata, 2)
            plot(tFold, squeeze(OOdata(:, j, inds(i))), '-',...
                'color', [1, 0, 0, 0.05]);
            plot(tFold, squeeze(DDdata(:, j, inds(i))), '-',...
                'color', [0, 0, 1, 0.05]);
        end
        plot(tFold, Odata(:, inds(i)), '-r');
        plot(tFold, Ddata(:, inds(i)), '-b');
        plot([1, 1]*dtAct, yl, ':k');
        ylim(yl);
        xlim([tFold(1), tFold(end)]);
        hold off;
        
        switch dataTyps{dTypInd}(1:2)
            case 'SD'
                prs=armt.SDprs(inds(i), :);
                
                title(sprintf('%d%c',...
                    prs(1), char(prs(2)+64)));
            case 'SS'
                prs=armt.SSprs(:, :, inds(i));
                prsSrc=unique(prs(:, 1));
                prsDet=unique(prs(:, 2));
                
                if length(prsDet)>1
                    title(sprintf('%d%c%c',...
                        prsSrc, char(prsDet(1)+64), char(prsDet(2)+64)));
                else
                    title(sprintf('%d%c%d',...
                        prsSrc(1), char(prsDet+64), prsSrc(2)));
                end
            case 'DS'
                prs=armt.DSprs(:, :, :, inds(i));
                tmp=prs(:, 1, :);
                prsSrc=unique(tmp(:));
                tmp=prs(:, 2, :);
                prsDet=unique(tmp(:));
                
                title(sprintf('%d%c%c%d',...
                    prsSrc(1), char(prsDet(1)+64),...
                    char(prsDet(2)+64), prsSrc(2)));
            otherwise
        end
        
        if i>(prod(sp)-sp(2))
            xlabel('t (sec)');
        else
            ax=gca;
            ax.XTickLabels={};
        end
        if mod(i, sp(2))==1
            ylabel('\DeltaO & \DeltaD (\muM)');
        else
            ax=gca;
            ax.YTickLabels={};
        end
    end
    
    switch dataTyps{dTypInd}(3)
        case 'I'
            titStr='I';
        case 'P'
            titStr='\phi';
        otherwise
            titStr=[];
    end
    sgtitle(sprintf('%s%s', dataTyps{dTypInd}(1:2), titStr));
    
    saveFig;
end