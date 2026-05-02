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

E=makeE('OD', lambda);

DSn_closeOptPropAvg=2;

%% Find DPF, dmua, and dblood SD
for SDind=1:size(SD.rho, 2)
    DSlocDist=vecnorm(DS.loc-SD.loc(:, SDind), 2, 1);
    [~, DSinds]=sort(DSlocDist);
    
    muaTmp=nanmean(DS.mua(DSinds(1:DSn_closeOptPropAvg), :), 1);
    muspTmp=nanmean(DS.musp(DSinds(1:DSn_closeOptPropAvg), :), 1);
    if any(isnan(muaTmp)) || any(isnan(muspTmp))
        muaTmp=nanmean(DS.mua, 1);
        muspTmp=nanmean(DS.musp, 1);
        
        badLam_mua=isnan(muaTmp);
        if sum(badLam_mua)==2
            muaTmp(:)=0.01; %1/mm
        elseif sum(badLam_mua)==1
            muaTmp(badLam_mua)=muaTmp(~badLam_mua);
        end
        
        badLam_musp=isnan(muspTmp);
        if sum(badLam_musp)==2
            muspTmp(:)=1; %1/mm
        elseif sum(badLam_musp)==1
            muspTmp(badLam_musp)=muspTmp(~badLam_musp);
        end
    end
    
    SD.DPF_I(SDind, :)=DPF_DSF_calc(...
        [0, 0, 0], [SD.rho(SDind), 0, 0], muaTmp, muspTmp, 'DPF_I');
    SD.DPF_P(SDind, :)=DPF_DSF_calc(...
        [0, 0, 0], [SD.rho(SDind), 0, 0], muaTmp, muspTmp, 'DPF_Ph');
    
    for lInd=1:length(lambda)
        SD.I0(SDind, lInd)=nanmean(SD.I(BLinds, SDind, lInd));
        
        inds=~isnan(SD.P(:, SDind, lInd));
        SD.P0(SDind, lInd)=wrapToPi(circ_mean(...
            SD.P(and(inds, BLinds), SDind, lInd)));
        
        SD.dmua_I(:, SDind, lInd)=...
            -(SD.I(:, SDind, lInd)-SD.I0(SDind, lInd))/...
            (SD.I0(SDind, lInd)*SD.rho(SDind)*SD.DPF_I(SDind, lInd));
        SD.dmua_P(:, SDind, lInd)=...
            -wrapToPi(SD.P(:, SDind, lInd)-SD.P0(SDind, lInd))/...
            (SD.rho(SDind)*SD.DPF_P(SDind, lInd));
    end
    
    muaTmp=squeeze(SD.dmua_I(:, SDind, :)).';
    tmp=(E\muaTmp).';
    SD.dO_I(:, SDind)=tmp(:, 1);
    SD.dD_I(:, SDind)=tmp(:, 2);
    
    muaTmp=squeeze(SD.dmua_P(:, SDind, :)).';
    tmp=(E\muaTmp).';
    SD.dO_P(:, SDind)=tmp(:, 1);
    SD.dD_P(:, SDind)=tmp(:, 2);
end

%% Find DSF, dmua, and dblood SS
for SSind=1:size(SS.rhos, 2)
    DSlocDist=vecnorm(DS.loc-SS.loc(:, SSind), 2, 1);
    [~, DSinds]=sort(DSlocDist);
    
    muaTmp=nanmean(DS.mua(DSinds(1:DSn_closeOptPropAvg), :), 1);
    muspTmp=nanmean(DS.musp(DSinds(1:DSn_closeOptPropAvg), :), 1);
    if any(isnan(muaTmp)) || any(isnan(muspTmp))
        muaTmp=nanmean(DS.mua, 1);
        muspTmp=nanmean(DS.musp, 1);
    end
    
    rdTmp=[SS.rhos(1:2, SSind), zeros(2, 2)];
    SS.DSF_I(SSind, :)=DPF_DSF_calc(...
        [0, 0, 0; 0, 0, 0], rdTmp, muaTmp, muspTmp, 'DSF_I');
    
    rdTmp=[SS.rhos(1:2, SSind), zeros(2, 2)];
    SS.DSF_P(SSind, :)=DPF_DSF_calc(...
        [0, 0, 0; 0, 0, 0], rdTmp, muaTmp, muspTmp, 'DSF_Ph');
    
    for lInd=1:length(lambda)
        SS.I0(SSind, lInd)=nanmean(SS.I(BLinds, SSind, lInd));
        SS.P0(SSind, lInd)=nanmean(SS.P(BLinds, SSind, lInd));
        
        SS.dmua_I(:, SSind, lInd)=...
            -(SS.I(:, SSind, lInd)-SS.I0(SSind, lInd))/...
            (SS.DSF_I(SSind, lInd));
        SS.dmua_P(:, SSind, lInd)=...
            -(SS.P(:, SSind, lInd)-SS.P0(SSind, lInd))/...
            (SS.DSF_P(SSind, lInd));
    end
    
    muaTmp=squeeze(SS.dmua_I(:, SSind, :)).';
    tmp=(E\muaTmp).';
    SS.dO_I(:, SSind)=tmp(:, 1);
    SS.dD_I(:, SSind)=tmp(:, 2);
    
    muaTmp=squeeze(SS.dmua_P(:, SSind, :)).';
    tmp=(E\muaTmp).';
    SS.dO_P(:, SSind)=tmp(:, 1);
    SS.dD_P(:, SSind)=tmp(:, 2);
end

%% Find DSF, dmua, and dblood DS
for DSind=1:size(DS.rhos, 2)
    DSlocDist=vecnorm(DS.loc-DS.loc(:, DSind), 2, 1);
    [~, DSinds]=sort(DSlocDist);
    
    muaTmp=nanmean(DS.mua(DSinds(1:DSn_closeOptPropAvg), :), 1);
    muspTmp=nanmean(DS.musp(DSinds(1:DSn_closeOptPropAvg), :), 1);
    if any(isnan(muaTmp)) || any(isnan(muspTmp))
        muaTmp=nanmean(DS.mua, 1);
        muspTmp=nanmean(DS.musp, 1);
    end
    
    DS.mua_avg(DSind, :)=muaTmp;
    DS.musp_avg(DSind, :)=muspTmp;
    
    rdTmp=[DS.rhos(1:2, DSind), zeros(2, 2)];
    DSF_tmp(1, :)=DPF_DSF_calc(...
        [0, 0, 0; 0, 0, 0], rdTmp, muaTmp, muspTmp, 'DSF_I');
    rdTmp=[DS.rhos(3:4, DSind), zeros(2, 2)];
    DSF_tmp(2, :)=DPF_DSF_calc(...
        [0, 0, 0; 0, 0, 0], rdTmp, muaTmp, muspTmp, 'DSF_I');
    DS.DSF_I(DSind, :)=mean(DSF_tmp, 1);
    
    rdTmp=[DS.rhos(1:2, DSind), zeros(2, 2)];
    DSF_tmp(1, :)=DPF_DSF_calc(...
        [0, 0, 0; 0, 0, 0], rdTmp, muaTmp, muspTmp, 'DSF_Ph');
    rdTmp=[DS.rhos(3:4, DSind), zeros(2, 2)];
    DSF_tmp(2, :)=DPF_DSF_calc(...
        [0, 0, 0; 0, 0, 0], rdTmp, muaTmp, muspTmp, 'DSF_Ph');
    DS.DSF_P(DSind, :)=mean(DSF_tmp, 1);
    
    for lInd=1:length(lambda)
        DS.I0(DSind, lInd)=nanmean(DS.I(BLinds, DSind, lInd));
        DS.P0(DSind, lInd)=nanmean(DS.P(BLinds, DSind, lInd));
        
        DS.dmua_I(:, DSind, lInd)=...
            -(DS.I(:, DSind, lInd)-DS.I0(DSind, lInd))/...
            (DS.DSF_I(DSind, lInd));
        DS.dmua_P(:, DSind, lInd)=...
            -(DS.P(:, DSind, lInd)-DS.P0(DSind, lInd))/...
            (DS.DSF_P(DSind, lInd));
    end
    
    muaTmp=squeeze(DS.dmua_I(:, DSind, :)).';
    tmp=(E\muaTmp).';
    DS.dO_I(:, DSind)=tmp(:, 1);
    DS.dD_I(:, DSind)=tmp(:, 2);
    
    muaTmp=squeeze(DS.dmua_P(:, DSind, :)).';
    tmp=(E\muaTmp).';
    DS.dO_P(:, DSind)=tmp(:, 1);
    DS.dD_P(:, DSind)=tmp(:, 2);
end

%% Save
save([filename '_analOutputA.mat']);

%% Plot SD
figure(10); clf;
subplot(2, 1, 1);
cols=turbo(length(SD.rho)*2);
for i=1:length(SD.rho)
    plot(t/(60*60), SD.dO_I(:, i), '-', 'color', cols(i+length(SD.rho), :));
    hold on;
end
xlim([t(1), t(end)]/(60*60));
xlabel('t (hr)');
ylabel('\DeltaO (\muM)');

subplot(2, 1, 2);
cols=turbo(length(SD.rho)*2);
for i=1:length(SD.rho)
    plot(t/(60*60), SD.dD_I(:, i), '-', 'color', cols(i, :)); hold on;
end
xlim([t(1), t(end)]/(60*60));
xlabel('t (hr)');
ylabel('\DeltaD (\muM)');

sgtitle('All SDI Sets');
saveFig;


figure(11); clf;
subplot(2, 1, 1);
cols=turbo(length(SD.rho)*2);
for i=1:length(SD.rho)
    plot(t/(60*60), SD.dO_P(:, i), '-', 'color', cols(i+length(SD.rho), :));
    hold on;
end
xlim([t(1), t(end)]/(60*60));
xlabel('t (hr)');
ylabel('\DeltaO (\muM)');

subplot(2, 1, 2);
cols=turbo(length(SD.rho)*2);
for i=1:length(SD.rho)
    plot(t/(60*60), SD.dD_P(:, i), '-', 'color', cols(i, :)); hold on;
end
xlim([t(1), t(end)]/(60*60));
xlabel('t (hr)');
ylabel('\DeltaD (\muM)');

sgtitle('All SD\phi Sets');
saveFig;


%% Plot DS
figure(12); clf;
subplot(2, 1, 1);
cols=turbo(size(DS.rhos, 2)*2);
for i=1:size(DS.rhos, 2)
    plot(t/(60*60), DS.dO_I(:, i), '-', 'color', cols(i+size(DS.rhos, 2), :));
    hold on;
end
xlim([t(1), t(end)]/(60*60));
xlabel('t (hr)');
ylabel('\DeltaO (\muM)');

subplot(2, 1, 2);
cols=turbo(size(DS.rhos, 2)*2);
for i=1:size(DS.rhos, 2)
    plot(t/(60*60), DS.dD_I(:, i), '-', 'color', cols(i, :)); hold on;
end
xlim([t(1), t(end)]/(60*60));
xlabel('t (hr)');
ylabel('\DeltaD (\muM)');

sgtitle('All DSI Sets');
saveFig;


figure(13); clf;
subplot(2, 1, 1);
cols=turbo(size(DS.rhos, 2)*2);
for i=1:size(DS.rhos, 2)
    plot(t/(60*60), DS.dO_P(:, i), '-', 'color', cols(i+size(DS.rhos, 2), :));
    hold on;
end
xlim([t(1), t(end)]/(60*60));
xlabel('t (hr)');
ylabel('\DeltaO (\muM)');

subplot(2, 1, 2);
cols=turbo(size(DS.rhos, 2)*2);
for i=1:size(DS.rhos, 2)
    plot(t/(60*60), DS.dD_P(:, i), '-', 'color', cols(i, :)); hold on;
end
xlim([t(1), t(end)]/(60*60));
xlabel('t (hr)');
ylabel('\DeltaD (\muM)');

sgtitle('All DS\phi Sets');
saveFig;