%% Setup
clear; home;

%% Find File
filesTMP=dir('*.set');
if length(filesTMP)>1
    error(['More than one .set file found, '...
        'place only one dataset in same folder']);
end

filename=filesTMP.name(1:(end-4));

load([filename '_analOutputB.mat']);

% a=0.01;
a=1;

%% Recon
for DTind=1:length(dataTyps)
    MTnm=sprintf('%s', dataTyps{DTind}(1:2));
    DTnm=sprintf('%s', dataTyps{DTind}(3));
    
    eval(sprintf(...
        'useSet=%s.useSetFold_%s;',...
        MTnm, DTnm));
    eval(sprintf(...
        'r_xy=sen.%s.rxy_%s;',...
        MTnm, DTnm));
    
    if ~isempty(sen.(MTnm).(DTnm))
        for Lind=1:length(lambda)
            S=[squeeze(sen.(MTnm).(DTnm)(:, 1, useSet, Lind)).',...
                squeeze(sen.(MTnm).(DTnm)(:, 2, useSet, Lind)).'];

            alpha=a*max(diag(S*S'));
            Sp=S'*inv(S*S'+alpha*eye(size(S, 1)));

            eval(sprintf(...
                'dmua_all=%s.dmuaFold_%s(:, useSet, Lind);',...
                MTnm, DTnm));
            eval(sprintf(...
                'dmuadmua_all=%s.dmuadmuaFold_%s(:, :, useSet, Lind);',...
                MTnm, DTnm));

            dmua_recon=NaN(size(r_xy, 1), 2, size(dmua_all, 1));
            dmuadmua_recon=...
                NaN(size(r_xy, 1), 2,...
                size(dmuadmua_all, 1), size(dmuadmua_all, 2));
            for i=1:size(dmua_all, 1)
                dmua=dmua_all(i, :).';

                if sum(isnan(dmua))==0
                    dmua_tmp1=Sp*dmua;
                    dmua_tmp2=[dmua_tmp1(1:(size(Sp, 1)/2)),...
                        dmua_tmp1((size(Sp, 1)/2+1):end)];

                    dmua_recon(:, :, i)=dmua_tmp2;
                elseif sum(isnan(dmua))<(length(dmua)-1)
                    DSinds_noNaN=~isnan(dmua);

                    S_noNaN=S(DSinds_noNaN, :);
                    alpha_noNaN=a*max(diag(S_noNaN*S_noNaN'));
                    Sp_noNaN=S_noNaN'*...
                        inv(S_noNaN*S_noNaN'+alpha*eye(size(S_noNaN, 1)));

                    dmua_tmp1=Sp_noNaN*dmua(DSinds_noNaN);
                    dmua_tmp2=[dmua_tmp1(1:(size(Sp_noNaN, 1)/2)),...
                        dmua_tmp1((size(Sp_noNaN, 1)/2+1):end)];

                    dmua_recon(:, :, i)=dmua_tmp2;
                else
                    dmua_recon(:, :, i)=NaN(size(r_xy, 1), 2);
                end
                
                for j=1:size(dmuadmua_all, 2)
                    dmua=squeeze(dmuadmua_all(i, j, :));
                    
                    if sum(isnan(dmua))==0
                        dmua_tmp1=Sp*dmua;
                        dmua_tmp2=[dmua_tmp1(1:(size(Sp, 1)/2)),...
                            dmua_tmp1((size(Sp, 1)/2+1):end)];

                        dmuadmua_recon(:, :, i, j)=dmua_tmp2;
                    elseif sum(isnan(dmua))<(length(dmua)-1)
                        DSinds_noNaN=~isnan(dmua);

                        S_noNaN=S(DSinds_noNaN, :);
                        alpha_noNaN=a*max(diag(S_noNaN*S_noNaN'));
                        Sp_noNaN=S_noNaN'*...
                            inv(S_noNaN*S_noNaN'+alpha*eye(size(S_noNaN, 1)));

                        dmua_tmp1=Sp_noNaN*dmua(DSinds_noNaN);
                        dmua_tmp2=[dmua_tmp1(1:(size(Sp_noNaN, 1)/2)),...
                            dmua_tmp1((size(Sp_noNaN, 1)/2+1):end)];

                        dmuadmua_recon(:, :, i, j)=dmua_tmp2;
                    else
                        dmuadmua_recon(:, :, i, j)=NaN(size(r_xy, 1), 2);
                    end
                end
                
            end

            recon.(MTnm).(DTnm).r_xy=r_xy;
            % [lambda, r_xy, layer, t]
            recon.(MTnm).(DTnm).dmua(Lind, :, :, :)=dmua_recon;
            % [lambda, r_xy, layer, t, nEvent]
            recon.(MTnm).(DTnm).dmuadmua(Lind, :, :, :, :)=dmuadmua_recon;
        end

        dmuaTMP=recon.(MTnm).(DTnm).dmua;
        dmuadmuaTMP=recon.(MTnm).(DTnm).dmuadmua;
        clear tmp tmptmp;
        for layInd=1:size(dmuaTMP, 3)
            for i=1:size(dmuaTMP, 4)
                tmp(:, :, layInd, i)=E\dmuaTMP(:, :, layInd, i);
                
                for j=1:size(dmuadmuaTMP, 5)
                    tmptmp(:, :, layInd, i, j)=E\dmuaTMP(:, :, layInd, i);
                end
            end
        end

        % [r_xy, layer, t]
        recon.(MTnm).(DTnm).dO=squeeze(tmp(1, :, :, :));
        recon.(MTnm).(DTnm).dD=squeeze(tmp(2, :, :, :));
        recon.(MTnm).(DTnm).dOdO=squeeze(tmptmp(1, :, :, :, :));
        recon.(MTnm).(DTnm).dDdD=squeeze(tmptmp(2, :, :, :, :));
        clear dmuaTMP dmuadmuaTMP tmp tmptmp;
    else
        recon.(MTnm).(DTnm)=[];
    end
end

%% Save
save([filename '_analOutputB.mat'], '-v7.3');