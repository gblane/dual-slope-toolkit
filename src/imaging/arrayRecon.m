function [recn] = arrayRecon(dstr, sens, NVA)
% arrayRecon Performs image reconstruction for NIRS array data.
%
% [recn] = arrayRecon(dstr, sens, NVA)
%
% Written by Giles Blaney, Ph.D. (Winter 2023)
%
% Inputs:
%   dstr - Struct containing NIRS data (dmua, etc.)
%   sens - Struct containing sensitivity maps (Jacobians)
%
% Optional Name-Value-Arguments (NVA):
%   datTyps - String array of data types to reconstruct (Default: ["C", "I", "P"])
%   dmuaNm  - Name of the dmua field to use ('Fold', or '') (Default: 'Fold')
%   a       - Regularization parameter multiplier (Default: 1)
%
% Outputs:
%   recn - Struct containing the reconstructed 2D maps of dmua and hemodynamics.
%
% Shared-repo dependencies:
%   makeE is provided by ../dos-inverse-models.

    %% Parse Input
    arguments
        dstr struct;
        sens struct;

        NVA.datTyps = ["C", "I", "P"];

        NVA.dmuaNm string ...
            {mustBeMember(NVA.dmuaNm,{'Fold',''})} = 'Fold';
        NVA.a = 1;
    end
    
    E=makeE('OD', dstr.lambda); % ../dos-inverse-models
    
    %% dmua Recon
    for Lind=1:length(dstr.lambda)
        for mTypsCom=NVA.datTyps
            try
                mTyps=split(mTypsCom, "_").';
                
                % Pull Meas
                tmpNm=join(['dmua' NVA.dmuaNm '_' mTyps(1)], '');
                tmpNum=size(dstr.(tmpNm), 2);
                dmuaM=NaN(tmpNum*length(mTyps), size(dstr.(tmpNm), 1));
                useSets=false(tmpNum*length(mTyps), 1);
                for i=1:length(mTyps)
                    i0=i-1;
                    tmpNm=join(['dmua' NVA.dmuaNm '_' mTyps(i)], '');
                    
                    dmuaM((i0*tmpNum+1):(i*tmpNum), :)=...
                        dstr.(tmpNm)(:, :, Lind).';
                    
                    tmpNm=join(['useSet_' mTyps(i)], '');
                    useSets((i0*tmpNum+1):(i*tmpNum))=...
                        dstr.(tmpNm)(:, Lind);
                end
                useSets(any(isnan(dmuaM), 2))=false;
    
                % Pull rxy
                rxy=[];
                for mTyp=mTyps
                    tmpNm=join(['rxy_' mTyp], '');
                    if strcmp(mTyp, mTyps(1))
                        rxy=sens.(tmpNm);
                    else
                        kp=false(size(rxy, 1), 1);
                        for i=1:size(sens.(tmpNm), 1)
                            kp=or(kp, all(sens.(tmpNm)(i, :)==rxy, 2));
                        end
                        rxy(~kp, :)=[];
                    end
                end
                tmpNm=join(['rxy_' mTypsCom], '');
                recn.(tmpNm)=rxy;
    
                % Pull sens
                tmpNm=join(['dmua' NVA.dmuaNm '_' mTyps(1)], '');
                tmpNum=size(dstr.(tmpNm), 2);
                S=NaN(size(dmuaM, 1), size(rxy, 1)*2);
                for i=1:length(mTyps)
                    i0=i-1;
                    mTyp=mTyps(i);
                    
                    tmpNm=join(['rxy_' mTyp], '');
                    inds=NaN(size(rxy, 1), 1);
                    for j=1:size(rxy, 1)
                        inds(j)=find(all(sens.(tmpNm)==rxy(j, :), 2));
                    end
                    
                    tmpNm=join(['S_' mTyp], '');
                    S((i0*tmpNum+1):(i*tmpNum), :)=...
                        [squeeze(sens.(tmpNm)(inds, 1, :, Lind)).',...
                        squeeze(sens.(tmpNm)(inds, 2, :, Lind)).'];
                end
    
                % Recon
                Suse=S(useSets, :);
                alpha=NVA.a*max(diag(Suse*Suse'));
                dmuaV=Suse'*((Suse*Suse' + alpha*eye(size(Suse, 1)))\...
                    dmuaM(useSets, :));
                
                tmpNm=join(['dmua' NVA.dmuaNm '_' mTypsCom], '');
                recn.(tmpNm)(:, :, 1, Lind)=...
                    dmuaV(1:(size(dmuaV, 1)/2), :).';
                recn.(tmpNm)(:, :, 2, Lind)=...
                    dmuaV((size(dmuaV, 1)/2+1):end, :).';

            catch ME
                warning('arrayRecon:DatatypeSkipped', ...
                    'Skipping reconstruction for %s: %s', ...
                    char(mTypsCom), ME.message);
                tmpNm=join(['rxy_' mTypsCom], '');
                recn.(tmpNm)=[];
                tmpNm=join(['dmua' NVA.dmuaNm '_' mTypsCom], '');
                recn.(tmpNm)=[];
            end
        end
        
    end

    %% Convert to Blood
    for mTypsCom=NVA.datTyps
        tmpNm=join(['dmua' NVA.dmuaNm '_' mTypsCom], '');
        tmpNmO=join(['dO' lower(NVA.dmuaNm) '_' mTypsCom], '');
        tmpNmD=join(['dD' lower(NVA.dmuaNm) '_' mTypsCom], '');

        if ~isempty(recn.(tmpNm))
            for tInd=1:size(recn.(tmpNm), 1)
                for lInd=1:size(recn.(tmpNm), 3)
                    dmuaTmp=squeeze(recn.(tmpNm)(tInd, :, lInd, :)).';
                    dCtmp=E\dmuaTmp;
    
                    recn.(tmpNmO)(tInd, :, lInd)=dCtmp(1, :);
                    recn.(tmpNmD)(tInd, :, lInd)=dCtmp(2, :);
                end
            end
        else
            recn.(tmpNmO)=[];
            recn.(tmpNmD)=[];
        end
    end
end
