function [armt] = findGoodOpts(datStruct, armt)
% findGoodOpts Updates the array struct with usable source/detector pairs.
%
% [armt] = findGoodOpts(datStruct, armt)
%
% Written by Giles Blaney, Ph.D. (Winter 2022)
% Modified by Cristianne Fernandez (August 2023)
%
% Inputs:
%   datStruct - Struct containing NIRS data (SD, SS, or DS format).
%               Must contain 'useSet' variables from rmBadChans.m.
%   armt      - Struct containing the array arrangement and optode positions.
%
% Outputs:
%   armt      - Updated array struct with useSrc/useDet logical maps.

    %% Parse Input
    arguments
        datStruct struct;
        armt struct;
    end

    datTyp=datStruct.typ;
    lamN=size(datStruct.I, 3);

    %% Run
    switch datTyp
        case 'SD'
            armt.SDuseSrc_C=false(size(armt.rSrc, 1), lamN);
            armt.SDuseSrc_I=false(size(armt.rSrc, 1), lamN);
            armt.SDuseSrc_P=false(size(armt.rSrc, 1), lamN);
            armt.SDuseDet_C=false(size(armt.rDet, 1), lamN);
            armt.SDuseDet_I=false(size(armt.rDet, 1), lamN);
            armt.SDuseDet_P=false(size(armt.rDet, 1), lamN);
            for SDind=1:size(armt.SDprs, 1)
                SrcInd=armt.SDprs(SDind, 1);
                DetInd=armt.SDprs(SDind, 2);
                
                for Lind=1:lamN
                    if datStruct.useSet_C(SDind, Lind)
                        armt.SDuseSrc_C(SrcInd, Lind)=true;
                        armt.SDuseDet_C(DetInd, Lind)=true;
                    end
                    if datStruct.useSet_I(SDind, Lind)
                        armt.SDuseSrc_I(SrcInd, Lind)=true;
                        armt.SDuseDet_I(DetInd, Lind)=true;
                    end
                    if datStruct.useSet_P(SDind, Lind)
                        armt.SDuseSrc_P(SrcInd, Lind)=true;
                        armt.SDuseDet_P(DetInd, Lind)=true;
                    end
                end
            end
        case 'SS'
            armt.SSuseSrc_C=false(size(armt.rSrc, 1), lamN);
            armt.SSuseSrc_I=false(size(armt.rSrc, 1), lamN);
            armt.SSuseSrc_P=false(size(armt.rSrc, 1), lamN);
            armt.SSuseDet_C=false(size(armt.rDet, 1), lamN);
            armt.SSuseDet_I=false(size(armt.rDet, 1), lamN);
            armt.SSuseDet_P=false(size(armt.rDet, 1), lamN);
            for SSind=1:size(armt.SSprs, 3)
                SrcInds=armt.SSprs(:, 1, SSind);
                DetInds=armt.SSprs(:, 2, SSind);
                
                for Lind=1:lamN
                    if datStruct.useSet_C(SSind, Lind)
                        armt.SSuseSrc_C(SrcInds, Lind)=true;
                        armt.SSuseDet_C(DetInds, Lind)=true;
                    end
                    if datStruct.useSet_I(SSind, Lind)
                        armt.SSuseSrc_I(SrcInds, Lind)=true;
                        armt.SSuseDet_I(DetInds, Lind)=true;
                    end
                    if datStruct.useSet_P(SSind, Lind)
                        armt.SSuseSrc_P(SrcInds, Lind)=true;
                        armt.SSuseDet_P(DetInds, Lind)=true;
                    end
                end
            end
        case 'DS'
            armt.DSuseSrc_C=false(size(armt.rSrc, 1), lamN);
            armt.DSuseSrc_I=false(size(armt.rSrc, 1), lamN);
            armt.DSuseSrc_P=false(size(armt.rSrc, 1), lamN);
            armt.DSuseDet_C=false(size(armt.rDet, 1), lamN);
            armt.DSuseDet_I=false(size(armt.rDet, 1), lamN);
            armt.DSuseDet_P=false(size(armt.rDet, 1), lamN);
            for DSind=1:size(armt.DSprs, 4)
                SSprsTmp=armt.DSprs(:, :, :, DSind);
                
                for Lind=1:lamN
                    SrcIndsTmp=SSprsTmp(:, 1, :);
                    DetIndsTmp=SSprsTmp(:, 2, :);
                    if datStruct.useSet_C(DSind, Lind)
                        armt.DSuseSrc_C(SrcIndsTmp(:), Lind)=true;
                        armt.DSuseDet_C(DetIndsTmp(:), Lind)=true;
                    end
                    if datStruct.useSet_I(DSind, Lind)
                        armt.DSuseSrc_I(SrcIndsTmp(:), Lind)=true;
                        armt.DSuseDet_I(DetIndsTmp(:), Lind)=true;
                    end
                    if datStruct.useSet_P(DSind, Lind)
                        armt.DSuseSrc_P(SrcIndsTmp(:), Lind)=true;
                        armt.DSuseDet_P(DetIndsTmp(:), Lind)=true;
                    end
                end
            end
        otherwise
            error('Unknown data type');
    end
end