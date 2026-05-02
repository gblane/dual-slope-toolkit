function [SDinds] = findDS_SDinds(DSprs, SDprs)
% findDS_SDinds Maps Dual-Slope (DS) pairs to their constituent Single-Distance (SD) indices.
%
% [SDinds] = findDS_SDinds(DSprs, SDprs)
%
% Written by Giles Blaney, Ph.D.
%
% Inputs:
%   DSprs - 4D matrix of Dual-Slope pairs.
%   SDprs - 2D matrix of all Single-Distance pairs.
%
% Outputs:
%   SDinds - Matrix mapping each DS set to its SD constituent indices.
    SDinds=[];
    for i=1:size(DSprs, 4)
        SDindsTmp=[];
        for j=1:size(DSprs, 1)
            for k=1:size(DSprs, 3)
                SDpr=DSprs(j, :, k, i);
                SDindsTmp=[SDindsTmp, find(all(SDpr==SDprs, 2))];
            end
        end
        SDinds=[SDinds; SDindsTmp];
    end
end