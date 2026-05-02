function [SDinds] = findSS_SDinds(SSprs, SDprs)
% findSS_SDinds Maps Single-Slope (SS) pairs to their constituent Single-Distance (SD) indices.
%
% [SDinds] = findSS_SDinds(SSprs, SDprs)
%
% Written by Giles Blaney, Ph.D.
%
% Inputs:
%   SSprs - 3D matrix of Single-Slope pairs.
%   SDprs - 2D matrix of all Single-Distance pairs.
%
% Outputs:
%   SDinds - Matrix mapping each SS set to its SD constituent indices.
    SDinds=[];
    for i=1:size(SSprs, 3)
        SDindsTmp=[];
        for j=1:size(SSprs, 1)
            SDpr=SSprs(j, :, i);
            SDindsTmp=[SDindsTmp, find(all(SDpr==SDprs, 2))];
        end
        SDinds=[SDinds; SDindsTmp];
    end
end