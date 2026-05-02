function [SDinds] = findSS_SDinds(SSprs, SDprs)
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