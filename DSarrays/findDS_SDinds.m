function [SDinds] = findDS_SDinds(DSprs, SDprs)
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