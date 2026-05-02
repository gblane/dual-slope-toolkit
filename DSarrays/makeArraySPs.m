function [SP, SPinds] = makeArraySPs(locs)
    arguments
        locs (3,:) double;
    end
    
    N=size(locs, 2);

    SP=numSubplots(N);
    
    %% Make Coor Grid
    xGrd=linspace(min(locs(1, :)), max(locs(1, :)), SP(2));
    yGrd=fliplr(linspace(min(locs(2, :)), max(locs(2, :)), SP(1)));
    
    [XXgrd, YYgrd]=meshgrid(xGrd, yGrd);
    
    tmp=XXgrd'; 
    xVect=tmp(:);
    
    tmp=YYgrd'; 
    yVect=tmp(:);
    
    %% Make Gird Order
    ix=1:ceil(SP(2)/2);
    jx=SP(2)-(0:1:floor(SP(2)/2-1));
    kx=NaN(SP(2), 1);
    kx(1:2:end)=ix;
    kx(2:2:end)=jx;
    
    iy=(1:ceil(SP(1)/2))-1;
    jy=SP(1)-(0:1:floor(SP(1)/2-1))-1;
    ky=NaN(SP(1), 1);
    ky(1:2:end)=iy;
    ky(2:2:end)=jy;
    
    kVect=repmat(kx, [SP(1), 1])+repelem(ky, SP(2))*SP(2);
    
    %% Loop
    SPinds=NaN(1, N);
    for n=1:N
        k=kVect(n);
        x=xVect(k);
        y=yVect(k);
    
        [~, ind]=min(abs(sqrt((locs(1, :)-x).^2+(locs(2, :)-y).^2)));
        
        SPinds(ind)=k;
        locs(:, ind)=Inf;
    end
end