function h = plotVectorizedMap(r_xy, map, NVA)
    arguments
        r_xy (:,2) double;
        map (:,1) double;

        NVA.doSmooth logical = false;
        NVA.smSz (1,1) = 25; %mm

        NVA.doComplex string ...
            {mustBeMember(NVA.doComplex,{'no', 'amp', 'ang'})} =...
            'no';
        NVA.complexAmpRng (1,2) double = [0, inf]; 

        NVA.angWrap (1,1) double = pi;
    end
    
    if ~isempty(r_xy)
        dx=min(diff(unique(r_xy(:, 1))));
        dy=min(diff(unique(r_xy(:, 2))));
        [XX, YY]=meshgrid(min(r_xy(:, 1)):dx:max(r_xy(:, 1)),...
            min(r_xy(:, 2)):dy:max(r_xy(:, 2)));
    
    %     ZZ=griddata(r_xy(:, 1), r_xy(:, 2), map, XX, YY);
        ZZ=NaN(size(XX));
        for i=1:length(map)
            [~, i1]=min(abs(XX(1, :)-r_xy(i, 1)));
            [~, i2]=min(abs(YY(:, 1)-r_xy(i, 2)));
            ZZ(i2, i1)=mean([ZZ(i2, i1), map(i)], 'omitnan');
        end
        
        if ~strcmp(NVA.doComplex, 'no')
            minAmp=NVA.complexAmpRng(1);
            maxAmp=NVA.complexAmpRng(2);
    
            ZZ(abs(ZZ)>maxAmp)=maxAmp.*exp(1i*angle(ZZ(abs(ZZ)>maxAmp)));
            ZZ(abs(ZZ)<minAmp)=minAmp.*exp(1i*angle(ZZ(abs(ZZ)<minAmp)));
        end
        
        if NVA.doSmooth
            HH=zeros(2*NVA.smSz/dx, 2*NVA.smSz/dy);
            [II, JJ]=meshgrid(-size(HH, 1)/2+0.5:size(HH, 1)/2-0.5,...
                -size(HH, 2)/2+0.5:size(HH, 2)/2-0.5);
            HH=exp(-(II.^2+JJ.^2)/(2*((NVA.smSz/2).^2)));
            HH=HH-min(HH(round(size(HH, 1)/2), :));
            HH(HH<0)=0;
            HH=HH/sum(HH(:));
    
            ZZ=nanconv(ZZ, HH, 'nanout');
        end
    
        switch NVA.doComplex
            case 'amp'
                ZZ=abs(ZZ);
            case 'ang'
                ZZ=wrapTo(angle(ZZ), NVA.angWrap);
            otherwise
        end
    
        h=pcolor(XX, YY, ZZ);
        shading flat;
        axis equal;
    else
        h=pcolor(NaN, NaN, NaN);
        shaing flat;
        axis equal;
    end
end