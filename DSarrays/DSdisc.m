function [armt, fh] = DSdisc(coors, nameValArgs)
    %% Parse Input
    arguments
        coors struct;
        
        nameValArgs.scl double = 1;
        nameValArgs.rRng (1,2) double = [0, 45]; %mm
        nameValArgs.lRng (1,2) double =[5-1e-3, 25]; %mm
        nameValArgs.lTol double = 1; %mm
        nameValArgs.plot ...
            string {mustBeMember(nameValArgs.plot,{'none', 'summary', 'full'})}...
            = 'none';
    end
    armt=[];
    fh=[];

    AllSrcs=coors.AllSrcs;
    AllDets=coors.AllDets;
    
    armt.rSrc=AllSrcs*nameValArgs.scl;
    armt.rDet=AllDets*nameValArgs.scl;
    
    if isfield(coors, 'name')
        armt.name=coors.name;
    else
        armt.name='Unnamed_Array';
    end

    switch nameValArgs.plot
        case 'none'
            doPlot=false;
            doFullPlot=false;
        case 'summary'
            doPlot=true;
            doFullPlot=false;
        case 'full'
            doPlot=true;
            doFullPlot=true;
        otherwise
            error('Unknown plot option');
    end

    
    %% Find All SD Pairs
    armt.SDprs=[];
    r_all=NaN(size(armt.rSrc, 1), size(armt.rDet, 1));
    for sInd=1:size(armt.rSrc, 1)
        for dInd=1:size(armt.rDet, 1)
            % Check distance range requirement (nameValArgs.rRng)
            r=norm(armt.rSrc(sInd, :)-armt.rDet(dInd, :));
            r_all(sInd, dInd)=r;
            if r<nameValArgs.rRng(1) || r>nameValArgs.rRng(2)
                continue;
            end
            
            armt.SDprs(end+1, :)=[sInd, dInd];
        end
    end
    
    %% Find All SS Pairs
    armt.SSprs=[];
    n=0;
    for sd1Ind=1:size(armt.SDprs, 1)
        for sd2Ind=1:size(armt.SDprs, 1)
            if sd1Ind>=sd2Ind
                continue;
            end
            
            s1Ind=armt.SDprs(sd1Ind, 1);
            d1Ind=armt.SDprs(sd1Ind, 2);
            s2Ind=armt.SDprs(sd2Ind, 1);
            d2Ind=armt.SDprs(sd2Ind, 2);
            
            if xor(s1Ind==s2Ind, d1Ind==d2Ind)
                r1=norm(armt.rSrc(s1Ind, :)-armt.rDet(d1Ind, :));
                r2=norm(armt.rSrc(s2Ind, :)-armt.rDet(d2Ind, :));
                
                % Check lever arm requirement (nameValArgs.lRng)
                if r1~=r2 && and(abs(r1-r2)>=nameValArgs.lRng(1),...
                        abs(r1-r2)<=nameValArgs.lRng(2))
                    n=n+1;
                    
                    if r1<r2
                        armt.SSprs(:, :, n)=[...
                            s1Ind, d1Ind;...
                            s2Ind, d2Ind];
                    elseif r2<r1
                        armt.SSprs(:, :, n)=[...
                            s2Ind, d2Ind;...
                            s1Ind, d1Ind];
                    end
                end
            end
        end
    end
    
    %% Find All DS Pairs
    armt.DSprs=[];
    sdAllHist={};
    n=0;
    for ss1Ind=1:size(armt.SSprs, 3)
        for ss2Ind=1:size(armt.SSprs, 3)
            if ss1Ind>=ss2Ind
                continue;
            end
            
            s11Ind=armt.SSprs(1, 1, ss1Ind);
            d11Ind=armt.SSprs(1, 2, ss1Ind);
            s21Ind=armt.SSprs(1, 1, ss2Ind);
            d21Ind=armt.SSprs(1, 2, ss2Ind);
            s12Ind=armt.SSprs(2, 1, ss1Ind);
            d12Ind=armt.SSprs(2, 2, ss1Ind);
            s22Ind=armt.SSprs(2, 1, ss2Ind);
            d22Ind=armt.SSprs(2, 2, ss2Ind);
            
            if length(unique([s11Ind, s21Ind, s12Ind, s22Ind]))==2 && ...
                length(unique([d11Ind, d21Ind, d12Ind, d22Ind]))==2
                
                sd11=[s11Ind, d11Ind];
                sd21=[s21Ind, d21Ind];
                sd12=[s12Ind, d12Ind];
                sd22=[s22Ind, d22Ind];
                
                sdAll=[sd11; sd21; sd12; sd22];
                
                if sum(size(sdAll)==size(unique(sdAll, 'rows')))==2
                    r11=norm(armt.rSrc(s11Ind, :)-armt.rDet(d11Ind, :));
                    r21=norm(armt.rSrc(s21Ind, :)-armt.rDet(d21Ind, :));
                    r12=norm(armt.rSrc(s12Ind, :)-armt.rDet(d12Ind, :));
                    r22=norm(armt.rSrc(s22Ind, :)-armt.rDet(d22Ind, :));
                    
                    [~, ord1]=sort([r11, r12]); %Order of slope 1
                    [~, ord2]=sort([r21, r22]); %Order of slope 2
    
                    s1Inds=[s11Ind, s12Ind];
                    s1Inds=s1Inds(ord1); %Ordered slope 1 sources
                    s2Inds=[s21Ind, s22Ind];
                    s2Inds=s2Inds(ord2); %Ordered slope 2 sources
    
                    d1Inds=[d11Ind, d12Ind];
                    d1Inds=d1Inds(ord1); %Ordered slope 1 detectors
                    d2Inds=[d21Ind, d22Ind];
                    d2Inds=d2Inds(ord2); %Ordered slope 2 detectors
    
                    if sum([s1Inds==s2Inds, d1Inds==d2Inds])>0
                        continue;
                    end
                    
                    if abs(abs(r12-r11)-abs(r22-r21))<=nameValArgs.lTol
                        newSet=true;
                        for i=1:length(sdAllHist)
                            if size(unique([sdAllHist{i}; sdAll], 'rows'))<=4
                                newSet=false;
                                break;
                            end
                        end
    
                        if newSet
                            n=n+1;
                            sdAllHist{n}=sdAll;
    
                            armt.DSprs(:, :, 1, n)=armt.SSprs(:, :, ss1Ind);
                            armt.DSprs(:, :, 2, n)=armt.SSprs(:, :, ss2Ind);
                        end
                    end
                end
            end
        end
    end
    
    %% Plotting
    % Plot Setup
    if doPlot || doFullPlot
        col=hsv(size(armt.DSprs, 4));
        yl=[min([armt.rSrc(:, 2); armt.rDet(:, 2)])-10,...
            max([armt.rSrc(:, 2); armt.rDet(:, 2)])+10];
        xl=[min([armt.rSrc(:, 1); armt.rDet(:, 1)])-10,...
            max([armt.rSrc(:, 1); armt.rDet(:, 1)])+10];
    
        rhos=NaN(size(armt.DSprs, 1)*2, size(armt.DSprs, 4));
        drhos=NaN(size(armt.DSprs, 1), size(armt.DSprs, 4));
        for i=1:size(armt.DSprs, 4)
            SSpr1=armt.DSprs(:, :, 1, i);
            SSpr2=armt.DSprs(:, :, 2, i);
            SDprs=[SSpr1; SSpr2];
            
            rhos(:, i)=...
                vecnorm(armt.rSrc(SDprs(:, 1), :)-...
                armt.rDet(SDprs(:, 2), :), 2, 2);
            
            SS1_rhos=...
                vecnorm(armt.rSrc(SSpr1(:, 1), :)-...
                armt.rDet(SSpr1(:, 2), :), 2, 2);
            SS2_rhos=...
                vecnorm(armt.rSrc(SSpr2(:, 1), :)-...
                armt.rDet(SSpr2(:, 2), :), 2, 2);
            drhos(:, i)=abs([diff(SS1_rhos); diff(SS2_rhos)]);
        end
    end
    
    if doPlot
        % Plot Summary
        fh(1)=figure('Name', [armt.name '_summary']); clf;
        subplot(2, 2, [1, 3]);
        for i=1:size(armt.DSprs, 4)
            SSpr1=armt.DSprs(:, :, 1, i);
            SSpr2=armt.DSprs(:, :, 2, i);
            SDprs=[SSpr1; SSpr2];
        
            for k=1:size(SDprs, 1)
                SDprs=[SSpr1; SSpr2];
        
                sInd=SDprs(k, 1);
                dInd=SDprs(k, 2);
        
                plot([armt.rSrc(sInd, 1), armt.rDet(dInd, 1)],...
                    [armt.rSrc(sInd, 2), armt.rDet(dInd, 2)],...
                    '-', 'color', col(i, :));
                hold on;
            end
        end
        
        for sInd=1:size(armt.rSrc, 1)
            plot(armt.rSrc(sInd, 1), armt.rSrc(sInd, 2), 'sr');
            text(armt.rSrc(sInd, 1)+3, armt.rSrc(sInd, 2)+3,...
                sprintf('%d', sInd),...
                'Color', 'r');
        end
        for dInd=1:size(armt.rDet, 1)
            plot(armt.rDet(dInd, 1), armt.rDet(dInd, 2), 'ob');
            text(armt.rDet(dInd, 1)+3, armt.rDet(dInd, 2)+3,...
                sprintf('%s', char(64+dInd)),...
                'Color', 'b');
        end
        
        hold off;
        xlabel('x (mm)', 'interpreter', 'latex');
        ylabel('y (mm)', 'interpreter', 'latex');
        axis equal tight;
        xlim(xl);
        ylim(yl);
        title(sprintf('%d Sources, %d Detectors\n%d Dual-Slope sets', ...
            size(armt.rSrc, 1), size(armt.rDet, 1), size(armt.DSprs, 4)),...
            'interpreter', 'latex');
        
        subplot(2, 2, 2);
        histogram(rhos(:));
        ylabel('counts', 'interpreter', 'latex');
        xlabel('$\rho$ (mm)', 'interpreter', 'latex');
        title(sprintf(['$\\rho_{SD,min}=%.1f$ mm, $\\rho_{SD,max}=%.1f$ mm\n'...
            '$\\rho_{DS,min}=%.1f$ mm, $\\rho_{DS,max}=%.1f$ mm'],...
            min(r_all(:)), max(r_all(:)), min(rhos(:)), max(rhos(:))),...
            'interpreter', 'latex');
        
        subplot(2, 2, 4);
        histogram(drhos(:));
        ylabel('counts', 'interpreter', 'latex');
        xlabel('$\Delta\rho$ (mm)', 'interpreter', 'latex');
        title(sprintf(...
            '$\\Delta\\rho_{DS,min}=%.1f$ mm, $\\Delta\\rho_{DS,max}=%.1f$ mm',...
            min(drhos(:)), max(drhos(:))),...
            'interpreter', 'latex');
        
        sgtitle(armt.name, 'interpreter', 'none');
    end
    
    % Plot Full
    if doFullPlot
        sp=numSubplots(size(armt.DSprs, 4));
        fh(2)=figure('Name', [armt.name '_full']); clf;
        
        for i=1:size(armt.DSprs, 4)
            SSpr1=armt.DSprs(:, :, 1, i);
            SSpr2=armt.DSprs(:, :, 2, i);
            SDprs=[SSpr1; SSpr2];
            
            subplot(sp(1), sp(2), i);
            for sInd=1:size(armt.rSrc, 1)
                plot(armt.rSrc(sInd, 1), armt.rSrc(sInd, 2), 'sr'); hold on;
                text(armt.rSrc(sInd, 1)+3, armt.rSrc(sInd, 2)+3,...
                    sprintf('%d', sInd),...
                    'Color', 'r');
            end
            for dInd=1:size(armt.rDet, 1)
                plot(armt.rDet(dInd, 1), armt.rDet(dInd, 2), 'ob');
                text(armt.rDet(dInd, 1)+3, armt.rDet(dInd, 2)+3,...
                    sprintf('%s', char(64+dInd)),...
                    'Color', 'b');
            end
            for k=1:size(SDprs, 1)
                SDprs=[SSpr1; SSpr2];
        
                sInd=SDprs(k, 1);
                dInd=SDprs(k, 2);
        
                plot([armt.rSrc(sInd, 1), armt.rDet(dInd, 1)],...
                    [armt.rSrc(sInd, 2), armt.rDet(dInd, 2)],...
                    '-', 'color', col(i, :));
                
                plot([armt.rSrc(sInd, 1), armt.rDet(dInd, 1)],...
                    [armt.rSrc(sInd, 2), armt.rDet(dInd, 2)],...
                    'xk');
            end; hold off;
            if i>(prod(sp)-sp(2))
                xlabel('x (mm)', 'interpreter', 'latex');
            else
                set(gca, 'XTickLabel', {});
            end
            if mod(i, sp(2))==1
                ylabel('y (mm)', 'interpreter', 'latex');
            else
                set(gca, 'YTickLabel', {});
            end
            axis equal tight;
            xlim(xl);
            ylim(yl);

            title(sprintf(['%d/%d\n'...
                '$\\rho=[%.0f, %.0f, %.0f, %.0f]$ mm\n'...
                '$\\Delta\\rho=[%.0f, %.0f]$ mm'],...
                i, size(armt.DSprs, 4),...
                rhos(:, i), drhos(:, i)), 'interpreter', 'latex');
        end
        sgtitle(armt.name, 'interpreter', 'none');
    end
end