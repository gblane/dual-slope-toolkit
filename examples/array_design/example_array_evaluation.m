%% Setup
clear; home;
% Giles Blaney Summer 2020

% G. Blaney, A. Sassaroli, and S. Fantini, “Design of a source-detector 
% array for dual-slope diffuse optical imaging,” Review of Scientific 
% Instruments, https://doi.org/10.1063/5.0015512.

load('S_HEXarray.mat');
m=size(S, 1);
n=size(S, 2);
sz_grid=size(X);
clear X Y Z
dr=[median(diff(y)), median(diff(x)), median(diff(z))];

doPlot=true;

%% Calc SpS
SpS=pinv(S)*S;
clear S;

%% Calc R and C Maps
% Reconstruction
R=diag(SpS)/m;

% Crosstalk
C=(sum(SpS, 2)-diag(SpS))/(n-1);

% Vec to Grid
R_grid=reshape(R, sz_grid);
C_grid=reshape(C, sz_grid);
clear R C;

%% Gamma and Delta Maps
% Gamma (Resolution) and Delta (Localization)
[Gamma_grid, Delta_grid]=GammaDelta(SpS, sz_grid, dr);

%% Plot
if doPlot
    col=jet(100);
    col(1, :)=[0, 0, 0];
    col(end, :)=[1, 1, 1];
    
    Delta_xy=sqrt(Delta_grid(:, :, 2, 1).^2+Delta_grid(:, :, 2, 2).^2);
    Gamma_xy=sqrt(Gamma_grid(:, :, 2, 1).^2+Gamma_grid(:, :, 2, 2).^2);
    
    temp=R_grid(:, :, 2);
    clR=[quantile(temp(:), 0.15), quantile(temp(:), 0.99)];
    temp=C_grid(:, :, 2);
    clC=[quantile(temp(:), 0.15), quantile(temp(:), 0.99)];
    clG=[quantile(Gamma_xy(:), 0.01), quantile(Gamma_xy(:), 0.99)];
    clD=[quantile(Delta_xy(:), 0.001), quantile(Delta_xy(:), 0.95)];

    figure(3); clf; colormap(col);
    subaxis(2, 2, 1);
    imagesc(x, y, R_grid(:, :, 2));
    caxis(clR);
    cb=colorbar;
    ylabel(cb, 'diag(S^+S)/m');
    set(gca, 'YDir', 'Normal');
    set(gca, 'XTickLabels', {});
    axis equal tight;
    ylabel('y (mm)');
    title('Reconstruction (R)');
    
    subaxis(2, 2, 2);
    imagesc(x, y, C_grid(:, :, 2));
    caxis(clC);
    cb=colorbar;
    ylabel(cb, '\Sigma_{offdiag}S^+S/(n-1)');
    set(gca, 'YDir', 'Normal');
    set(gca, 'XTickLabels', {});
    set(gca, 'YTickLabels', {});
    axis equal tight;
    title('Crosstalk (C)');
    
    subaxis(2, 2, 3);
    imagesc(x, y, Gamma_xy);
    caxis(clG);
    cb=colorbar;
    ylabel(cb, '(\Gamma_x^2+\Gamma_y^2)^{1/2} (mm)');
    set(gca, 'YDir', 'Normal');
    xlabel('x (mm)');
    ylabel('y (mm)');
    axis equal tight;
    title('Resolution (\Gamma)');
    
    subaxis(2, 2, 4);
    imagesc(x, y, Delta_xy);
    caxis(clD);
    cb=colorbar;
    ylabel(cb, '(\Delta_x^2+\Delta_y^2)^{1/2} (mm)');
    set(gca, 'YDir', 'Normal');
    xlabel('x (mm)');
    set(gca, 'YTickLabels', {});
    axis equal tight;
    title('Localization (\Delta)');
end