function [Gamma_grid, Delta_grid] = GammaDelta(SpS, sz_grid, dr)
% GammaDelta Calculates crosstalk (Gamma) and localization error (Delta) maps.
%
% [Gamma_grid, Delta_grid] = GammaDelta(SpS, sz_grid, dr)
%
% Written by Giles Blaney, Ph.D. (Summer 2020)
%
% Reference:
%   Blaney, G., Sassaroli, A., & Fantini, S. (2020). Design of a 
%   source-detector array for dual-slope diffuse optical imaging. 
%   Review of Scientific Instruments. https://doi.org/10.1063/5.0015512
%
% Inputs:
%   SpS     - (n x n) Resolution matrix [pinv(S) * S].
%   sz_grid - (1 x 3) Grid size [ny, nx, nz].
%   dr      - (1 x 3) Grid spacing [dy, dx, dz] [mm].
%
% Outputs:
%   Gamma_grid - (ny x nx x nz x 3) FWHM resolution maps.
%   Delta_grid - (ny x nx x nz x 3) Localization error maps.

    initVar=NaN([sz_grid, 3]);
    Gamma_grid=initVar;
    Delta_grid=initVar;
    clear initVar;
    for i=1:3
        yesInd=i+3;
        noInds=4:6;
        noInds(noInds==yesInd)=[];

        switch i
            case 1
                [~, R, ~]=meshgrid(1:sz_grid(1), 1:sz_grid(2), 1:sz_grid(3));
            case 2
                [R, ~, ~]=meshgrid(1:sz_grid(1), 1:sz_grid(2), 1:sz_grid(3));
            case 3
                [~, ~, R]=meshgrid(1:sz_grid(1), 1:sz_grid(2), 1:sz_grid(3));
            otherwise
                R=[];
        end

        Xm_grid=sum(sum(reshape(SpS, [sz_grid, sz_grid]),...
            noInds(1)), noInds(2));

        Gamma_grid(:, :, :, i)=sum(Xm_grid>=(max(Xm_grid, [], yesInd)/2),...
            yesInd)*dr(i);
        
        [~, maxInds]=max(Xm_grid, [], yesInd);
        Delta_grid(:, :, :, i)=(maxInds-R)*dr(i);

    end
end

