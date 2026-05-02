function [outp] = DPF_DSF_calc(rs, rd, muaLamb, muspLamb, calcmode,...
    fmod, nin, opt2L)
% [outp] = DPF_DSF_calc(rs, rd, muaLamb, muspLamb, calcmode, fmod, nin)
% Updated by Giles Blaney Spring 2022
% Inputs:
%   rs          - Pencil beam source positions in [x, y, z] with direction 
%                   vector of [0, 0, 1].
%                   Size: [numSD, 3]; Units: mm
%   rd          - Detector positions in [x, y, z].
%                   Size: [numSD, 3]; Units: mm
%   muaLamb     - Vector of absorption coefficients.
%                   Size: [numOptProp, 1] OR [1, numOptProp]; Units: 1/mm
%   muspLamb    - Vector of reduced scattering coefficients.
%                   Size: [numOptProp, 1] OR [1, numOptProp]; Units: 1/mm
%   calcmode    - String with four possible values:
%                   - 'DPF_I': To calculate the intensity DPF
%                   - 'DPF_Ph' or 'DPF_P': To calculate the phase DPF
%                   - 'DSF_I': To calculate the intensity DSF
%                   - 'DSF_Ph' or 'DSF_P': To calculate the phase DSF
%                   - ...
%   fmod        - OPTIONAL, default=140.625e6 Hz
%                   Modulation frequency
%                   Size: scalar; Units: Hz
%   nin         - OPTIONAL, default=1.4
%                   Index of refraction inside the medium (outside assumed
%                   to be 1).
%                   Size: scalar; Units: unitless
%               %%% Added by Cristy on 03/08/2023 %%%%
%  opts2L         - OPTIONAL; struct - when you want to do the 2-Layer you
%                   need to have these inputs
%                   totL   - String controlling the total pathlength to be used
%                        with values:
%                        '2L'      - Use true two-layer total pathlengths
%                        'homoEff' - [DO NOT USE] Use homogenous total pathlengths from
%                        effective homogenous optical properties, as-if
%                        effective homogenous optical properties were used
%                        to find DPF and DSF in a real-world measurment
%                   thk     - Layer thickness. (mm)
%                   nin  - (default=[1.4, 1.4]) Index of refraction inside. (-)
%                   nout - (default=1) Index of refraction outside. (-)
%                   musp - Size lambda x layer (default=[1.20, 0.25] 1/mm) Reduced scattering.
%                           (1/mm)
%                   mua  - Size lambda x layer (default=[0.008, 0.020] 1/mm) Absorption. (1/mm) 
%                   fmod   - (default=1.40625 Hz) Modulation frequency {Hz}
%                   h_end  - (default=2000) Number of Bessel function zeros
%                   B      - (default=150 mm) Radius of cylindrical boundary {mm}
% Outputs:
%   outd        - Value requested by calcmode input.
%                   Units: unitless
%                   - For DSF or DPF if numSD=1
%                       Size: same as muaLamb input
%                   - For DPF if numSD>1
%                       Size: [numOptProp, numSD]

    %% Setup
    if exist('opt2L', 'var')
        if strcmp(opt2L.totL, 'homoEff')
            error('Do not use opts2L.totL=homoEff')
        end
    end
    
    if nargin<6 % only works if you dont have any of these
        fmod=1.40625e8; %Hz
        nin=1.4;
        opt2L.totL='homoEff';
    end
    nout=1;
    % Added for scenario where nargin<6 isnt true bc you input fmod and nin
    if ~exist('opt2L', 'var')
       opt2L.totL='homoEff';
    end 

    
        
    switch opt2L.totL
        case '2L'
            % deal with mua0Lam and muspLam so they dont actually need to be
            % input 
            muaLamb = ones(size(opt2L.mua,1), 1);
            muspLamb = ones(size(opt2L.mua,1), 1);
    end
    

    rho=(sqrt(sum((rs-rd).^2, 2))).';
    rhoAvg=mean(rho);

    tmp=split(calcmode, '_');
    coefNm=tmp{1};
    mTyp=tmp{2};
    if strcmp(coefNm, 'DPF') && length(rho)~=1
        outp=NaN(length(muaLamb), length(rho));
    else
        outp=NaN(size(muaLamb));
    end

    if strcmp(mTyp, 'Ph')
        mTyp="P";
    end  

    load('zeroOrdBesselRoots.mat', 'en');
    %% Calc DPF or DSF
    for lamInd=1:length(muaLamb)     
        
        L=NaN(size(rho));
        Lcw=NaN(size(rho));
        R=NaN(size(rho));
        Rcw=NaN(size(rho));
        
        % Change to a switch 
switch opt2L.totL
    case 'homoEff'
        optProp.nin=nin;
        optProp.nout=nout;
        optProp.mua=muaLamb(lamInd);
        optProp.musp=muspLamb(lamInd);
        rs_iso=rs+[0, 0, 1./optProp.musp];
        for i=1:length(rho)
            [L(i), R(i)]=complexTotPathLen(rs_iso(i, :), rd(i, :),...
                2*pi*fmod, optProp);
            [Lcw(i), Rcw(i)]=complexTotPathLen(rs_iso(i, :), rd(i, :),...
               0, optProp);
        end
    case '2L'
        thk = opt2L.thk; 
        [optProp, opts] = addNeededInput2L(opt2L, lamInd); 
        rs_iso=rs+[0, 0, 1./optProp.musp(1)];
        for i=1:length(rho)
        [L(i), ~, R(i)] =...
            complexTotPathLen2L(rs_iso(i, :), rd(i, :),...
            thk, en, optProp, opts);
        optsCW = opts; 
        optsCW.fmod = 0; 
        [Lcw(i), ~, Rcw(i)] =...
            complexTotPathLen2L(rs_iso(i, :), rd(i, :),...
            thk, en, optProp, optsCW);
        end
end 
        LY=calcPathLen_datTyp(R, L, Rcw, Lcw, mTyp);
        
        switch coefNm
            case 'DPF'
                if length(rho)~=1
                    outp(lamInd, :)=LY./rho;
                else
                    outp(lamInd)=LY./rho;
                end

            case 'DSF'
                outp(lamInd)=sum((rho-rhoAvg).*LY, 2)/...
                    (length(rho)*var(rho, 1));

            otherwise
                error('Unknown calcmode');
        end
    end
end


function [optProp, opts] = addNeededInput2L(opt2L, lamInd)
% Add inputs needed to calculate the complex total path length for a 2
% -layer medium. This allows for some variables to be assumed

% mua
if isfield(opt2L, 'mua')
    optProp.mua = opt2L.mua(lamInd, :);
else
    optProp.mua = [0.008, 0.020];
end
% musp
if isfield(opt2L, 'musp')
    optProp.musp = opt2L.musp(lamInd, :);
else
    optProp.musp = [1.2, 0.25];
end
% n in 
if isfield(opt2L, 'nin')
    optProp.nin = opt2L.nin;
else
    optProp.nin = [1.4, 1.4];
end
% n out
if isfield(opt2L, 'nout')
    optProp.nout = opt2L.nout;
else
    optProp.nout = 1;
end
% Modulation frequency
if isfield(opt2L, 'fmod')
    opts.fmod = opt2L.fmod;
else
    opts.fmod = 1.40625e8;
end
% Number of Bessel function zeros
if isfield(opt2L, 'h_end')
    opts.h_end = opt2L.h_end;
else
    opts.h_end = 2000;
end
% Radius of cylindrical boundary {mm}
if isfield(opt2L, 'B')
    opts.B = opt2L.B;
else
    opts.B = 150;
end

end



