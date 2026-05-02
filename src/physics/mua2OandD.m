function [O, D] = mua2OandD(mua, lambda)
% mua2OandD Calculates hemodynamics (HbO, HbR) from absorption spectra.
%
% [O, D] = mua2OandD(mua, lambda)
%
% Written by Giles Blaney, Ph.D. (Fall 2019)
%
% Reference:
%   Bigio, I., & Fantini, S. (2016). Quantitative Biomedical Optics:
%   Theory, Methods, and Applications (Cambridge Texts in Biomedical
%   Engineering). Cambridge: Cambridge University Press.
%   doi:10.1017/CBO9781139029797
%
% Inputs:
%   mua    - (nTime x nLambda) Matrix of absorption coefficients [1/cm].
%   lambda - (1 x nLambda) Vector of wavelengths [nm].
%
% Outputs:
%   O      - Oxyhemoglobin concentration [uM].
%   D      - Deoxyhemoglobin concentration [uM].

    %% Parse Input
    if size(mua, 1)~=length(lambda)
        mua=mua.';
    end
    if size(lambda, 1)~=1
        lambda=lambda.';
    end
    
    %% Interpolate Extinction Spectra
    spectra=load('OandDspect.mat');
    Oext=interp1(spectra.lambda, spectra.Oext, lambda); %1/(mM cm)
    Dext=interp1(spectra.lambda, spectra.Dext, lambda); %1/(mM cm)
    
    %% Slove for O and D
    X=linsolve([Oext.', Dext.'], mua);
    O=X(1, :).'*1000; %uM
    D=X(2, :).'*1000; %uM
    
end

