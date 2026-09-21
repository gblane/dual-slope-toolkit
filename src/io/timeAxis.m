function [ t ] = timeAxis( fs, l )
% timeAxis Column vector of sample times for a run of length l at rate fs.
%
%   t = timeAxis(fs, l) returns an l-by-1 column of times in seconds, starting
%   at 0. Used by load_Imagent to build data.timemat.

    t=(0:(1/fs):((l-1)/fs))';
end

