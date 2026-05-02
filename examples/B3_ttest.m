%% Setup
clear; home;

%% Find File
filesTMP=dir('*.set');
if length(filesTMP)>1
    error(['More than one .set file found, '...
        'place only one dataset in same folder']);
end

filename=filesTMP.name(1:(end-4));

load([filename '_analOutputB.mat']);

dt_buffAct=1; %sec
dt_buffRst=10; %sec
dt_ttest=10; %sec

%% Do Tests
tIndsAct=and(tFold<(dtAct-dt_buffAct), tFold>(dtAct-dt_ttest-dt_buffAct));
tIndsRst=and(tFold<(dtAct+dt_ttest+dt_buffRst), tFold>(dtAct+dt_buffRst));
nRstnAct=floor([sum(tIndsAct), sum(tIndsRst)]*(fLP/fs));

for DTind=1:length(dataTyps)
    MTnm=sprintf('%s', dataTyps{DTind}(1:2));
    DTnm=sprintf('%s', dataTyps{DTind}(3));
    
    % On fold avg channels
    eval(sprintf([...
        '[%s.hOfold_%s, %s.pOfold_%s, %s.RmA_Ofold_%s]='...
        'ttest2_nxny(%s.Ofold_%s(tIndsRst, :), %s.Ofold_%s(tIndsAct, :),'...
        '''Tail'', ''left'', ''Dim'', 1, ''Alpha'', 0.05,'...
        '''Vartype'', ''unequal'', ''n'', nRstnAct);'],...
        MTnm, DTnm, MTnm, DTnm, MTnm, DTnm, MTnm, DTnm, MTnm, DTnm));
    eval(sprintf([...
        '[%s.hDfold_%s, %s.pDfold_%s, %s.RmA_Dfold_%s]='...
        'ttest2_nxny(%s.Dfold_%s(tIndsRst, :), %s.Dfold_%s(tIndsAct, :),'...
        '''Tail'', ''right'', ''Dim'', 1, ''Alpha'', 0.05,'...
        '''Vartype'', ''unequal'', ''n'', nRstnAct);'],...
        MTnm, DTnm, MTnm, DTnm, MTnm, DTnm, MTnm, DTnm, MTnm, DTnm));
    
    % On folds maps
    tmpOrst=[];
    tmpDrst=[];
    tmpOact=[];
    tmpDact=[];
    for j=1:eval(sprintf('size(%s.OOfold_%s, 2)', MTnm, DTnm))
        if j==1
            tmpOrst=eval(sprintf(...
                '%s.OOfold_%s(tIndsRst, j, :)',...
                MTnm, DTnm));
            tmpDrst=eval(sprintf(...
                '%s.DDfold_%s(tIndsRst, j, :)',...
                MTnm, DTnm));
            
            tmpOact=eval(sprintf(...
                '%s.OOfold_%s(tIndsAct, j, :)',...
                MTnm, DTnm));
            tmpDact=eval(sprintf(...
                '%s.DDfold_%s(tIndsAct, j, :)',...
                MTnm, DTnm));
        else
            tmpOrst((end+1):(end+sum(tIndsRst)), :)=...
                eval(sprintf(...
                '%s.OOfold_%s(tIndsRst, j, :)',...
                MTnm, DTnm));
            tmpDrst((end+1):(end+sum(tIndsRst)), :)=...
                eval(sprintf(...
                '%s.DDfold_%s(tIndsRst, j, :)',...
                MTnm, DTnm));

            tmpOact((end+1):(end+sum(tIndsAct)), :)=...
                eval(sprintf(...
                '%s.OOfold_%s(tIndsAct, j, :)',...
                MTnm, DTnm));
            tmpDact((end+1):(end+sum(tIndsAct)), :)=...
                eval(sprintf(...
                '%s.DDfold_%s(tIndsAct, j, :)',...
                MTnm, DTnm));
        end
    end

    [hOOtmp, pOOtmp, RmA_OOtmp]=...
        ttest2_nxny(tmpOrst, tmpOact,...
        'Tail', 'left', 'Dim', 1, 'Alpha', 0.05,...
            'Vartype', 'unequal', 'n',...
            nRstnAct*eval(sprintf('size(%s.OOfold_%s, 2)', MTnm, DTnm)));
    [hDDtmp, pDDtmp, RmA_DDtmp]=...
        ttest2_nxny(tmpDrst, tmpDact,...
        'Tail', 'right', 'Dim', 1, 'Alpha', 0.05,...
            'Vartype', 'unequal', 'n',...
            nRstnAct*eval(sprintf('size(%s.DDfold_%s, 2)', MTnm, DTnm)));
    eval(sprintf(...
        '%s.hOOfold_%s=squeeze(hOOtmp);',...
        MTnm, DTnm));
    eval(sprintf(...
        '%s.pOOfold_%s=squeeze(pOOtmp);',...
        MTnm, DTnm));
    eval(sprintf(...
        '%s.RmA_OOfold_%s=squeeze(RmA_OOtmp);',...
        MTnm, DTnm));
    eval(sprintf(...
        '%s.hDDfold_%s=squeeze(hDDtmp);',...
        MTnm, DTnm));
    eval(sprintf(...
        '%s.pDDfold_%s=squeeze(pDDtmp);',...
        MTnm, DTnm));
    eval(sprintf(...
        '%s.RmA_DDfold_%s=squeeze(RmA_DDtmp);',...
        MTnm, DTnm));
    
    % On fold avg maps
    if ~isempty(recon.(MTnm).(DTnm))
        [recon.(MTnm).(DTnm).hO, recon.(MTnm).(DTnm).pO,...
            recon.(MTnm).(DTnm).RmA_O]=...
            ttest2_nxny(...
            recon.(MTnm).(DTnm).dO(:, :, tIndsRst),...
            recon.(MTnm).(DTnm).dO(:, :, tIndsAct),...
            'Tail', 'left', 'Dim', 3, 'Alpha', 0.05,...
            'Vartype', 'unequal', 'n', nRstnAct);
        [recon.(MTnm).(DTnm).hD, recon.(MTnm).(DTnm).pD,...
            recon.(MTnm).(DTnm).RmA_D]=...
            ttest2_nxny(...
            recon.(MTnm).(DTnm).dD(:, :, tIndsRst),...
            recon.(MTnm).(DTnm).dD(:, :, tIndsAct),...
            'Tail', 'right', 'Dim', 3, 'Alpha', 0.05,...
            'Vartype', 'unequal', 'n', nRstnAct);
    end
    
    % On folds maps
    tmpOrst=[];
    tmpDrst=[];
    tmpOact=[];
    tmpDact=[];
    if ~isempty(recon.(MTnm).(DTnm))
        for j=1:size(recon.(MTnm).(DTnm).dOdO, 4)
            if j==1
                tmpOrst=recon.(MTnm).(DTnm).dOdO(:, :, tIndsRst, j);
                tmpDrst=recon.(MTnm).(DTnm).dDdD(:, :, tIndsRst, j);
                
                tmpOact=recon.(MTnm).(DTnm).dOdO(:, :, tIndsAct, j);
                tmpDact=recon.(MTnm).(DTnm).dDdD(:, :, tIndsAct, j);
            else
                tmpOrst(:, :, (end+1):(end+sum(tIndsRst)))=...
                    recon.(MTnm).(DTnm).dOdO(:, :, tIndsRst, j);
                tmpDrst(:, :, (end+1):(end+sum(tIndsRst)))=...
                    recon.(MTnm).(DTnm).dDdD(:, :, tIndsRst, j);
                
                tmpOact(:, :, (end+1):(end+sum(tIndsAct)))=...
                    recon.(MTnm).(DTnm).dOdO(:, :, tIndsAct, j);
                tmpDact(:, :, (end+1):(end+sum(tIndsAct)))=...
                    recon.(MTnm).(DTnm).dDdD(:, :, tIndsAct, j);
            end
        end
        
        [recon.(MTnm).(DTnm).hOO, recon.(MTnm).(DTnm).pOO,...
            recon.(MTnm).(DTnm).RmA_OO]=...
            ttest2_nxny(tmpOrst, tmpOact,...
            'Tail', 'left', 'Dim', 3, 'Alpha', 0.05,...
            'Vartype', 'unequal', 'n',...
            nRstnAct*size(recon.(MTnm).(DTnm).dOdO, 4));
        [recon.(MTnm).(DTnm).hDD, recon.(MTnm).(DTnm).pDD,...
            recon.(MTnm).(DTnm).RmA_DD]=...
            ttest2_nxny(tmpDrst, tmpDact,...
            'Tail', 'right', 'Dim', 3, 'Alpha', 0.05,...
            'Vartype', 'unequal', 'n',...
            nRstnAct*size(recon.(MTnm).(DTnm).dDdD, 4));
    end
end

%% Save
save([filename '_analOutputB.mat'], '-v7.3');

%% Functions
function [h,p,difference,ci,stats] = ttest2_nxny(x,y,varargin)
% Giles Blaney Summer 2021
% Modified to accept ..., 'n', [nx, ny],...
% and output difference (the differnce of the means x-y)

%TTEST2 Two-sample t-test with pooled or unpooled variance estimate.
%   H = TTEST2(X,Y) performs a t-test of the hypothesis that two
%   independent samples, in the vectors X and Y, come from distributions
%   with equal means, and returns the result of the test in H.  H=0
%   indicates that the null hypothesis ("means are equal") cannot be
%   rejected at the 5% significance level.  H=1 indicates that the null
%   hypothesis can be rejected at the 5% level.  The data are assumed to
%   come from normal distributions with unknown, but equal, variances.  X
%   and Y can have different lengths.
%
%   This function performs an unpaired two-sample t-test. For a paired
%   test, use the TTEST function.
%
%   X and Y can also be matrices or N-D arrays.  For matrices, TTEST2
%   performs separate t-tests along each column, and returns a vector of
%   results.  X and Y must have the same number of columns.  For N-D
%   arrays, TTEST2 works along the first non-singleton dimension.  X and Y
%   must have the same size along all the remaining dimensions.
%
%   TTEST2 treats NaNs as missing values, and ignores them.
%
%   [H,P] = TTEST2(...) returns the p-value, i.e., the probability of
%   observing the given result, or one more extreme, by chance if the null
%   hypothesis is true.  Small values of P cast doubt on the validity of
%   the null hypothesis.
%
%   [H,P,CI] = TTEST2(...) returns a 100*(1-ALPHA)% confidence interval for
%   the true difference of population means.
%
%   [H,P,CI,STATS] = TTEST2(...) returns a structure with the following fields:
%      'tstat' -- the value of the test statistic
%      'df'    -- the degrees of freedom of the test
%      'sd'    -- the pooled estimate of the population standard deviation
%                 (for the equal variance case) or a vector containing the
%                 unpooled estimates of the population standard deviations
%                 (for the unequal variance case)
%
%   [...] = TTEST2(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies one or
%   more of the following name/value pairs:
%
%       Parameter       Value
%       'alpha'         A value ALPHA between 0 and 1 specifying the
%                       significance level as (100*ALPHA)%. Default is
%                       0.05 for 5% significance.
%       'dim'           Dimension DIM to work along. For example, specifying
%                       'dim' as 1 tests the column means. Default is the
%                       first non-singleton dimension.
%       'tail'          A string specifying the alternative hypothesis:
%           'both'  "means are not equal" (two-tailed test)
%           'right' "mean of X is greater than mean of Y" (right-tailed test)
%           'left'  "mean of X is less than mean of Y" (left-tailed test)
%       'vartype'       'equal' to perform the default test assuming equal
%                       variances, or 'unequal', to perform the test
%                       assuming that the two samples come from normal
%                       distributions with unknown and unequal variances.
%                       This is known as the Behrens-Fisher problem. TTEST2
%                       uses Satterthwaite's approximation for the
%                       effective degrees of freedom.
%
%   See also TTEST, RANKSUM, VARTEST2, ANSARIBRADLEY.

%   References:
%      [1] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, section 13.4. (Table 13.4.1 on page 210)

%   Copyright 1993-2017 The MathWorks, Inc.


if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 2
    error(message('stats:ttest2:TooFewInputs'));
end

% Process remaining arguments
alpha = 0.05;
tail = 0;    % code for two-sided;
vartype = '';
dim = '';

if nargin>=3
    if isnumeric(varargin{1})
        % Old syntax
        %    TTEST2(X,Y,ALPHA,TAIL,VARTYPE,DIM)
        alpha = varargin{1};
        if nargin>=4
            tail = varargin{2};
              if nargin>=5
                  vartype =  varargin{3};
                  if nargin>=6
                      dim = varargin{4};
                  end
             end
        end
        
    elseif nargin==3
            error(message('stats:ttest2:BadAlpha'));
    
    else
        % Calling sequence with named arguments
        okargs =   {'alpha' 'tail' 'vartype' 'dim' 'n'};
        defaults = {0.05    'both'    ''      ''   []};
        [alpha, tail, vartype, dim, n] = ...
                         internal.stats.parseArgs(okargs,defaults,varargin{:});
    end
end

if isempty(alpha)
    alpha = 0.05;
elseif ~isscalar(alpha) || alpha <= 0 || alpha >= 1
    error(message('stats:ttest2:BadAlpha'));
end

if isempty(tail)
    tail = 0;
elseif isnumeric(tail) && isscalar(tail) && ismember(tail,[-1 0 1])
    % OK, grandfathered
else
    [~,tail] = internal.stats.getParamVal(tail,{'left','both','right'},'''tail''');
    tail = tail - 2;
end

if isempty(vartype)
    vartype = 1;
elseif isnumeric(vartype) && isscalar(vartype) && ismember(vartype,[1 2])
    % OK, grandfathered
else
    [~,vartype] = internal.stats.getParamVal(vartype,{'equal','unequal'},'''vartype''');
end

if isempty(dim)
    % Figure out which dimension mean will work along by looking at x.  y
    % will have be compatible. If x is a scalar, look at y.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = find(size(y) ~= 1, 1); end
    if isempty(dim), dim = 1; end
    
    % If we haven't been given an explicit dimension, and we have two
    % vectors, then make y the same orientation as x.
    if isvector(x) && isvector(y)
        if dim == 2
            y = y(:)';
        else % dim == 1
            y = y(:);
        end
    end
end

% Make sure all of x's and y's non-working dimensions are identical.
sizex = size(x); sizex(dim) = 1;
sizey = size(y); sizey(dim) = 1;
if ~isequal(sizex,sizey)
    error(message('stats:ttest2:InputSizeMismatch'));
end

xnans = isnan(x);
ynans = isnan(y);
if isempty(n)
    if any(xnans(:))
        nx = sum(~xnans,dim);
    else
        nx = size(x,dim); % a scalar, => a scalar call to tinv
    end
    if any(ynans(:))
        ny = sum(~ynans,dim);
    else
        ny = size(y,dim); % a scalar, => a scalar call to tinv
    end
else
    nx=n(1);
    ny=n(2);
end


s2x = nanvar(x,[],dim);
s2y = nanvar(y,[],dim);
xmean = nanmean(x,dim);
ymean = nanmean(y,dim);
difference = xmean - ymean;

% Check for rounding issues causing spurious differences                                                                                  -
sqrtn = sqrt(nx)+sqrt(ny);
fix = (difference~=0) & ...                                     % non-zero
    (abs(difference) < sqrtn.*100.*max(eps(xmean),eps(ymean))); % but small                                                                                 -
if any(fix(:))
    % Fix any columns that are constant, even if computed difference is
    % non-zero but small
    constvalue = min(x,[],dim);
    fix = fix & all(x==constvalue | isnan(x),dim) ...
              & all(y==constvalue | isnan(y),dim);
    difference(fix) = 0;
end

if vartype == 1 % equal variances
    dfe = nx + ny - 2;
    sPooled = sqrt(((nx-1) .* s2x + (ny-1) .* s2y) ./ dfe);
    sPooled(fix) = 0;
    
    se = sPooled .* sqrt(1./nx + 1./ny);
    ratio = difference ./ se;

    if (nargout>3)
        stats = struct('tstat', ratio, 'df', cast(dfe,'like',ratio), ...
                       'sd', sPooled);
        if isscalar(dfe) && ~isscalar(ratio)
            stats.df = repmat(stats.df,size(ratio));
        end
    end
elseif vartype == 2 % unequal variances
    s2xbar = s2x ./ nx;
    s2ybar = s2y ./ ny;
    dfe = (s2xbar + s2ybar) .^2 ./ (s2xbar.^2 ./ (nx-1) + s2ybar.^2 ./ (ny-1));
    se = sqrt(s2xbar + s2ybar);
    se(fix) = 0;
    ratio = difference ./ se;

    if (nargout>3)
        stats = struct('tstat', ratio, 'df', cast(dfe,'like',ratio), ...
                       'sd', sqrt(cat(dim, s2x, s2y)));
        if isscalar(dfe) && ~isscalar(ratio)
            stats.df = repmat(stats.df,size(ratio));
        end
    end
    
    % Satterthwaite's approximation breaks down when both samples have zero
    % variance, so we may have gotten a NaN dfe.  But if the difference in
    % means is non-zero, the hypothesis test can still reasonable results,
    % that don't depend on the dfe, so give dfe a dummy value.  If difference
    % in means is zero, the hypothesis test returns NaN.  The CI can be
    % computed ok in either case.
    if all(se(:) == 0), dfe = 1; end
end

% Compute the correct p-value for the test, and confidence intervals
% if requested.
if tail == 0 % two-tailed test
    p = 2 * tcdf(-abs(ratio),dfe);
    if nargout > 2
        spread = tinv(1 - alpha ./ 2, dfe) .* se;
        ci = cat(dim, difference-spread, difference+spread);
    end
elseif tail == 1 % right one-tailed test
    p = tcdf(-ratio,dfe);
    if nargout > 2
        spread = tinv(1 - alpha, dfe) .* se;
        ci = cat(dim, difference-spread, Inf(size(p)));
    end
elseif tail == -1 % left one-tailed test
    p = tcdf(ratio,dfe);
    if nargout > 2
        spread = tinv(1 - alpha, dfe) .* se;
        ci = cat(dim, -Inf(size(p)), difference+spread);
    end
end

% Determine if the actual significance exceeds the desired significance
h = cast(p <= alpha, 'like', p);
h(isnan(p)) = NaN; % p==NaN => neither <= alpha nor > alpha
end