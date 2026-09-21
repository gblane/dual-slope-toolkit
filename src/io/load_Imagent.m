function [data,setups] = load_Imagent(filename,rowstart,swNum,doMultCorr)
% load_Imagent Read an ISS Imagent BOXY ASCII record into detector structs.
%
%   Promoted into this repository on 2026-09-21 from a loose copy. Behavior is
%   unchanged. BOXY 0.86.2 reports phase with the opposite sign convention, and
%   this reader corrects for it; removing that correction silently flips the
%   sign of every phase-derived quantity.
%
%   parseArrayData consumes the struct this returns.
%
% data: Main and auxillary data from text file
% -data.Add: including flag and marking
% -data.A: data from channel A (from 1 to 4 is 830 nm, from 5 to 8 is 690 nm)
% -data.AUX: auxillary data
% -data.timemat: time trace of the measurement

% setups: set-ups of the recorded data
% -setups.Det: number of detectors
% -setups.MUXCh: number of NIRS channels
% -setups.AuxCh: number of Auxillary channels
% -setups.Fs: sampling frequency

%% Set defaults
if nargin<=1
    rowstart=149;
end
if nargin<=2
    swNum=NaN;
    doMultCorr=false;
end

%% open file
rowstart = rowstart - 1;
idx = find(filename=='\',1,'last');
% fprintf(strcat('loading --', filename(idx+1:end), '...'));
fid=fopen(filename);
C = textscan(fid,'%s');
fclose(fid);

data = {};
data0=dlmread(filename,'\t',rowstart,1);
[r c] = size(data0);
data.Add = struct('group',data0(:,1),'step',data0(:,2),'mark',data0(:,3),'flag',data0(:,4));
data0 = data0(1:end-1,5:end);
data.data0 = data0;
% fprintf('(flag(1) = %i)\n',data.Add.flag(1));

if data.Add.flag(1)~=512 % Giles Blaney 181114
    warning('Rowstart may not be correct (flag(1)~=512)');
end

% timemat
% data.timemat = [];
% x = find(strcmp(C{:},'digaux'))+1;
% for i = 1:r-1
%     val = C{1,1}(x+(i-1)*c);
%     data.timemat(i,1) = str2num(val{1,1});
% end

%% setups
x = find(strcmp(C{:},'Version'));
setups.Ver = C{1,1}{x(1)+1,1};

x = find(strcmp(C{:},'Detector'));
setups.Det = str2double(C{1,1}{x(1)-1,1});
x = find(strcmp(C{:},'External'));
setups.MUXCh = str2double(C{1,1}{x(1)-1,1});
x = find(strcmp(C{:},'Auxiliary'));
setups.AuxCh = str2double(C{1,1}{x(1)-1,1});

x = find(strcmp(C{:},'Waveform'));
setups.CCF = str2double(C{1,1}{x(1)-1,1});
x = find(strcmp(C{:},'Waveforms'));
setups.wavSkip = str2double(C{1,1}{x(1)-1,1});
setups.wavAvg = str2double(C{1,1}{x(2)-1,1});
x = find(strcmp(C{:},'Cycles'));
setups.cycAvg = str2double(C{1,1}{x(1)-1,1});
x = find(strcmp(C{:},'Acquisitions'));
setups.aquPerWav = str2double(C{1,1}{x(1)-1,1});

x = find(strcmp(C{:},'Updata'));
if isempty(x)
    x = find(strcmp(C{:},'Update'));
    setups.Fs = str2double(C{1,1}{x(1)+2,1});
else
    setups.Fs = str2double(C{1,1}{x(1)-1,1});
end

if round(setups.Fs)==11 % if there is error in Fs reported by ISS
    setups.Fs = 9.9306;
    warning('Fs assumed to be 9.93 Hz');
end

%% Data
% x = find(data.data0 == 255);
% if isempty(x)
% x = find(data.data0 == 31);
% end 
%col = x(end)/size(data.data0,1);
col = size(data.data0,2); 

if data0(1,end) == 0
data.AUX = data0(:,col-setups.AuxCh-1:col-2);
else
data.AUX = data0(:,col-setups.AuxCh:col-1);
end 


data.fs = setups.Fs;
data.timemat = timeAxis(data.fs, size(data0, 1));

mark_start = find(strcmp(C{:},'flag'))+1;
mark_end = find(strcmp(C{:},'digaux'))-setups.AuxCh-1;

ix = 0;
store_det = C{1,1}{mark_start,1}(1:4);    
comp_det = [];
MUX0 = [];
for i = mark_start:mark_end
    ix = ix+1;
    det = C{1,1}{i,1}(1);  
    type = C{1,1}{i,1}(3:4);
    comp_det = [comp_det,strcmp(store_det,C{1,1}{i,1}(1:4))];
    
    Iraw = data0(:,ix);
    if strcmp(setups.Ver, '0.86.2') && strcmp(type, 'Ph')
        Iraw=-Iraw;
    end
    timemat = [0:length(Iraw)-1]'/setups.Fs;
    sidx = str2num(C{1,1}{i,1}(5:end)); % source index
    
    if doMultCorr
        % Correction for multiplexing
        Tcycle=1/(setups.Fs*setups.cycAvg); %Sec
        del = mod(sidx,swNum)*Tcycle/swNum;
        I = interp1(timemat+del,Iraw,timemat);
        I(1)=I(2);
    else
        I=Iraw;
    end
    
    % Store data
    if strcmp(store_det,C{1,1}{i,1}(1:4))
        MUX0 = [MUX0,I];
        data = setfield(data,det,type,MUX0);
    else
        MUX0 = [];
        MUX0 = [MUX0,I];
        store_det = C{1,1}{i,1}(1:4);
    end    
end
%% End function


