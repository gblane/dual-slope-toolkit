%% Setup
% Giles Blaney Summer 2021
clear; home;

% Name of experiment
dataPrefix='testLJ'; % <---------------------------------------------------
% Your name
operatorName='TestPerson';
% Data will be saved in:
% ~/Desktop/DRSdata/*operatorName*/*dataPrefix*_YYMMDD_hhmmss.mat

% The faster (less time) of the two determines run time
% Number of DS samples
numSamp=1e6; % <----------------------------------------------------------
% Total run time                           OR
runTime=60; %sec <---------------------------------------------------------

% Number of repeated SD integrations for averaging
meas_num=1;
% Length of SD integration
int_time=10; %ms <--------------------------------------------------------

% Assumed scattering values for real-time plot of absorption
musp0=0.7; %1/mm
lam0=830; %nm
b=1;
% Source-detector distnaces in [S1D1, S1D2, S2D2, S2D1]
rhos=[25, 35, 25, 35]; %mm
% Wavelength range to fit for chromophores in real-time plot
lamRng=[650, 1064]; %nm

% Plot data in real time
doPlot=true;
% Sollect dark spectrum
doDark=true; 

%% Init
fprintf('Setting Up\n');

addpath(genpath('C:\AvaSpecX64-DLL_9.9.3.0'));
addpath(genpath('C:\Users\DOITuser\Documents\GitHub\DOIT-internDoc\LydiaSameeRyanHan_senior2020'));

% Save filename and folder
dataFolder=['C:/Users/DOITuser/Desktop/DRSdata/' operatorName];
if ~exist(dataFolder, 'dir')
    mkdir(dataFolder);
end
dataFilename=[dataFolder '/' dataPrefix '_' datestr(now,'yymmdd_HHMMSS')];

% Begin instrument communication
try
    %% Configure IO
    fprintf('Configuring IO\n');

    % Switch
    s1=serial('COM10','BaudRate',9600);
%     s1=serial('COM9','BaudRate',9600);
%     s1.Terminator = "CR/LF";
    fopen(s1);
    dtPauseSwitch=0.5; %sec
    
    % LabJack
    [ljudObj, ljhandle, e_NoMoreData]=setupLabJackU3_AI();

    % Light Sources
    a=arduino('COM4');
    dtPauseLight=0.1; %sec
    % Source 1
    configurePin(a, 'D2', 'DigitalOutput')
    writeDigitalPin(a, 'D2', 0)
    % Source 2
    configurePin(a, 'D7', 'DigitalOutput')
    writeDigitalPin(a, 'D7', 0)

    %% Run
    % Start on channel 1 for the detector
%     fprintf(s1,'C&D&C1&');
    fwrite(s1, [0x01 0x12 0x00 0x00]);
%     fscanf(s1);
    D1=true;
    fprintf('Press Any Key to Start\n');
    pause;
    tic;
    fprintf('\n');

    % Taking the initial dark spectrum measurement from the spectrometer
    
    if ~doDark
        tempSpect=take_meas(1, 10, 0);
        dark=tempSpect(:, 2);
        dark=zeros(size(dark));
    else
        tempSpect=take_meas(1, int_time, 0); 
        dark=tempSpect(:, 2);
    end
    
    if doPlot
        figure(1); clf;
    end
    
    t_dark=toc;
    n=0;
    n_AUX=0;
    
    if ~exist('runTime', 'var')
        runTime=inf;
    elseif ~exist('numSamp', 'var')
        numSamp=inf;
    end
    
    while toc<runTime && n<numSamp
        n=n+1;
        fprintf('Sample %d\n', n);
        fprintf('\tCollecting...');
        t_startSamp(n) = toc;
        
        % AUX
        n_AUX=n_AUX+1;
        [AUX(n_AUX, :), ljudObj, ljhandle]=...
            getLabJackU3_AI(ljudObj, ljhandle, e_NoMoreData);
        t_AUX(n_AUX)=toc;
        
        % Turning on source 1
        writeDigitalPin(a, 'D2', 1);
        writeDigitalPin(a, 'D7', 0);
        pause(dtPauseLight);
        % Taking a measurement
        if D1
            fprintf('S1D1...');
            tempSpect=take_meas(meas_num, int_time, dark);
            t_S1D1(n)=toc;
            S1D1(n, :, :)=tempSpect(:, 2:end);
        else
            fprintf('S1D2...');
            tempSpect=take_meas(meas_num, int_time, dark);
            t_S1D2(n)=toc;
            S1D2(n, :, :)=tempSpect(:, 2:end);
        end
        
        % AUX
        n_AUX=n_AUX+1;
        [AUX(n_AUX, :), ljudObj, ljhandle]=...
            getLabJackU3_AI(ljudObj, ljhandle, e_NoMoreData);
        t_AUX(n_AUX)=toc;
        
        % Changing to source 2
        writeDigitalPin(a, 'D2', 0)
        writeDigitalPin(a, 'D7', 1)
        pause(dtPauseLight)
        % Taking a measurement
        if D1
            fprintf('S2D1...');
            tempSpect=take_meas(meas_num, int_time, dark);
            t_S2D1(n)=toc;
            S2D1(n, :, :)=tempSpect(:, 2:end);
        else
            fprintf('S2D2...');
            tempSpect=take_meas(meas_num, int_time, dark);
            t_S2D2(n)=toc;
            S2D2(n, :, :)=tempSpect(:, 2:end);
        end
        fprintf('\n');
        
        % AUX
        n_AUX=n_AUX+1;
        [AUX(n_AUX, :), ljudObj, ljhandle]=...
            getLabJackU3_AI(ljudObj, ljhandle, e_NoMoreData);
        t_AUX(n_AUX)=toc;
        
        % Switching detectors
        if D1
        % Switch to channel 2
            fprintf('\tSwitching D1 -> D2\n');
%             fprintf(s1,'C&D&C2&');
            fwrite(s1, [0x01 0x12 0x00 0x01]);
%             fscanf(s1);
        else
        % Switch to channel 1
            fprintf('\tSwitching D2 -> D1\n');
%             fprintf(s1,'C&D&C1&');
            fwrite(s1, [0x01 0x12 0x00 0x00]);
%             fscanf(s1);
        end
        
        % Plot while switching 
        swSt=toc;
        if doPlot && n>1
            lambda=tempSpect(:, 1);
            t_midSamp=(t_stopSamp+t_startSamp(1:(end-1)))/2;

            figure(1);
            subplot(2, 2, 1);
            plot(lambda, mean(S1D1(n-1, :, :), 3), 'r-'); hold on;
            plot(lambda, mean(S1D2(n-1, :, :), 3), 'b-');
            plot(lambda, mean(S2D1(n-1, :, :), 3), 'b--');
            plot(lambda, mean(S2D2(n-1, :, :), 3), 'r--'); hold off;
            set(gca, 'YScale', 'log');
            xlim([lamRng(1), lamRng(2)]);
            xlabel('\lambda (nm)');
            ylabel('SDI (counts)');
            legend('S1D1', 'S1D2', 'S2D1', 'S1D2', 'location', 'best');
            drawnow;
            
            [~, i1]=min(abs(lambda-lamRng(1)));
            [~, i2]=min(abs(lambda-lamRng(2)));
            CWdata.lambda=lambda(i1:i2);
            CWdata.rhos=rhos;
            CWdata.II=[...
                mean(S1D1(n-1, i1:i2, :), 3).',...
                mean(S1D2(n-1, i1:i2, :), 3).',...
                mean(S2D2(n-1, i1:i2, :), 3).',...
                mean(S2D1(n-1, i1:i2, :), 3).'];
            mua=analyzeCWdata(CWdata, musp0*(CWdata.lambda/lam0).^-b);
            E=makeE('ODWL', CWdata.lambda);
            C(:, n-1)=E\mua;
            muaHat=E*C(:, n-1);
            
            subplot(2, 2, 2);
            yyaxis left;
            plot(CWdata.lambda, mua, 'k-'); hold on;
            plot(CWdata.lambda, muaHat, 'k:'); hold off;
            ylabel('\mu_a (1/mm)');
            yyaxis right;
            plot(CWdata.lambda, mua-muaHat, 'm-');
            ylabel('\mu_a-\mu_{a,fit} (1/mm)');
            xlim([lamRng(1), lamRng(2)]);
            xlabel('\lambda (nm)');
            legend('Meas', 'Fit', 'Res', 'location', 'best');
            ax=gca;
            ax.YAxis(1).Color='k';
            ax.YAxis(2).Color='m';
            drawnow;
            
            subplot(2, 2, 3);
            yyaxis left;
            plot(t_midSamp, C(1, :), '-or'); hold on;
            plot(t_midSamp, C(2, :), '-ob'); hold off;
            ylabel('O & D (\muM)');
            yyaxis right;
            plot(t_midSamp, C(3, :), '-oc'); hold on;
            plot(t_midSamp, C(4, :), '-og'); hold off;
            ylabel('W & L (v/v)');
            xlabel('t (sec)');
            legend('O', 'D', 'W', 'L', 'location', 'southwest');
            ax=gca;
            ax.YAxis(1).Color='k';
            ax.YAxis(2).Color='k';
            drawnow;
            
            subplot(2, 2, 4);
            plot(t_AUX, AUX);
            xlabel('t (sec)');
            ylabel('V');
            legend('AI0', 'AI1', 'AI2', 'AI3', 'location', 'southwest');
            drawnow;
            
            sgtitle(sprintf('t=%.1f min, n=%d', t_midSamp(end)/60, n-1));
        end
        
        % AUX
        n_AUX=n_AUX+1;
        [AUX(n_AUX, :), ljudObj, ljhandle]=...
            getLabJackU3_AI(ljudObj, ljhandle, e_NoMoreData);
        t_AUX(n_AUX)=toc;
        
        while (toc-swSt)<dtPauseSwitch
            pause(0.001);
        end
        D1=~D1;

        fprintf('\tCollecting...')
        % Taking a measurement
        if D1
            fprintf('S2D1...');
            tempSpect=take_meas(meas_num, int_time, dark);
            t_S2D1(n)=toc;
            S2D1(n, :, :)=tempSpect(:, 2:end);
        else
            fprintf('S2D2...');
            tempSpect=take_meas(meas_num, int_time, dark);
            t_S2D2(n)=toc;
            S2D2(n, :, :)=tempSpect(:, 2:end);
        end
        
        % AUX
        n_AUX=n_AUX+1;
        [AUX(n_AUX, :), ljudObj, ljhandle]=...
            getLabJackU3_AI(ljudObj, ljhandle, e_NoMoreData);
        t_AUX(n_AUX)=toc;
        
        % Changing to source 1
        writeDigitalPin(a, 'D2', 1);
        writeDigitalPin(a, 'D7', 0);
        pause(dtPauseLight)
        % Taking a measurement
        if D1
            fprintf('S1D1...');
            tempSpect=take_meas(meas_num, int_time, dark);
            t_S1D1(n)=toc;
            S1D1(n, :, :)=tempSpect(:, 2:end);
        else
            fprintf('S1D2...');
            tempSpect=take_meas(meas_num, int_time, dark);
            t_S1D2(n)=toc;
            S1D2(n, :, :)=tempSpect(:, 2:end);
        end
        t_stopSamp(n)=toc;
        fprintf('\n');
    end

    lambda=tempSpect(:, 1);

    dt_samp=t_stopSamp-t_startSamp;
    t_midSamp=(t_stopSamp+t_startSamp)/2;

    %% Save
    save([dataFilename '.mat'],...
        'S1D1', 'S1D2', 'S2D1', 'S2D2',...
        't_S1D1', 't_S1D2', 't_S2D1', 't_S2D2',...
        't_startSamp', 't_stopSamp', 't_midSamp', 'dt_samp',...
        'lambda', 'dark',...
        't_AUX', 'AUX');
    errOcc=false;
catch ERR
    fprintf('An error occurred\n');
    errOcc=true;
end

%% Cleanup
if ~errOcc
    writeDigitalPin(a, 'D2', 0)
    writeDigitalPin(a, 'D7', 0)
end
fprintf('End of Data Collection\n')
fclose(s1);

if errOcc
    rethrow(ERR)
end

%% Functions
function [mua, CWdata_out]=analyzeCWdata(CWdata_in, musp)
% [mua, CWdata_out]=analyzeCWdata(CWdata_in, musp)
% Giles Blaney Spring 2021
% Inputs:   - CWdata_in: Struct with the following fields:
%               - lambda: Vector of wavelengths
%               - rhos: Vector source detector distances (4 long)
%               - II: Vector of real reflectance (4 long, same order as
%                     rhos)
%           - musp: Vector of assumed reduced scattering coefficients
% 
% Outputs:  - mua: Vector of recovered absorption coefficients
%           - CWdata_out: Struct with same fields as CWdata_in with the
%                         following added:
%               - mua: Same as mua output
    
    lambda=CWdata_in.lambda;
    rhos=CWdata_in.rhos;
    II=CWdata_in.II;
    
    mua=NaN(size(lambda));
    for lamInd=1:length(lambda)
        [mua(lamInd), iter]=...
            DSI2muaEB_iterRecov(rhos, II(lamInd, :), musp(lamInd));
        if iter.n>=10
            mua(lamInd)=NaN;
        end
    end
    
    CWdata_out=CWdata_in;
    CWdata_out.mua=mua;
    
end