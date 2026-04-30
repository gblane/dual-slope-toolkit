%% Setup
% Giles Blaney Summer 2021
clear; home; 
% Name of experiment
dataPrefix='legTest'; % <---------------------------------------------------
% Your name
operatorName='Giles'; % <---------------------------------------------------
% Data will be saved in:
% ~/Desktop/DRSdata/*operatorName*/*dataPrefix*_YYMMDD_hhmmss.mat

% The faster (less time) of the two determines run time
% Number of DS samples
numSamp=1; % <----------------------------------------------------------
% Total run time                           OR
runTime=3000; %sec <---------------------------------------------------------

% Number of repeated SD integrations for averaging
meas_num=1;
% Length of SD integration
int_time=3000; %ms <--------------------------------------------------------
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
    
    t_dark=toc;
    n=0;
    
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
        pause(dtPauseSwitch);
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

        if doPlot
            lambda=tempSpect(:, 1);

            figure(1); clf;
            plot(lambda, mean(S1D1(n, :, :), 3), 'r-'); hold on;
            plot(lambda, mean(S1D2(n, :, :), 3), 'b-');
            plot(lambda, mean(S2D1(n, :, :), 3), 'b--');
            plot(lambda, mean(S2D2(n, :, :), 3), 'r--'); hold off;
            set(gca, 'YScale', 'log');
            xlim([lambda(1), lambda(end)]);
            xlabel('\lambda (nm)');
            ylabel('SDI (counts)');
            legend('S1D1', 'S1D2', 'S2D1', 'S1D2', 'location', 'northwest');
            title(sprintf('t=%.1f min, n=%d', toc/60, n));
            drawnow;
        end
    end

    lambda=tempSpect(:, 1);

    dt_samp=t_stopSamp-t_startSamp;
    t_midSamp=(t_stopSamp+t_startSamp)/2;

    %% Save
    save([dataFilename '.mat'],...
        'S1D1', 'S1D2', 'S2D1', 'S2D2',...
        't_S1D1', 't_S1D2', 't_S2D1', 't_S2D2',...
        't_startSamp', 't_stopSamp', 't_midSamp', 'dt_samp',...
        'lambda', 'dark');
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