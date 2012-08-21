%% Initial

clear all
clc

% This file loads waveforms (Time + Amplitude) saved by the WavePro7100
% The variable 'path' below should be the absolute path to a folder
% containing measurements in the form 'C?*.dat', for example
% 'C1mcp00012.dat'.

%% Settings
plotOffSets = 0;
plotFourierTransform = 1;
plotCharges = 1;
plotSignals = 1;
plotPositions = 1;
importSavedData = 1;
saveData = 1;

%set(0,'DefaultFigureNumberTitle', 'off')

%% Load data

if importSavedData
    disp('Loading saved data...')
    load mcpData
else
    path = '/home/thorleif/mcp/tests/lotsofmeas/';
    %path = '/home/thorleif/mcp/tests/gooddata/';
    channels = 4;

    C1files = dir([path 'C1*.dat']);
    nbrOfFiles = length(C1files);

    dummyData = importdata([path C1files(1).name]);
    measPerFile = length(dummyData);
    nbrOfMeas = length(C1files);
    %nbrOfMeas = 1000;

    %This will contain all the measurements. The four channels will be on top
    %of each other.
    
    data = zeros(measPerFile*channels, nbrOfMeas);

    i = 1;
    if nbrOfMeas > 100
        modCheck = floor(nbrOfMeas/100);
    else
        modCheck = 1;
    end

    disp(['Loading ' num2str(nbrOfMeas*channels) ' files...'])
    %for measFile = C1files'
    for i = 1:nbrOfMeas
        measFile = C1files(i);
        fileName = measFile.name;
        for j = 1:channels
            fp = [path 'C' int2str(j) fileName(3:end)];
            importedData = importdata(fp);
            data((1:measPerFile) + (measPerFile * (j - 1)), i) = importedData(:, 2);
        end
        if mod(i, modCheck) == 0
            percentProgress = ceil(i/nbrOfMeas*100);
            disp([num2str(percentProgress) '% done'])
        end
        %i = i + 1;
    end
    if saveData
        save mcpData
    end
end

%% Post Loading

channelPairs = [1 2 3 4]; %Real
channelGroups = [channelPairs(1:2); channelPairs(3:4)];
freqCut = 0.2e9;

inputImpedance = 50; %Impedance of the oscilloscope in Ohms

T = dummyData(:, 1); %Time vector
t = T(2) - T(1); %Sampling time
Fs = 1/t; %Sampling frequency
L = measPerFile; %Length of signal
f = Fs/2*linspace(0,1,L/2+1); %Frequency vector
fCutLength = length(f) - find(f > freqCut, 1, 'first');
fZeroMask = zeros(2*fCutLength, 1);

riseTime = 1e-8;
nRiseTime = floor(riseTime/t);

%% New Remove Offsets and filter out bad measurements

disp('Removing offsets and finding bad signals...')

if plotOffSets
    figure(22)
    clf(22)
    set(gcf, 'Name', 'Signal Offsets')
    suptitle('Before and after removing the offset')
end

good = ones(nbrOfMeas, 1);

for i = 1:nbrOfMeas
    for j = 1:channels
        meas = data((1:measPerFile) + (measPerFile * (j - 1)), i);
        [minVal minIndex] = min(meas);
        potentialStart = minIndex - nRiseTime;
        if potentialStart < measPerFile/15 || measPerFile - potentialStart < 2*nRiseTime + 1
            good(i) = 0;
        else
            measStd = std(meas(1:potentialStart));
            measMean = mean(meas(1:potentialStart));
            upperlimit = 4*measStd + measMean;
            lowerlimit = -4*measStd + measMean;
            if length(find(meas < lowerlimit)) < nRiseTime/2
                good(i) = 0;
            end
        end
        meanCut = find(meas < lowerlimit, 1, 'first');
        data((1:measPerFile) + (measPerFile * (j - 1)), i) = meas - mean(meas(1:meanCut));
        %subplot(2, 1, 2)
        %plot(T, meas)
        %line([T(1) T(end)], [0 0])
        
        if plotOffSets %&& i == 100 && j == 1
            subplot(2, 1, 1)
            hold off
            title('With offset')
            xlabel('Time [s]')
            ylabel('Voltage [V]')
            plot(T, meas)
            line([T(1); T(2)], [measStd; measStd])
            line([T(1) T(end)], [upperlimit upperlimit])
            line([T(1) T(end)], [lowerlimit lowerlimit])
            line([T(1) T(end)], [measMean measMean], 'Color', 'g')
            subplot(2, 1, 2)
            hold off
            title('Without offset')
            xlabel('Time [s]')
            ylabel('Voltage [V]')
            plot(T, data((1:measPerFile) + (measPerFile * (j - 1)), i))
            hold on
            if good(i)
                disp('Good signal')
            else
                disp('Bad signal')
            end
            %pause
        end
    end
end

disp(['Found ' num2str(length(find(good == 0))) ' bad signals'])
disp('Removing bad signals...')
nbrOfGoods = length(find(good == 1));
goodData = zeros(measPerFile*4, nbrOfGoods);
goodsExtracted = 0;
for i = 1:nbrOfMeas
    if good(i)
        goodData(:, goodsExtracted + 1) = data(:, i);
        goodsExtracted = goodsExtracted + 1;
    end
end

data = goodData;
nbrOfMeas = size(data, 2);

%% Clean signals from noise using the Fourier Transform

disp('Cleaning with Fourier Transform...')

if plotFourierTransform
    figure(11)
    clf(11)
    set(gcf, 'Name', 'Signal Fourier Transform')
    suptitle(['Before and after frequency cut at ' num2str(freqCut, '%.2e')])
    hold on
end

fprintf(1, '  0%% done\n')
totalIter = channels*nbrOfMeas;
if channels*nbrOfMeas > 100
    modCheck = floor(totalIter/100);
else
    modCheck = 1;
end
loopCounter = 1;
for i = 1:nbrOfMeas
    for j = 1:channels
        meas = data((1:measPerFile) + (measPerFile * (channelPairs(j) - 1)), i);
        MEAS = fft(meas)/L;
        if plotFourierTransform && i == 1 && j == 1
            subplot(2, 2, 1)
            hold on
            title('Before low pass')
            xlabel('Time [s]')
            ylabel('Voltage [V]')
            plot(T, meas)
            
            subplot(2, 2, 2)
            hold on
            title('Uncut Fourier spectrum')
            xlabel('Frequency (Hz)')
            ylabel('Fourier transform [Vs]')
            plot(f, 2*abs(MEAS(1:L/2+1))) 
        end
        
        MEAS(L/2 - fCutLength + 1:L/2 + fCutLength) = fZeroMask;
        cleanedMeas = real(ifft(MEAS));
        
        if plotFourierTransform && i == 1 && j == 1
            subplot(2, 2, 3)
            hold on
            title('After low pass')
            xlabel('Time [s]')
            ylabel('Voltage [V]')
            plot(T, cleanedMeas);
            subplot(2, 2, 4)
            hold on
            title('Cut Fourier spectrum')
            xlabel('Frequency (Hz)')
            ylabel('Fourier transform [Vs]')
            plot(f, 2*abs(MEAS(1:L/2+1))) 
        end
        
        data((1:measPerFile) + (measPerFile * (channelPairs(j) - 1)), i) = cleanedMeas;
        if mod(loopCounter, modCheck) == 0
            percentProgress = ceil(loopCounter/totalIter*100);
            fprintf(1, '\b\b\b\b\b\b\b\b\b\b%3d%% done\n', percentProgress)
        end
        loopCounter = loopCounter + 1;
    end
end

%% Calculate charge

%FIXME: Units seem not to be right. Expecting something like 1e7 elementary
%charges per event. At least I hope so!!

bins = 50;
disp('Calculating total charge...')
ePerCoulomb = 1/1.602e-19;
charge = zeros(channels, nbrOfMeas);
for j = 1:channels
    meas = data((1:measPerFile) + (measPerFile * (channelPairs(j) - 1)), :);
    charge(j, :) = trapz(meas(:, :), 1)*t / inputImpedance;
end
totalCharge = zeros(channels/2, nbrOfMeas);
totalCharge(1, :) = charge(channelPairs(1), :) + charge(channelPairs(2), :);
totalCharge(2, :) = charge(channelPairs(3), :) + charge(channelPairs(4), :);

if plotCharges
    figure(30)
    clf(30)
    set(gcf, 'Name', 'Individual charge histograms')
    suptitle('Histograms of charges for the different channels')
    for j = 1:channels
        subplot(4, 1, j)
        hold on
        title(['Charge deposited on channel ' num2str(channelPairs(j))])
        xlabel('Charge [e]')
        ylabel('Counts')
        hist(charge(j, :)*ePerCoulomb, bins)
    end

    figure(31)
    clf(31)
    set(gcf, 'Name', 'Total charge histograms')
    suptitle('Histograms of total charge for the two delay lines')
    for k = 1:2
        subplot(2, 1, k)
        hold on
        title(['Total charge deposited on channels ' num2str(channelGroups(k, 1)) ' and ' num2str(channelGroups(k, 2))])
        xlabel('Charge [e]')
        ylabel('Counts')
        hist(totalCharge(k, :)*ePerCoulomb, bins)
    end
end

%% Calculate pulse shape by averaging

%Keep working here:
meas = data(1:measPerFile, :);
[tjo mins] = min(meas);
figure(40)
clf(40)
hold on
title('Pulses overlaid')
pulseShaper = zeros(2*nRiseTime + 1, size(data, 2));
colors = ['b', 'r', 'y', 'g']
for i = 1:length(meas)
%for i = 1:100
    nRange = (mins(i) - nRiseTime):(mins(i) + nRiseTime);
    %plot(T(nRange), meas(nRange, i))
    i
    plot(meas(nRange, i), colors(mod(i, 4) + 1))
    pulseShaper(:, i) = meas(nRange, i);
end
pulse = mean(pulseShaper, 2);
figure(41)
clf(41)
hold on
title('Mean pulse')
plot(pulse)

%% Locate peaks

disp('Locating peaks...')
signalIndices = zeros(4, nbrOfMeas);

for i = 1:nbrOfMeas
    for j = 1:channels
        meas = data((1:measPerFile) + (measPerFile * (j - 1)), i);
        [minValue, minIndex] = min(meas);
        signalIndices(j, i) = minIndex;
    end
end




%% Plot signals

%Look into correlation between signal heights and delays

if plotSignals
    disp('Plotting signals...')
    figure(1)
    clf(1)
    set(gcf, 'Name', 'Signal plots')
    suptitle('Delay Line signals')
    colors = ['r', 'g', 'b', 'y'];
    pic = 1;
    for i = 1:1
        i
        for j = 1:channels
        color = colors(j);
            meas = data((1:measPerFile) + (measPerFile * (channelPairs(j) - 1)), i);
            subplot(2, 1, ceil(j/2));
            hold on
            title('Delay Line signals')
            xlabel('Time [s]')
            ylabel('Voltage [V]')
            %plot(T, data((1:measPerFile) + (measPerFile * (channelPairs(j) - 1)), i), color)
            plot(T, meas, color)
            plot(T(signalIndices(channelPairs(j), i)), meas(signalIndices(channelPairs(j), i)), 'o')
        end
        %pause
        %clf(1)
    end
end

%% Calculate total time

disp('Calculating sums of times...')
timeSum = zeros(2, nbrOfMeas);

for i = 1:nbrOfMeas
    for k = 1:2
        timeSum(k, i) = T(signalIndices(channelGroups(k, 1), i)) + T(signalIndices(channelGroups(k, 2), i));
    end
end

disp('Plotting results in histogram and x-y plot...')
bins = 500;

if plotPositions
    figure(200)
    clf(200)
    set(gcf, 'Name', 'Normalized histograms of time sums')
    suptitle('Normalized histograms of time sums for the two delay lines')
    subplot(2, 1, 1)
    hold on
    title(['Delayline for channels ' num2str(channelPairs(1)) ' and ' num2str(channelPairs(2))])
    xlabel('$t_1 + t_2$ [s]', 'Interpreter', 'LaTeX')
    ylabel('Normalized counts')
    histnorm(timeSum(1, :), bins)
    subplot(2, 1, 2)
    hold on
    title(['Delayline for channels ' num2str(channelPairs(3)) ' and ' num2str(channelPairs(4))])
    xlabel('$t_1 + t_2$ [s]', 'Interpreter', 'LaTeX')
    ylabel('Normalized counts')
    histnorm(timeSum(2, :), bins)
end



%% Calculate positions

disp('Calculating spatial coordinates...')
timeDiff = zeros(2, nbrOfMeas);

for i = 1:nbrOfMeas
    for k = 1:2
        timeDiff(k, i) = T(signalIndices(channelGroups(k, 1), i)) - T(signalIndices(channelGroups(k, 2), i));
    end
end

disp('Plotting results in histogram and x-y plot...')
bins = 500;

if plotPositions
    figure(2)
    clf(2)
    set(gcf, 'Name', 'Histograms of time differences')
    suptitle('Histograms of time differences for the two delay lines')
    subplot(2, 1, 1)
    hold on
    title(['Delayline for channels ' num2str(channelPairs(1)) ' and ' num2str(channelPairs(2))])
    xlabel('$\Delta t$ [s]', 'Interpreter', 'LaTeX')
    ylabel('Counts')
    hist(timeDiff(1, :), bins)
    subplot(2, 1, 2)
    hold on
    title(['Delayline for channels ' num2str(channelPairs(3)) ' and ' num2str(channelPairs(4))])
    xlabel('$\Delta t$ [s]', 'Interpreter', 'LaTeX')
    ylabel('Counts')
    hist(timeDiff(2, :), bins)

    figure(4)
    clf(4)
    set(gcf, 'Name', 'MCP 2D-plot')
    hold on
    suptitle('Reconstruction of particle hits on the MCP')
    plot(timeDiff(1, :), timeDiff(2, :), '.')
    xlabel('$x\propto \Delta t_x$', 'Interpreter', 'LaTeX');
    ylabel('$y\propto \Delta t_y$', 'Interpreter', 'LaTeX');
    axis square
end


%%
cut1 = 9.12e-8;
less1 = find(timeSum(1, :) < cut1);
more1 = find(timeSum(1, :) > cut1);
cut2 = 9.6e-8;
less2 = find(timeSum(2, :) < cut2);
more2 = find(timeSum(2, :) > cut2);

figure(301)
clf(301)
set(gcf, 'Name', 'MCP 2D-plot divided events, low times')
title('Events with high low sums')
hold on
plot(timeDiff(1, less1), timeDiff(2, less1), '.b')
plot(timeDiff(1, less2), timeDiff(2, less2), '.r')
xlabel('$x\propto \Delta t_x$', 'Interpreter', 'LaTeX');
ylabel('$y\propto \Delta t_y$', 'Interpreter', 'LaTeX');
axis square

figure(302)
clf(302)
set(gcf, 'Name', 'MCP 2D-plot divided events, high times')
title('Events with high time sums')
hold on
plot(timeDiff(1, more1), timeDiff(2, more1), '.b')
plot(timeDiff(1, more2), timeDiff(2, more2), '.r')
xlabel('$x\propto \Delta t_x$', 'Interpreter', 'LaTeX');
ylabel('$y\propto \Delta t_y$', 'Interpreter', 'LaTeX');
axis square

return
%% Run this cell to plot a selected event

figureNbr = 4;
figure(figureNbr)
dcmObj = datacursormode(figureNbr);
returnData = dcmObj.getCursorInfo.Position;
returnData = [1 2];
xT = timeDiff(1, :);
yT = timeDiff(2, :);
[foo xevent] = min(abs(xT - returnData(1)));
[foo yevent] = min(abs(yT - returnData(2)));

if xevent == yevent
    i = xevent;
    disp(['Found event! Number: ' num2str(i) '. Plotting...'])
    figure(1000)
    clf(1000)
    set(gcf, 'Name', 'Single event plot')
    suptitle('Delay Line signals')
    colors = ['r', 'g', 'b', 'y'];
    for j = 1:channels
    color = colors(j);
        meas = data((1:measPerFile) + (measPerFile * (channelPairs(j) - 1)), i);
        subplot(2, 1, ceil(j/2));
        hold on
        title('Delay Line signals')
        xlabel('Time [s]')
        ylabel('Voltage [V]')
        plot(T, meas, color)
        plot(T(signalIndices(channelPairs(j), i)), meas(signalIndices(channelPairs(j), i)), 'o')
    end
else
    disp('Could not find event.')
end