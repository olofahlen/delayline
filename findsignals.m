%% Initial

%NOTE TO SELF: Work now on optimizing the code by doing away with all these
%horrible for loops, now that data is a 3-dim array.

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
    
    data = zeros(measPerFile, nbrOfMeas, channels);
    disp(['Loading ' num2str(nbrOfMeas*channels) ' files...'])
    
    modCheck = max(floor(nbrOfMeas/100), 1);
    fprintf(1, '  0%% done\n')
    for i = 1:nbrOfMeas
        measFile = C1files(i);
        fileName = measFile.name;
        for j = 1:channels
            fp = [path 'C' int2str(j) fileName(3:end)];
            importedData = importdata(fp);
            data(:, i, j) = importedData(:, 2);
        end
        if mod(i, modCheck) == 0
            percentProgress = ceil(i/nbrOfMeas*100);
            fprintf(1, '\b\b\b\b\b\b\b\b\b\b%3d%% done\n', percentProgress)
        end
    end
    if saveData
        disp('Saving data...')
        save mcpData
    end
end

%% Post Loading

channelPairs = [1 2 3 4]; %1 2 3 4 is the correct configuration
channelGroups = [channelPairs(1:2); channelPairs(3:4)];
freqCut = 0.2e9;

inputImpedance = 50; %Impedance of the oscilloscope in Ohms

T = dummyData(:, 1); %Time vector
t = T(2) - T(1); %Sampling time
Fs = 1/t; %Sampling frequency
L = measPerFile; %Length of signal
f = Fs/2*linspace(0,1,L/2+1); %Frequency vector
fCutLength = length(f) - find(f > freqCut, 1, 'first');
fZeroMask = zeros(2*fCutLength, size(data, 2), size(data, 3));

riseTime = 1e-8;
nRiseTime = floor(riseTime/t);

%% Clean signals from noise using the Fourier Transform

disp('Cleaning with Fourier Transform...')

if plotFourierTransform
    i = 1;
    j = 1;
    figure(11)
    clf(11)
    set(gcf, 'Name', 'Signal Fourier Transform')
    suptitle(['Before and after frequency cut at ' num2str(freqCut, '%.2e')])
    hold on
    subplot(2, 2, 1)
    hold on
    title('Before low pass')
    xlabel('Time [s]')
    ylabel('Voltage [V]')
    plot(T, data(:, i, j))
end

%Memory may be saved here by overwriting the data-variable instead of
%creating a new variable (DATA).
DATA = fft(data)/L;
if plotFourierTransform
    subplot(2, 2, 2)
    hold on
    title('Uncut Fourier spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Fourier transform [Vs]')
    plot(f, 2*abs(DATA(1:L/2+1, i, j))) 
end
DATA(L/2 - fCutLength + 1:L/2 + fCutLength, :, :) = fZeroMask;
data = real(ifft(DATA));
if plotFourierTransform
    subplot(2, 2, 3)
    hold on
    title('After low pass')
    xlabel('Time [s]')
    ylabel('Voltage [V]')
    plot(T, data(:, i, j));
    subplot(2, 2, 4)
    hold on
    title('Cut Fourier spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Fourier transform [Vs]')
    plot(f, 2*abs(DATA(1:L/2+1, i, j))) 
end

%% Remove Offsets and filter out bad measurements

disp('Removing offsets and finding bad signals...')

if plotOffSets
    figure(22)
    clf(22)
    set(gcf, 'Name', 'Signal Offsets')
    suptitle('Before and after removing the offset')
end

good = ones(size(data, 2), 1);

[minVal minIndex] = min(data);;
potentialStart = squeeze(minIndex - nRiseTime);
[row col] = find(potentialStart < measPerFile/15);
good(row) = 0;
[row col] = find(potentialStart > measPerFile - (2*nRiseTime + 1));
good(row) = 0;
%Performance can be improved here by removin the bad measurements before
%continuing

for i = 1:nbrOfMeas
    for j = 1:channels
        if good(i)
            meas = data(:, i, j);
            measStd = std(meas(1:potentialStart(i, j)));
            measMean = mean(meas(1:potentialStart(i, j)));
            upperlimit = 4*measStd + measMean;
            lowerlimit = -4*measStd + measMean;
            if length(find(meas(potentialStart(i, j):potentialStart(i, j) + 2*nRiseTime) < lowerlimit)) < nRiseTime/2
                good(i) = 0;
            end
        end
    end
end

goods = find(good == 1);
nbrOfGoods = length(goods);
disp(['Found ' num2str(nbrOfMeas - nbrOfGoods) ' bad signals from signal shape. Removing...'])
data = data(:, goods, :);
nbrOfMeas = size(data, 2);
good = ones(size(data, 2), 1);

[minValues minIndices] = min(data);
minIndices = squeeze(minIndices);
tTot = T(minIndices(:, [channelGroups(1, 1) channelGroups(2, 1)])) + T(minIndices(:, [channelGroups(1, 2) channelGroups(2, 2)]));
tMean = mean(tTot);
tStd = std(tTot);

for k = 1:channels/2
    [row col] = find(abs([tTot(:, k) - tMean(k)]) > 3*tStd(k));
    good(row) = 0;
end

nbrOfGoods = length(find(good == 1));
disp(['Found ' num2str(nbrOfMeas - nbrOfGoods) ' bad signals from time sum. Removing...'])
data = data(:, find(good == 1), :);
nbrOfMeas = size(data, 2);

%% Calculate charge

%FIXME: Units seem not to be right. Expecting something like 1e7 elementary
%charges per event. At least I hope so!!

bins = 50;
disp('Calculating total charge...')
ePerCoulomb = 1/1.602e-19;
%charge = zeros(nbrOfMeas, channels);
charge = squeeze(sum(data))*t / inputImpedance;
totalCharge = [sum(charge(:, [channelGroups(1, :)]), 2) sum(charge(:, [channelGroups(2, :)]), 2)];

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
        hist(charge(:, j)*ePerCoulomb, bins)
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
        hist(totalCharge(:, k)*ePerCoulomb, bins)
    end
end

%% Calculate pulse shape by averaging. This needs some more work

disp('Calculatin mean pulse shape...')
[foo mins] = min(data);
mins = squeeze(mins);
figure(40)
clf(40)
hold on
title('Pulses overlaid')
pulseShaper = zeros(2*nRiseTime + 1, size(data, 2), size(data, 3));
colors = ['b', 'r', 'y', 'g'];
for i = 1:nbrOfMeas
    for j = 1:channels
        nRange = (mins(i, j) - nRiseTime):(mins(i, j) + nRiseTime);
        %plot(T(nRange), meas(nRange, i))
        pulseShaper(:, i, j) = data(nRange, i, j);
        if j == 1
            plot(pulseShaper(:, i, j), colors(mod(i, 4) + 1))
        end
    end
end

pulse = squeeze(mean(pulseShaper, 2));
figure(41)
clf(41)
hold on
title('Mean pulse for channel 1')
plot(pulse(:, 1))

%% Locate peaks

disp('Locating peaks...')
signalIndices = zeros(nbrOfMeas, 4);
signals = zeros(nbrOfMeas, 4);

for i = 1:nbrOfMeas
    for j = 1:channels
        meas = data(:, i, j);
        %figure(123321)
        %clf(123321)
        %hold on
        [minValue minIndex] = min(meas);
        interval = [minIndex - 2:minIndex + 2];
        %plot(T(interval), meas(interval))
        [p, S, mu] = polyfit(T(interval), meas(interval), 2);
        %fineT = linspace(T(interval(1)), T(interval(end)), 100);
        %[fittedData delta] = polyval(p, fineT, S, mu);
        %plot(fineT, fittedData, 'r')
        minT = -p(2)/(2*p(1)) * mu(2) + mu(1);
        %minV = polyval(p, minT, S, mu);
        %plot(minT, polyval(p, minT, S, mu), 'go')
        signalIndices(i, j) = minIndex;
        signals(i, j) = minT;
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
        for j = 1:channels
            color = colors(j);
            meas = data(:, i, channelPairs(j));
            subplot(2, 1, ceil(j/2));
            hold on
            title('Delay Line signals')
            xlabel('Time [s]')
            ylabel('Voltage [V]')
            plot(T, meas, color)
            plot(T(signalIndices(i, channelPairs(j))), meas(signalIndices(i, channelPairs(j))), 'o')
            %The y-value in the followin plot is not exact
            plot(signals(i, j), meas(signalIndices(i, channelPairs(j))), '*')
        end
        %pause
        %clf(1)
    end
end
return

%% Calculate total time

disp('Calculating sums of times...')
timeSum = zeros(2, nbrOfMeas);

for i = 1:nbrOfMeas
    for k = 1:channels/2
        %timeSum(k, i) = T(signalIndices(channelGroups(k, 1), i)) + T(signalIndices(channelGroups(k, 2), i));
        timeSum(k, i) = signals(channelGroups(k, 1), i) + signals(channelGroups(k, 2), i);
    end
end

disp('Plotting results in histogram and x-y plot...')
bins = 500;

if plotPositions
    figure(200)
    clf(200)
    set(gcf, 'Name', 'Normalized histograms of time sums')
    suptitle('Normalized histograms of time sums for the two delay lines')
    for k = 1:2
        subplot(2, 1, k)
        hold on
        title(['Delayline for channels ' num2str(channelGroups(k, 1)) ' and ' num2str(channelGroups(k, 2))])
        xlabel('$t_1 + t_2$ [s]', 'Interpreter', 'LaTeX')
        ylabel('Normalized counts')
        [xo, no] = histnorm(timeSum(k, :), bins, 'plot')
        tMean = mean(timeSum(k, :));
        tStd = std(timeSum(k, :));
        plot(no, normpdf(no, tMean, tStd), 'r')
    end
end



%% Calculate positions

disp('Calculating spatial coordinates...')
timeDiff = zeros(2, nbrOfMeas);

for i = 1:nbrOfMeas
    for k = 1:2
        timeDiff(k, i) = T(signalIndices(channelGroups(k, 1), i)) - T(signalIndices(channelGroups(k, 2), i));
        %timeDiff(k, i) = signals(channelGroups(k, 1), i) - signals(channelGroups(k, 2), i);
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
    scatter(timeDiff(1, :), timeDiff(2, :), 20, mean(totalCharge), 'filled')
    %surf(timeDiff(1, :), timeDiff(2, :), mean(totalCharge(:, :)))
    xlabel('$x\propto \Delta t_x$', 'Interpreter', 'LaTeX');
    ylabel('$y\propto \Delta t_y$', 'Interpreter', 'LaTeX');
    axis square
end


%% Select the events corresponding to the left and right peaks of the time sums
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
% To use this, select a point in figure <figureNbr> using the data cursor,
% then run this cell

figureNbr = 4;
figure(figureNbr)
datacursormode on
dcmObj = datacursormode(figureNbr);
datacursormode off
i = dcmObj.getCursorInfo.DataIndex;

if i >= 1
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