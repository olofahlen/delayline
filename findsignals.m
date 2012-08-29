%% Initial

%Toolboxes: DSP, SIgnal Processing, Wavelet
clear all
clc
tic

% This file loads waveforms (Time + Amplitude) saved by the WavePro7100
% The variable 'path' below should be the absolute path to a folder
% containing measurements in the form 'C?*.dat', for example
% 'C1mcp00012.dat'.

%% Load data

importSavedData = 1;
saveData = 1;

if importSavedData
    disp('Loading saved data...')
    load mcpData3
else
    path = '/home/thorleif/mcp/tests/lotsofmeas/';
    %path = '/home/thorleif/mcp/tests/gooddata/';
    channels = 4;

    C1files = dir([path 'C1*.dat']);
    nbrOfFiles = length(C1files);

    dummyData = importdata([path C1files(1).name]);
    measPerFile = length(dummyData);
    nbrOfMeas = length(C1files);
    %nbrOfMeas = 50;

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
        save('mcpData3')
    end
end

%% Settings
plotOffsets = 1;
plotFourierTransform = 1;
plotCharges = 1;
plotManyPulses = 1;
plotFittedPeaks = 1;
plotSignals = 1;
plotPositions = 1;

chosenChannel = 1;
chosenSignal = 1;

%% Post Loading

channelPairs = [1 2 3 4]; %1 2 3 4 is the correct configuration
channelGroups = [channelPairs(1:2); channelPairs(3:4)];
freqCut = 0.2e9;

colors = ['y', 'r', 'b', 'g'];

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

%% Remove offsets

if plotOffsets
    offsetPlot = figure(22);
    clf(22)
    set(gcf, 'Name', 'Signal Offsets')
    subplot(2, 1, 1)
    hold on
    title('Before removing offset')
    for j = 1:channels
        plot(T, data(:, 1, channelPairs(j)), colors(channelPairs(j)))
    end
    xlabel('Time [s]')
    ylabel('Voltage [V]')
end

pedestal = mean(data(1:floor(measPerFile/15), :, :));
data = bsxfun(@minus, data, pedestal);

if plotOffsets
    subplot(2, 1, 2)
    hold on
    title('After removing offset')
    i = 1;
    j = 1;
    for j = 1:channels
        plot(T, data(:, 1, channelPairs(j)), colors(channelPairs(j)))
    end
    xlabel('Time [s]')
    ylabel('Voltage [V]')
end

%% Clean signals from noise using the Fourier Transform

disp('Cleaning with Fourier Transform...')

if plotFourierTransform
    fourierPlot = figure(11);
    clf(11)
    set(gcf, 'Name', 'Signal Fourier Transform')
    hold on
    subplot(2, 2, 1)
    hold on
    title('Before low pass')
    xlabel('Time [s]')
    ylabel('Voltage [V]')
    plot(T, data(:, chosenSignal, chosenChannel))
end

%Memory may be saved here by overwriting the data-variable instead of
%creating a new variable (DATA).
DATA = fft(data)/L;
if plotFourierTransform
    subplot(2, 2, 2)
    semilogy(f, 2*abs(DATA(1:L/2+1, i, j))) 
    hold on
    title('Uncut Fourier spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Fourier transform [Vs]')
end
%DATA(L/2 - fCutLength + 1:L/2 + fCutLength, :, :) = fZeroMask;
%data = real(ifft(DATA))*L;
sc = t/4e-9;
data = filter(sc, [1 sc-1], data);
DATA = fft(data)/L;
if plotFourierTransform
    subplot(2, 2, 3)
    hold on
    title('After low pass')
    xlabel('Time [s]')
    ylabel('Voltage [V]')
    plot(T, data(:, chosenSignal, chosenChannel));
    subplot(2, 2, 4)
    semilogy(f, 2*abs(DATA(1:L/2+1, i, j))) 
    hold on
    title('Cut Fourier spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Fourier transform [Vs]')
    suptitle(['Before and after frequency cut at ' num2str(freqCut, '%.2e')])
end

%% Filter out bad measurements

disp('Removing offsets and finding bad signals...')

good = ones(size(data, 2), 1);

[minVal minIndex] = min(data);;
potentialStart = squeeze(minIndex - nRiseTime);
[row col] = find(potentialStart < measPerFile/15);
good(row) = 0;
[row col] = find(potentialStart > measPerFile - (2*nRiseTime + 1));
good(row) = 0;
%Performance can be improved here by removing the bad measurements before
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

if plotOffsets
    suptitle('Before and after removing the offset')
end

%% Calculate charge

disp('Calculating total charge...')

[minValues minIndices] = min(data);
minIndices = squeeze(minIndices);

eCharge = -1.602e-19;
charge = zeros(nbrOfMeas, channels);

for i = 1:nbrOfMeas
    for j = 1:channels
        interval = minIndices(i, j) - nRiseTime : minIndices(i, j) + floor(1.1*nRiseTime);
        charge(i, j) = sum(data(interval, i, j));
    end
end

charge = charge * t / (inputImpedance * eCharge);
totalCharge = [sum(charge(:, [channelGroups(1, :)]), 2) sum(charge(:, [channelGroups(2, :)]), 2)];

if plotCharges
    bins = 100;
    interval = linspace(0, 14e6, 100);
    individualChargePlot = figure(30);
    clf(30)
    set(gcf, 'Name', 'Individual charge histograms')
    for j = 1:channels
        subplot(4, 1, j)
        hold on
        title(['Charge deposited on channel ' num2str(channelPairs(j))])
        xlabel('Charge [e]')
        ylabel('Counts')
        hist(charge(:, j), interval)
        %hist(charge(:, j), bins)
    end

    interval = linspace(0, 3e7, 100);
    totalChargePlot = figure(31);
    clf(31)
    set(gcf, 'Name', 'Total charge histograms')
    suptitle('Histograms of total charge for the two delay lines')
    for k = 1:2
        subplot(2, 1, k)
        hold on
        title(['Total charge deposited on channels ' num2str(channelGroups(k, 1)) ' and ' num2str(channelGroups(k, 2))])
        xlabel('Charge [e]')
        ylabel('Counts')
        hist(totalCharge(:, k), interval)
        %hist(totalCharge(:, k), bins)
    end
end
if plotCharges
    suptitle('Histograms of charges for the different channels')
    chargeScatterPlot = figure(32);
    clf(32)
    set(gcf, 'Name', 'Charge scatter plot')
    scatter(totalCharge(:, 1), totalCharge(:, 2))
    xlabel(['Total charge deposited on channels ' num2str(channelGroups(1, 1)) ' and ' num2str(channelGroups(1, 2)) ' [e]'])
    ylabel(['Total charge deposited on channels ' num2str(channelGroups(2, 1)) ' and ' num2str(channelGroups(2, 2)) ' [e]'])
end

%% Calculate FDHM of pulses. Good idea to check the code below...

[minVals minIndices] = min(data);
minVals = squeeze(minVals);
threshold = minVals / 2;
longer = zeros(nbrOfMeas, channels);
longerX = longer;
loopCounter = 0;
for j = 1:channels
    for i = 1:nbrOfMeas
        longerX(i, j) = find(data(:, i, j) < threshold(i, j), 1, 'first');
        longer(i, j) = longerX(i, j) + measPerFile*loopCounter;
        loopCounter = loopCounter + 1;
    end
end
shorter = longer - 1;
t1 = longerX + 1./(data(longer) - data(shorter)).*(threshold - data(longer));

loopCounter = 0;
for j = 1:channels
    for i = 1:nbrOfMeas
        longerX(i, j) = find(data(:, i, j) < threshold(i, j), 1, 'last') + 1;
        longer(i, j) = longerX(i, j) + measPerFile*loopCounter;
        loopCounter = loopCounter + 1;
    end
end
shorter = longer - 1;
t2 = longerX + 1./(data(longer) - data(shorter)).*(threshold - data(longer));
fdhm = (t2 - t1)*t;

fdhmHistPlot = figure(35);
clf(fdhmHistPlot)
set(gcf, 'Name', 'FDHM of the four channels')
hold on
cutFdhm = 1.15e-8;
for j = 1:channels
    subplot(4, 1, j)
    hold on
    title(['FDHM for channel ' num2str(channelPairs(j))])
    xlabel('FDHM [s]')
    ylabel('Counts')
    hist(fdhm(find(fdhm(:, j) < cutFdhm), j), 200)
end
suptitle(['Distribution of FDHM. This distribution is cut at ' num2str(cutFdhm) ' s.'])

%% Calculate pulse shape by averaging. This needs some more work, like cutting from the histogram fdhmHistPlot

disp('Calculating mean pulse shape...')
[foo mins] = min(data);
mins = squeeze(mins);
meanPulsePlot = figure(40);
clf(40)
set(gcf, 'Name', 'Mean pulse calculation')
subplot(1, 2, 1)
xlabel('Shifted time indices')
ylabel('Voltage [V]')
hold on
title(['Pulses from channel ' num2str(chosenChannel) ' overlaid'])
pulseShaper = zeros(4*nRiseTime + 1, size(data, 2), size(data, 3));
for i = 1:nbrOfMeas
    for j = 1:channels
        nRange = (mins(i, j) - 2*nRiseTime):(mins(i, j) + 2*nRiseTime);
        %plot(T(nRange), meas(nRange, i))
        pulseShaper(:, i, j) = data(nRange, i, j);
        if j == chosenChannel
            plot(pulseShaper(:, i, j), colors(mod(i, 4) + 1))
        end
    end
end

pulse = squeeze(mean(pulseShaper, 2));
subplot(1, 2, 2)
hold on
title(['Mean pulse for channel ' num2str(chosenChannel)])
xlabel('Shifted time indices')
ylabel('Voltage [V]')
plot(pulse(:, 1))
suptitle('Calculation of the average pulse')

%% Locate peaks

disp('Locating peaks...')
signalIndices = zeros(nbrOfMeas, 4);
signals = zeros(nbrOfMeas, 4);

for i = 1:nbrOfMeas
    for j = 1:channels
        meas = data(:, i, j);
        [minValue minIndex] = min(meas);
        interval = [minIndex - 2:minIndex + 2];
        [p, S, mu] = polyfit(T(interval), meas(interval), 2);
        minT = -p(2)/(2*p(1)) * mu(2) + mu(1);
        signalIndices(i, j) = minIndex;
        signals(i, j) = minT;
    end
end

if plotFittedPeaks
    fittedPeakPlot = figure(400);
    clf(400)
    set(gcf, 'Name', 'Fitting parabola')
    hold on
    title('Fitting of a parbola to find the true minimum')
    plot(T(interval), meas(interval))
    fineT = linspace(T(interval(1)), T(interval(end)), 100);
    [fittedData delta] = polyval(p, fineT, S, mu);
    plot(fineT, fittedData, 'r')
    minV = polyval(p, minT, S, mu);
    plot(minT, polyval(p, minT, S, mu), 'g*')
    xlabel('Time [s]')
    ylabel('Voltage [V]')
end


%% Plot signals

%Look into correlation between signal heights and delays

if plotSignals
    disp('Plotting signals...')
    signalPlot = figure(1);
    clf(1)
    set(gcf, 'Name', 'Signal plots')
    pic = 1;
    for i = chosenSignal:chosenSignal
        for j = 1:channels
            color = colors(j);
            meas = data(:, i, channelPairs(j));j
            subplot(2, 1, ceil(j/2));
            hold on
            title(['Channels ' num2str(channelGroups(ceil(j/2), 1)) ' and ' num2str(channelGroups(ceil(j/2), 2))])
            xlabel('Time [s]')
            ylabel('Voltage [V]')
            plot(T, meas, color)
            plot(T(signalIndices(i, channelPairs(j))), meas(signalIndices(i, channelPairs(j))), 'o')
            %The y-value in the following plot is not exact
            plot(signals(i, j), meas(signalIndices(i, channelPairs(j))), '*')
        end
        %pause
        %clf(1)
    end
    suptitle('Delay Line signals')
end

%% Calculate total time

disp('Calculating sums of times...')
timeSum = zeros(nbrOfMeas, 2);

%timeSum = [sum(signalIndices(:, channelGroups(1, :)), 2) sum(signalIndices(:, channelGroups(2, :)), 2)];
timeSum = [sum(signals(:, channelGroups(1, :)), 2) sum(signals(:, channelGroups(2, :)), 2)];

disp('Plotting results in histogram and x-y plot...')
bins = 500;

if plotPositions
    timeSumHistPlot = figure(200);
    clf(200)
    set(gcf, 'Name', 'Normalized histograms of time sums')
    for k = 1:2
        subplot(2, 1, k)
        hold on
        title(['Delayline for channels ' num2str(channelGroups(k, 1)) ' and ' num2str(channelGroups(k, 2))])
        xlabel('$t_1 + t_2$ [s]', 'Interpreter', 'LaTeX')
        ylabel('Normalized counts')
        [xo, no] = histnorm(timeSum(:, k), bins, 'plot');
        tMean = mean(timeSum(:, k));
        tStd = std(timeSum(:, k));
        plot(no, normpdf(no, tMean, tStd), 'r')
    end
    suptitle('Normalized histograms of time sums for the two delay lines')
end

%% Calculate positions

disp('Calculating spatial coordinates...')

%timeDiff = [diff(signalIndices(:, channelGroups(1, :)), 2) diff(signalIndices(:, channelGroups(2, :)), 2)];
timeDiff = [diff(signals(:, channelGroups(1, :)), 1, 2) diff(signals(:, channelGroups(2, :)), 1, 2)];

%This minus sign is arbitrary, only mirrors the image in the origin.
timeDiff = -timeDiff;

disp('Plotting results in histogram and x-y plot...')
bins = 500;

if plotPositions
    timeDiffHistPlot = figure(2);
    clf(2)
    set(gcf, 'Name', 'Histograms of time differences')
    subplot(2, 1, 1)
    hold on
    title(['Delayline for channels ' num2str(channelPairs(1)) ' and ' num2str(channelPairs(2))])
    xlabel('$\Delta t$ [s]', 'Interpreter', 'LaTeX')
    ylabel('Counts')
    hist(timeDiff(:, 1), bins)
    subplot(2, 1, 2)
    hold on
    title(['Delayline for channels ' num2str(channelPairs(3)) ' and ' num2str(channelPairs(4))])
    xlabel('$\Delta t$ [s]', 'Interpreter', 'LaTeX')
    ylabel('Counts')
    hist(timeDiff(:, 2), bins)
    suptitle('Histograms of time differences for the two delay lines')

    mcpHitmapPlot = figure(4);
    clf(4)
    set(gcf, 'Name', 'MCP 2D-plot')
    hold on
    suptitle('Reconstruction of particle hits on the MCP')
    scatter(timeDiff(:, 1), timeDiff(:, 2), 20, sum(totalCharge,2 ), 'filled')
    xlabel('$x\propto \Delta t_x$', 'Interpreter', 'LaTeX');
    ylabel('$y\propto \Delta t_y$', 'Interpreter', 'LaTeX');
    axis square
    cb = colorbar;
    ylabel(cb, 'Charge [e]')
end


%% Select the events corresponding to the left and right peaks of the time sums
cut1 = 9.16e-8;
less1 = find(timeSum(:, 1) < cut1);
more1 = find(timeSum(:, 1) > cut1);

cut2 = 9.6e-8;
less2 = find(timeSum(:, 2) < cut2);
more2 = find(timeSum(:, 2) > cut2);

timeCutHitmapPlot = figure(301);
clf(301)
hold on
set(gcf, 'Name', 'MCP 2D-plot with time cuts')

subplot(1, 2, 1)
hold on
title(['Cut at ' num2str(cut1) 's for channels ' num2str(channelGroups(1, 1)) ' and ' num2str(channelGroups(1, 2))])
scatter(timeDiff(less1, 1), timeDiff(less1, 2), 'b')
scatter(timeDiff(more1, 1), timeDiff(more1, 2), 'r')
xlabel('$x\propto \Delta t_x$', 'Interpreter', 'LaTeX');
ylabel('$y\propto \Delta t_y$', 'Interpreter', 'LaTeX');
legend('Short times', 'Long times')
axis square

subplot(1, 2, 2)
hold on
title(['Cut at ' num2str(cut2) 's for channels ' num2str(channelGroups(2, 1)) ' and ' num2str(channelGroups(2, 2))])
scatter(timeDiff(less2, 1), timeDiff(less2, 2), 'b')
scatter(timeDiff(more2, 1), timeDiff(more2, 2), 'r')
xlabel('$x\propto \Delta t_x$', 'Interpreter', 'LaTeX');
ylabel('$y\propto \Delta t_y$', 'Interpreter', 'LaTeX');
legend('Short times', 'Long times')
axis square
suptitle('Events cut in time histograms')

%%
toc
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
    suptitle('Delay Line signals')
else
    disp('Could not find event.')
end

%% Save figures

figures = {'fourierPlot', 'offsetPlot', 'individualChargePlot', 'totalChargePlot', 'chargeScatterPlot', 'fdhmHistPlot', 'meanPulsePlot', 'fittedPeakPlot', 'signalPlot', 'timeSumHistPlot', 'timeDiffHistPlot', 'mcpHitmapPlot', 'timeCutHitmapPlot'};
formats = {'png', 'epsc', 'fig', 'pdf'};
path = '/home/thorleif/mcp/tests/pics/matlab/';

for k = 1:length(figures)
    figureName = figures{k};
    for l = 1:length(formats)
        format = formats{l};
        disp(['Saving ' figureName ' as ' format '...'])
        saveas(eval(figureName), [path figureName], format)
    end
end