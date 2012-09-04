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
    load mcpData
else
    path = '/home/thorleif/mcp/tests/hugeshortmeas/';
    %path = '/home/thorleif/mcp/tests/gooddata/';
    channels = 4;

    C1files = dir([path 'C1*.dat']);
    nbrOfFiles = length(C1files);

    timeAndAmplitudeMeas = 'timeandamplitudehugeshortmeas.dat';
    dummyData = importdata([path timeAndAmplitudeMeas]);
    measPerFile = length(dummyData);
    nbrOfMeas = length(C1files);
    nbrOfMeas = 5000;

    %This will contain all the measurements. The four channels will be on top
    %of each other.
    
    data = zeros(measPerFile, nbrOfMeas, channels);
    disp(['Loading ' num2str(nbrOfMeas*channels) ' files...'])
    
    modCheck = max(floor(nbrOfMeas/100), 1);
    fprintf(1, '  0%% done')
    for i = 1:nbrOfMeas
        measFile = C1files(i);
        fileName = measFile.name;
        for j = 1:channels
            fp = [path 'C' int2str(j) fileName(3:end)];
            importedData = importdata(fp);
            data(:, i, j) = importedData;
        end
        if mod(i, modCheck) == 0
            percentProgress = ceil(i/nbrOfMeas*100);
            fprintf(1, '\b\b\b\b\b\b\b\b\b\b%3d%% done', percentProgress)
        end
    end
    fprintf(1, '\n')
    if saveData
        disp('Saving data...')
        save('mcpData')
    end
end

%% Settings
plotOffsets = 1;
plotFourierTransform = 1;
plotSignals = 1;
plotMeanPulse = 1;
plotFittedPeaks = 1;

chosenChannel = 1;
chosenSignal = 1;

%% Post Loading

channelPairs = [1 2 3 4]; %1 2 3 4 is the correct configuration
channelGroups = [channelPairs(1:2); channelPairs(3:4)];

colors = ['y', 'r', 'b', 'g'];

inputImpedance = 50; %Impedance of the oscilloscope in Ohms

T = dummyData(:, 1); %Time vector
t = T(2) - T(1); %Sampling time
Fs = 1/t; %Sampling frequency
L = measPerFile; %Length of signal
f = Fs/2*linspace(0,1,L/2+1); %Frequency vector

riseTime = 1e-8;
nRiseTime = floor(riseTime/t);

%% Remove offsets

disp('Removing offsets...')
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
    suptitle('Before and after removing the offset')
end

%% Clean signals from noise with low pass filters

disp('Cleaning with low pass filters...')

if plotFourierTransform
    meas = data(:, chosenSignal, chosenChannel);
    fourierPlot = figure(11);
    clf(11)
    set(gcf, 'Name', 'Signal Fourier Transform')
    hold on
    subplot(2, 2, 1)
    hold on
    title('Before low pass')
    xlabel('Time [s]')
    ylabel('Voltage [V]')
    plot(T, meas)
    MEAS = fft(meas)/L;
    subplot(2, 2, 2)
    semilogy(f, 2*abs(MEAS(1:L/2+1))) 
    hold on
    title('Uncut Fourier spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Fourier transform [Vs]')
end
sc = t/1e-9;
for n = 1:4
    data = filter(sc, [1 sc-1], data);
end
if plotFourierTransform
    meas = data(:, chosenSignal, chosenChannel);
    MEAS = fft(meas)/L;
    subplot(2, 2, 3)
    hold on
    title('After low pass')
    xlabel('Time [s]')
    ylabel('Voltage [V]')
    plot(T, meas);
    subplot(2, 2, 4)
    semilogy(f, 2*abs(MEAS(1:L/2+1)))
    hold on
    title('Cut Fourier spectrum')
    xlabel('Frequency (Hz)')
    ylabel('Fourier transform [Vs]')
    suptitle('Before and after filter')
end

%% Remove bad signals by shape

disp('Removing bad signals by shape...')

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

%% Locate peaks

disp('Locating peaks by minimum of fitted quadratic...')
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
        if plotFittedPeaks && i == chosenSignal && j == chosenChannel
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
    end
end

%% Remove bad signals from time sum

disp('Removing bad signals from time sum...')
timeSum = [sum(signals(:, channelGroups(1, :)), 2) sum(signals(:, channelGroups(2, :)), 2)];
tMean = mean(timeSum);
tStd = std(timeSum);

for k = 1:channels/2
    [row col] = find(abs([timeSum(:, k) - tMean(k)]) > 3*tStd(k));
    good(row) = 0;
end

nbrOfGoods = length(find(good == 1));
disp(['Found ' num2str(nbrOfMeas - nbrOfGoods) ' bad signals from time sum. Removing...'])
data = data(:, find(good == 1), :);
nbrOfMeas = size(data, 2);
signals = signals(find(good == 1), :);
signalIndices = signalIndices(find(good == 1), :);

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
            meas = data(:, i, channelPairs(j));
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
grandTotalCharge = sum(totalCharge, 2);

%% Calculate FDHM of pulses. Good idea to check the code below...

disp('Calculating FDHM for the pulses...')

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

%% Calculate pulse shape by averaging. This needs some more work, like cutting from the histogram fdhmHistPlot

disp('Calculating mean pulse shape...')
[foo mins] = min(data);
mins = squeeze(mins);
pulseShaper = zeros(4*nRiseTime + 1, size(data, 2), size(data, 3));
if plotMeanPulse
    meanPulsePlot = figure(40);
    clf(40)
    set(gcf, 'Name', 'Mean pulse calculation')
    subplot(1, 2, 1)
    xlabel('Shifted time indices')
    ylabel('Voltage [V]')
    hold on
    title(['Pulses from channel ' num2str(chosenChannel) ' overlaid'])
end
for i = 1:nbrOfMeas
    for j = 1:channels
        nRange = (mins(i, j) - 2*nRiseTime):(mins(i, j) + 2*nRiseTime);
        pulseShaper(:, i, j) = data(nRange, i, j);
        if plotMeanPulse && i < 150 && j == chosenChannel
            plot(pulseShaper(:, i, j), colors(mod(i, 4) + 1))
        end
    end
end

pulse = squeeze(mean(pulseShaper, 2));
if plotMeanPulse
    subplot(1, 2, 2)
    hold on
    title(['Mean pulse for channel ' num2str(chosenChannel)])
    xlabel('Shifted time indices')
    ylabel('Voltage [V]')
    plot(pulse(:, 1))
    suptitle('Calculation of the average pulse')
end

%% Calculate total time and fit double Gaussian

disp('Calculating sums of times and fitting double Gaussian...')

timeMinSum = [sum(signalIndices(:, channelGroups(1, :)), 2) sum(signalIndices(:, channelGroups(2, :)), 2)];
timeSum = [sum(signals(:, channelGroups(1, :)), 2) sum(signals(:, channelGroups(2, :)), 2)];

disp('Plotting results in histogram and x-y plot...')

for k = 1:2
    tMean = mean(timeSum(:, k));
    tStd = std(timeSum(:, k));
    interval = tMean - 3*tStd : t/2 : tMean + 3*tStd;
    [N, x] = hist(timeSum(:, k), interval);
    x = x(2:end-1)';
    N = N(2:end-1)';
    timeSumX{k} = x;
    timeSumN{k} = N;

    f = fittype('gauss2');
    options = fitoptions('gauss2');
    options.Lower = [0 -Inf 0 0 -Inf 0];
    [fittedGaussians gof output] = fit(x, N, f, options);%, 'Weight', sqrt(N));
    gaussianFits{k} = fittedGaussians;
end

%% Calculate positions

disp('Calculating spatial coordinates...')

timeMinDiff = [diff(signalIndices(:, channelGroups(1, :)), 2) diff(signalIndices(:, channelGroups(2, :)), 2)];
timeDiff = [diff(signals(:, channelGroups(1, :)), 1, 2) diff(signals(:, channelGroups(2, :)), 1, 2)];

%This minus sign is arbitrary, only mirrors the image in the origin.
timeDiff = -timeDiff;

%% End timing for the data processing
toc

%% Run the data analysis script to produce figures
analysis

%%
toc
