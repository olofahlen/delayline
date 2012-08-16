%% Initial

clear all
clc

% This file loads waveforms (Time + Amplitude) saved by the WavePro7100
% The variable 'path' below should be the absolute path to a folder
% containing measurements in the form 'C?*.dat', for example
% 'C1mcp00012.dat'.

%% Settings
plotFourierTransform = 1;

%% Load data

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
%of each other
data = zeros(measPerFile*4, nbrOfMeas);

i = 1;
if nbrOfMeas > 100
    modCheck = floor(nbrOfMeas/100);
else
    modCheck = 1;
end

disp('Loading files...')
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

%% Post Loading

channelPairs = [1 2 3 4]; %Real
freqCut = 0.2e9;

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
figure(22)
clf(22)

good = ones(nbrOfMeas, 1);
for i = 1:nbrOfMeas
    for j = 1:channels
        meas = data((1:measPerFile) + (measPerFile * (j - 1)), i);
        [minVal minIndex] = min(meas);
        potentialStart = minIndex - nRiseTime;
        if potentialStart < measPerFile/15
            good(i) = 0;
        else
            nicestd = std(meas(1:potentialStart));
            nicemean = mean(meas(1:potentialStart));
            upperlimit = 4*nicestd + nicemean;
            lowerlimit = -4*nicestd + nicemean;
            if length(find(meas < lowerlimit)) < nRiseTime/2
                good(i) = 0;
            end
        end
        if j == 1 && i < 21
            subplot(2, 1, 1)
            plot(T, meas)
            hold on
        end
        %subplot(2, 1, 1)
        %plot(T, meas)
        %line([T(1); T(2)], [nicestd; nicestd])
        %line([T(1) T(end)], [upperlimit upperlimit])
        %line([T(1) T(end)], [lowerlimit lowerlimit])
        %line([T(1) T(end)], [nicemean nicemean], 'Color', 'g')
        meanCut = find(meas < lowerlimit, 1, 'first');
        data((1:measPerFile) + (measPerFile * (j - 1)), i) = meas - mean(meas(1:meanCut));
        %subplot(2, 1, 2)
        %plot(T, meas)
        %line([T(1) T(end)], [0 0])
        if j == 1 && i == 1
            subplot(2, 1, 2)
            plot(T, meas)
            hold on
        end
       
%         if good(i)
%             xlabel('Good signal')
%         else
%             xlabel('Bad signal')
%         end

        %pause
    end
end

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

nbrOfMeas = nbrOfGoods;
data = goodData;
    


%% Old method to remove offsets

% disp('Removing offsets...')
% for j = 1:channels
%     offset = mean2(data((1:measPerFile) + (measPerFile * (j - 1)), 1:end/10));
%     data((1:measPerFile) + (measPerFile * (j - 1)), :) = data((1:measPerFile) + (measPerFile * (j - 1)), :) - offset;
% end

%% Clean signals from noise using the Fourier Transform

disp('Cleaning with Fourier Transform...')

if plotFourierTransform
    figure(11)
    clf(11)
    hold on
end

totalIter = channels*nbrOfMeas;
if channels*nbrOfMeas > 100
    modCheck = floor(totalIter/100);
else
    modCheck = 1;
end
loopCounter = 1;
for i = 1:nbrOfMeas
    for j = 1:channels
        if j == 1 && i == 46
            plotFourierTransform = 1;
        else
            plotFourierTransform = 0;
        end
        meas = data((1:measPerFile) + (measPerFile * (channelPairs(j) - 1)), i);
        MEAS = fft(meas)/L;
        if plotFourierTransform
            subplot(2, 2, 1)

            hold off
            plot(T, meas)

            hold on
            plot(T(signalIndices(channelPairs(j), i)), meas(signalIndices(channelPairs(j), i)), 'o')
            subplot(2, 2, 2)

            plot(f, 2*abs(MEAS(1:L/2+1))) 
            xlabel('Frequency (Hz)')

            subplot(2, 2, 4)
            hold off
        end
        
        MEAS(L/2 - fCutLength + 1:L/2 + fCutLength) = fZeroMask;
        cleanedMeas = real(ifft(MEAS));
        
        if plotFourierTransform
            plot(f, 2*abs(MEAS(1:L/2+1))) 
            subplot(2, 2, 3)
            plot(T, cleanedMeas);
        end
        
        data((1:measPerFile) + (measPerFile * (channelPairs(j) - 1)), i) = cleanedMeas;
        if mod(loopCounter, modCheck) == 0
            percentProgress = ceil(loopCounter/totalIter*100);
            disp([num2str(percentProgress) '% done'])
        end
        loopCounter = loopCounter + 1;
    end
end

%% Locate peaks

disp('Locating peaks...')
%signals = zeros(4, nbrOfMeas);
signalIndices = zeros(4, nbrOfMeas);

for i = 1:nbrOfMeas
    for j = 1:channels
        meas = data((1:measPerFile) + (measPerFile * (j - 1)), i);
        [minValue, minIndex] = min(meas);
        %signals(j, i) = T(minIndex);
        signalIndices(j, i) = minIndex;
    end
end

%% Plot signals

disp('Plotting signals...')
figure(1)
clf(1)
hold on
colors = ['r', 'g', 'b', 'y'];
for j = 1:channels
    color = colors(j);
    for i = 1:1
        meas = data((1:measPerFile) + (measPerFile * (channelPairs(j) - 1)), i);
        subplot(2, 1, ceil(j/2));
        hold on
        %plot(T, data((1:measPerFile) + (measPerFile * (channelPairs(j) - 1)), i), color)
        plot(T, meas, color)
        plot(T(signalIndices(channelPairs(j), i)), meas(signalIndices(channelPairs(j), i)), 'o')
    end
end

%% Calculate positions

disp('Calculating spatial coordinates...')
channelPairs = [1 2 3 4];
channelGroups = [channelPairs(1:2); channelPairs(3:4)];
timeDiff = zeros(2, nbrOfMeas);

for i = 1:nbrOfMeas
    for k = 1:2
        timeDiff(k, i) = T(signalIndices(channelGroups(k, 1), i)) - T(signalIndices(channelGroups(k, 2), i));
    end
end

disp('Plotting results in histogram and x-y plot...')
bins = 500;

figure(2)
clf(2)
subplot(2, 1, 1)
hist(timeDiff(1, :), bins)
subplot(2, 1, 2)
hist(timeDiff(2, :), bins)

figure(4)
clf(4)
plot(timeDiff(1, :), timeDiff(2, :), '.')
axis square

%% Plot times
%figure
%plot(timeDiff(1, :))

%% Plot Histograms

%figure(3)
%clf(3)
%subplot(2, 1, 1)
%hist(timeDiff(1, :))
%subplot(2, 1, 2)
%hist(timeDiff(2, :))