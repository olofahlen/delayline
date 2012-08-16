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
nbrOfMeas = 1000;

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
freqCut = 1e9;

T = dummyData(:, 1); %Time vector
t = T(2) - T(1); %Sampling time
Fs = 1/t; %Sampling frequency
L = measPerFile; %Length of signal
f = Fs/2*linspace(0,1,L/2+1); %Frequency vector
fCutLength = length(f) - find(f > freqCut, 1, 'first');
fZeroMask = zeros(2*fCutLength, 1);

%% Remove offsets

disp('Removing offsets...')
for j = 1:channels
    offset = mean2(data((1:measPerFile) + (measPerFile * (j - 1)), 1:end/10));
    data((1:measPerFile) + (measPerFile * (j - 1)), :) = data((1:measPerFile) + (measPerFile * (j - 1)), :) - offset;
end

%% Filter out bad measurements -- UNDER CONSTRUCTION!

bad = zeros(nbrOfMeas, 1);

meas = data(1:measPerFile, 46)
%meas = meas - mean(meas)
nicestd = std(meas)
limit = nicestd;
figure(21)
clf(21)
plot(T, meas)
hold on
%line([T(1); T(2)], [nicestd; nicestd])
line([T(1) T(end)], [limit limit])
line([T(1) T(end)], [-limit -limit])



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
        if j == 2 && i == 418
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
            plot(f, 2*abs(Y(1:L/2+1))) 
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

%% Plot signals

figure(1)
clf(1)
hold on
colors = ['r', 'g', 'b', 'y'];
for j = 1:channels
    color = colors(j);
    for i = 10:10
        meas = data((1:measPerFile) + (measPerFile * (channelPairs(j) - 1)), i);
        subplot(2, 1, ceil(j/2));
        hold on
        %plot(T, data((1:measPerFile) + (measPerFile * (channelPairs(j) - 1)), i), color)
        plot(T, meas, color)
        plot(T(signalIndices(channelPairs(j), i)), meas(signalIndices(channelPairs(j), i)), 'o')
    end
end

%% Calculate positions
channelPairs = [1 2 3 4];
channelGroups = [channelPairs(1:2); channelPairs(3:4)];
timeDiff = zeros(2, nbrOfMeas);

for i = 1:nbrOfMeas
    for k = 1:2
        timeDiff(k, i) = T(signalIndices(channelGroups(k, 1), i)) - T(signalIndices(channelGroups(k, 2), i));
    end
end

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