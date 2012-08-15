clear all
clc
close all

% This file loads waveforms (Time + Amplitude) saved by the WavePro7100
% The variable 'path' below should be the absolute path to a folder
% containing measurements in the form 'C?*.dat', for example
% 'C1mcp00012.dat'.

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
T = dummyData(:, 1);

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

%% Remove offsets

for j = 1:channels
    offset = mean2(data((1:measPerFile) + (measPerFile * (j - 1)), 1:end/10));
    data((1:measPerFile) + (measPerFile * (j - 1)), :) = data((1:measPerFile) + (measPerFile * (j - 1)), :) - offset;
end

%% Locate peaks

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
        plot(T(signalIndices(channelPairs(j), i)), meas(signalIndices(channelPairs(j))), 'o')
    end
end

%% Calculate positions
channelPairs = [1 2 3 4];
channelGroups = [channelPairs(1:2); channelPairs(3:4)];
timeDiff = zeros(2, nbrOfMeas);

for i = 1:nbrOfMeas
    for k = 1:2
        timeDiff(k, i) = T(signalIndices(channelGroups(k, 1), i)) + T(signalIndices(channelGroups(k, 2), i));
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
plot(timeDiff(1, :), timeDiff(2, :), '*')

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