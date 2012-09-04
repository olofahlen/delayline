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

%figures = {'fourierPlot', 'offsetPlot', 'individualChargePlot', 'totalChargePlot', 'chargeScatterPlot', 'fdhmHistPlot', 'meanPulsePlot', 'fittedPeakPlot', 'signalPlot', 'timeSumHistPlot', 'timeDiffHistPlot', 'mcpHitmapPlot', 'timeCutHitmapPlot'};
figureNames = fieldnames(figures);
formats = {'png', 'epsc', 'fig', 'pdf'};
path = '/home/thorleif/mcp/tests/pics/';

for figureName = figureNames
    for l = 1:length(formats)
        format = formats{l};
        disp(['Saving ' figureName ' as ' format '...'])
        saveas(eval(figureName), [path figureName], format)
    end
end
