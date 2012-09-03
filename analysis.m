plotCharges = true;
plotFdhm = true;
plotSignals = true;
plotTimeSums = true;
plotPositions = true;
plotTimeCutHitmap = true;

%% Plot histograms with charges
if plotCharges
    disp('Plotting charges...')
    bins = 100;
    interval = linspace(0, 14e6, bins);
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
        subplot(3, 1, k)
        hold on
        title(['Total charge deposited on channels ' num2str(channelGroups(k, 1)) ' and ' num2str(channelGroups(k, 2))])
        xlabel('Charge [e]')
        ylabel('Counts')
        hist(totalCharge(:, k), interval)
        %hist(totalCharge(:, k), bins)
    end
    subplot(3, 1, 3)
    hist(grandTotalCharge, bins)
    title(['Total charge deposited both delay lines'])
    xlabel('Charge [e]')
    ylabel('Counts')
    suptitle('Histograms of charges for the different channels')
    chargeScatterPlot = figure(32);
    clf(32)
    hold on
    title('Correlation of charge deposited on the two delay lines')
    set(gcf, 'Name', 'Charge scatter plot')
    scatter(totalCharge(:, 1), totalCharge(:, 2))
    axis square
    xlabel(['Total charge deposited on channels ' num2str(channelGroups(1, 1)) ' and ' num2str(channelGroups(1, 2)) ' [e]'])
    ylabel(['Total charge deposited on channels ' num2str(channelGroups(2, 1)) ' and ' num2str(channelGroups(2, 2)) ' [e]'])
end

%% Plot histograms for FDHM
if plotFdhm
    disp('Plptting FDHM histograms...')
    fdhmHistPlot = figure(35);
    clf(fdhmHistPlot)
    set(gcf, 'Name', 'FDHM of the four channels')
    hold on

    fdhmMeans = mean(fdhm);
    fdhmStds = std(fdhm);
    fdhmLimits = [fdhmMeans - 3*fdhmStds; fdhmMeans + 3*fdhmStds];
    for j = 1:channels
        binX = fdhmLimits(1, j):t/2:fdhmLimits(2, j);
        subplot(4, 1, j)
        hold on
        title(['FDHM for channel ' num2str(channelPairs(j))])
        xlabel('FDHM [s]')
        ylabel('Counts')
        hist(fdhm(:, j), binX)
        axis([min(fdhmLimits(1, :)) max(fdhmLimits(2, :)) 0 1])
        axis 'auto y'
    end
    suptitle(['Distribution of FDHM'])
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

%%  Plot histograms for time sums and the fitted double Gaussian
if plotTimeSums
    disp('Plotting time sum histograms with double Guassian fit...')
    timeSumHistPlot = figure(200);
    clf(200)
    set(gcf, 'Name', 'Normalized histograms of time sums')
    for k = 1:2
        subplot(2, 1, k)
        hold on
        xlabel('$t_1 + t_2$ [s]', 'Interpreter', 'LaTeX')
        ylabel('Normalized counts')
        bar(timeSumX{k}, timeSumN{k})

        fittedGaussians = gaussianFits{k};
        fitPlot = plot(fittedGaussians, 'r');
        h = legend(fitPlot, '$a\cdot N(\mu_1, \sigma_1) + (1-a)\cdot N(\mu_2, \sigma_2)$');
        set(h, 'Interpreter', 'LaTeX')
        title(['Delayline for channels ' num2str(channelGroups(k, 1)) ' and ' num2str(channelGroups(k, 2)) sprintf('. a = %.4f, mu_1 = %.4f ns, mu_2 = %.4f ns, sigma_1 = %.4f ns, sigma_2 = %f ns', fittedGaussians.a1, 1e9*fittedGaussians.b1, 1e9*fittedGaussians.b2, 1e9*fittedGaussians.c1, 1e9*fittedGaussians.c2)])
    end
    suptitle('Normalized histograms of time sums for the two delay lines')
end

%% Plot histograms for time differences and MCP hitmap
if plotPositions
    disp('Plotting time difference histograms and MCP hitmap...')
    bins = 500;
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

if plotTimeCutHitmap
    disp('Plotting time cut hitmap...')
    timeCutHitmapPlot = figure(301);
    clf(301)
    hold on
    set(gcf, 'Name', 'MCP 2D-plot with time cuts')

    cuts = [9.68e-8 1.011e-7];

    for k = 1:2
        less = find(timeSum(:, k) < cuts(k));
        more = find(timeSum(:, k) > cuts(k));

        subplot(1, 2, k)
        hold on
        title(['Cut at ' num2str(cuts(k)) 's for channels ' num2str(channelGroups(k, 1)) ' and ' num2str(channelGroups(k, 2))])
        scatter(timeDiff(less, 1), timeDiff(less, 2), 'b')
        scatter(timeDiff(more, 1), timeDiff(more, 2), 'r')
        xlabel('$x\propto \Delta t_x$', 'Interpreter', 'LaTeX');
        ylabel('$y\propto \Delta t_y$', 'Interpreter', 'LaTeX');
        legend('Short times', 'Long times')
        axis square
    end

    suptitle('Events cut in time histograms')
end
