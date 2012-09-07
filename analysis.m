%%Initial
disp('Starting data analysis...')

tic
importProcessedData = false;

if importProcessedData
    disp('Importing saved processed data...')
    load('processed96705')
end

plotCharges = true;
plotFdhm = true;
plotTimeSums = true;
plotPositions = true;
plotCutHitmap = true;
plotSkewness = true;

markerSize = 4;
bins = 200;
nbrOfStd = 4;

allMask = [1:nbrOfMeas]';
smallLeftSquareMask = find(-6e-8 < timeDiff(:, 1) & timeDiff(:, 1) < -4e-8 & 0 < timeDiff(:, 2) & timeDiff(:, 2) < 2e-8);
secondQuadrantMask = find(timeDiff(:, 1) < 0 & timeDiff(:, 2) > 0);
shortFdhmMask = find(fdhm(:, 1) < 7.76e-9);
skewBlobMask = find(-0.65 < signalSkewness(:, 1) & signalSkewness(:, 1) < -0.45 & -0.65 < signalSkewness(:, 2) & signalSkewness(:, 2) < -0.5);

%skewLineMask = zeros(1, 1);
upperCut = @(sn) (-0.446 - -0.2)/(-0.232 - -0.5)*(sn - -0.5) + -0.2;
lowerCut = @(sn) (-0.529 - -0.3013)/(-0.395 - -0.581)*(sn - -0.581) + -0.3013;
l = 1;
for i = 1:nbrOfMeas
    snx = signalSkewness(i, 1);
    sny = signalSkewness(i, 2);
    if sny < upperCut(snx) && lowerCut(snx) < sny && -0.55 < snx && snx < -0.25
        skewLineMask(l) = i;
        l = l + 1;
    end
end

voltageMean = mean(signalVoltages(: ,1));
voltageStd = std(signalVoltages(: ,1));
lowPulseHeightMask = find(signalVoltages(:, 1) < voltageMean-voltageStd);

mask = allMask;

maskedCharge = charge(mask, :);
maskedSignalSkewness = signalSkewness(mask, :);
maskedTotalCharge = totalCharge(mask, :);
maskedGrandTotalCharge = grandTotalCharge(mask, :);
maskedFdhm = fdhm(mask, :);
maskedTimeSum = timeSum(mask, :);
maskedTimeDiff = timeDiff(mask, :);

compare = false;
compareMask = shortFdhmMask;
intersectMask = intersect(mask, compareMask);

compareCharge = charge(compareMask, :);
compareSignalSkewness = signalSkewness(compareMask, :);
compareTotalCharge = totalCharge(compareMask, :);
compareGrandTotalCharge = grandTotalCharge(compareMask, :);
compareFdhm = fdhm(compareMask, :);
compareTimeSum = timeSum(compareMask, :);
compareTimeDiff = timeDiff(compareMask, :);

intersectCharge = charge(intersectMask, :);
intersectSignalSkewness = signalSkewness(intersectMask, :);
intersectTotalCharge = totalCharge(intersectMask, :);
intersectGrandTotalCharge = grandTotalCharge(intersectMask, :);
intersectFdhm = fdhm(intersectMask, :);
intersectTimeSum = timeSum(intersectMask, :);
intersectTimeDiff = timeDiff(intersectMask, :);


if compare
    disp('Compare mode is on')
else
    disp('Compare mode is off')
end

%% Create scatter plot for skewness

if plotSkewness
    disp('Creating skewness scatter plot...')
    figures.skewnessPlot = figure(210);
    clf(figures.skewnessPlot)
    hold on
    set(gcf, 'Name', 'Signal Skewness')
    for k = 1:2
        subplot(1, 2, k)
        hold on
        xlabel(num2str(channelGroups(k, 1), 'Signal skewness for channel %d'))
        ylabel(num2str(channelGroups(k, 2), 'Signal skewness for channel %d'))
        title(num2str(k, 'Correlations for delay line %d'))
        scatter(maskedSignalSkewness(:, channelGroups(k, 1)), maskedSignalSkewness(:, channelGroups(k, 2)), markerSize)
        if compare
            scatter(compareSignalSkewness(:, channelGroups(k, 1)), compareSignalSkewness(:, channelGroups(k, 2)), markerSize, 'r')
            scatter(intersectSignalSkewness(:, channelGroups(k, 1)), intersectSignalSkewness(:, channelGroups(k, 2)), markerSize, 'm')
        end
    end
    suptitle('Correlation for signal skewness')
end


%% Plot histograms with charges
if plotCharges
    disp('Plotting charges...')
    interval = linspace(0, 14e6, bins);
    figures.individualChargePlot = figure(21);
    clf(figures.individualChargePlot)
    set(gcf, 'Name', 'Individual charge histograms')
    intMean = mean(charge);
    intStd = std(charge);
    limits = [intMean - nbrOfStd*intStd; intMean + nbrOfStd*intStd];
    for j = 1:channels
        interval = linspace(limits(1, j), limits(2, j), bins);
        subplot(4, 1, j)
        hold on
        title(['Charge deposited on channel ' num2str(channelPairs(j))])
        xlabel('Charge [e]')
        ylabel('Counts')
        hist(maskedCharge(:, j), interval)
        if compare
            h = findobj(gca, 'Type', 'patch');
            set(h, 'FaceColor', 'b', 'EdgeColor', 'w', 'facealpha', 0.75)
            hist(compareCharge(:, j), interval)
            newH = findobj(gca, 'Type', 'patch');
            newH = newH(find(newH ~= h));
            set(newH, 'FaceColor', 'r', 'EdgeColor', 'w', 'facealpha', 0.75)
        end
    end
    suptitle('Charge arrived at each channel')

    figures.totalChargePlot = figure(22);
    clf(figures.totalChargePlot)
    set(gcf, 'Name', 'Total charge histograms')
    intMean = mean(totalCharge);
    intStd = std(totalCharge);
    limits = [intMean - nbrOfStd*intStd; intMean + nbrOfStd*intStd];
    for k = 1:2
        interval = linspace(limits(1, k), limits(2, k), bins);
        subplot(3, 1, k)
        hold on
        title(['Total charge deposited on channels ' num2str(channelGroups(k, 1)) ' and ' num2str(channelGroups(k, 2))])
        xlabel('Charge [e]')
        ylabel('Counts')
        hist(maskedTotalCharge(:, k), interval)
        if compare
            h = findobj(gca, 'Type', 'patch');
            set(h, 'FaceColor', 'b', 'EdgeColor', 'w', 'facealpha', 0.75)
            hist(compareTotalCharge(:, k), interval)
            newH = findobj(gca, 'Type', 'patch');
            newH = newH(find(newH ~= h));
            set(newH, 'FaceColor', 'r', 'EdgeColor', 'w', 'facealpha', 0.75)
        end
    end
    subplot(3, 1, 3)
    intMean = mean(grandTotalCharge);
    intStd = std(grandTotalCharge);
    limits = [intMean - nbrOfStd*intStd; intMean + nbrOfStd*intStd];
    interval = linspace(limits(1), limits(2), bins);
    hist(maskedGrandTotalCharge, interval)
    hold on
    if compare
        h = findobj(gca, 'Type', 'patch');
        set(h, 'FaceColor', 'b', 'EdgeColor', 'w', 'facealpha', 0.75)
        hist(compareGrandTotalCharge, interval)
        newH = findobj(gca, 'Type', 'patch');
        newH = newH(find(newH ~= h));
        set(newH, 'FaceColor', 'r', 'EdgeColor', 'w', 'facealpha', 0.75)
    end
    title(['Total charge deposited both delay lines'])
    xlabel('Charge [e]')
    ylabel('Counts')
    suptitle('Charge deposited on each delay line and total charge deposited')

    figures.chargeScatterPlot = figure(23);
    clf(figures.chargeScatterPlot)
    hold on
    title('Correlation of charge deposited on the two delay lines')
    set(gcf, 'Name', 'Charge scatter plot')
    scatter(maskedTotalCharge(:, 1), maskedTotalCharge(:, 2), markerSize)
    if compare
        scatter(compareTotalCharge(:, 1), compareTotalCharge(:, 2), markerSize, 'r')
        scatter(intersectTotalCharge(:, 1), intersectTotalCharge(:, 2), markerSize, 'm')
    end
    axis square
    maxCharge = max(max(totalCharge));
    axis([0 maxCharge 0 maxCharge])
    xlabel(['Total charge deposited on channels ' num2str(channelGroups(1, 1)) ' and ' num2str(channelGroups(1, 2)) ' [e]'])
    ylabel(['Total charge deposited on channels ' num2str(channelGroups(2, 1)) ' and ' num2str(channelGroups(2, 2)) ' [e]'])
end

%% Plot histograms for FDHM
if plotFdhm
    disp('Plotting FDHM histograms...')
    figures.fdhmHistPlot = figure(24);
    clf(figures.fdhmHistPlot)
    set(gcf, 'Name', 'FDHM of the four channels')
    hold on

    intMean = mean(fdhm);
    intStd = std(fdhm);
    limits = [min(fdhm); intMean + nbrOfStd/4*intStd];
    %limits(1, :) = max(limits(1, :), zeros(1, size(limits, 2)));
    interval = linspace(min(limits(1, :)), max(limits(2, :)), bins);
    for j = 1:channels
        %interval = linspace(limits(1, j), limits(2, j), bins);
        subplot(4, 1, j)
        hold on
        title(['FDHM for channel ' num2str(channelPairs(j))])
        xlabel('FDHM [ns]')
        ylabel('Counts')
        hist(1e9*maskedFdhm(:, j), 1e9*interval)
        if compare
            h = findobj(gca, 'Type', 'patch');
            set(h, 'FaceColor', 'b', 'EdgeColor', 'w', 'facealpha', 0.75)
            hist(1e9*compareFdhm(:, j), 1e9*interval)
            newH = findobj(gca, 'Type', 'patch');
            newH = newH(find(newH ~= h));
            set(newH, 'FaceColor', 'r', 'EdgeColor', 'w', 'facealpha', 0.75)
        end
        xlim(1e9*[min(limits(1, :)) 1.03*max(limits(2, :))])
        axis 'auto y'
    end
    suptitle('Distribution of FDHM')
end

%%  Plot histograms for time sums and fit a double Gaussian
disp('Fitting double Gaussian...')
intMean = mean(timeSum);
intStd = std(timeSum);
limits = [intMean - nbrOfStd*intStd; intMean + nbrOfStd*intStd];
for l = 1:1 + compare
    if l == 1
        calcTimeSum = maskedTimeSum;
    else
        calcTimeSum = compareTimeSum;
    end
    for k = 1:2
        interval = linspace(limits(1, k), limits(2, k), bins);
        [N, x] = hist(calcTimeSum(:, k), interval);
        timeSumX{k, l} = x;
        timeSumN{k, l} = N;

        %This part cuts away the tails to better fit the Gaussian
        x = x(2:end-1)';
        N = N(2:end-1)';

        fitObj = fittype('gauss2');
        options = fitoptions('gauss2');
        options.Lower = [0 -Inf 0 0 -Inf 0];
        [fittedGaussians gof output] = fit(x, N, fitObj, options);
        gaussianFits{k, l} = fittedGaussians;
    end
end

if plotTimeSums
    disp('Plotting time sum histograms with double Guassian fit...')
    figures.timeSumHistPlot = figure(25);
    clf(figures.timeSumHistPlot)
    set(gcf, 'Name', 'Histograms of time sums')
    for k = 1:2
        subplot(2, 1, k)
        hold on
        xlabel('$t_1 + t_2$ [ns]', 'Interpreter', 'LaTeX')
        ylabel('Counts')
        bar(1e9*timeSumX{k, 1}, timeSumN{k, 1})
        if compare
            h = findobj(gca, 'Type', 'patch');
            set(h, 'FaceColor', 'b', 'EdgeColor', 'w', 'facealpha', 0.75)
            bar(1e9*timeSumX{k, 2}, timeSumN{k, 2})
            newH = findobj(gca, 'Type', 'patch');
            newH = newH(find(newH ~= h));
            set(newH, 'FaceColor', 'r', 'EdgeColor', 'w', 'facealpha', 0.75)
        end

        fittedGaussians = gaussianFits{k, 1};
        fitPlot = plot(1e9*timeSumX{k, 1}, fittedGaussians(timeSumX{k, 1}), 'k');
        title(['Delayline for channels ' num2str(channelGroups(k, 1)) ' and ' num2str(channelGroups(k, 2)) sprintf('. a_1 = %.2f, a_2 = %.2f, mu_1 = %.2f ns, mu_2 = %.2f ns, sigma_1 = %.2f ns, sigma_2 = %.2f ns', fittedGaussians.a1, fittedGaussians.a2, 1e9*fittedGaussians.b1, 1e9*fittedGaussians.b2, 1e9*fittedGaussians.c1, 1e9*fittedGaussians.c2)])
        if compare
            fittedGaussians = gaussianFits{k, 2};
            fitPlot = plot(1e9*timeSumX{k, 2}, fittedGaussians(timeSumX{k, 2}), 'k');
            oldTitle = get(get(gca, 'Title'), 'string');
            title([oldTitle sprintf('\na_1 = %.2f, a_2 = %.2f, mu_1 = %.2f ns, mu_2 = %.2f ns, sigma_1 = %.2f ns, sigma_2 = %.2f ns', fittedGaussians.a1, fittedGaussians.a2, 1e9*fittedGaussians.b1, 1e9*fittedGaussians.b2, 1e9*fittedGaussians.c1, 1e9*fittedGaussians.c2)])
        end
    end
    suptitle('Time sums for the delay lines')
    figures.timeSumScatterPlot = figure(29);
    clf(figures.timeSumScatterPlot)
    set(gcf, 'Name', 'Scatter plot of time sums')
    title('Time sums for the two delay lines')
    hold on
    xlabel('Time sum [ns]')
    ylabel('Time sum [ns]')
    scatter(1e9*maskedTimeSum(:, 1), 1e9*maskedTimeSum(:, 2), markerSize)
    if compare
        scatter(1e9*compareTimeSum(:, 1), 1e9*compareTimeSum(:, 2), markerSize, 'r')
        scatter(1e9*intersectTimeSum(:, 1), 1e9*intersectTimeSum(:, 2), markerSize, 'm')
    end
    axis square
end

%% Plot histograms for time differences and MCP hitmap
if plotPositions
    disp('Plotting time difference histograms and MCP hitmap...')
    figures.timeDiffHistPlot = figure(26);
    clf(figures.timeDiffHistPlot)
    set(gcf, 'Name', 'Histograms of time differences')
    intMean = mean(timeDiff);
    intStd = std(timeDiff);
    limits = [intMean - nbrOfStd*intStd; intMean + nbrOfStd*intStd];
    for k = 1:2
        interval = linspace(limits(1, k), limits(2, k), bins);
        subplot(2, 1, k)
        hold on
        title(['Time difference for channels ' num2str(channelGroups(1, k)) ' and ' num2str(channelGroups(2, k))])
        xlabel('$\Delta t$ [ns]', 'Interpreter', 'LaTeX')
        ylabel('Counts')
        hist(1e9*maskedTimeDiff(:, k), 1e9*interval)
        if compare
            h = findobj(gca, 'Type', 'patch');
            set(h, 'FaceColor', 'b', 'EdgeColor', 'w', 'facealpha', 0.75)
            hist(1e9*compareTimeDiff(:, k), 1e9*interval)
            newH = findobj(gca, 'Type', 'patch');
            newH = newH(find(newH ~= h));
            set(newH, 'FaceColor', 'r', 'EdgeColor', 'w', 'facealpha', 0.75)
        end
    end
    suptitle('Histograms of time differences for the two delay lines')
    figures.mcpHitmapPlot = figure(27);
    clf(figures.mcpHitmapPlot)
    set(gcf, 'Name', 'MCP 2D-plot')
    hold on
    for l = 1:1+compare
        if compare
            subplot(1, 2, l)
            hold on
        end
        if l == 1
            title('MCP hitmap')
            calcTimeDiff = maskedTimeDiff;
            calcTotalCharge = maskedTotalCharge;
        else
            title('MCP hitmap Compare')
            calcTimeDiff = compareTimeDiff;
            calcTotalCharge = compareTotalCharge;
        end
        xlabel('$x\propto \Delta t_x$', 'Interpreter', 'LaTeX');
        ylabel('$y\propto \Delta t_y$', 'Interpreter', 'LaTeX');
        scatter(calcTimeDiff(:, 1), calcTimeDiff(:, 2), markerSize, sum(calcTotalCharge, 2 ))
        axis square
        cb = colorbar;
        ylabel(cb, 'Charge [e]')
    end
    if compare
        suptitle('Reconstruction of particle hits on the MCP')
    end
end

%% Select the events corresponding to the left and right peaks of the time sums

if plotCutHitmap
    disp('Plotting time cut hitmap...')
    figures.timeCutHitmapPlot = figure(28);
    clf(figures.timeCutHitmapPlot)
    hold on
    set(gcf, 'Name', 'MCP 2D-plot with time cuts')

    gaussians = gaussianFits{1};
    cut = (gaussians.b2 - gaussians.b1)/(gaussians.c2 + gaussians.c1) * gaussians.c1 + gaussians.b1;

    for l = 1:1+compare
        if compare
            subplot(1, 2, l)
            hold on
        end
        if l == 1
            title('Cut MCP hitmap')
            calcTimeDiff = maskedTimeDiff;
            calcTimeSum = maskedTimeSum;
        else
            title('Cut MCP hitmap Compare')
            calcTimeDiff = compareTimeDiff;
            calcTimeSum = compareTimeSum;
        end
        less = find(calcTimeSum(:, 1) < cut);
        more = find(calcTimeSum(:, 1) > cut);
        title(['Cut at ' num2str(cut) 's for channels ' num2str(channelGroups(1, 1)) ' and ' num2str(channelGroups(1, 2))])
        scatter(calcTimeDiff(less, 1), calcTimeDiff(less, 2), markerSize, 'b')
        scatter(calcTimeDiff(more, 1), calcTimeDiff(more, 2), markerSize, 'r')
        xlabel('$x\propto \Delta t_x$', 'Interpreter', 'LaTeX');
        ylabel('$y\propto \Delta t_y$', 'Interpreter', 'LaTeX');
        legend('Short times', 'Long times')
        axis square
    end
    if compare
        suptitle('Events cut in time histograms')
    end
end

%% End timing for the analysis
toc

%Set urgency hint for the MATLAB terminal after execution
!sh seturgent.sh
