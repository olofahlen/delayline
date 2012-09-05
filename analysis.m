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
plotTimeCutHitmap = true;

allMask = [1:nbrOfMeas]';
smallLeftSquareMask = find(-6e-8 < timeDiff(:, 1) & timeDiff(:, 1) < -4e-8 & 0 < timeDiff(:, 2) & timeDiff(:, 2) < 2e-8);
secondQuadrantMask = find(timeDiff(:, 1) < 0 & timeDiff(:, 2) > 0);
shortFdhmMask = find(fdhm(:, 1) < 7.76e-9);

mask = allMask;

maskedCharge = charge(mask, :);
maskedTotalCharge = totalCharge(mask, :);
maskedGrandTotalCharge = grandTotalCharge(mask, :);
maskedFdhm = fdhm(mask, :);
maskedTimeSum = timeSum(mask, :);
maskedTimeDiff = timeDiff(mask, :);

compare = false;
compareMask = secondQuadrantMask;

compareCharge = charge(compareMask, :);
compareTotalCharge = totalCharge(compareMask, :);
compareGrandTotalCharge = grandTotalCharge(compareMask, :);
compareFdhm = fdhm(compareMask, :);
compareTimeSum = timeSum(compareMask, :);
compareTimeDiff = timeDiff(compareMask, :);

%% Plot histograms with charges
if plotCharges
    disp('Plotting charges...')
    bins = 100;
    interval = linspace(0, 14e6, bins);
    figures.individualChargePlot = figure(21);
    clf(figures.individualChargePlot)
    set(gcf, 'Name', 'Individual charge histograms')
    for j = 1:channels
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

    interval = linspace(0, 3e7, 100);
    figures.totalChargePlot = figure(22);
    clf(figures.totalChargePlot)
    set(gcf, 'Name', 'Total charge histograms')
    for k = 1:2
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
    interval = linspace(min(maskedGrandTotalCharge), max(maskedGrandTotalCharge), 100);
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
    scatter(maskedTotalCharge(:, 1), maskedTotalCharge(:, 2), 4)
    if compare
        scatter(compareTotalCharge(:, 1), compareTotalCharge(:, 2), 4, 'r')
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

    fdhmMeans = mean(maskedFdhm);
    fdhmStds = std(maskedFdhm);
    if compare
        fdhmMeansCompare = mean(compareFdhm);
        fdhmStdsCompare = std(compareFdhm);
        fdhmLimits = [max(0, min(fdhmMeans - 3*fdhmStds, fdhmMeansCompare - 3*fdhmStdsCompare)); max(fdhmMeans + 3*fdhmStds, fdhmMeansCompare + 3*fdhmStdsCompare)];
    else
        fdhmLimits = [max(0, fdhmMeans - 3*fdhmStds); fdhmMeans + 3*fdhmStds];
    end
    for j = 1:channels
        binX = fdhmLimits(1, j):t/2:fdhmLimits(2, j);
        subplot(4, 1, j)
        hold on
        title(['FDHM for channel ' num2str(channelPairs(j))])
        xlabel('FDHM [ns]')
        ylabel('Counts')
        hist(1e9*maskedFdhm(:, j), 1e9*binX)
        if compare
            h = findobj(gca, 'Type', 'patch');
            set(h, 'FaceColor', 'b', 'EdgeColor', 'w', 'facealpha', 0.75)
            hist(1e9*compareFdhm(:, j), 1e9*binX)
            newH = findobj(gca, 'Type', 'patch');
            newH = newH(find(newH ~= h));
            set(newH, 'FaceColor', 'r', 'EdgeColor', 'w', 'facealpha', 0.75)
        end
        axis(1e9*[min(fdhmLimits(1, :)) max(fdhmLimits(2, :)) 0 1])
        axis 'auto y'
    end
    suptitle('Distribution of FDHM')
end

%%  Plot histograms for time sums and fit a double Gaussian
disp('Fitting double Gaussian...')
for l = 1:1 + compare
    if l == 1
        calcTimeSum = maskedTimeSum;
    else
        calcTimeSum = compareTimeSum;
    end
    for k = 1:2
        tMean = mean(calcTimeSum(:, k));
        tStd = std(calcTimeSum(:, k));
        interval = tMean - 3*tStd : t/2 : tMean + 3*tStd;
        [N, x] = hist(calcTimeSum(:, k), interval);
        x = x(2:end-1)';
        N = N(2:end-1)';
        timeSumX{(l-1)*2 + k} = x;
        timeSumN{(l-1)*2 + k} = N;

        fitObj = fittype('gauss2');
        options = fitoptions('gauss2');
        options.Lower = [0 -Inf 0 0 -Inf 0];
        [fittedGaussians gof output] = fit(x, N, fitObj, options);
        gaussianFits{(l-1)*2 + k} = fittedGaussians;
    end
end

if plotTimeSums
    disp('Plotting time sum histograms with double Guassian fit...')
    figures.timeSumHistPlot = figure(25);
    clf(figures.timeSumHistPlot)
    set(gcf, 'Name', 'Normalized histograms of time sums')
    for k = 1:2
        subplot(2, 1, k)
        hold on
        xlabel('$t_1 + t_2$ [ns]', 'Interpreter', 'LaTeX')
        ylabel('Counts')
        bar(1e9*timeSumX{k}, timeSumN{k})
        if compare
            h = findobj(gca, 'Type', 'patch');
            set(h, 'FaceColor', 'b', 'EdgeColor', 'w', 'facealpha', 0.75)
            bar(1e9*timeSumX{2 + k}, timeSumN{2 + k})
            newH = findobj(gca, 'Type', 'patch');
            newH = newH(find(newH ~= h));
            set(newH, 'FaceColor', 'r', 'EdgeColor', 'w', 'facealpha', 0.75)
        end
        %bar(timeSumX{k}, timeSumN{k})

        fittedGaussians = gaussianFits{k};
        fitPlot = plot(1e9*timeSumX{k}, fittedGaussians(timeSumX{k}), 'r');
        %fitPlot = plot(fittedGaussians, 'r');
        title(['Delayline for channels ' num2str(channelGroups(k, 1)) ' and ' num2str(channelGroups(k, 2)) sprintf('. a = %.4f, b = %.4f, mu_1 = %.4f ns, mu_2 = %.4f ns, sigma_1 = %.4f ns, sigma_2 = %f ns', fittedGaussians.a1, fittedGaussians.a2, 1e9*fittedGaussians.b1, 1e9*fittedGaussians.b2, 1e9*fittedGaussians.c1, 1e9*fittedGaussians.c2)])
        if compare
            fittedGaussians = gaussianFits{2 + k};
            fitPlot = plot(1e9*timeSumX{k}, fittedGaussians(timeSumX{k}), 'r');
        end
    end
    suptitle('Time sums for the two delay lines')
end

%% Plot histograms for time differences and MCP hitmap
if plotPositions
    disp('Plotting time difference histograms and MCP hitmap...')
    bins = 500;
    figures.timeDiffHistPlot = figure(26);
    clf(figures.timeDiffHistPlot)
    set(gcf, 'Name', 'Histograms of time differences')
    subplot(2, 1, 1)
    hold on
    title(['Delayline for channels ' num2str(channelPairs(1)) ' and ' num2str(channelPairs(2))])
    xlabel('$\Delta t$ [ns]', 'Interpreter', 'LaTeX')
    ylabel('Counts')
    hist(1e9*maskedTimeDiff(:, 1), bins)
    subplot(2, 1, 2)
    hold on
    title(['Delayline for channels ' num2str(channelPairs(3)) ' and ' num2str(channelPairs(4))])
    xlabel('$\Delta t$ [ns]', 'Interpreter', 'LaTeX')
    ylabel('Counts')
    hist(1e9*maskedTimeDiff(:, 2), bins)
    suptitle('Histograms of time differences for the two delay lines')

    figures.mcpHitmapPlot = figure(27);
    clf(figures.mcpHitmapPlot)
    set(gcf, 'Name', 'MCP 2D-plot')
    hold on
    for l = 1:2
        subplot(1, 2, l)
        if l == 1
            calcTimeDiff = maskedTimeDiff;
            calcTotalCharge = maskedTotalCharge;
        else
            calcTimeDiff = compareTimeDiff;
            calcTotalCharge = compareTotalCharge;
        end
        scatter(calcTimeDiff(:, 1), calcTimeDiff(:, 2), 4, sum(calcTotalCharge, 2 ))
        xlabel('$x\propto \Delta t_x$', 'Interpreter', 'LaTeX');
        ylabel('$y\propto \Delta t_y$', 'Interpreter', 'LaTeX');
        axis square
        cb = colorbar;
        ylabel(cb, 'Charge [e]')
    end
    suptitle('Reconstruction of particle hits on the MCP')
end

%% Select the events corresponding to the left and right peaks of the time sums

if plotTimeCutHitmap
    disp('Plotting time cut hitmap...')
    figures.timeCutHitmapPlot = figure(28);
    clf(figures.timeCutHitmapPlot)
    hold on
    set(gcf, 'Name', 'MCP 2D-plot with time cuts')

    for k = 1:2
        gaussians = gaussianFits{k};
        cut = (gaussians.b2 - gaussians.b1)/(gaussians.c2 + gaussians.c1) * gaussians.c1 + gaussians.b1;
        less = find(maskedTimeSum(:, k) < cut);
        more = find(maskedTimeSum(:, k) > cut);

        subplot(1, 2, k)
        hold on
        title(['Cut at ' num2str(cut) 's for channels ' num2str(channelGroups(k, 1)) ' and ' num2str(channelGroups(k, 2))])
        scatter(maskedTimeDiff(less, 1), maskedTimeDiff(less, 2), 4, 'b')
        scatter(maskedTimeDiff(more, 1), maskedTimeDiff(more, 2), 4, 'r')
        xlabel('$x\propto \Delta t_x$', 'Interpreter', 'LaTeX');
        ylabel('$y\propto \Delta t_y$', 'Interpreter', 'LaTeX');
        legend('Short times', 'Long times')
        axis square
    end

    suptitle('Events cut in time histograms')
end

%% End timing for the analysis
toc

%Set urgency hint for the MATLAB terminal after execution
!sh seturgent.sh
