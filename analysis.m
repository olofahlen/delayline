%%Initial
disp('Starting data analysis...')

tic

plotCharges = true;
plotFdhm = true;
plotTimeSums = true;
plotPositions = true;
plotTimeCutHitmap = true;

allMask = [1:nbrOfMeas]';
secondQuadrantMask = find(timeDiff(:, 1) < 0 & timeDiff(:, 2) > 0);

mask = secondQuadrantMask;

maskedCharge = charge(mask, :);
maskedTotalCharge = totalCharge(mask, :);
maskedGrandTotalCharge = grandTotalCharge(mask, :);
maskedFdhm = fdhm(mask, :);
maskedTimeSum = timeSum(mask, :);
maskedTimeDiff = timeDiff(mask, :);

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
        hist(maskedCharge(:, j), interval)
        %hist(maskedCharge(:, j), bins)
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
        hist(maskedTotalCharge(:, k), interval)
        %hist(maskedTotalCharge(:, k), bins)
    end
    subplot(3, 1, 3)
    hist(maskedGrandTotalCharge, bins)
    title(['Total charge deposited both delay lines'])
    xlabel('Charge [e]')
    ylabel('Counts')
    suptitle('Histograms of charges for the different channels')
    chargeScatterPlot = figure(32);
    clf(32)
    hold on
    title('Correlation of charge deposited on the two delay lines')
    set(gcf, 'Name', 'Charge scatter plot')
    scatter(maskedTotalCharge(:, 1), maskedTotalCharge(:, 2))
    axis square
    xlabel(['Total charge deposited on channels ' num2str(channelGroups(1, 1)) ' and ' num2str(channelGroups(1, 2)) ' [e]'])
    ylabel(['Total charge deposited on channels ' num2str(channelGroups(2, 1)) ' and ' num2str(channelGroups(2, 2)) ' [e]'])
end

%% Plot histograms for FDHM
if plotFdhm
    disp('Plotting FDHM histograms...')
    fdhmHistPlot = figure(35);
    clf(fdhmHistPlot)
    set(gcf, 'Name', 'FDHM of the four channels')
    hold on

    fdhmMeans = mean(maskedFdhm);
    fdhmStds = std(maskedFdhm);
    fdhmLimits = [fdhmMeans - 3*fdhmStds; fdhmMeans + 3*fdhmStds];
    for j = 1:channels
        binX = fdhmLimits(1, j):t/2:fdhmLimits(2, j);
        subplot(4, 1, j)
        hold on
        title(['FDHM for channel ' num2str(channelPairs(j))])
        xlabel('FDHM [s]')
        ylabel('Counts')
        hist(maskedFdhm(:, j), binX)
        axis([min(fdhmLimits(1, :)) max(fdhmLimits(2, :)) 0 1])
        axis 'auto y'
    end
    suptitle(['Distribution of FDHM'])
end

%%  Plot histograms for time sums and fit a double Gaussian
disp('Fitting double Gaussian...')
for k = 1:2
    tMean = mean(maskedTimeSum(:, k));
    tStd = std(maskedTimeSum(:, k));
    interval = tMean - 3*tStd : t/2 : tMean + 3*tStd;
    [N, x] = hist(maskedTimeSum(:, k), interval);
    x = x(2:end-1)';
    N = N(2:end-1)';
    timeSumX{k} = x;
    timeSumN{k} = N;

    fitObj = fittype('gauss2');
    options = fitoptions('gauss2');
    options.Lower = [0 -Inf 0 0 -Inf 0];
    [fittedGaussians gof output] = fit(x, N, fitObj, options);%, 'Weight', sqrt(N));
    gaussianFits{k} = fittedGaussians;
end

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
        title(['Delayline for channels ' num2str(channelGroups(k, 1)) ' and ' num2str(channelGroups(k, 2)) sprintf('. a = %.4f, b = %.4f, mu_1 = %.4f ns, mu_2 = %.4f ns, sigma_1 = %.4f ns, sigma_2 = %f ns', fittedGaussians.a1, fittedGaussians.a2, 1e9*fittedGaussians.b1, 1e9*fittedGaussians.b2, 1e9*fittedGaussians.c1, 1e9*fittedGaussians.c2)])
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
    hist(maskedTimeDiff(:, 1), bins)
    subplot(2, 1, 2)
    hold on
    title(['Delayline for channels ' num2str(channelPairs(3)) ' and ' num2str(channelPairs(4))])
    xlabel('$\Delta t$ [s]', 'Interpreter', 'LaTeX')
    ylabel('Counts')
    hist(maskedTimeDiff(:, 2), bins)
    suptitle('Histograms of time differences for the two delay lines')

    mcpHitmapPlot = figure(4);
    clf(4)
    set(gcf, 'Name', 'MCP 2D-plot')
    hold on
    suptitle('Reconstruction of particle hits on the MCP')
    scatter(maskedTimeDiff(:, 1), maskedTimeDiff(:, 2), 20, sum(maskedTotalCharge, 2 ), 'filled')
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
