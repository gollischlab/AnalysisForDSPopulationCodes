path = ['sample data', filesep];
spikeFiles = {'13_SP_C807.txt', '13_SP_C1803.txt'}; % files contain spike times in seconds
stimFile = 'motionSteps.mat';
frameFile = '13_OMB_bg4x4corr8_C150_Gsteps3_frametimings.mat';
frameOffset = 25;  % offset in ms
samplingRate = 30; % in Hz
stimLen = 15*60*samplingRate; % 15 min trajectory
colors = lines(4);

% parameters for linear filters and CCA analysis
filterWindow = .8; % 800 ms filters
filterWindow2 = 2; % 2 second segments for CCA
filterLen = ceil(filterWindow*samplingRate);
filterLen2 = ceil(2*samplingRate);

% parameters for nonlinearity
nBins = 25;
NLfitType = 'u';
tbin = 1/samplingRate;

load([path, frameFile], 'ftimes');   % frametimes in milliseconds
ftimes = (ftimes(1:stimLen)-frameOffset)/1000;  

load([path, stimFile], 'stimulus');
stimFrames = stimulus(:, 1:stimLen) - .5;
nDims = size(stimFrames, 1);

nCells = numel(spikeFiles);
spikeCounts = zeros(2, stimLen);
for cellIdx = 1:nCells
	spikes = load([path, spikeFiles{cellIdx}], '-ascii');
	spikeCounts(cellIdx, :) = histc(spikes, ftimes);
end

% stimulus reconstruction and mutual information (Fig. 1)
[reconstrFrames, reconstrStimBins, populationFilter] = getLinearPopulationReadout(spikeCounts, stimFrames, filterLen);
[totalInfo, infoDensity, stimDensity, reconstrDensity, errorDensity] = calcMutualInformation(stimFrames(:, reconstrStimBins), reconstrFrames, filterLen);

% cell filters and nonlinearities (Fig. 5)
filters = zeros(nCells, nDims*filterLen);
NLrates = zeros(nCells, nBins);
NLbins = zeros(nCells, nBins);
for cellIdx = 1:nCells
	filters(cellIdx, :) = getFilter(spikeCounts(cellIdx, :), stimFrames, filterLen);
    [NLrates(cellIdx, :), NLbins(cellIdx, :)] = getNL(spikeCounts(cellIdx, :), stimFrames, filters(cellIdx, :), tbin, nBins);
end

% LN model (Fig. 5)
PoissonSpikeCounts = LNmodel(stimFrames(1, :), filters(:, 1:filterLen), NLbins(1, :), NLrates(1, :), NLfitType, max(spikeCounts(:)));
simFilters = zeros(nCells, filterLen);
simNLrates = zeros(nCells, nBins);
simNLbins = zeros(nCells, nBins);
for cellIdx = 1:nCells
	simFilters(cellIdx, :) = getFilter(PoissonSpikeCounts(cellIdx, :), stimFrames(1, :), filterLen);
    [simNLrates(cellIdx, :), simNLbins(cellIdx, :)] = getNL(PoissonSpikeCounts(cellIdx, :), stimFrames(1, :), simFilters(cellIdx, :), tbin, nBins);
end
[simReconstrFrames, simReconstrStimBins, simPopulationFilter] = getLinearPopulationReadout(PoissonSpikeCounts, stimFrames(1, :), filterLen);
[simTotalInfo, simInfoDensity, simStimDensity, simReconstrDensity, simErrorDensity] = calcMutualInformation(stimFrames(1, simReconstrStimBins), simReconstrFrames, filterLen);

% canonical correlation analysis (Fig. 7)
[respComp, stimComp, corrCoeff] = getCCAcomponents(spikeCounts, stimFrames, filterLen2);

%% plot filters and nonlinearities
figure;
for cellIdx = 1:nCells
    subplot(2, 2, 2*(cellIdx-1)+1);
    plot(-filterWindow+tbin:tbin:0, filters(cellIdx, 1:filterLen), 'color', colors(3, :), 'linewidth', 1.5);
    hold on
    plot(-filterWindow+tbin:tbin:0, filters(cellIdx, filterLen+1:end), 'color', colors(4, :), 'linewidth', 1.5);
    xlim([-filterWindow, 0]);
    title(['Filter - cell ', int2str(cellIdx)]);
    if cellIdx == 2
        xlabel('Time (s)');
        legend('x-dir.', 'y-dir.', 'Location', 'northwest');
    end
    
    subplot(2, 2, 2*cellIdx);
    plot(NLbins(cellIdx, :), NLrates(cellIdx, :), 'k', 'linewidth', 1.5);
    ylim([0 10]);
    ylabel('Average response (Hz)');
    title(['Nonlinearity - cell ', int2str(cellIdx)]);
    if cellIdx == 2
        xlabel('Stimulus projection');
    end
end

%% plot filters and nonlinearities of simulated data
figure;
for cellIdx = 1:nCells
    subplot(2, 2, 2*(cellIdx-1)+1);
    plot(-filterWindow+tbin:tbin:0, filters(cellIdx, 1:filterLen), 'color', [.7, .7, .7]);
    hold on
    plot(-filterWindow+tbin:tbin:0, simFilters(cellIdx, :), 'k', 'linewidth', 1.5);
    xlim([-filterWindow, 0]);
    title(['Simulated x-filter - cell ', int2str(cellIdx)]);
    if cellIdx == 2
        xlabel('Time (s)');
        legend('original', 'simulated', 'Location', 'northwest');
    end
    
    subplot(2, 2, 2*cellIdx);
    plot(NLbins(cellIdx, :), NLrates(cellIdx, :), 'color', [.7, .7, .7]);
    hold on
    plot(simNLbins(cellIdx, :), simNLrates(cellIdx, :), 'k', 'linewidth', 1.5);
    ylim([0 10]);
    ylabel('Average response (Hz)');
    title(['Nonlinearity - cell ', int2str(cellIdx)]);
    if cellIdx == 2
        xlabel('Stimulus projection');
    end
end

%% plot first five CCA components
figure;
for comp = 1:5
    subplot(5, 3, 3*(comp-1)+1);
    plot(0:tbin:filterWindow2-tbin, stimComp(1:filterLen2, comp), 'color', colors(3, :), 'linewidth', 1.5);
    hold on
    plot(0:tbin:filterWindow2-tbin, stimComp(filterLen2+1:end, comp), 'color', colors(4, :), 'linewidth', 1.5);
    ylabel(['Comp. ', int2str(comp)]);
    if comp == 1
        title('Stimulus motion');
    end
    if comp == 5
        xlabel('Time (s)');
    end
    
    subplot(5, 3, 3*(comp-1)+2);
    plot(0:tbin:filterWindow2-tbin, respComp(1:filterLen2, comp), 'color', colors(1, :), 'linewidth', 1.5);
    hold on
    plot(0:tbin:filterWindow2-tbin, respComp(filterLen2+1:end, comp), 'color', colors(2, :), 'linewidth', 1.5);
    if comp == 1
        title('Cell activity');
    end
    if comp == 5
        xlabel('Time (s)');
    end
end

subplot(2, 3, 3);
plot(corrCoeff, '.k');
title('Corr. coeff.');
xlim([0, numel(corrCoeff)]);
xlabel('Components');

%% plot stimulus reconstruction and Gaussian smoothed stimulus
figure;
subplot(2, 2, 1);
plot(0:tbin:2*filterWindow-tbin, populationFilter, 'linewidth', 1.5);
legend('cell 1', 'cell 2');
title('Population filter');

subplot(2, 2, 2);
plot(linspace(0, samplingRate/2, size(infoDensity, 2)), sum(infoDensity, 1), 'k', 'linewidth', 1.5);
xlim([0, 10]);
xlabel('Frequency (Hz)');
title('Information density');

subplot(2, 1, 2);
takeFrames = 1:1000;
gaussKernel = exp(-(-5:5).^2/9);
gaussKernel = gaussKernel/sum(gaussKernel);
plot(ftimes(reconstrStimBins(takeFrames)), reconstrFrames(1, takeFrames), 'k', 'linewidth', 1.5);
hold on
plot(ftimes(reconstrStimBins(takeFrames)), conv(stimFrames(1, reconstrStimBins(takeFrames)), gaussKernel, 'same'), 'color', [.7, .7, .7]);
legend('reconstr.', 'smoothed stim');
title('Stimulus reconstruction');
xlabel('Time (s)');
ylabel('Motion steps (x-dir.)');

%% plot stimulus reconstruction from simulated responses and Gaussian smoothed stimulus
figure;
subplot(2, 2, 1);
plot(0:tbin:filterWindow-tbin, reshape(simPopulationFilter, [], 2), 'linewidth', 1.5);
legend('cell 1', 'cell 2');
title('Population filter - simulated data');

subplot(2, 2, 2);
plot(linspace(0, samplingRate/2, numel(simInfoDensity)), simInfoDensity, 'k', 'linewidth', 1.5);
xlim([0, 10]);
xlabel('Frequency (Hz)');
title('Information density');

subplot(2, 1, 2);
takeFrames = 1:1000;
gaussKernel = exp(-(-5:5).^2/9);
gaussKernel = gaussKernel/sum(gaussKernel);
plot(ftimes(simReconstrStimBins(takeFrames)), simReconstrFrames(1, takeFrames), 'k', 'linewidth', 1.5);
hold on
plot(ftimes(simReconstrStimBins(takeFrames)), conv(stimFrames(1, simReconstrStimBins(takeFrames)), gaussKernel, 'same'), 'color', [.7, .7, .7]);
legend('reconstr.', 'smoothed stim');
title('Stimulus reconstruction');
xlabel('Time (s)');
ylabel('Motion steps (x-dir.)');