function [PoissonSpikeCounts, fitParams] = LNmodel(stimFrames, filters, NLbins, NLrates, NLfitType, maxBinCount)
% LNmodel - simulated responses to 1D motion trajectory using LN model with exponential or u-shaped nonlinearities
% see also Chichilnisky, 2001, "A simple white noise analysis of neuronal light responses"
%   spikeCounts: N x T array of spike counts of N cells
%   stimFrames: 1 x T array of motion stimulus in one direction
%   filters: N x L array of one-dimensional motion filters from N cells of lenght L
%   NLbins: 1 x M array of M nonlinearity bins
%   NLrates: 1 x M array of average firing rates per bin
%   NLfitType: type of applied nonlinearity function, can either be exponential ('exp') 
%           or an exponential multiplied with a quadratic ('u-shaped')  
%   maxBinCount: maximum number of spikes within a bin, for normalization

if nargin < 5
	NLfitType = 'u-shaped';
end
if nargin < 6
	maxBinCount = 4;
end

% assign nonlinearity functions
switch NLfitType
	case 'exp'
		NLfunc = @(A, x) A(1)*exp(A(2)*x);
	case {'u', 'u-shaped', 'ushaped'}
		NLfunc = @(A, x) A(1)*exp(A(2)*x).*x.^2;
	otherwise
		error("NLfitType must either be 'exp' or 'u-shaped'.");
end

% fit function to nonlinearity of first cell
OLS = @(A) sum((NLfunc(A, NLbins) - NLrates).^2); % Ordinary Least Squares cost function
opts = optimset('MaxFunEvals', 50000, 'MaxIter', 10000);
fitParams = fminsearch(OLS, rand(2,1), opts); % Use ‘fminsearch’ to minimise the ‘OLS’ function

% convolve stimulus with filters
% model spike output of all N cells, using determined nonlinearity function of first cell
nCells = size(filters, 1);
stimLen = size(stimFrames, 2);
filterLen = size(filters, 2);

effectiveStimLen = stimLen - filterLen;

% reshape stimulus into matrix of stimulus fragments of size 'stimFragmentLen'
stimMatrix = zeros(effectiveStimLen, filterLen);
for k = 1:effectiveStimLen
    stimMatrix(k, :) = reshape(stimFrames(:, k:k+filterLen-1)', [], 1)';
end

stimFiltered = zeros(nCells, effectiveStimLen);
for cell = 1:nCells
    filterNorm = sqrt(sum(filters(cell, :).*filters(cell, :)));
    stimFiltered(cell, :) = stimMatrix*filters(cell, :)'/filterNorm;
end

% obtain mean spike count of cells by applying nonlinearity,
% scale spike rate to spike count in tbin
meanSpikeCounts = NLfunc(fitParams, stimFiltered);

% adjust mean spike count to maximum spike count in 'spikeCounts'
meanSpikeCounts = meanSpikeCounts/max(meanSpikeCounts(:))*maxBinCount;

% plot nonlinearity
figure;
NLrange = NLbins(1):.1:NLbins(end);
plot(NLrange, NLfunc(fitParams, NLrange)/max(meanSpikeCounts(:))*maxBinCount, 'k', 'linewidth', 1.5);
hold on
plot(NLbins, NLrates, 'color', [.7 .7 .7]);
title('Nonlinearity for simulation');
legend('fitted', 'original');
ylabel('Average response (Hz)');
xlabel('Stimulus projection');

% calculate Poisson distributed spike counts
% with mean and standard deviation of 'meanSpikeCounts'
PoissonSpikeCounts = [zeros(nCells, filterLen), poissrnd(meanSpikeCounts)];
end