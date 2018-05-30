function [respComp, stimComp, corrCoeff] = getCCAcomponents(spikeCounts, stimFrames, filterLen)
% getCCAcomponents - get stimulus and response components from Hotelling's Canonical Correlation Analysis 
% adopted from Macke, Zeck & Bethge, 2008, "Receptive fields without spike-triggering"
%   spikeCounts: 2 x T array of spike counts of two cells
%   stimFrames: 2 x T array of motion steps in x- and y-direction
%   filterLen: filter length in bins

if nargin == 0
	[fileName, pathName] = uigetfile('*.mat', 'Spike counts missing, get data from file...');
	load([pathName, fileName], 'spikeCounts');
end
if nargin < 2
	[fileName, pathName] = uigetfile('*.mat', 'Stimulus missing, get data from file...');
	load([pathName, fileName], 'stimFrames');
end
if nargin < 3
		disp('No filter length given, using default of 25 bins');
		filterLen = 25;
end

trainFrac = 0.7; % fraction used for training
stimLen = size(spikeCounts, 2);
totalFilterLen = filterLen*2;

% assign training and test bins
trainBins = 1:floor(stimLen*trainFrac);
testBins = (floor(stimLen*trainFrac)+1):stimLen;

trainLen = numel(trainBins)-filterLen;
testLen = numel(testBins)-filterLen;

trainStim = stimFrames(:, trainBins);
testStim = stimFrames(:, testBins);

trainSpikeCounts = spikeCounts(:, trainBins);
testSpikeCounts = spikeCounts(:, testBins);

% reshape spike count matrices into response fragments of filter length
trainR = zeros(trainLen, totalFilterLen);
testR = zeros(testLen, totalFilterLen);
for k = 1:trainLen
    trainR(k, :) = reshape(trainSpikeCounts(:, k:k+filterLen-1)', 1, []);
end
for k = 1:testLen
    testR(k, :) = reshape(testSpikeCounts(:, k:k+filterLen-1)', 1, []);
end

% reshape stimulus into matrix of stimulus fragments of size 'stimFragmentLen'
trainS = zeros(size(trainR));
testS = zeros(size(testR));
for k = 1:trainLen
    trainS(k, :) = reshape(trainStim(:, k:k+filterLen-1)', [], 1)';
end
for k = 1:testLen
    testS(k, :) = reshape(testStim(:, k:k+filterLen-1)', [], 1)';
end

% calculate covariances for training stimulus and response
covTrainStim = cov(trainS);
covTrainResp = cov(trainR);
covTrainStimResp = trainS'*trainR/size(trainS, 1)-(mean(trainS, 1)'*mean(trainR, 1));
whiteXcov = covTrainStim^(-0.5)*covTrainStimResp*covTrainResp^(-0.5);

% calculate SVD of whitened cross-covariance
[stimBase, D, respBase] = svd(whiteXcov);

% get stimulus and response components from SVD bases with test stimulus and response matrices
covTestStim = cov(testS);
covTestResp = cov(testR);
stimComp = covTestStim^(-0.5)*stimBase;
respComp = covTestResp^(-0.5)*respBase;

% get correlation coefficients between stimulus and response components
corrCoeff = diag(D);

end