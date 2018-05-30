function [totalInfo, infoDensity, stimDensity, reconstrDensity, errorDensity] = calcMutualInformation(stimFrames, reconstrFrames, filterLen)
% calcMutualInformation - calculate mutual information between stimulus and reconstruction
% adopted from Warland, Reinagel & Meister, 1997, "Decoding visual information from a population of retinal ganglion cells"
%   stimFrames: 2 x T array of motion steps in x- and y-direction
%   reconstrFrames: array of reconstructed motion steps of same size as 'stimFrames'
%   filterLen: filter length in bins

% Determine stimulus, reconstruction and reconstruction error in frequency domain
% by calculating the Fourier transform of segments of length 'filterLen'
nDims = size(stimFrames, 1);
nSegments = floor(size(stimFrames, 2)/filterLen);
NFFT = power(2, nextpow2(filterLen));

stimFFT = zeros(nDims, nSegments, NFFT);
reconstrFFT = zeros(nDims, nSegments, NFFT);
errorFFT = zeros(nDims, nSegments, NFFT);
for k = 1:nSegments
    stimFFT(:, k, :) = fft(stimFrames(:, (k-1)*filterLen+1:k*filterLen), NFFT, 2);
	reconstrFFT(:, k, :) = fft(reconstrFrames(:, (k-1)*filterLen+1:k*filterLen), NFFT, 2);
	errorFFT(:, k, :) = fft(reconstrFrames(:, (k-1)*filterLen+1:k*filterLen) - stimFrames(:, (k-1)*filterLen+1:k*filterLen), NFFT, 2);
end

% take average across segments
stimDensity = mean(2*power(abs(stimFFT(:, :, 1:NFFT/2+1)), 2), 2);
reconstrDensity = mean(2*power(abs(reconstrFFT(:, :, 1:NFFT/2+1)), 2), 2);
errorDensity = mean(2*power(abs(errorFFT(:, :, 1:NFFT/2+1)), 2), 2);

stimDensity = squeeze(stimDensity)/power(filterLen, 2);
reconstrDensity = squeeze(reconstrDensity)/power(filterLen, 2);
errorDensity = squeeze(errorDensity)/power(filterLen, 2);

infoDensity = log2(stimDensity./errorDensity);
totalInfo = sum(infoDensity(:));
end
