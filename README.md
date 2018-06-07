# AnalysisForDSPopulationCodes
Analysis methods for populations of direction-selective ganglion cells. Accompanying a manuscript by Kühn and Gollisch.

Contact: Tim Gollisch, tim.gollisch@med.uni-goettingen.de

Website: http://www.retina.uni-goettingen.de/

Description:
-----------
Functions are written for time-binned data, here, with 30 Hz sampling rate.
A minimal example for the use of each function can be found in 'MAIN.m'.

Inputs:
-----------
spikeCounts: activity of "N" neurons over a time of "T" time bins in N x T matrices counting spikes per time bin

stimFrames: motion steps in um per frame (time bin)

filterLen: length "L" of the filter in time bins

tbin: time bin size in seconds

nBins: number of nonlinearity bins, NLrates

NLfitType: type of applied nonlinearity function, either exponential ('exp') or exponential multiplied with a quadratic ('u-shaped')

maxBinCount: maximum number of spikes within a time bin for normalization of LN model parameters


Functions:
-----------
[reconstrFrames, reconstrStimBins, filter] = getLinearPopulationReadout(spikeCounts, stimFrames, filterLen)

[totalInfo, freqBins, infoDensity, stimDensity, reconstrDensity, errorDensity] = calcMutualInformation(stimFrames, reconstrFrames, filterLen, tbin)

filter = getFilter(spikeCounts, stimFrames, filterLen)

[NLrates, NLbins] = getNL(spikeCounts, stimFrames, filter, tbin, nBins)

[PoissonSpikeCounts, fitParams] = LNmodel(stimFrames, filters, NLbins, NLrates, NLfitType, maxBinCount)

[respComp, stimComp, corrCoeff] = getCCAcomponents(spikeCounts, stimFrames, filterLen)
