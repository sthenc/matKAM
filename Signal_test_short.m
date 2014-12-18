% stjepan.henc@fer.hr

infile = 'data\britney-short0.wav';
outfile = 'output\britney-short0-inverse.wav';
stftfile = 'test-orig\britney-short0-stft.mat';
istftfile = 'test-orig\britney-short0-istft.mat';


addpath('includes');

%Signal and STFT
windowLength = 50; % window length in ms
overlap = 0.85;    % overlap ratio between 0 and 1

Lmax = 180; %length to consider in seconds -- wut?

tic;

sig = Signal(infile);

sig.singlePrecision = 0;

len = size(sig.s,1);
range = max(1,round(len/2-Lmax*sig.fs/2)):min(len, round(len/2+Lmax*sig.fs/2));
sig.s = sig.s(range,:);

sig.windowLength = windowLength;
sig.overlapRatio = overlap;
fprintf('Computing STFT\n');
sig.STFT; % transform

stft = sig.S; % output STFT data


% test inverse transform
extractedSig = sig.iSTFT; 
istft = extractedSig;


wavwrite(extractedSig, sig.fs, outfile);

toc;

save(stftfile, 'stft');
save(istftfile, 'istft');
