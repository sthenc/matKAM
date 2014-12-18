%*************************************************************************
% Source separation with kernel additive modelling.
%*************************************************************************
% MATLAB implementation of the vocal/music extraction method presented in
% "Kernel additive models for source separation", submitted at IEEE TSP
% ------------------------------------------------------------------------
%
% This script performs separation of the poly-repeating part of music
% (usually the accompaniement) from the lead signal (usually the voice).
%
% Separation is performed through spatial kernel backfitting of a given
% number of repeating patterns whose period is estimated and of a source
% whose kernel is a simple local selection. An additional stable source is
% considered to account for slowly varying, stable and non-repeating 
% instrumental parts (assumed to be non vocals).
%
% As such, this program is only for reviewing purpose. Please mention the
% following paper when referencing this program:
%
% @article{KAM2014, 
%  AUTHOR = {A. Liutkus and D. Fitzgerald and Z. Rafii and B. Pardo and L. Daudet}, 
%  TITLE = {{Kernel Additive Models for Source Separation}}, 
%  JOURNAL = {{IEEE Transactions on Signal Processing}}, 
%  YEAR = {2014}, 
%  NOTE = {submitted}}
%
%*************************************************************************
% Copyright (C) 2014, Inria, Antoine Liutkus
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Affero General Public License as
%    published by the Free Software Foundation, either version 3 of the
%    License, or (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Affero General Public License for more details.
%
%    You should have received a copy of the GNU Affero General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.



%**************************************************************************
%Parameters
%**************************************************************************

%Where to find the files
% directory = '/mnt/data/Fer/diplomski/test_data/dom_govor_muzika';
directory = '/mnt/data/Fer/diplomski/test_data2/glazba_glasnija_od_vokala';
%directory = './data';
files=dir(sprintf('%s/*.wav', directory));
files =  {files(:).name};
Lmax = 180; %length to consider in seconds
mixdownRepet = 1; %if nonzero, will mix down all non vocal sources

%what prefix to append to the file for extracted signals
processing_name = 'KAM';

%where to save separated wav files
%output_dir = '/mnt/data/Fer/diplomski/test_data/dom_govor_muzika_output';
output_dir = '/mnt/data/Fer/diplomski/test_data2/glazba_glasnija_od_vokala_output';

%Signal and STFT
windowLength = 50; % window length in ms
overlap = 0.85;    % overlap ratio between 0 and 1


%Parameters for the kernels
%--------------------------
%1) Repeating patterns (music)
fMaxTempo = 18000; %max frequency to use for period selection
nPeriods = 5;      %number of periods to separate 
minPeriod = 0.7;     %in seconds
maxPeriod = 20;    %in seconds

%definition of the different frequency bands. Each one may have a different
%stability parameter
fstop_bands = [3000  , 8000];
%The stability P indicates for how many periods each pattern is stable. A
%high value will lead to the use of more points during the median fitlering
%and thus will eliminate voice more. However, it will lead to a separated
%voice which contains more music.
P           = [3,     2]; %stability for each band

%2) Stable harmonic source
%stability in seconds for each band. If it is two small, it will catch the
%voice. If it is too large, it will be useless and the voice will catch the
%slowly evolving harmonic parts
includeStableSourceModel = 1;
harmonicStability = [1.0, 0.8]; 

%3) Vocals (lead)
includeVoiceModel = 1;
voiceKernel = [50, 0.04];% [Hertz, seconds], size of the voice kernel
fMinVoice=70;      %min frequency of voice


%Parameters for the backfitting algorithm
%----------------------------------------
niter= 5;          %number of iterations
parallel = 0;       %parallel version. may be faster, but uses a lot of RAM

%Display parameters
%------------------
displaySpectrograms = 0;
displayPeriods = 0;

%**************************************************************************
% END OF PARAMETERS.
%**************************************************************************


%Including helper files in directory 'includes'
addpath('includes')

%Closing parallel workers in case parallelization is not wanted
if ~parallel&&exist('matlabpool')&&matlabpool('size')
	matlabpool close
end

%Loop over the files to separate
for ifile = 1:length(files)
	file = sprintf('%s/%s',directory,files{ifile});
	[pathstr, name, ext] = fileparts(file);
	fprintf('Handling file %s\n',files{ifile});
	
	%Loading data 
	%------------
	fprintf('Loading signal\n');
	sig = Signal(file);

	
	%Truncating audio if needs be
	%----------------------------
	len = size(sig.s,1);
	range = max(1,round(len/2-Lmax*sig.fs/2)):min(len, round(len/2+Lmax*sig.fs/2));
	sig.s = sig.s(range,:);
	
	%Compute STFT
	%------------
	sig.windowLength = windowLength;
	sig.overlapRatio = overlap;
	fprintf('Computing STFT\n');
	sig.STFT;
	hopSize = sig.windowLength*(1-sig.overlapRatio)/1000;
	
	%Getting power spectrogram for periods estimation
	%------------------------------------------------
	X = single(sig.S(1:sig.nfftUtil,:,:));
	V = single(max(eps,abs(X).^0.3)); %to the 0.3th power for a reason to make clear in a further study
	[F,T,I] = size(V);
	
	%Normalizing spectrograms so that all frequencies have same mean energy
	%this is for the automatic period detection
	%----------------------------------------------------------------------
	weights = mean(V,2);
	V = single(bsxfun(@times,V,1./weights));
	
	%compute min and max periods in terms of number of frames   
	%--------------------------------------------------------
	nminPeriod = ceil(minPeriod/hopSize);
	nmaxPeriod= min(floor(T/max(P)),ceil(maxPeriod/hopSize));
	Imax = floor(min(fMaxTempo,sig.fs)/sig.fs*F);
	if ~isempty(fstop_bands)
	    poscuts=floor(min(fstop_bands,sig.fs)/sig.fs*sig.nfft);
	else
	    poscuts = floor(sig.nfft/2);
	end
	
	%compute beat data, simple sum of autocorrelations of the lines of the 
	%spectrogram: xcorr(x)=ifft(abs(fft(x)))
	%---------------------------------------------------------------------
	tempo = sum(ifft(abs(fft((sum(V(1:Imax,:,:),3)'))).^2),2);
	tempo(nminPeriod:nmaxPeriod) = ...
	    tempo(nminPeriod:nmaxPeriod).* linspace(1,0.5,nmaxPeriod-nminPeriod+1)';
	
	%Pick periods using magic peak peaking algorithm =)
	%--------------------------------------------------
	periods = pickpeaks(tempo(nminPeriod:nmaxPeriod),nPeriods,0)+nminPeriod-2;
	if displayPeriods
	    figure('units','normalized','outerposition',[0 0 1 1], 1);
	    clf
	    timePeriods= (0:length(tempo)-1)*hopSize;
	    plot(timePeriods(nminPeriod:nmaxPeriod),tempo(nminPeriod:nmaxPeriod),'LineWidth',2);
	    hold on;
	    plot(timePeriods(periods+1),tempo(periods+1),'ro','LineWidth',2,'MarkerSize',15);
	    grid on
	    xlabel('period (s)');
	    ylabel('spectral autocorrelation');
	    legend('average auto-correlation ov','selected periods')
	    title('Periods selection');
	    drawnow;
	end

	%Creating kernels
	%----------------
	%1) for repeating part
	kernels = {};
	for j = 1:length(periods)
	    %For all repeating patterns, kernel is a periodic line
	    kernels{end+1} = {};
	    for band = 1:length(poscuts)
	        kernels{end}{band} = periodicKernel(periods(j),P(band),0,0);
	    end
	end
	
	if includeStableSourceModel
	    %2) Harmonic part
	    kernels{end+1} = {};
	    for band = 1:length(poscuts)
	        nharmT = max(1,round(harmonicStability(band)/hopSize));
	        kernels{end}{band}=ones(1,nharmT);
	    end
	end
	
	if includeVoiceModel
	    %3) For voice part, it's a 'cross' like
	    %  0 1 0
	    %  1 1 1
	    %  0 1 0
	    nvoiceF = max(1,round(voiceKernel(1)*sig.nfft/sig.fs));
	    nvoiceT = max(1,round(voiceKernel(2)/hopSize));
	    kernels{end+1} = zeros(2*nvoiceF+1,2*nvoiceT+1);
	    kernels{end}(nvoiceF+1,:)=1;
	    kernels{end}(:,nvoiceT+1)=1;
	    %Computing min voice frequency
	    IminVoice = floor(min(fMinVoice,sig.fs)/sig.fs*sig.nfft);
	end
	
	%J is the number of sources. 
	J = length(kernels);
	
	tic 
	
	% Initializing model
	%-------------------
	%PSD are initialized simply as mixtures PSD /J
	S = single(repmat(sum(V,3),[1,1,J]))/J; 
	%All spatial covariance matrices as initially identity
	R = zeros(1,I,I,J); 
	for j = 1:J
	    R(1,:,:,j) = eye(I);
	end
	R = repmat(R,F,[1,1,1]);
	
	%Performing estimation through kernel backfitting
	%--------------------------------------------------------
	for it=1:niter
	    %Separating sources with current model
	    out = posterior(X,S,'R',R,'Ymu',1,'parallel',parallel);
	    Ymu = out.Ymu;
	    
	    %Backfitting each source
	    parfor j=1:J 
	        fprintf('Backfitting, iteration %d / %d, source %d / %d \n',it,niter,j,J);
	        tempS = zeros(F,T);
	        
	        %learning spatial covariance and PSD from separated image
	        [Z,Rj] = learnDSP(Ymu(:,:,:,j),1);
	        
	        %nan or inf will occur if all PSD estimates are 0, typical for voice below
	        %minimal frequency
	        Rj(isnan(Rj)) = 0;
	        Rj(isinf(Rj)) = 0;
	        R(:,:,:,j)=Rj;
	        
	        %median filter the estimated PSD Z
	        if iscell(kernels{j})
	            %we handle the case were we have several frequency bands
	            pos = 1;
	            for k = 1:length(poscuts)
	                order= round(length(find(kernels{j}{k}))/2);
	                tempS(pos:poscuts(k),:) = ordfilt2(Z(pos:poscuts(k),:),order,kernels{j}{k},'symmetric');
	                pos = poscuts(k)+1;
	            end
	        else
	            order= round(length(find(kernels{j}))/2);
	            tempS = ordfilt2(Z,order,kernels{j},'symmetric');
	        end
	            

	        if includeVoiceModel && fMinVoice && (j==J)
	            %if source is voice, then set to 0 its low frequencies if
	            %needed
	            tempS(1:IminVoice,:) = 0;                                
	        end 
	        S(:,:,j) = tempS;
	    end
	    

	    
	    %display if needed
	    if displaySpectrograms
	        figure(2)
	        clf
	        for j=1:J
	            prows =  min(factor(J));
	            if prows == J
	                subplot(1,J,j);
	            else
	                subplot(prows,J/prows,j);
	            end
	            %imagesc(log(max(1E-10,S(:,:,j))))
	            imagesc(real(S(:,:,j).^(0.3/2)));
	            colorbar
	        end
	        drawnow
	        pause(0.02)
	    end
	end
	
	%Finished. Now doing the final separation
	%----------------------------------------
	out = posterior(X,S,'R',R,'Ymu',1);
	disp('Rendering final wave files...')
	extractedSig = zeros(size(sig.s));
	for j = 1:J
	    %Rendering each each source, and adding it to the background if a
	    %mixdown is asked for, or exporting its wave in the contrary case
	    sig.S = sig.buildComplete(out.Ymu(:,:,:,j));
	    if mixdownRepet
	        extractedSig = extractedSig + sig.iSTFT();
	        if j<J, sourceName = 'background'; else sourceName='lead';end
	    else
	        extractedSig = sig.iSTFT();
	        sourceName = num2str(j);
	    end
	    if (j>=J-1) || (~mixdownRepet)
	        estimatedFilename=fullfile(output_dir,sprintf('%s_%s_%s.wav',processing_name,name,sourceName));
	        wavwrite(extractedSig,sig.fs,estimatedFilename);
	        extractedSig = zeros(size(sig.s));
	    end
	end  
	temps(ifile) = toc;
	duree(ifile) = sig.sLength/sig.fs;
	fprintf('Done. %0.1fs : computing speed %0.2f s/s \n',temps(ifile),duree(ifile)/temps(ifile))
end

