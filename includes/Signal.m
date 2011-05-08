classdef Signal < handle
%-------------------------------------------------------------------------
% Class name : Signal
%-------------------------------------------------------------------------
% Description : handles basic Signal operations and transforms. 
%               * wave loading
%               * STFT, with any overlap ratio, any frames length,
%                 any weighting function
%               * MDCT and its inverse
%               * Constant Q Transform (by Jacques Prado)
%               * splitting into frames
%               * pitch detection (for harmonic sounds)
%               * onset detection
%
% Main properties that are read/write 
%   * s : signal
%   * windowLength (ms)
%   * nfft (samples)
%   * overlapRatio (>=0 and <1)
%   * S : stft data
%   * Smdct : MDCT data
%   * CQTBinsPerOctave : number of bins per octave for cQt
%   * CQTMinFreq : min freq for the cQt
%   * CQTMaxFreq : max freq for the cQt
%   * CQTAlign : where to center cQt lobe weights, either 'l' for left or
%               'c' for center
%   * weightingFunction: handle to a function taking the number of points
%     N as an argument and returning the Nx1 weighting window
%
% Main properties that are read only : 
%   * sLength : signal length
%   * nChans : number of channels
%   * nfftUtil : number of bins in the positive frequency domain
%   * framesPositions, nFrames : positions and number of frames
%   * sWin, sWeights : windowed data
%   * Sq : cQt data
%
% examples : 
%   sig = Signal('myfile.wav'); %creates signal for a wave file
%   sig.windowLength = 70;      %sets windows of 70ms
%   sig.overlapRatio = 0.8;     %sets overlap ratio 0 < overlapRatio < 1
%   sig.STFT;                   %performs STFT
%   sig.CQT;                    %same for CQT
%   sig.MDCT;                   %same for MDCT (note that overlap is set to
%                               %               50 percents)
%   mySTFT = sig.S;
%   myMDCT = sig.S;
%   myCQT = sig.Sq
%
%   %modifying s.S for example
%   ...
%
%   waveform = sig.iSTFT();     % gets inverse transform
%
% Note that :
%   * all properties are automatically set relevantly in case of
%     modifications. for example, when nfft is set, windowLength is changed
%     accordingly
%   * STFT, MDCT and CQT produce exactly aligned data
%-------------------------------------------------------------------------
%
% Copyright (c) 2014, Antoine Liutkus, Inria
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions 
% are met:
%
%    *  Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%    *  Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the 
%       distribution.
%    *  Neither the name of Inria nor the names of its 
%       contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%-------------------------------------------------------------------------
    methods 
        %Constructor
        function self = Signal(varargin)
            %signal properties
            self.singlePrecision = 0;
            self.s = [];
            self.fs = 44100;
            self.sLength = 0;
            self.nChans = 0;
            self.weightingFunction = @hamming;
            
            %STFT properties
            self.S = [];
            self.windowLength = 60;
            self.nfft = 0;
            self.nfftUtil = 0;
            self.overlapRatio = 0.5;
            self.framesPositions = [];
            self.nFrames = 0;
            self.weightingWindow = [];
            self.overlap = 0;
            

            %Windowing properties
            self.sWin = [];
            self.sWeights = [];
            
            %Properties listeners
            self.propertiesListeners = {};
            self.propertiesListeners{end+1} = addlistener(self,'s','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'singlePrecision','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'fs','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'weightingFunction','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'windowLength' ,'PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'nfft' ,'PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'overlapRatio','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'S','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));    
            
            
            switch nargin
                case 0
                    return;
                case 1
                    switch class(varargin{1})
                        case 'char'
                            %varargin{1} is a file
                            self.LoadFromFile(varargin{1});
                        case 'Signal'
                            self = varargin{1}.copy;
                    end
                case 2
                    %varargin{1} is a signal and varargin{2} is the
                    %sampling frequency
                    self.LoadWF(varargin{1},varargin{2});

            end
        end
        
        function SigOut = getOnsets(self,fMin, fMax)
            %This function returns a vector of length nFrames, containing
            % an onset detection in the signal between frequencies fMin and
            % fMax. It uses the complex spectral difference method.
            if isempty(self.S)
                self.STFT
            end
            
            IFmin = max(1,round(fMin/self.fs*self.nfft));
            IFmax = round(min(fMax/2,fMax/self.fs*self.nfft));
            
            SigFrame = self.S(IFmin:IFmax,:,:);
            
            LenWindow = size(SigFrame,1);
            nCanaux = size(SigFrame,3);
            Nbre = size(SigFrame,2);
            SigOut = zeros(size(SigFrame,2),nCanaux);
            
            for canal = 1:nCanaux
                SigOutInst = zeros(LenWindow,1);
                CurrPhi = angle(SigFrame(:,1,canal));
                SigOut(1) = 0;
                Phi = (unwrap(angle(SigFrame(:,:,canal))));
                DevPhi = [zeros(LenWindow,2) (Phi(:,3:end) - 2*Phi(:,2:end-1) +...
                        Phi(:,1:end-2))];
                    for n = 2:Nbre
                            SigOutInst = sqrt(abs(SigFrame(:,n-1)).^2 + abs(SigFrame(:,n)).^2 -2*abs(SigFrame(:,n)).*abs(SigFrame(:,n-1)).*cos(DevPhi(:, n)));
                            SigOut(n,canal) = sum(SigOutInst);
                    end
                SigOut(:,canal) = 1/max(abs(SigOut(:,canal))).*SigOut(:,canal);
            end
        end  
        %Data loaders
        function LoadFromFile(self, file)
            [self.s, self.fs] = wavread(file);
            [self.sLength, self.nChans] = size(self.s);
        end 
        function LoadWF(self, waveform, fs)
            self.s = waveform;
            self.fs = fs;
            [self.sLength, self.nChans] = size(self.s);
        end
        
        %Transforms
        function S = STFT(self)
            %Builds Signal STFT, puts it in self.S and returns it.
            weighted_frames = self.split(0);
            self.S = fft(weighted_frames);  
            if self.singlePrecision
                self.S = single(self.S);
            end
            if nargout
                S = self.S;
            end
        end
        
        function V = specgram(self,nSteps)
            %Builds Signal spectrogram, 
            % through averaging nSteps delayed spectrograms
            if nargin == 1
                nSteps = 20;
            end
            weighted_frames = self.split(0);
            V = abs(fft(weighted_frames)).^2;  
            for index = 1:nSteps
                disp(sprintf('specgram, pass %d / %d',index,nSteps));drawnow
                weighted_frames = self.split(index);
                V = V*(1- 1/(index+1)) +1/(index+1)* abs(fft(weighted_frames)).^2;  
            end
            if self.singlePrecision
                V = single(self.S);
            end
        end
        
        function s = iSTFT(self)
            %Computes inverse STFT, puts it in self.s and returns it.
            s = [];
            if isempty(self.S)
                return;
            end
            tempSignal = real(ifft(self.S));
    
            tempLength = self.framesPositions(end)+self.nfft - 1;     
            oldLength = self.sLength;       
            
            s_temp=zeros(tempLength,self.nChans);
            W = zeros(tempLength,self.nChans);

            for index_frame=1:self.nFrames
                    for index_chan=1:self.nChans
                        s_temp(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan)= ...
                            s_temp(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan) + tempSignal(:,index_frame,index_chan).*self.weightingWindow;
                        W(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan) = ...
                            W(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan) + self.weightingWindow.^2;
                    end
            end
            s_temp = (s_temp./W);      
            self.s = s_temp(1:oldLength,:);

            if nargout
                s = self.s;
            end
        end
        
        function s = unsplit(self,tempSignal)
            tempLength = self.framesPositions(end)+self.nfft - 1;     
            oldLength = self.sLength;       
            
            s_temp=zeros(tempLength,self.nChans);
            W = zeros(tempLength,self.nChans);

            for index_frame=1:self.nFrames
                    for index_chan=1:self.nChans
                        s_temp(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan)= ...
                            s_temp(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan) + tempSignal(:,index_frame,index_chan).*self.weightingWindow;
                        W(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan) = ...
                            W(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan) + self.weightingWindow.^2;
                    end
            end
            s_temp = (s_temp./W);            
            s = s_temp(1:oldLength,:);
        end
        
        function [sWin, sWeights] = split(self,store,initialPos)
            %default value
            if (nargin <= 1)
                store = 1;
            end
            if nargin <= 2
                initialPos = 1;
            end
            
            %computing splitting parameters
            step = self.nfft-self.overlap;
            self.framesPositions = initialPos:step:self.sLength;
            self.framesPositions((self.framesPositions + self.nfft -1) > self.sLength ) = [];
            self.framesPositions(end+1) = self.framesPositions(end) + step;
            self.nFrames = length(self.framesPositions);
            
            %initializing splitted signal matrix
            sWin = zeros(self.nfft,self.nFrames, self.nChans);
            
            %Initializing weighting window
            win = self.weightingWindow(:)*ones(1,self.nChans);

            %Initializing weighting signal
            tempLength = self.framesPositions(end)+self.nfft - 1;     
            sWeights = zeros(tempLength,1);
            
            %For each frame, compute signal
            for index=1:(self.nFrames - 1)
                sWin(:,index,:) = self.s(self.framesPositions(index):self.framesPositions(index)+self.nfft-1,:).*win;
                sWeights(self.framesPositions(index):(self.framesPositions(index)+self.nfft-1)) = ...
                            sWeights(self.framesPositions(index):(self.framesPositions(index)+self.nfft-1)) + self.weightingWindow;                
            end
            
            %Handle last frame on which zeropadding may be necessary
            lEnd = (self.framesPositions(end) + self.nfft - 1) - self.sLength;
            finalZeroPadding = zeros(lEnd, self.nChans);
            sWin(:,end,:) = [self.s(self.framesPositions(end):end,:) ; finalZeroPadding].*win;
            
            %handle last frame for weighting signal
            sWeights(self.framesPositions(end):(self.framesPositions(end)+self.nfft-1)) = ...
                        sWeights(self.framesPositions(end):(self.framesPositions(end)+self.nfft-1)) + self.weightingWindow;                

            %stores result if necessary
            if store 
                self.sWin = sWin;
                self.sWeights = sWeights;
            end
            if nargout == 0
                clear sWin
                clear sWeights
            end
                
        end
        
        function Scomp = buildComplete(self, S)
            %builds complete spectrogram/STFT/mask using a version up to nfft/2
            Scomp = zeros(size(self.S));
            for c = 1:self.nChans
                Scomp(1:self.nfftUtil,:,c) = S(:,:,c);
                Scomp(end:-1:(end-self.nfftUtil+2),:,c) = conj(S(2:end,:,c));
            end
        end
        
        %copy
        function copied = copy(self)
            copied = Signal;
            copied.suspendListeners
            Meta = ?Signal;
            for index = 1:length(Meta.Properties)
                if strcmp(Meta.Properties{index}.Name, 'propertiesListeners')
                    continue
                end
                eval(['copied.' (Meta.Properties{index}.Name) '= eval([''self.'' (Meta.Properties{index}.Name)]);']);
            end
            copied.activateListeners
        end
        
        function stripCanal(self, channel)
            self.suspendListeners;
            toKeep = setdiff(1:self.nChans,channel);
            self.s = self.s(:,toKeep);
            if ~isempty(self.S)
                self.S = self.S(:,:,toKeep);
            end
            if ~isempty(self.Sq)
                self.Sq = self.Sq(:,:,toKeep);
            end
            if ~isempty(self.sWin)
                self.sWin = self.sWin(:,:,toKeep);
            end
            self.nChans = length(toKeep);
            self.activateListeners;
        end
    end
        
    properties (Access=public, SetObservable=true)
        %Waveform, sampling frequency
        s, fs, singlePrecision
        
        %weighting parameters
        weightingFunction
        
        %STFT parameters and data
        windowLength, nfft, overlapRatio
        S
        

    end
    
    properties (SetAccess = private, GetAccess = public, SetObservable=false)
        %waveform properties
        sLength, nChans, nfftUtil

        %STFT parameters and data
        framesPositions, nFrames
        
        %Windowed data
        sWin, sWeights

    end
    properties (Access = private)
        %Properties listeners
        propertiesListeners
        
        %STFT 
        weightingWindow, overlap
        

    end
    methods (Static)
        function handlePropertyEvents(self,src,evnt)
            switch src.Name 
                case 'singlePrecision'
                    self.S = single(self.S);
                    self.Sq = single(self.Sq);
                case 's' 
                    [self.sLength, self.nChans] = size(self.s);
                case 'fs' 
                    self.nfft = round(self.fs*self.windowLength/1000);
                    self.nfftUtil = round(self.nfft/2);
                    self.overlap = round(self.nfft*self.overlapRatio); 
                case 'windowLength' 
                    self.nfft = round(self.fs*self.windowLength/1000);
                    self.nfftUtil = round(self.nfft/2);
                    self.overlap = round(self.nfft*self.overlapRatio); 
                case 'nfft' 
                    self.windowLength = (self.nfft*1000/self.fs);
                    self.overlap = round(self.nfft*self.overlapRatio); 
                case 'overlapRatio'
                    self.overlap = round(self.nfft*self.overlapRatio); 
                case 'S'
                    self.nFrames = size(self.S,2);
                    self.nfft = size(self.S,1);
                    self.nfftUtil = round(self.nfft/2);
                    self.nChans = size(self.S,3);
            end
        
            if self.nfft
                nArgWeightingFunc = nargin(self.weightingFunction);
                if (nArgWeightingFunc<0)||(nArgWeightingFunc==1)
                    self.weightingWindow = feval(self.weightingFunction, self.nfft);
                elseif nArgWeightingFunc==2
                    self.weightingWindow = feval(self.weightingFunction, self.nfft,self.overlapRatio);
                else
                    error('weighting function badly defined.');
                end
            end
        end  
        function w = hann_nonzero(N)
            n = 0:(N-1);
            phi0 = pi / N;
            Delta = 2 * phi0;
            w = 0.5 * (1 - cos( n*Delta + phi0));
        end
    end
    methods (Access = private)
        function suspendListeners(self)
            for index = 1:length(self.propertiesListeners)
                self.propertiesListeners{index}.Enabled = false;
            end
        end
        function activateListeners(self)
            for index = 1:length(self.propertiesListeners)
                self.propertiesListeners{index}.Enabled = true;
            end
        end     
    end
end
