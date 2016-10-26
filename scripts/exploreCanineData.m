%% HELP
%
%		This script loads data from a canine dataset and estimates
%		activation times and recovery times directly from the measured EGM
%
%


%% PARAMS
	% QRS segmentation params
	startQRS = -30;
	lengthQRS = 60;
	
	% T-wave segmentation params
	startT = 80; % start recording T-wave wrt Qwave
	lengthT = 100; % end recording Twave wrt Qwave
	
	% baseline removal params
	prevMargin = 80; % select before Qwave
	intLength = 10; % averagie interval to take reference


%% LOAD DATA
	load data_bspiral/TWA/realData/rsm-12-10-03/mat/Run0234-sock.mat
	potvals = ts.potvals;
	
	M = size(potvals,1);
	
%% CLEAN DATA
	% remove badleads
	potvals(130,:) = zeros(1,size(potvals,2));
	potvals(247,:) = zeros(1,size(potvals,2));
	
	% baseline removal
	[filt_potvals, baseline, refPoints, refIntervals] = baselineCorrection_Splines_auto(potvals, prevMargin, intLength);
	

%% SEGMENT HEARTBEATS
		[QRSseg, qrs_peaks] = segmentHeartBeats_bySections_fixedLength(filt_potvals, startQRS, lengthQRS, 1);
		[Tseg, qrs_peaks] = segmentHeartBeats_bySections_fixedLength(filt_potvals, startT, lengthT, 1);
		NBT = numel(Tseg);
		NBQRS = numel(QRSseg);

%% ESTIMATE ACTIVATION TIMES
	activationTimes = zeros(M,NBQRS);
	for bb = 1:NBQRS
		[~, dy] = findMinDVDT(QRSseg{bb}, 11, 2);
		[~, activationTimes(:,bb)] = min(dy,[],2);
	end
	
%% ESTIMATE RCOVERY TIMES
	recoveryTimes = zeros(M,NBT);
	for bb = 1:NBT
		[~, dy] = findMinDVDT(Tseg{bb}, 19, 2);
		[~, recoveryTimes(:,bb)] = max(dy(:,20:end),[],2);
		recoveryTimes(:,bb) = recoveryTimes(:,bb) +20;
	end

%% PLOT
	
	geomFile = sprintf(' data_bspiral/geometry_modification/inputGeometry/rsm10-27-2014_sock.mat');
	geomCommand = {geomFile, geomFile};
	
	recoveryTimes = [recoveryTimes; max(recoveryTimes(:))*ones(337-M,NBT)];
	save('tmp/heart2.mat','recoveryTimes');
	activationTimes = [activationTimes; max(activationTimes(:))*ones(337-M, NBQRS)];
	save('tmp/heart1.mat','activationTimes');
	potentialCommand = { sprintf(' tmp/heart%d.mat',1),sprintf(' tmp/heart%d.mat',2) };

	positionCoords =    {sprintf(' -as	%d  %d  %d	%d',1 , 270, 10, 300 ), ...
                          sprintf(' -as	%d  %d  %d	%d',1 , 270, 311, 600 ) };
					  
	optionsCommand = {strcat(positionCoords{1},'  -bg 255 255 255 -fg 0 0 0 -sc 1 0 0 -sm 1 -rm 0 -el 1 '), ...
                      strcat(positionCoords{2},'  -bg 255 255 255 -fg 0 0 0 -sc 1 0 0 -sm 1 -rm 0 -el 1 ')...
                              };
	
	plot_map3d(	geomCommand, potentialCommand, optionsCommand);
	