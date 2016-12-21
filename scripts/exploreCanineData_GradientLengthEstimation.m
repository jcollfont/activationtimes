%% HELP
%
%		This script loads data from a canine dataset and estimates
%		activation times and recovery times directly from the measured EGM
%
%

clc;clear all;close all;
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
	
	% Gradient estimator
	pathLength = [2, 3, 4, 5, 6, 7];	% number of jumps on graph to determine neighborhood
	
	% activation times estimator
	windowQRS = 11;
	windowT = 19;
	
	% inverses
	vec_lambda = 10.^linspace(-8,-4,1000);
	SNR = 30;


%% LOAD DATA
	load data_bspiral/TWA/realData/rsm-12-10-03/mat/Run0234-sock.mat;
	potvals = ts.potvals;
	load data_bspiral/geometry_modification/inputGeometry/rsm10-27-2014_sock.mat;
	
	load data_bspiral/TWA/geometries/utahtankandsockgeometries/cappedheart.mat
	load data_bspiral/TWA/geometries/utahtankandsockgeometries/tankclosed.mat
	body = struct('surface',[]);
	body.surface{1} = struct('pts',tankclosed.node,'fac',tankclosed.face,'sigma', [0 1]);
	body.surface{2} = struct('pts',cappedheart.node,'fac',cappedheart.face,'sigma', [1 0]);
	A = real( bemMatrixPP2(body) );
	A = A(1:192,:);
	
	M = size(potvals,1);
	Mh = max(size(heart.node));
	
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
		
%% COMPUTE BSPM
	QRSbspm = cell(1,NBQRS);
	Tbspm = cell(1,NBT);
	for bb = 1:NBQRS
		QRSbspm{bb} = awgn(A(:,1:247)*QRSseg{bb},SNR,'measured');
	end
	for bb = 1:NBT
		Tbspm{bb} = awgn(A(:,1:247)*Tseg{bb},SNR,'measured');
	end
	
%% Compute Inverse Solutions
	QRSinv = cell(1,NBQRS);
	Tinv = cell(1,NBT);
	for bb = 1:NBQRS
		QRSinv{bb} = tikhonov_jcf(A,eye(Mh),[],QRSbspm{bb},vec_lambda, false);
	end
	for bb = 1:NBT
		Tinv{bb} = tikhonov_jcf(A,eye(Mh),[],Tbspm{bb},vec_lambda, false);
	end
		
%% DETERMINE GRADIENT MATRICES
	AdjMtrx = cell(1,numel(pathLength));
	D = cell(1,numel(pathLength));
	Ltan = cell(1,numel(pathLength));
	for pp = 1:numel(pathLength)
		
		% compute adjacency matrix
		[AdjMtrx{pp}] = computeAdjacencyMatrix(heart, pathLength(pp));

		% compute gradient estimator
		wghFcn = @(indx) AdjMtrx{pp}(indx,:);
		[D{pp}, H] = meshVolDiffHessMatrix(heart,wghFcn);	
		Ltan{pp}=LaplacianMatrixFromHessianMatrix(H);
		
	end
	
%% COMPUTE ACTIVATION GRADIENTS
	activationGradient = cell(numel(pathLength),NBQRS);
	recoveryGradient = cell(numel(pathLength),NBT);
	activationGradientNorm = cell(numel(pathLength),NBQRS);
	recoveryGradientNorm = cell(numel(pathLength),NBT);
	for pp = 1:numel(pathLength)
		for bb = 1:NBQRS
			[activationGradientNorm{pp,bb}, activationGradient{pp,bb}] = gradientEstimation(QRSinv{bb},D{pp});
		end
		for bb = 1:NBT
			[recoveryGradientNorm{pp,bb}, recoveryGradient{pp,bb}] = gradientEstimation(Tinv{bb},D{pp});
		end
	end
	
%% Valiadate gradient estimation
	SS = [eye(Mh), eye(Mh), eye(Mh)];
	
	activationAngleVal = cell([numel(pathLength),numel(pathLength),Mh]);
	for p1 = 1:numel(pathLength)
		for p2 = 1:numel(pathLength)
			NN = kron([1,1,1], AdjMtrx{p2});
			tmp =  cell2mat(activationGradient(p1,:));
			L = size(tmp,2);
			tmp = permute(reshape(tmp',[L,Mh,3]),[3,2,1]);
			
			for m = 1:Mh
				neighbors = find(NN(:,m));
				tmp2 = tmp(:,neighbors,:);
				
				ix = find(tril(ones(numel(neighbors))));
				tmp4 = tmp(:,neighbors,:);
				activationAngleVal{p1,p2,m} = 0;
				for l = 1:L
					tmp3 = squeeze(tmp4(:,:,l)'*tmp4(:,:,l));
					activationAngleVal{p1,p2,m} = activationAngleVal{p1,p2,m} + mean(tmp3(ix))/(L);
				end
			end
		end
	end
	
%% PLOT
	averageActivationAngleVal = zeros(numel(pathLength));
	for p1 = 1:numel(pathLength)
		for p2 = 1:numel(pathLength)
			averageActivationAngleVal(p1,p2) = mean(cell2mat(activationAngleVal(p1,p2,:)));
		end
	end
	figure;
	plot(averageActivationAngleVal);
	xlabel('Neighborhood for D');
	ylabel('Average angle deviation');








%% MAP3d
	geomFile = sprintf(' data_bspiral/geometry_modification/inputGeometry/rsm10-27-2014_sock.mat');
	geomCommand = {geomFile, geomFile};
	
% 	recoveryTimes = [recoveryTimes; max(recoveryTimes(:))*ones(337-M,NBT)];
	save('tmp/heart2.mat','recoveryTimes');
% 	activationTimes = [activationTimes; max(activationTimes(:))*ones(337-M, NBQRS)];
	save('tmp/heart1.mat','activationTimes');
	potentialCommand = { sprintf(' tmp/heart%d.mat',1),sprintf(' tmp/heart%d.mat',2) };

	positionCoords =    {sprintf(' -as	%d  %d  %d	%d',1 , 270, 10, 300 ), ...
                          sprintf(' -as	%d  %d  %d	%d',1 , 270, 311, 600 ) };
					  
	optionsCommand = {strcat(positionCoords{1},'  -bg 255 255 255 -fg 0 0 0 -sc 1 0 0 -sm 1 -rm 0 -el 1 '), ...
                      strcat(positionCoords{2},'  -bg 255 255 255 -fg 0 0 0 -sc 1 0 0 -sm 1 -rm 0 -el 1 ')...
                              };
	
	plot_map3d(	geomCommand, potentialCommand, optionsCommand);
	