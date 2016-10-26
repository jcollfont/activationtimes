%% HELP:
%
%	This function computes the re-entry vulnerability index (RVI) as
%	defined in the paper:
%
%			"An activation-repolarization time metric to predict localized
%			regions of high susceptibility to re-entry", Child et. al.,
%			Heart Rhythm 2015
%
%	The RVI consists on calculating the difference between the
%	repolarization time and the activation times 
%
%	The RVI maps are computed as....
%
%	INPUT:
%		-repolTimes - <M,1>double - set of repolarization times on the
%									heart.
%		-actTimes - <M,1>double - set of activation times on the heart.
%		- geom - struct - heart geometry defined as a triangular mesh.
%
%	OUTPUT:
%		- RVI - <M,maxDeg>double - all pairs of RVI indices defined at each
%								node. The nodes with less neighbors than
%								the maxDeg have entries set to NaN.
%		- RVImax - <M,1>double - RVI map computed as the max of the RVI pairs on a node.
%		- RVImean - <M,1>double - RVI map computed as the mean RVI on each
%								node.
%
%	AUTHOR:
%		Jaume Coll-Font <jcollfont@gmail.com>
%
%

function [RVI, RVImax, RVImean] = computeRVI(repolTimes, actTimes, geom)

	%% define
		[d,M] = size(geom.node);
		[numFac] = max(size(geom.face));

	%% compute adjacency matrix

		if d ~= 3
			M = d;
			geom.face = geom.face';
			geom.node = geom.node';
		end

		% calculate distances on graph
		distance = PairwiseDistance(geom.node);

		% create adjacency matrix
		adjMtrx = zeros(M);
		for ii = 1:numFac;
			adjMtrx(geom.face(1,ii),geom.face(2,ii)) = 1;%distance(geom.face(1,ii),geom.face(2,ii));
			adjMtrx(geom.face(2,ii),geom.face(3,ii)) = 1;%distance(geom.face(2,ii),geom.face(3,ii));
			adjMtrx(geom.face(3,ii),geom.face(1,ii)) = 1;%distance(geom.face(3,ii),geom.face(1,ii));
			adjMtrx(geom.face(2,ii),geom.face(1,ii)) = 1;%distance(geom.face(2,ii),geom.face(1,ii));
			adjMtrx(geom.face(3,ii),geom.face(2,ii)) = 1;%distance(geom.face(3,ii),geom.face(2,ii));
			adjMtrx(geom.face(1,ii),geom.face(3,ii)) = 1;%distance(geom.face(1,ii),geom.face(3,ii));
		end
	
		% add another hop in the neighborhood
		tempAdj = zeros(M);
		for ii = 1:M
			indx = find(adjMtrx(ii,:));
			tempAdj(ii,:)  = adjMtrx(ii,:);
			for jj = 1:numel(indx)
				tempAdj(ii,:) =  tempAdj(ii,:) + adjMtrx(indx(jj),:);
			end
			tempAdj(ii,ii) = 0;
			tempAdj = min(tempAdj,1);
		end
		
	%% compute "downstream" order
		downstreamMtrx = zeros(M);
		
		for ii =1:M
			for jj = (ii+1):M
				downstreamMtrx(ii,jj) = actTimes(ii) >= actTimes(jj);
				downstreamMtrx(jj,ii) = ~downstreamMtrx(ii,jj);
			end
		end
		
		downstreamMtrx = downstreamMtrx.*tempAdj;
		
	
	%% compute all RVI of each node
		maxDeg = max(sum(downstreamMtrx,2));
		
		RVI = zeros(M,maxDeg);
		for ii = 1:M
			indx = find(downstreamMtrx(ii,:));
			for jj = 1:numel(indx)
				RVI(ii,jj) = repolTimes(ii) - actTimes(indx(jj));
			end
		end
	
	%% compute RVI map
		RVImax = max(RVI,[],2);
		RVImean = mean(RVI,2);
	
end