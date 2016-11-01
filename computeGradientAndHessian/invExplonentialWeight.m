function [wghFcn] = invExplonentialWeight(geom, sigma, trunc)
%% HELP:
%		[wghFcn] = invExplonentialWeight(geom, sigma, trunc)
%			This function returns a function handle that calculates the
%			corresponding squared root of weight of all nodes with respect 
%			to the reference node.
%
%		INPUT:
%			- geom - struct - struct containing the geometry where to
%			measure the gradient and hessian. The structure of it is:
%			geom.node for the nodes and geom.face for the faces.
%			- sigma - double / <N,N>double - this parameter controls the
%			area of effect of the inverse exponential. can be defined as a
%			single value for all nodes or as a matrix that has the sigma
%			relative to all pairwise  points. In the matrix option, the
%			directionality is defined so that each row corresponds to the
%			reference node.
%			- trunc - double - minimum value to be considered part of the
%			support. if below thi threshold, the weight will automatically
%			be 0.
%
%		OUTPUT:
%			- wghFcn - function handle - returns the weightening matrix diagonal to
%			use for the reference node. This weightening function is equal
%			to the squared root of the inverse exponential of the distance between nodes.
%				w(i,j) = exp( -dist(i,j)/sigma)
%			The function should be called as:
%				w = wghFcn(ref);
%			where ref is the reference node being used and w is the
%			corresponding weightening matrix diagonal vector.
%
%		PROCESS:
%			- calculate the distance between all nodes.
%			- calculates the inverse exponential of the distances.
%			- truncates if necessary and sets the outside values to 0.
%			- build function handle.
%
%		DEPENDENCES:
%			-PairwiseDistance.m
%
%		AUTHOR:
%			Jaume Coll-Font <jcollfont@ece.neu.edu>
%
%
	%% DEFINE
	
	% check if the geometry is oriented correctly:
		if(size(geom.node,1)<size(geom.node,2))
			geom.node=geom.node.';
		end
		if(size(geom.face,1)<size(geom.face,2))
			geom.face=geom.face.';
		end
		
	%% calculate the distance between all nodes.
		[D] = PairwiseDistance(geom.node');
		
	%% calculates the squared root inverse exponential of the distances.
		iE = exp( -D./(2*sigma));
		
	%% truncates if necessary and sets the outside values to 0.
		iE(iE < trunc) = 0;
	
	%% build function handle
		wghFcn = @(ref) (iE(ref,:));

end