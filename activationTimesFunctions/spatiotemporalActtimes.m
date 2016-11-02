%% HELP:
%
%	This function computes the activation times in a cardiac sequence using
%	a spatiotemporal approach.
%	In this spatiotemporal approach, the objective is to find the minimum
%	dVdt weighted by the norm of the spatial gradient at each time
%	instance.
%
%	INPUT:
%		- X - <N,T>double - sequence of potentials on the heart.
%		- D - <3N,N>double - estimator of spatial gradient of potentials.
%		It should be organized such that the first N columns correspond to
%		the gradient in X direction of each node, the following N columns
%		are Y and the last Z.
%		- window - int - (OPTIONAL) window of time over which the temporal
%		derivative is calculated.
%		- alpha - double [0,1] - weight to be used for the gradient. If 1
%		it is direct multiplication, using 0 implies no weighting.
%
%


function [tau,objfun,DX,dXdt] = spatiotemporalActtimes(X,D,varargin)

	%% parse inputs
		if(size(varargin,2)>0),window=varargin{1};else window=1;end;
% 		if(size(varargin,2)>1),alpha=varargin{2};else alpha=1;end;

		[N,T] = size(X);

	%% Compute gradient
		SS = [eye(N), eye(N), eye(N)];		% add up X,Y,Z entries of gradient
		
		DX = sqrt( SS*(D(:,1:N)*X).^2 );		% compute gradient norm
		DX = DX./repmat(max(DX,[],2),[1,T]);	% normalize to 1 (over time)
		DX(isnan(DX)) = 0;						% eliminate outliers

		
	%% Compute min dVdT
		[~, dXdt] = findMinDVDT( X, window, 2);
	
	%% Weight dVdt with spatial gradient norm
		objfun= DX.*dXdt; % objfun now contains the gradient norms weighted by negative temporal slope

	%% FIND MINIMUM OF THE OBJECTIVE FUNCTION
		[val,tau]=min(objfun,[],2);

end
