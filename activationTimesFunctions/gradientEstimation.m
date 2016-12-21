%% HELP:
%
%	This function computes the gradient of a cardiac sequence.
%
%	INPUT:
%		- X - <N,T>double - sequence of potentials on the heart.
%		- D - <3N,N>double - estimator of spatial gradient of potentials.
%		It should be organized such that the first N columns correspond to
%		the gradient in X direction of each node, the following N columns
%		are Y and the last Z.
%
%


function [DX, gradX] = gradientEstimation(X,D)

	%% parse inputs
		[N,T] = size(X);

	%% Compute gradient
		SS = [eye(N), eye(N), eye(N)];		% add up X,Y,Z entries of gradient
		
		gradX = D(:,1:N)*X;
		
		DX = sqrt( SS*(gradX).^2 );		% compute gradient norm
		DX = DX;%./repmat(max(DX,[],2),[1,T]);	% normalize to 1 (over time)
		DX(isnan(DX)) = 0;						% eliminate outliers

		gradX = gradX./kron( [1;1;1], DX);
		
end
