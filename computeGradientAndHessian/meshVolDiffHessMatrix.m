function [cDf cHf] = meshVolDiffHessMatrix(geom,wghFcn)
%% HELP
%		[cDf cHf] = meshVolDiffHessMatrix(geom,wghFcn)
%			This function calculates the partial derivatives of first and
%			second order (Df and Hf) for the canonical basis of functions
%			on the geometry.
%			This means that it calculates the partial derivtives for all
%			possible combinations of one node set to 1 and the rest to 0
%			for a given geometry.
%
%		INPUT:
%			-  geom - struct - struct containing the geometry where to
%			measure the gradient and hessian. The structure of it is:
%			geom.node for the nodes and geom.face for the faces.
%			- wghFcn - function handle -  returns the squared root of the
%			weightening matrix diagonal to use for the reference node. The function 
%			should be called as:
%				w = wghFcn(ref);
%			where ref is the reference node being used and w is the
%			corresponding weightening matrix diagonal vector.
%
%		OUTPUT:
%			- cDf - <N,dim*N>double - This is the first order derivative values for
%			the canonical basis of the function. Each row has the gradients
%			of each node stacked. The input function equals 1 for the
%			node with the same index and 0 for the rest.
%			- cHf - <N,N*halfVecDim>double - This is the second order derivative values for
%			the canonical basis of the function. Each row has the gradients
%			of each node stacked. The input function equals 1 for the
%			node with the same index and 0 for the rest.
%
%		PROCESS:
%			- for all nodes
%				- set function value of current node to 1, the rest to 0.
%				- calculate the first and second order derivatives.
%				- stack the result.
%
%		DEPENDENCES:
%			meshVolHessian.m
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
	
		[N dim] = size(geom.node);
		
		%prealocate variables
		cDf = zeros(N*dim,N);
		halfVecDim = sum(sum(tril(ones(dim,dim))));
		cHf = zeros(N*halfVecDim,N);
	
	%% FOR ALL NODES
	for ii = 1:N
		
		%% set function value of current node to 1, the rest to 0.
			fcn = zeros(N,1);
			fcn(ii) = 1;
		
		%% calculate the first and second order derivatives.
			 [Df Hf] = meshVolHessian(geom, fcn, wghFcn);
		
		%% stack the result
			% OJO! the matrices in Df and Hf are being vectorized by (:)
			cDf(:,ii) = Df(:)';
			cHf(:,ii) = Hf(:)';

	end
end