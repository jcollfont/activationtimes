function [Df Hf] = meshVolHessian(geom, fcn, wghFcn)
%% HELP:
%		[Df Hf] = meshVolHessian(geom, fcn, wghFcn)
%			This function returns the approximated gradient and hessian
%			in each node of a geometry.
%
%		INPUT:
%			- geom - struct - struct containing the geometry where to
%			measure the gradient and hessian. The structure of it is:
%			geom.node for the nodes and geom.face for the faces.
%			- fcn - <N,1>double - value of the function at each node. 
%			the order must be the same as the nodes in geom.
%			- wghFcn - function handle - returns the squared root of the
%			weightening matrix diagonal to use for the reference node. The function 
%			should be called as:
%				w = wghFcn(ref);
%			where ref is the reference node being used and w is the
%			corresponding weightening matrix diagonal vector.
%
%		OUTPUT:
%			- Df - <N,dim>double - value of the gradient at each node. dim
%			is the dimension of the space.
%			- Hf - <N,halfVecDim>double - value of the half vectorized hessian at
%			each node.
%
%		PROCESS:
%			- for all nodes
%				- calculate weights.
%				- search for neighbors.
%				- set current node as origin; build P.
%				- build matrix P2.
%				- build matrix PP = [P P2];
%				- build vector df = (fi - f0);
%				- solve least squares.
%
%		DEPENDENCES:
%
%		AUTHOR:
%			Jaume Coll-font <jcollfont@ece.neu.edu>
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
	
		epsilon = 0; % minimum weight (but not including) to be considered neighbor. it is expected to be 0.
		[N dim] = size(geom.node);
		halfVecDim = sum(sum(tril(ones(dim,dim))));
		Df = zeros(N,dim);
		Hf = zeros(N, halfVecDim);
	
	%% FOR ALL NODES THAT MATTER
	matterIX = find(wghFcn(find(fcn)));
	for ref = matterIX
		
		%% CALCULATE WEIGHTS
			w = wghFcn(ref);
		
		%% SEARCH FOR NEIGHBORS
			neigh = false(1,N);
			neigh(w > epsilon) = true;
			neigh(ref) = false;
			NN = numel(find(neigh));
			
			% select the corresponding weight only for the neighbors
			W = diag(w(neigh));
		
		%% SET CURRENT NODE AS ORIGIN AND BUILD P
			P = geom.node(neigh,:);
			P = P - repmat(geom.node(ref,:),NN,1);
		
		%% BUILD P2
			P2=zeros(NN,halfVecDim);
			for n=1:NN % for all neighbors
				temp=P(n,:)'*P(n,:); % this is the quantity we want, but to enforce symmetry we use some tricks to merge redundant unknowns
				temp=2*temp-diag(diag(temp)); % multiply off-diagonal terms by 2, keep diags as-is
				temp=temp(tril(ones(dim,dim))>0); % keep only the lower triangular part of the matrix %%% half vectorization
				P2(n,:)=temp(:).';%%% fi de la half vectorization
			end
		
		%% BUILD PP =[P P2]
			PP=[P,0.5*P2];
		
		%% BUILD DF VECTOR df = (fi - f0)
			df = fcn(neigh) - fcn(ref);
		
		%% SOLVE LEAST  df = PP*[Df Hf]'
			temp=pinv(W*PP)*W*df(:);
			Df(ref,:)=temp(1:dim);
			Hf(ref,:)=temp(dim+1:end);
	
	end
	
end