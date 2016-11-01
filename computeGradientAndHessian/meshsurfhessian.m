function [Df1,Tf1,Df2,Tf2] = meshsurfhessian(geom,fcn,varargin)

if(~isempty(varargin))
    normals=varargin{1};
    if(size(normals,1)<size(normals,2))
        normals=normals.';
    end
end

% Orient the data structures in 'geom' correctly:
if(size(geom.node,1)<size(geom.node,2))
    geom.node=geom.node.';
end
if(size(geom.face,1)<size(geom.face,2))
    geom.face=geom.face.';
end

dim=size(geom.node,2);

Vh2V=vech2vec(dim);

Df1=zeros(length(fcn),size(geom.node,2));
Tf1=zeros(length(fcn),size(geom.node,2));
Df2=zeros(length(fcn),sum(sum(tril(ones(dim,dim)))));
Tf2=zeros(length(fcn),sum(sum(tril(ones(dim,dim)))));

for i=1:size(geom.node,1)
    % get all the nodes that share triangles with node i
    [m,n]=find(geom.face==i);
    neighbface=geom.face(m,:);
    neighbface=neighbface(neighbface~=i);
    neighb=unique(neighbface(:));
    
    % form the difference vector matrix
    P=geom.node(neighb,:);
    P=P-repmat(geom.node(i,:),size(P,1),1);
    P2=zeros(size(P,1),sum(sum(tril(ones(dim,dim)))));
    for n=1:size(P,1)
        temp=P(n,:)'*P(n,:); % this is the quantity we want, but to enforce symmetry we use some tricks to merge redundant unknowns
        temp=2*temp-diag(diag(temp)); % multiply off-diagonal terms by 2, keep diags as-is
        temp=temp(tril(ones(dim,dim))>0); % keep only the lower triangular part of the matrix
        P2(n,:)=temp(:).';
    end
    d=fcn(neighb)-fcn(i);
    
    % the projection matrix to the tangent space
    if(~isempty(varargin))
        TanProj=eye(dim)-normals(i,:)'*normals(i,:)/(norm(normals(i,:))^2);
    else
        [U,S,V]=svd(P);
        TanProj=V(:,1:dim-1)*V(:,1:dim-1)';
    end
    PP=[P,0.5*P2];
    
%   We want to solve for the gradient and Hessian subject to the constraint
%   that they reside in the null space of Pdiag
    NorProj=eye(dim)-TanProj;
    Pdiag=blkdiag(NorProj,kron(NorProj,NorProj)*Vh2V);
    B=null(Pdiag); % Let B be a basis for the null space of Pdiag
    b=pinv(PP*B)*d(:);
    temp=B*b;
    Tf1(i,:)=temp(1:dim);
    Tf2(i,:)=temp(dim+1:end);

%     temp=([P,0.5*P2]\d(:));
    temp=pinv(PP)*d(:);
    Df1(i,:)=temp(1:dim);
    Df2(i,:)=temp(dim+1:end);

% % Enforce the constraints on our solution above
%     Tf1(i,:)=(TanProj*Df1(i,:).').';
%     
%     temp=zeros(dim,dim);
%     temp(tril(ones(dim,dim))>0)=Df2(i,:).'; % half vectorization to lower triangular
%     temp=temp+temp.'-diag(diag(temp)); % complete the matrix from its lower triangular form
%     temp=TanProj*temp*TanProj'; % project out any compotent in the normal direction
%     
%     Tf2(i,:)=temp(tril(ones(dim,dim))>0); % again, keep only the lower triangular part, vectorized    
end
