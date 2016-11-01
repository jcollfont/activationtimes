function [D,Dtan,H,Htan] = meshsurfdiffhessmatrix(geom,varargin)

if(~isempty(varargin)) % Normals are provided
    normals=varargin{1};
    if(size(normals,1)<size(normals,2)),normals=normals.';end;
else % Precalculate the normals
    % Orient the data structures in 'geom' correctly:
    if(size(geom.node,1)<size(geom.node,2))
        geom.node=geom.node.';
    end
    if(size(geom.face,1)<size(geom.face,2))
        geom.face=geom.face.';
    end
    dim=size(geom.node,2);
    
    normals=zeros(size(geom.node));
    
    for i=1:size(geom.node,1)
        % get all the nodes that share triangles with node i
        [m,n]=find(geom.face==i);
        neighbface=geom.face(m,:);
        neighbface=neighbface(neighbface~=i);
        neighb=unique(neighbface(:));

        % form the difference vector matrix
        P=geom.node(neighb,:);
        P=P-repmat(geom.node(i,:),size(P,1),1);

    % the following approximates the tangent space using the svd
            [U,S,V]=svd(P);
            normals(i,:)=V(:,end)';
    end
end

D=zeros(numel(geom.node),max(size(geom.node)));
Dtan=zeros(numel(geom.node),max(size(geom.node)));

for i=1:max(size(geom.node))
    fcn=zeros(max(size(geom.node)),1);
    fcn(i)=1;
    if(~isempty(varargin))
        [Df1,Tf1,Df2,Tf2]=meshsurfhessian(geom,fcn,normals);
    else
        [Df1,Tf1,Df2,Tf2]=meshsurfhessian(geom,fcn);
    end
    
    D(:,i)=Df1(:);
    Dtan(:,i)=Tf1(:);
    H(:,i)=Df2(:);
    Htan(:,i)=Tf2(:);
end
