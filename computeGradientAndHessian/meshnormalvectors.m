function normalvecs = meshnormalvectors(geom)
% calculates area weighted normal vectors to each node in a surface
% triangle mesh
% Author: Burak Erem

% Orient the data structures in 'geom' correctly:
if(size(geom.node,1)<size(geom.node,2))
    geom.node=geom.node.';
end
if(size(geom.face,1)<size(geom.face,2))
    geom.face=geom.face.';
end

dim=size(geom.node,2);
normalvecs=zeros(size(geom.node));

for i=1:size(geom.node,1)
    % get all the nodes that share triangles with node i
    [m,n]=find(geom.face==i);
    neighbface=geom.face(m,:);
    
    % calculate triangle normals weighted by areas (<ABC)
    for k=1:size(neighbface,1)
        AB=geom.node(neighbface(k,2),:)-geom.node(neighbface(k,1),:);
        AC=geom.node(neighbface(k,3),:)-geom.node(neighbface(k,1),:);
        normalvecs(i,:)=normalvecs(i,:)+0.5*cross(AB,AC);
    end
    
    normalvecs(i,:)=normalvecs(i,:)/norm(normalvecs(i,:)); % normalize
end
