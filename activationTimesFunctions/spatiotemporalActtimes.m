function [tau,objfun,DX,dXdt] = spatiotemporalActtimes(X,D,varargin)
% Each column of X is a vector of spatial potentials
% Each row of X are the samples of a temporal waveform of a node in space
% The matrix D computes the spatial gradient and is organized as
% D=[Dx;Dy;Dz] in block matrix form

if(size(varargin,2)>0),window=varargin{1};else window=1;end;

N=size(X,1);

DX=D*X;
% DX=abs(DX(1:N,1:end-1)).^2+abs(DX(N+1:2*N,1:end-1)).^2+abs(DX(2*N+1:3*N,1:end-1)).^2;
if(window==1)
    dXdt=diff(X,1,2);
else
    for i=1:size(X,1)
        dXdt(i,:)=polydiff(X(i,:),window,2);
    end
    dXdt=dXdt(:,1:end-1);
end
% DX=sqrt(DX); % DX now contains the gradient norms over time

% objfun=DX.*dXdt; % objfun now contains the gradient norms weighted by negative temporal slope
objfun=dXdt;
% the times with the greatest spatial gradient weighted by negative
% temporal slope are the activation times we choose
[val,tau]=min(objfun,[],2);

% need to shift back the taus
if(window==1)
tau=tau-1;
end

end

function dy = polydiff(sig,win,deg)
cen = ceil(win/2);
X = zeros(win,(deg+1));
L = [-(cen-1):(cen-1)]'; for p=1:(deg+1), X(:,p) = L.^((deg+1)-p); end
E = inv(X'*X)*X';
sig = [sig sig(end)*ones(1,cen-1)];
a = filter(E(deg,[win:-1:1]),1,sig);
dy = a(cen:end);
end