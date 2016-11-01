function L = LaplacianMatrixFromHessianMatrix(H)
% Author: Burak Erem
% This function takes a Hessian approximation matrix (e.g. as created by
% meshVolDiffHessMatrix or meshsurfdiffhessmatrix) and creates a Laplacian
% approximation matrix from it
% 

N=size(H,2); % the number of nodes in the geometry on which the Hessian is approximated

L=H(1:N,:)+H(3*N+1:4*N,:)+H(5*N+1:6*N,:);
