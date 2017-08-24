function Dupmatrix = vech2vec(n)
% Creates a duplication matrix that converts a half-vectorization of a 
% (n x n) symmetric matrix to full vectorization by left multiplication

N=n*n;
Nhalf=sum(reshape(tril(ones(n,n)),N,1)); % keep only the lower triangular part of the matrix

Dupmatrix=zeros(N,Nhalf);

for i=1:Nhalf
    temp=zeros(Nhalf,1);
    temp2=zeros(n,n);
    
    temp(i)=1;
    
    temp2(tril(ones(n,n))>0)=temp; % half vectorization to lower triangular
    temp2=temp2+temp2.'-diag(diag(temp2)); % complete the matrix from its lower triangular form
    
    Dupmatrix(:,i)=temp2(:);
end
