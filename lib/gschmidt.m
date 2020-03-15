function [Q,R]=gschmidt(A)
% Input: V is an m by n matrix of full rank m<=n
% Output: an m-by-n upper triangular matrix R
% and an m-by-m unitary matrix Q so that A = Q*R.
[m,n]=size(A);
Q=zeros(m,n);
R=zeros(n,n);

for j=1:n
    v=A(:,j);
    for i=1:j-1
        R(i,j)=Q(:,i)'*A(:,j);
        v=v-R(i,j)*Q(:,i);    
    end
    R(j,j)=norm(v);
    Q(:,j)=v/R(j,j);
end

end
