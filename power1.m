function [lamda,eigenvector]=power1(A,X,eps,max1)
% function to compute the maximum eigenvalue (lambda) of a Matrix A and
% its corresponding eigenvector of unit length (eigenvector)
% X is some base vector row matrix, input
% eps is the tolerance you want the eigenvalue to be computed to
% max1 is the maximum number of iteractions allowed
lamda=0;
cnt=1;
xnew=X;
xold=0*X;
err=norm(X);
xnew = xnew/max(abs(xnew));
while err>eps & cnt < max1
    xold=xnew;
    ck = max(abs(xnew));
    xnew=A*xnew/ck;
    cnt = cnt+1;
    err = norm(xold-xnew);
end
% if (cnt >=max1)
%     error('error in power2, max number of iterations exceeded')
% end
eigenvector = xnew/norm(xnew);
lamda = ck;