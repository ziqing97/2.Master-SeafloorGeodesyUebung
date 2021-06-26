x = [1,2,3,4,5,6,7];
A = x' * x;
% A = ones(7)*5.0;
[Q,R] = qr(A);
eig(A)