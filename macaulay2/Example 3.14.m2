--The diagonal matrix diag(1,2,3) is not a product of two Toeplitz matrices

--Settings
n=3;--the size of matrices
m=2*n-1; --the number of independent variables for n by n Toeplitz matrices
R=QQ[x_(0,0)..x_(n-1,n-1),t_0..t_(m-1),s_0..s_(m-1)];

--Two general Teoplitz matrices and their multiplication
A=matrix table(n,n,(i,j) -> t_(j-i+n-1));
B=matrix table(n,n,(i,j) -> s_(j-i+n-1));
C=A*B;

--The diagonal matrix
D=diagonalMatrix{1,2,3};

--Parametrization
P=D-C;
l={};
for i from 0 to n-1 do(
    for j from 0 to n-1 do(
	l=append(l,P_(i,j));
	)
    )

I=ideal(l);

--Check the existence of solution
time Gb = gens gb I

--Result: (1), so there is no solution.
