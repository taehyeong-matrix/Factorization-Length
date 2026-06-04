--The degree of the Toeplitz facotrization variety \mu_2(\mathcal{T}_4)

--Package
needsPackage "NumericalAlgebraicGeometry"

--Settings
kk=CC_100 --the base field; put C_100 for more precise calculations
n=4;--the size of matrices
m=2*n-1; --the number of independent variables for n by n Toeplitz matrices
R=kk[x_(0,0)..x_(n-1,n-1)];
S=kk[t_0..t_(m-1),s_0..s_(m-1)];

--Two general Teoplitz matrices and their multiplication
A=matrix table(n,n,(i,j) -> t_(j-i+n-1));
B=matrix table(n,n,(i,j) -> s_(j-i+n-1));
C=A*B;

--Parametrization
L={};
for i from 0 to n-1 do(
    for j from 0 to n-1 do(
	L=append(L,C_(i,j));
	)
    )

LL={};
for i from 0 to n-1 do(
    for j from 0 to n-1 do(
	LL=append(LL,x_(i,j)=>L_(n*i+j));
    )
)

--Dimension
r=2; --the number of Toeplitz matrices to multiply
D=min(n^2,2*r*(n-1)+1); --D is the dimension of the Toeplitz factorization variety

--Generic affine subspace
v=matrix {{x_(0,0)}..{x_(n-1,n-1)}};
RM=random(kk^D,kk^(n^2));
w=random(kk^D,kk^1); 
RE=RM*v-w;
M= transpose RE; -- the list of generic affine linear equations (the number of the equations is the codimension)

LLLL={};
for i from 0 to D-1 do(
    LLLL=append(LLLL,M_(0,i));
    )

eqs1={};

for i from 0 to D-1 do(
    eqs1=append(eqs1,sub(LLLL_i,LL));
    ) 

--Resolve the system to let the number of solutions be finite
tv = transpose matrix{{t_0..t_(m-1)}};
gauge1 = (random(kk^1,kk^m) * tv)_(0,0) - 1;
system1=append(eqs1,gauge1);

--Degree
time sol=solveSystem(system1); 
length sol 

--Result
-- used 143.467s (cpu); 104.283s (thread); 0s (gc)
-- the number of solutions is 74, i.e. the degree is at least 74.
