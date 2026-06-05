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
tally apply(sol, p -> status p)

--Result
-- used 143.467s (cpu); 104.283s (thread); 0s (gc)
-- Tally{Regular => 74}
-- The number of solutions is 74 counted with multiplicity.


--The degree of generic fiber

-- Generic Toeplitz factorization
a = apply(m, i -> random(kk));
b = apply(m, i -> random(kk));
A0 = matrix table(n,n,(i,j) -> a_(j-i+n-1));
B0 = matrix table(n,n,(i,j) -> b_(j-i+n-1));
C0 = A0*B0;

eqs2 = flatten for i from 0 to n-1 list (
         for j from 0 to n-1 list (
             (A*B)_(i,j) - C0_(i,j)
         )
      );

tvv = transpose matrix{{t_0..t_(m-1)}};
gaugeCoeff = random(kk^1,kk^m);
gauge2 = (gaugeCoeff * tvv)_(0,0) - 1;

system2 = append(eqs2, gauge2);

time fsol = solveSystem(system2);
length fsol
tally apply(fsol, p -> status p)

--Result
-- used 468.342s (cpu); 195.721s (thread); 0s (gc)
-- Tally{Regular => 1}
-- The number of solutions is 1 counted with multiplicity.
