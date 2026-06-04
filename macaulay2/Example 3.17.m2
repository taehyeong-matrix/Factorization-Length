--The 5-diagonal matrix S is not a product of two tridiagonal matrices

--Settings
R=QQ[x_(0,0)..x_(3,3),t_(0,0)..t_(3,3),s_(0,0)..s_(3,3)];

--Two general tridiagonal matrices and their multiplication
A=matrix{{t_(0,0),t_(0,1),0,0},{t_(1,0),t_(1,1),t_(1,2),0},{0,t_(2,1),t_(2,2),t_(2,3)},{0,0,t_(3,2),t_(3,3)}};
B=matrix{{s_(0,0),s_(0,1),0,0},{s_(1,0),s_(1,1),s_(1,2),0},{0,s_(2,1),s_(2,2),s_(2,3)},{0,0,s_(3,2),s_(3,3)}};
C=A*B;

--The diagonal matrix
S=matrix{{0,0,1,0},{0,0,0,1},{1,0,0,0},{0,1,0,0}};

--Parametrization
P=S-C;
l={};
for i from 0 to 3 do(
    for j from 0 to 3 do(
	l=append(l,P_(i,j));
	)
    )

I=ideal(l);

--Check the existence of solution
Gb = gens gb I

--Result: (1), so there is no solution.
