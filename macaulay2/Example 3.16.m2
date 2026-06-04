--The upper 3-diagonal matrix U_3 is not a product of two upper bidiagonal matrices

--Settings
R=QQ[x_(0,0)..x_(4,4),t_(0,0)..t_(4,4),s_(0,0)..s_(4,4)];

--Two general tridiagonal matrices and their multiplication
A=matrix{{t_(0,0),t_(0,1),0,0,0},{0,t_(1,1),t_(1,2),0,0},{0,0,t_(2,2),t_(2,3),0},{0,0,0,t_(3,3),t_(3,4)},{0,0,0,0,t_(4,4)}};
B=matrix{{s_(0,0),s_(0,1),0,0,0},{0,s_(1,1),s_(1,2),0,0},{0,0,s_(2,2),s_(2,3),0},{0,0,0,s_(3,3),s_(3,4)},{0,0,0,0,s_(4,4)}};
C=A*B;

--The diagonal matrix
U3=matrix{{1,0,1,0,0},{0,1,0,0,0},{0,0,1,0,1},{0,0,0,1,0},{0,0,0,0,1}};

--Parametrization
P=U3-C;
l={};
for i from 0 to 4 do(
    for j from 0 to 4 do(
	l=append(l,P_(i,j));
	)
    )

I=ideal(l);

--Check the existence of solution
Gb = gens gb I

--Result: (1), so there is no solution.
