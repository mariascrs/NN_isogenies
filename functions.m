/*
	Useful functions to perform Fast Kummer arithmetic
*/


///////////////////////////////////////////
//////////// BASIC FUNCTIONS //////////////
///////////////////////////////////////////


function Hadamard(P);

    /////////////////////////////////////////////////////////////////
    // The four-way Hadamard transform 
    // INPUTS: a tuple of four projective field elements [x,y,z,t]
    // OUTPUTS: the Hadamard transform 
    //              [x+y+z+t,x+y-z-t,x-y+z-t,x-y-z+t]	
	// 
	// COST: 8a
    /////////////////////////////////////////////////////////////////

	T1:=P[1]+P[2]; T2:=P[3]+P[4]; T3:=P[1]-P[2]; T4:=P[3]-P[4];

	return [T1+T2,T1-T2,T3+T4,T3-T4];

end function;

function Squaring(P);

    /////////////////////////////////////////////////////////////////
    // The four-way squaring
    // INPUTS: a tuple of four projective field elements [x,y,z,t]
    // OUTPUTS: the coordinate-wise squarings [x^2,y^2,z^2,t^2]	
	// 
	// COST: 4 S
    /////////////////////////////////////////////////////////////////


	return [P[1]^2,P[2]^2,P[3]^2,P[4]^2];

end function;

function FourWayMult(V,W);

    /////////////////////////////////////////////////////////////////
    // The four-way coordinate-wise multiplication of 2 tuples
    // INPUTS: two tuples of four projective field elements
    // OUTPUTS: one 4-tuple of their coordinate-wise products
	// 
	// COST: 4M
    /////////////////////////////////////////////////////////////////

	return [V[1]*W[1],V[2]*W[2],V[3]*W[3],V[4]*W[4]];

end function;

function normalise(P);

    /////////////////////////////////////////////////////////////////
    // Normalizes Kummer points in P^3, mainly for equality and
    // correctness checking. Sets the 4th coordinate as 1, unless 
    // there is a zero coordinate, in which case it's a 2-torsion 
    // point, where we set the 1st coordinate as 1.
    // 
    // INPUTS: a point P in P^3, represented as a 4-tuple
    // OUTPUTS: normalised point equivalent to P in P^3
    /////////////////////////////////////////////////////////////////

	if &*P ne 0 then
		return [P[1]/P[4],P[2]/P[4],P[3]/P[4],P[4]/P[4]];
	else 
		if P[4] ne 0 then 
			return [P[1]/P[4],P[2]/P[4],P[3]/P[4],P[4]/P[4]];
		elif P[3] ne 0 then 
			return [P[1]/P[3],P[2]/P[3],P[3]/P[3],P[4]/P[3]];
		elif P[2] ne 0 then
			return [P[1]/P[2],P[2]/P[2],P[3]/P[2],P[4]/P[2]];
		else
			return [P[1]/P[1],P[2]/P[1],P[3]/P[1],P[4]/P[1]];
		end if;
	end if;

end function;

function Invert4Constants(P);

    /////////////////////////////////////////////////////////////////
    // The four-way inversion in projective 3-space
    // INPUTS: one tuple of 4-elements in P^3 
    // OUTPUTS: a tuple of elements projectively equivalent to 
    //              their inverses
	// 
	// COST: 6M
    /////////////////////////////////////////////////////////////////

	pi1:=P[3]*P[4];
	pi2:=pi1*P[1];
	pi1:=pi1*P[2];
	pi3:=P[1]*P[2];
	pi4:=pi3*P[3];
	pi3:=pi3*P[4];

	return [pi1,pi2,pi3,pi4];

end function;

function OnKummer(P,K);

    /////////////////////////////////////////////////////////////////
    // Checks if a point P is on a Fast Kummer surface K
    // INPUTS: one tuple of 4-elements in P^3 
    // OUTPUTS: true, if P is on K, and false otherwise
    /////////////////////////////////////////////////////////////////
 
	x:=P[1]; y:=P[2]; z:=P[3]; t:=P[4];
	E:=K[1][1]; F:=K[1][2]; G:=K[1][3]; H:=K[1][4];

	return ((x^4+y^4+z^4+t^4)+ 2*E*x*y*z*t - F*(x^2*t^2+y^2*z^2)
     -G*(x^2*z^2+y^2*t^2) - H*(x^2*y^2+t^2*z^2)) eq 0;

end function;

///////////////////////////////////////////
/////////// KUMMER ARITHMETIC /////////////
///////////////////////////////////////////


function DoubleKummer(P,K);


    /////////////////////////////////////////////////////////////////
    // (Pseudo-)doubling of a Kummer point
    // Inputs: a point P in P^3, represented as a 4-tuple, lying on
    //         Kummer surface K
    // Outputs: [2]P
    /////////////////////////////////////////////////////////////////

	X:=P[1]; Y:=P[2]; Z:=P[3]; T:=P[4]; 

	P:=Squaring(P);
	P:=Hadamard(P);
	P:=Squaring(P);
	P:=FourWayMult(P,K[4]);
	P:=Hadamard(P);
	P:=FourWayMult(P,K[3]);

	// assert OnKummer(P,K);
	// P := normalise(P);

	return P;

end function;

function DoubleAddKummer(P,Q,R,K);


    /////////////////////////////////////////////////////////////////
    // (Pseudo-)Double-Add of Kummer points
    // INPUTS: P,Q,R: Points in P^3, all represented as a 4-tuple, 
    //         where R=(Q+-P) lying on Kummer surface K
    // Outputs: (Q-+P), 2[P]
    /////////////////////////////////////////////////////////////////
	
	P:=Squaring(P);
	Q:=Squaring(Q);
	P:=Hadamard(P);
	Q:=Hadamard(Q);
	S:=FourWayMult(P,K[4]);
	Q:=FourWayMult(Q,S);
	P:=FourWayMult(P,S);
	P:=Hadamard(P);
	Q:=Hadamard(Q);
	P:=FourWayMult(P,K[3]);
	Q:=FourWayMult(Q,Invert4Constants(R));
	
	// assert OnKummer(P,K);
	// assert OnKummer(Q,K);

	return P,Q;

end function;

function PseudoAddKummer(P,Q,R,K);

    /////////////////////////////////////////////////////////////////
    // (Pseudo-)addition of Kummer points
    // Inputs: P,Q,R: Points in P^3, all represented as a 4-tuple, 
    //         where R=(Q+-P) lying on Kummer surface K
    // Outputs: (Q-+P)
    /////////////////////////////////////////////////////////////////

	P:=Squaring(P);
	Q:=Squaring(Q);
	P:=Hadamard(P);
	Q:=Hadamard(Q);
	S:=FourWayMult(P,K[4]);
	Q:=FourWayMult(Q,S);
	Q:=Hadamard(Q);
	Q:=FourWayMult(Q,Invert4Constants(R));
	
	// assert OnKummer(Q,K);

	return Q;

end function;

function ScalarMultKummer(P,n,K);
    
	/////////////////////////////////////////////////////////////////
    // Scalar multiplcation of Kummer points
    // INPUTS: P: A point in P^3, represented as a 4-tuple, 
	// 		   n: scalar 
	//         K: the Kummer surface P lies on
    //         
    // OUTPUTS: [n]P
    /////////////////////////////////////////////////////////////////

	bits:=Reverse(IntegerToSequence(n,2));
	
	Pm:=P; Pp:=DoubleKummer(P,K); R:=P;

	for i:=2 to #bits do
		if bits[i] eq 1 then
			Pp,Q:=DoubleAddKummer(Pp,Pm,R,K);
			Pm:=Q;
		else
			Pm,Q:=DoubleAddKummer(Pm,Pp,R,K);
			Pp:=Q;
		end if;
	end for;
	
	// assert OnKummer(Pm,K);

	return normalise(Pm);

end function;

function RandomKummerPoint(K);
	
	/////////////////////////////////////////////////////////////////
    // Sample a random Kummer point
    // INPUTS: K: the Kummer surface 
    // OUTPUTS: Random point on P
    /////////////////////////////////////////////////////////////////

	E:=K[1][1]; F:=K[1][2]; G:=K[1][3]; H:=K[1][4];
	FF:=Parent(E);

	repeat
		x := Random(FF);
		y := Random(FF);
		z := Random(FF);
		Poly<T>:=PolynomialRing(FF);
		Eq := (x^4+y^4+z^4+T^4) + 2*E*x*y*z*T - F*(x^2*T^2+y^2*z^2)
			-G*(x^2*z^2+y^2*T^2) - H*(x^2*y^2+T^2*z^2);
		test, t := HasRoot(Eq);
	until test;

	P:=[x,y,z,t];

	assert OnKummer(P,K);

	return P;

end function;

///////////////////////////////////////////
////////////////// SET-UP /////////////////
///////////////////////////////////////////


function random_richelot(C);
	// Chooses a (pseudo-)random genus-2 curve CC such that Jac(CC) is (2,2)-isogenous to Jac(C)
    isos:=RichelotIsogenousSurfaces(C);
	repeat
		CC:=isos[Random(1,#isos)];
	until Type(CC) eq CrvHyp;
	return CC;
end function;

function random_walk(C);
	// For a curve C, outputs a CC such that Jac(CC) is a random number of steps away from Jac(C) in the (2,2)-isogeny graph
	steps:=Random(5,12);
	CC:=C;
	for i:=1 to steps do
		CC:=random_richelot(CC);
	end for;
	return CC;
end function;

function all_rosenhains(C);
	// Computes all possible Rosenhains corresponding to the genus-2 curve C

    f:=HyperellipticPolynomials(C);
    poly<x>:=Parent(f);
    assert Degree(f) eq 6;
    rts:=Roots(f); rts:=[rts[i][1]: i in [1..6]];
    set:={};

    for r1 in rts do 
        others1:=Remove(rts,Index(rts,r1));
        for r2 in others1 do
            others2:=Remove(others1,Index(others1,r2));
            for r3 in others2 do
                others3:=Remove(others2,Index(others2,r3));
                g:=(x-r1)/(x-r2)*(r3-r2)/(r3-r1);
                for r4 in others3 do
                    others4:=Remove(others3,Index(others3,r4));
                    lambda:=Evaluate(g,r4);
                    for r5 in others4 do
                        mu:=Evaluate(g,r5);
                        others5:=Remove(others4,Index(others4,r5));
                        nu:=Evaluate(g,others5[1]);
                        Include(~set,[lambda,mu,nu]);
                    end for;
                end for;
            end for;
        end for;
    end for;

    return set;

end function;

function one_squared_thetas(C);	
	// Computes the sqaured theta constants associated to a genus-2 curve C

    rosen:=all_rosenhains(C);

    poly<x>:=Parent(HyperellipticPolynomials(C));

    thetas:={};

    repeat 

		lmn:=Random(rosen);
        lambda:=lmn[1]; mu:=lmn[2]; nu:=lmn[3];
        CC:=HyperellipticCurve(x*(x-1)*(x-lambda)*(x-mu)*(x-nu));
        assert IsSquare(lambda*mu/nu);
        assert IsSquare(mu*(mu-1)*(lambda-nu)/(nu*(nu-1)*(lambda-mu)));

        d2:=1;
		c2sign:=Random([-1,1]);
		b2sign:=Random([-1,1]);
		c2:=c2sign*Sqrt(lambda*mu/nu);
        b2:=b2sign*Sqrt(mu*(mu-1)*(lambda-nu)/(nu*(nu-1)*(lambda-mu)));
        a2:=b2*c2*nu/mu;

    until IsSquare(a2) and IsSquare(b2) and IsSquare(c2);

    return [a2,b2,c2,d2],[lambda,mu,nu];

end function;

function KummerFromCurve(C);

	// Given a genus-2 curve C, outputs the Fast Kummer surface associated to it

	sq_thetas,rosen:=one_squared_thetas(C);
	
	a:=Sqrt(sq_thetas[1]);
	b:=Sqrt(sq_thetas[2]);
	c:=Sqrt(sq_thetas[3]);
	d:=Sqrt(sq_thetas[4]);

	E := a*b*c*d*(a^2+b^2+c^2+d^2)*(d^2-a^2+b^2-c^2)*(d^2-a^2-b^2+c^2)*(d^2+a^2-b^2-c^2)/(a^2*d^2-b^2*c^2)/(b^2*d^2-a^2*c^2)/(c^2*d^2-a^2*b^2);
	F:=(a^4+d^4 - b^4-c^4)/(a^2*d^2-b^2*c^2);
	G:=(b^4+d^4 - a^4-c^4)/(b^2*d^2-a^2*c^2);
	H:=(c^4+d^4 - a^4-b^4)/(c^2*d^2-a^2*b^2);

	K1:=[E,F,G,H];
	K2:=[a,b,c,d];
	K3:=Invert4Constants(K2);
	K4:=Invert4Constants(Hadamard(sq_thetas));

	K:=[K1,K2,K3,K4];

	return K,rosen;

end function;

function KummerFromThetas(thetas);

	// Given fundamental theta constants 'thetas', outputs the Fast Kummer surface corresponding to it

	a,b,c,d:=Explode(thetas);

	E := a*b*c*d*(a^2+b^2+c^2+d^2)*(d^2-a^2+b^2-c^2)*(d^2-a^2-b^2+c^2)*(d^2+a^2-b^2-c^2)/(a^2*d^2-b^2*c^2)/(b^2*d^2-a^2*c^2)/(c^2*d^2-a^2*b^2);
	F:=(a^4+d^4 - b^4-c^4)/(a^2*d^2-b^2*c^2);
	G:=(b^4+d^4 - a^4-c^4)/(b^2*d^2-a^2*c^2);
	H:=(c^4+d^4 - a^4-b^4)/(c^2*d^2-a^2*b^2);

	K1:=[E,F,G,H];
	K2:=[a,b,c,d];
	K3:=Invert4Constants(K2);
	K4:=Invert4Constants(Hadamard(Squaring(thetas)));

	K:=[K1,K2,K3,K4];

	return K;

end function;


JtoK:=function(P,rosen,K)
	// Maps a point P on the Jacobian of the Rosenhain curve C (with Rosenhain invariants 'rosen')
	// to the corresponding Fast Kummer surface K

	// As this map is a (2,2)-isogeny, this will only work if P is a point of odd order

	a2:=K[2][1]^2; b2:=K[2][2]^2; c2:=K[2][3]^2; d2:=K[2][4]^2; 
	lambda:=rosen[1]; mu:=rosen[2]; nu:=rosen[3];

	q:=Coefficients(P[1])[2];
	r:=Coefficients(P[1])[1];
	t:=Coefficients(P[2])[1];

	X:=a2*(r*(mu-r)*(lambda+q+nu)-t^2);
	Y:=b2*(r*(nu*lambda-r)*(1+q+mu)-t^2);
	Z:=c2*(r*(nu-r)*(lambda+q+mu)-t^2);
	T:=d2*(r*(mu*lambda-r)*(1+q+nu)-t^2);

	P:=[X,Y,Z,T];

	P:=Hadamard(P);
	P:=Squaring(P);
	P:=FourWayMult(P,K[4]);
	P:=Hadamard(P);
	P:=FourWayMult(P,K[3]);
	assert OnKummer(P,K);

	return normalise(P);

end function;

