
/*
This file is used in isoNN.m to compute the final scaling that brings the image of the isogeny into the correct form

We stress that, for a fixed N, this can be done as a precomputation and so we don't count this in the total cost of the algorithm.
*/

load "partition.m";

///////////////////////////////////////////
//////////// SCALING FOR N = 5 ////////////
///////////////////////////////////////////

function Scaling_5(phi);
    /*
    Method of finding the final scaling for N = 5, as describe in Algorithm 3 of the paper.

    INPUT: phi, the (5,5)-isogeny (whose image is in incorrect form)

    OUTPUT: the scaled (5,5)-isogeny whose image is in fast Kummer form

    TOTAL COST: 6M
    */
    phi1, phi2, phi3, phi4 := Explode(phi);
    PP := Parent(phi1);
    X := PP.1; Y := PP.2; Z := PP.3; T := PP.4;

    lX := 1;

    cY1 := MonomialCoefficient(phi1, Y*Z*T^3);
    cY2 := MonomialCoefficient(phi2, X*Z^3*T);

    cZ1 := MonomialCoefficient(phi1, Y*Z*T^3);
    cZ2 := MonomialCoefficient(phi3, T*X*Y^3);

    cT1 := MonomialCoefficient(phi3, X*Y*T^3);
    cT2 := MonomialCoefficient(phi4, X*Y*Z^3);

    cZ2T2 := cZ2*cT2;
    cY2Z1 := cY2*cZ1;

    scals := [cY2*cZ2T2, cY1*cZ2T2, cY2Z1*cT2, cY2Z1*cT1];
    
    return FourWayMult(phi, scals); 

end function;


///////////////////////////////////////////
////////// METHOD 2 FOR SCALING ///////////
///////////////////////////////////////////

function Scaling_GE(phi, K);
    /*
    Method 2 of finding the final scaling for N >= 7, as describe in Algorithm 4 of the paper.

    INPUT: phi, the (N,N)-isogeny (whose image is in incorrect form)
           K, the domain fast Kummer surface  

    OUTPUT: the scaled (N,N)-isogeny whose image is in fast Kummer form
    */

    phi1, phi2, phi3, phi4 := Explode(phi);
    N := Degree(phi1);

    Fp2 := Parent(K[1][1]);
    n := #MonomialsOfDegree(PolynomialRing(Rationals(), 4), N-4);
    FF := FunctionField(Fp2, n+4);
    PP<X,Y,Z,T> := PolynomialRing(FF, 4);
    
    ms := find_partition(N-4);
    ms := [[PP!m : m in mons] : mons in ms];
    M := #ms[1];

    Keqn := PP!_computeKummer(K[1]);

    p1 := &+[FF.i*ms[1][i] : i in [1..M]];
    p2 := &+[FF.(M + i)*PP!ms[1][i] : i in [1..M]];
    p3 := &+[FF.(2*M + i)*PP!ms[1][i] : i in [1..M]];
    
    Phi1 := FF.(n+1)*PP!phi1;
    Phi2 := FF.(n+2)*PP!phi2;
    Phi3 := FF.(n+3)*PP!phi3;
    Phi4 := FF.(n+4)*PP!phi4;

    F1 := Evaluate(Phi2, [Y,X,T,Z]) - Evaluate(Phi1, [X,Y,Z,T]) + p1*Keqn;
    M1 := [MonomialCoefficient(F1, m) : m in Monomials(F1)];
    M1 := [[MonomialCoefficient(Numerator(m),Numerator(FF.i)) : i in [(n+2)] cat [1..M] cat [(n+1)]] : m in M1[1..(M+2)]];
    
    F2 := Evaluate(Phi3, [Z,T,X,Y]) - Evaluate(Phi1, [X,Y,Z,T]) + p2*Keqn;
    M2 := [MonomialCoefficient(F2, m) : m in Monomials(F2)];
    M2 := [[MonomialCoefficient(Numerator(m),Numerator(FF.i)) : i in [(n+3)] cat [(M+1)..2*M] cat [(n+1)]] : m in M2[1..(M+2)]];

    F3 := Evaluate(Phi4, [T,Z,Y,X]) - Evaluate(Phi1, [X,Y,Z,T]) + p3*Keqn;
    M3 := [MonomialCoefficient(F3, m) : m in Monomials(F3)];
    M3 := [[MonomialCoefficient(Numerator(m),Numerator(FF.i)) : i in [(n+4)] cat [(2*M+1)..3*M] cat [(n+1)]] : m in M3[1..(M+2)]];
    
    k1 := EchelonForm(Matrix(M1));
    c1 := NumberOfColumns(k1);
    lY := -k1[1][c1];
    k2 := EchelonForm(Matrix(M2));
    c2 := NumberOfColumns(k2);
    lZ := -k2[1][c2];
    k3 := EchelonForm(Matrix(M3));
    c3 := NumberOfColumns(k3);
    lT := -k3[1][c3];

    phi := FourWayMult(phi, [1,lY,lZ,lT]);

    return phi; 

end function;


///////////////////////////////////////////
////////// METHOD 3 FOR SCALING ///////////
///// 2 Sqrts needed to find the map ////// 
////// (but none to find the image) ///////
///////////////////////////////////////////

function _get_scals(phi, K)
    /*
    INPUT: phi, the (N,N)-isogeny (whose image is in incorrect form)
           K, the domain fast Kummer surface  

    OUTPUT: scaling data that will be used to find the correct scalings (using sqrts)

    */

    a,b,c,d := Explode(K[2]);

    lYX2 := Evaluate(phi[1], [a,b,c,d])*Evaluate(phi[1], [b,a,d,c])/Evaluate(phi[2], [a,b,c,d])/Evaluate(phi[2], [b,a,d,c]);
    b2a2 := Evaluate(phi[1], [b,a,d,c])*Evaluate(phi[2], [a,b,c,d])/Evaluate(phi[1], [a,b,c,d])/Evaluate(phi[2], [b,a,d,c]);
    lTZ2 := Evaluate(phi[3], [b,a,d,c])*Evaluate(phi[3], [a,b,c,d])/Evaluate(phi[4], [a,b,c,d])/Evaluate(phi[4], [b,a,d,c]);
    d2c2 := Evaluate(phi[3], [b,a,d,c])*Evaluate(phi[4], [a,b,c,d])/Evaluate(phi[3], [a,b,c,d])/Evaluate(phi[4], [b,a,d,c]);
    lZY2 := Evaluate(phi[2], [a,b,c,d])*Evaluate(phi[2], [d,c,b,a])/Evaluate(phi[3], [a,b,c,d])/Evaluate(phi[3], [d,c,b,a]);
    c2b2 := Evaluate(phi[2], [d,c,b,a])*Evaluate(phi[3], [a,b,c,d])/Evaluate(phi[2], [a,b,c,d])/Evaluate(phi[3], [d,c,b,a]);
    lXZlYT := Evaluate(phi[2], [a,b,c,d])*Evaluate(phi[4], [d,c,b,a])/Evaluate(phi[1], [a,b,c,d])/Evaluate(phi[3], [d,c,b,a]);
    acbd := Evaluate(phi[1], [a,b,c,d])*Evaluate(phi[2], [d,c,b,a])/Evaluate(phi[2], [a,b,c,d])/Evaluate(phi[1], [d,c,b,a]);
    c2a2 := c2b2*b2a2;
    d2a2 := d2c2*c2a2;

    scals := [lYX2, lTZ2, lZY2, lXZlYT, b2a2, d2c2, c2b2, acbd];

    return scals;
end function;


function _get_isogeny(phi, scals)
    /*
    INPUT: phi, the (N,N)-isogeny (whose image is in incorrect form)
           scals, scaling data that is used to find the correct scalings (using sqrts)

    OUTPUT: the scaled isogeny phi whose image is in fast Kummer form
    */

    Fp2 := Parent(scals[1]);
    poly := Parent(phi[1]);
    X := poly.1;
    Y := poly.2;
    Z := poly.3;
    T := poly.4;

    lYX2, lTZ2, lZY2, lXZlYT, b2a2, d2c2, c2b2, acbd := Explode(scals);

    try 
        s1 := Sqrt(Fp2!lYX2);
    catch e 
        print "Error in Scaling_sqrt as there is no square root over Fp2";
        error "Likely landing on a product of elliptic curves";
    end try;
    try 
        s2 := Sqrt(Fp2!lZY2);
    catch e 
        print "Error in Scaling_sqrt as there is no square root over Fp2";
        error "Likely landing on a product of elliptic curves";
    end try;
    
    
    lYXs := [s1, -s1];
    lZYs := [s2, -s2];
    
    for lYX in lYXs do 
        for lZY in lZYs do 
            lZX := lYX*lZY;
            lTX := lZY/lXZlYT;
            scals := [1, lYX, lZX, lTX];
            phi_scal := FourWayMult(phi, scals);
            phi_scal := [p/Coefficients(phi_scal[1])[1] : p in phi_scal];

            phi_scal1 := [Evaluate(p, [X,Y,-Z,-T]) : p in phi_scal];
            phi_scal1 := [p/Coefficients(phi_scal1[1])[1] : p in phi_scal1];
            phi_scal2 := [Evaluate(p, [X,-Y,Z,-T]) : p in phi_scal];
            phi_scal2 := [p/Coefficients(phi_scal2[1])[1] : p in phi_scal2];
            phi_scal3 := [Evaluate(p, [X,-Y,-Z,T]) : p in phi_scal];
            phi_scal3 := [p/Coefficients(phi_scal3[1])[1] : p in phi_scal3];

            if (FourWayMult(phi_scal, [1,1,-1,-1]) eq phi_scal1) and (FourWayMult(phi_scal, [1,-1,1,-1]) eq phi_scal2) and (FourWayMult(phi_scal, [1,-1,-1,1]) eq phi_scal3) then 
                return phi_scal;
            end if;
        end for;
    end for;
    "There's some issue with taking sqrts";
    assert false;

end function;


function Scaling_sqrt(phi, K)
    scals := _get_scals(phi, K);
    phi := _get_isogeny(phi, scals);
    return phi;
end function;