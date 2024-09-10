/*
This script contains the functions necessary to compute (N,N)-isogenies between Kummer surfaces in the Fast Model when N >= 5

Note that this can be script can be adapted for use with other Kummer surface models 
by simply changing the '_biquadratics' function

When relevant, we give the cost of each algorithm. Here S, M, I and a denote squarings, multiplications, inversions and additions(/subtractions) over the base field.
*/


function _computeKummer(consts)
    /*
    Input: Constants giving the defining equation of the fast Kummer surface K 
    Output: The equation defining the Kummer surface (in projective space)
    */
    Fp2 := Parent(consts[1]);
	P3<X,Y,Z,T> := ProjectiveSpace(Fp2, 3);
    E,F,G,H := Explode(consts);

	return ((X^4+Y^4+Z^4+T^4)+ 2*E*X*Y*Z*T - F*(X^2*T^2+Y^2*Z^2)-G*(X^2*Z^2+Y^2*T^2) - H*(X^2*Y^2+T^2*Z^2));
end function;

_biquadratics_full := function(P,Q,thetas)
    /* 
    Input: Points P, Q on fast Kummer surface defined by fundamental constants thetas
    Output: The biquadratic forms evaluted at P and Q. 

    This is the unoptimised function, when we need to output the biquadratic forms
    corresponding to all indices
    */

    a,b,c,d := Explode(thetas);
    XP:=P[1]; YP:=P[2]; ZP:=P[3]; TP:=P[4];
    XQ:=Q[1]; YQ:=Q[2]; ZQ:=Q[3]; TQ:=Q[4];  
    
    B11:=(TP^2+XP^2+YP^2+ZP^2)*(TQ^2+XQ^2+YQ^2+ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2))+(-TP^2+XP^2+YP^2-ZP^2)*(-TQ^2+XQ^2+YQ^2-ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2))+(-TP^2+XP^2-YP^2+ZP^2)*(-TQ^2+XQ^2-YQ^2+ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2))+(TP^2+XP^2-YP^2-ZP^2)*(TQ^2+XQ^2-YQ^2-ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B22:=(TP^2+XP^2+YP^2+ZP^2)*(TQ^2+XQ^2+YQ^2+ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2))+(-TP^2+XP^2+YP^2-ZP^2)*(-TQ^2+XQ^2+YQ^2-ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2))-(-TP^2+XP^2-YP^2+ZP^2)*(-TQ^2+XQ^2-YQ^2+ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2))-(TP^2+XP^2-YP^2-ZP^2)*(TQ^2+XQ^2-YQ^2-ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B33:=(TP^2+XP^2+YP^2+ZP^2)*(TQ^2+XQ^2+YQ^2+ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2))-(-TP^2+XP^2+YP^2-ZP^2)*(-TQ^2+XQ^2+YQ^2-ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2))+(-TP^2+XP^2-YP^2+ZP^2)*(-TQ^2+XQ^2-YQ^2+ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2))-(TP^2+XP^2-YP^2-ZP^2)*(TQ^2+XQ^2-YQ^2-ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B44:=(TP^2+XP^2+YP^2+ZP^2)*(TQ^2+XQ^2+YQ^2+ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2))-(-TP^2+XP^2+YP^2-ZP^2)*(-TQ^2+XQ^2+YQ^2-ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2))-(-TP^2+XP^2-YP^2+ZP^2)*(-TQ^2+XQ^2-YQ^2+ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2))+(TP^2+XP^2-YP^2-ZP^2)*(TQ^2+XQ^2-YQ^2-ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B12:=(2*(a*b*(TP*TQ*ZP*ZQ+XP*XQ*YP*YQ)-c*d*(TP*XQ*YQ*ZP+TQ*XP*YP*ZQ)))/(((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)-((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B13:=(2*(a*c*(TP*TQ*YP*YQ+XP*XQ*ZP*ZQ)-b*d*(TP*XQ*YP*ZQ+TQ*XP*YQ*ZP)))/(((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2)-((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B14:=(2*(a*d*(TP*TQ*XP*XQ+YP*YQ*ZP*ZQ)-b*c*(TP*XP*YQ*ZQ+TQ*XQ*YP*ZP)))/(((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)-((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2));
    B23:=(2*(b*c*(TP*TQ*XP*XQ+YP*YQ*ZP*ZQ)-a*d*(TP*XP*YQ*ZQ+TQ*XQ*YP*ZP)))/(((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2)-((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2));
    B24:=(2*(b*d*(TP*TQ*YP*YQ+XP*XQ*ZP*ZQ)-a*c*(TP*XQ*YP*ZQ+TQ*XP*YQ*ZP)))/(((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)-((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2));
    B34:=(2*(c*d*(TP*TQ*ZP*ZQ+XP*XQ*YP*YQ)-a*b*(TP*XQ*YQ*ZP+TQ*XP*YP*ZQ)))/(((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)-((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2));
    
    B21:=B12;
    B31:=B13;
    B41:=B14;
    B32:=B23;
    B42:=B24;
    B43:=B34;
    
    B:=[
        [B11,B12,B13,B14],
        [B21,B22,B23,B24],
        [B31,B32,B33,B34],
        [B41,B42,B43,B44]
       ];

    return B;
end function;

_biquadratics:=function(P,Q,thetas)
    /* 
    Input: Points P, Q on fast Kummer surface defined by fundamental constants thetas
    Output: The biquadratic forms evaluted at P and Q. 

            Note that here we only ouput the biquadratic forms with indices in 
                        \{(1,1), (2,2) (1,2), (2,1), (3,3), (4,4), (3,4), (4,3)\} 
            as this will be sufficient for our application (when using the basis from Conjecture 5.2).
            We optimise the formulae so that they are inversion-less and minimise the number of multiplications performed.
            
            TOTAL COST: 12S + 39M + 43a
    */
    a,b,c,d := Explode(thetas);
    XP:=P[1]; YP:=P[2]; ZP:=P[3]; TP:=P[4];
    XQ:=Q[1]; YQ:=Q[2]; ZQ:=Q[3]; TQ:=Q[4];  

    HSP := Hadamard(Squaring(P));
    HSQ := Hadamard(Squaring(Q));
    HSO := Hadamard(Squaring(thetas)); 

    HSO12 := HSO[1]*HSO[2];
    HSO13 := HSO[1]*HSO[3];
    HSO14 := HSO[1]*HSO[4];
    HSO23 := HSO[2]*HSO[3];
    HSO24 := HSO[2]*HSO[4];
    HSO34 := HSO[3]*HSO[4];
    
    HSO124 := HSO12*HSO[4];
    HSO1234 := HSO124*HSO[3];

    Hinv := [HSO23*HSO[4], HSO14*HSO[3], HSO124, HSO12*HSO[3]];

    H := Hadamard([HSP[i]*HSQ[i]*Hinv[i] : i in [1..4]]); 
    B11, B22, B33, B44 := Explode(H);

    HSO12m34 := HSO12-HSO34;
    HSO13m24 := HSO13-HSO24;
    HSO14m23 := HSO14-HSO23;

    ab := a*b;
    cd := c*d;

    XQYQ := XQ*YQ;
    XPYP := XP*YP;
    ZPTP := ZP*TP;
    ZQTQ := ZQ*TQ;
    a12 := XQYQ*XPYP + ZPTP*ZQTQ;
    b12 := XQYQ*ZPTP + XPYP*ZQTQ;

    lambda23 := HSO12m34*HSO13m24;
    lambda34 := HSO14m23*HSO13m24;
    lambda234 := lambda23*HSO14m23;

    lambda34 := lambda34*HSO1234;

    ab := lambda34*ab;
    cd := lambda34*cd;

    B12 := 16*(ab*a12-cd*b12);
    B34 := 16*(ab*b12-cd*a12);

    B21:=B12;
    B43:=B34;

    B11 := lambda234*B11;
    B22 := lambda234*B22;
    B33 := lambda234*B33;
    B44 := lambda234*B44;


    B:=[[B11, B12, 1, 1],
        [B21, B22, 1, 1],
        [1, 1, B33, B34],
        [1, 1, B43, B44]];

    return B;

end function;


_computeForm := function(P, R, B, RN, BN, inds, N)
    /*
    This function computes the form corresponding to a list of indices. 

    INPUT: P, a generic point on the Kummer surface
           R, an N-torsion point on the Kummer surface 
           B, biquadratic forms evaluated at P and R
           RN, a list of length (N-1)/2 containing the multiples of R, namely [2*R, 3*R, ..., (N-1)/2*R] 
           BN, a list of length (N-1)/2 containing the biquadratic forms evaluated at P and l*R, where l = 2, ..., (N-1)/2
           inds, a list of indices of length N with each values in the range {1,...,4}. 
           N, the order of the point R (we input this rather than computing it as this is costly)

    OUTPUT: the form corresponding to the list of indices given as input 
                            (we sum over these forms to get the invariant form as defined in Equation (5.2) of the paper)
    */

    assert #RN eq (N-1)/2 - 1;
    assert #BN eq (N-1)/2 - 1;

    if N eq 3 then 
        prod := 1;
    else
        I := [[inds[2*i], inds[2*i+1]] : i in [2..((N-1)/2)]];
        prod := &*[BN[i][I[i][1], I[i][2]] : i in [1..((N-1)/2-1)]];
    end if;
    return P[inds[1]] * B[inds[2], inds[3]] * prod;
end function;

_findOrbit := function(I)
    /*
    This function computes the orbit of an list of indices under the action of Sym(N), taking into account that 
    some of the lists will correspond to the same invariant form (as the biquadratic forms are 
    symmetric, i.e., B_{i,j} = B_{j,i}). See Proposition (5.2) in the paper. 
    We want to remove any superflouous computations because, the larger the orbit, the more time it will take to generate 
    the invariant forms in FindBasis (which is the bottleneck step of the algorithm). 

    INPUT: I, a list of indices of length N

    OUTPUT: The orbit of I as described above
    */

    N := #I;

    if #{i : i in I} eq 1 then 
        return [[I, [1]]], 0;
    end if;

    a := I[1];
    b := I[#I];
    Na := #[i : i in I | i eq a];
    Nb := #[i : i in I | i eq b];
    assert Na + Nb eq N;

    if Nb gt Ceiling(N/2) then 
        // swap the roles of a and b
        tmp := a;
        a := b;
        b := tmp;

        tmp := Na;
        Na := Nb;
        Nb := tmp;
    end if;

    L := {1} join {2*i : i in [1..Floor(N/2)]}; 
    L := Subsets(L, Nb);
    orbit := [];
    for list in L do
        tmp := [a : i in [1..N]];
        for l in list do 
            tmp[l] := b;
        end for;
        if tmp[1] eq b then 
            orbit cat:= [[tmp, [2^(Nb-1)]]];
        else 
            orbit cat:= [[tmp, [2^(Nb)]]];
        end if;
    end for;

    possible_pairs := {2*i : i in [1..Floor(N/2)]};
    max_num_pairs := Floor(Nb/2);
    
    for n in [1..max_num_pairs] do 
        pairs := Subsets(possible_pairs, n);
        for P in pairs do 
            tmp := [a : i in [1..N]];
            for l in P do 
                tmp[l] := b;
                tmp[l+1] := b;
            end for;
            if Nb-2*n ne 0 then 
                L := {1} join {2*i : i in [1..Floor(N/2)]} diff P;
                L := Subsets(L, Nb-2*n);
                for list in L do
                    tmp1 := tmp;
                    for l in list do 
                        tmp1[l] := b;
                    end for;
                    if tmp1[1] eq b then 
                        orbit cat:= [[tmp1, [2^(Nb-2*n-1)]]];
                    else
                        orbit cat:= [[tmp1, [2^(Nb-2*n)]]];
                    end if;
                end for;
            else
                orbit cat:= [[tmp, [1]]];
            end if;
        end for;
    end for;

    return orbit;

end function;

_invariantForm := function(P, R, B, RN, BN, inds, N)
    /*
    Computes the invariant form corresponding to some index list 

    INPUT: P, a generic point on the Kummer surface
           R, an N-torsion point on the Kummer surface 
           B, biquadratic forms evaluated at P and R
           RN, a list of length (N-1)/2 containing the multiples of R, namely [2*R, 3*R, ..., (N-1)/2*R] 
           BN, a list of length (N-1)/2 containing the biquadratic forms evaluated at P and l*R, where l = 2, ..., (N-1)/2
           inds, a list of indices of length N with each values in the range {1,...,4}. 
           N, the order of the point R (we input this rather than computing it as this is costly)
    
    OUTPUT: a form *invariant under translation by R* corresponding to the list of indices given as input (as described in Equation (5.2))

    */
    result := 0;

    orbit := _findOrbit(inds);
    for I in orbit do 
        g := _computeForm(P, R, B, RN, BN, I[1], N);
        result +:= I[2][1]*g;
    end for;  

    return result;
end function;


FindBasis := function(P, R, B, RN, BN, N, basis_inds)
    /*
    Computes the basis of forms invariant under translation by an N-torsion point

    INPUT: P, a generic point on the Kummer surface
        R, an N-torsion point on the Kummer surface 
        B, biquadratic forms evaluated at P and R
        RN, a list of length (N-1)/2 containing the multiples of R, namely [2*R, 3*R, ..., (N-1)/2*R] 
        BN, a list of length (N-1)/2 containing the biquadratic forms evaluated at P and l*R, where l = 2, ..., (N-1)/2
        N, the order of the point R (we input this rather than computing it as this is costly)
        basis_inds, an array of the form [I1, I2, I3, I_4] where I_k is an array containing the the indices used to 
                    construct the basis of the k-th part B^(k)_R. If basis_inds = [], we use the basis given
                    in the paper (from Conjecture 5.4).

    OUTPUT: the basis of invariant forms corresponding to R. Should always be of size (2N+2)

    */
    
    // "Generating the invariant forms..";
    if #basis_inds eq 0 then 
        F_ind := func< a,b,c | [a : i in [0..c-1]] cat [b : i in [c .. N-1]]>;
        inds12 := [F_ind(1,2,i) : i in [N..0 by -1]];
        inds34 := [F_ind(3,4,i) : i in [N..0 by -1]];
        inds := [[ inds12[i] : i in [1..#inds12] | i mod 2 eq 1 ]] cat [[ inds12[i] : i in [1..#inds12] | i mod 2 eq 0 ]] cat [[ inds34[i] : i in [1..#inds34] | i mod 2 eq 1 ]] cat [[ inds34[i] : i in [1..#inds34] | i mod 2 eq 0 ]];
    else 
        inds := basis_inds;
    end if;
    // "Number of inds is:", #inds;
    inv_forms := [[_invariantForm(P, R, B, RN, BN, i, N) : i in I] : I in inds];
    // "Found all the invariant forms..";
    return inv_forms;
    
end function;


FindIntersection:=function(XR, XS, N)
    /*
    Finds the intersection of two spaces 

    INPUT: XR, a list of forms generating the space of forms invariant under translation by an N-torsion point R
           XS, a list of forms generating the space of forms invariant under translation by an N-torsion point S
           N, the order of the points R and S

    OUTPUT: the basis of the intersection of the spaces generated by XR and XS. Should always be dimension 4
    */

    monomials := Monomials(XR[1]);
    Fp2 := Parent(Coefficients(XR[1])[1]);

    m := (N+1)/2;
    
    assert #XR eq m;
    assert #XS eq m;
    
	Rbasis := [[MonomialCoefficient(XR[i], j) : j in monomials[1..(N+1)]] : i in [1..m]]; //can be monomials[1..N] but let's just make it sqaure as in the paper
	Sbasis := [[MonomialCoefficient(-XS[i], j) : j in monomials[1..(N+1)]] : i in [1..m]]; //can be monomials[1..N] but let's just make it sqaure as in the paper
	RSbasis:= Rbasis cat Sbasis;
	RSmat := Matrix(Fp2, RSbasis);
	RSker := Basis(Nullspace(RSmat));
	basis := [];
	for i in [1..#RSker] do
		f:=0;
		for j in [1..#XR] do
			f +:= RSker[i, j]*XR[j];
		end for;
		Append(~basis, f);
	end for;
	return basis;
end function;

function quotient(f, g)
    /* Computes the image of f \in k[X,Y,Z,T] in k[X,Y,Z,T]/(g)
    // Ṭhis could probably be done without using these Magma functions, but we leave this for now
    */

    poly := Parent(f);
    I := ideal<poly|g>;
    QQ := poly/I;
    h := QQ!f;
    h := poly!h;

    return h;
end function;

load "scalings.m";

function GetIsogeny(P, R, S, K, N, method : timing := false, basis_inds := [])
    /*
    Computes the basis of forms invariant under translation by an N-torsion point

    INPUT: P, a generic point on the Kummer surface
           R and S, the N-torsion points generating the kernel of the isogeny 
           K, the domain Kummer surface
           N, the order of the point R (we input this rather than computing it as this is costly)
           method, the index of the method used for scaling. If method = 0 then we do not do the final scaling (useful if wanting to use GetImage)
           (optional) timing, for benchmarking reasons. Times each of the main subroutines in the algorithm separately
           (optional) basis_inds, an array of the form [I1, I2, I3, I_4] where I_k is an array containing the indices used to 
                                construct the basis of B^(k)_R and B^(k)_S (i.e., the k-th part of B_R and B_S). If this is not given as input, we use the indices given in the paper (from Conjecture 5.4).

    OUTPUT: phi, the equations for of the (N,N)-isogeny with kernel <R, S>. 
                If method = 2 or 3, the image of phi is in correct form. 
                If method = 0, the phi is unscaled and the image is not in the correct form

            (optional) the times if timings := true
    */

    // ṬODO: change so that if N = 5 then it calls iso55.m

    if timing then 
        times := [[], [], [], []];
        t1 := Cputime();
    end if;

    thetas := K[2];    
    RN := [ScalarMultKummer(R, j, K) : j in [2..((N-1)/2)]];
    SN := [ScalarMultKummer(S, j, K) : j in [2..((N-1)/2)]];

    if #basis_inds eq 0 then 
        // Use our optimised function when we use the basis in Conjecture 5.2
        BR := _biquadratics(P,R,thetas);
        BRN := [_biquadratics(P,RN[i],thetas) : i in [1..#RN]];
        BS := _biquadratics(P,S,thetas);
        BSN := [_biquadratics(P,SN[i],thetas) : i in [1..#SN]];
    else 
        BR := _biquadratics_full(P,R,thetas);
        BRN := [_biquadratics_full(P,RN[i],thetas) : i in [1..#RN]];
        BS := _biquadratics_full(P,S,thetas);
        BSN := [_biquadratics_full(P,SN[i],thetas) : i in [1..#SN]];
    end if;

    if timing then 
        t := Cputime();
    end if;
    Rbasis := FindBasis(P, R, BR, RN, BRN, N, basis_inds);
    Sbasis := FindBasis(P, S, BS, SN, BSN, N, basis_inds); 
    if timing then 
        Append(~times[1], Cputime(t));
    end if;

    // ̣quotient forms
    K_eqn:=_computeKummer(K[1]);
    Rset := [ [quotient(r, K_eqn) : r in RR] : RR in Rbasis];
    Sset := [ [quotient(s, K_eqn) : s in SS] : SS in Sbasis];
    if timing then 
        t := Cputime();
    end if;
    phi1:= FindIntersection(SetToIndexedSet(SequenceToSet(Rset[1])), SetToIndexedSet(SequenceToSet(Sset[1])), N)[1];
    phi2:= FindIntersection(SetToIndexedSet(SequenceToSet(Rset[2])), SetToIndexedSet(SequenceToSet(Sset[2])), N)[1]; 
    phi3:= FindIntersection(SetToIndexedSet(SequenceToSet(Rset[3])), SetToIndexedSet(SequenceToSet(Sset[3])), N)[1];
    phi4:= FindIntersection(SetToIndexedSet(SequenceToSet(Rset[4])), SetToIndexedSet(SequenceToSet(Sset[4])), N)[1];
    
    if timing then 
        Append(~times[2], Cputime(t));
    end if;
    
    phi := [phi1, phi2, phi3, phi4];

    if method eq 0 then
        // Here we return the isogeny with no final scaling 
        // The unscaled isogeny can then be used in GetImage 
        // To obtain the equation for the image Kummer
        // with no sqrts or linear algebra
        if timing then 
            Append(~times[4], Cputime(t1));
            return phi, times;
        else
            return phi;
        end if;
    elif N eq 5 then 
        // Use the scaling method for N=5
        if timing then 
            t := Cputime();
        end if;
        phi := Scaling_5(phi);
        if timing then
            Append(~times[3], Cputime(t));
            Append(~times[4], Cputime(t1));
            return phi, times;
        else 
            return phi;
        end if;
    elif method eq 2 then
        // METHOD 2: use linear algebra

        if timing then 
            t := Cputime();
        end if;
        phi := Scaling_GE(phi, K);
        if timing then
            Append(~times[3], Cputime(t));
            Append(~times[4], Cputime(t1));
            return phi, times;
        else 
            return phi;
        end if;
    elif method eq 3 then 
        // METHOD 3: take sqrts
        if timing then 
            t := Cputime();
        end if;
        phi := Scaling_sqrt(phi, K);
        if timing then 
            Append(~times[3], Cputime(t));
            Append(~times[4], Cputime(t1));
            return phi, times;
        else 
            return phi; 
        end if;
    else 
        "There is some problem, 'method' input incorrectly";
    return false;
    end if;

end function;

function GetImage(phi, K : timing := false)
    /*
    Computes the equation defining the image Kummer surface of an (N,N)-isogeny
    without taking square roots or using any linear algebra

    INPUT: phi, an (N,N)-isogeny whose image is not in the correct fast Kummer form 
           K, the domain Kummer surface
           (optional) timing, for benchmarking purposes
    
    OUTPUT: Kim, the equation defining the image Kummer surface in fast Kummer form
            (optional) t, the time taken to run the algorithm

    */

    // ṬODO: remove inversions if possible 
    if timing then 
        t := Cputime();
    end if;

    a,b,c,d := Explode(K[2]);
    b2a2 := Evaluate(phi[1], [b,a,d,c])*Evaluate(phi[2], [a,b,c,d])/Evaluate(phi[1], [a,b,c,d])/Evaluate(phi[2], [b,a,d,c]);
    d2c2 := Evaluate(phi[3], [b,a,d,c])*Evaluate(phi[4], [a,b,c,d])/Evaluate(phi[3], [a,b,c,d])/Evaluate(phi[4], [b,a,d,c]);
    c2b2 := Evaluate(phi[2], [d,c,b,a])*Evaluate(phi[3], [a,b,c,d])/Evaluate(phi[2], [a,b,c,d])/Evaluate(phi[3], [d,c,b,a]);
    acbd := Evaluate(phi[1], [a,b,c,d])*Evaluate(phi[2], [d,c,b,a])/Evaluate(phi[2], [a,b,c,d])/Evaluate(phi[1], [d,c,b,a]);
    c2a2 := c2b2*b2a2;
    d2a2 := d2c2*c2a2;
    a2, b2, c2, d2 := Explode([1, b2a2, c2a2, d2a2]);

    Ahat,Bhat,Chat,Dhat := Explode(Hadamard([a2,b2,c2,d2]));

    Ehat := acbd*b2*d2*Ahat*Bhat*Chat*Dhat/(a2*d2-b2*c2)/(a2*c2-b2*d2)/(a2*b2-c2*d2);
    Fhat := (a2^2-b2^2-c2^2+d2^2)/(a2*d2-b2*c2);
    Ghat := (a2^2-b2^2+c2^2-d2^2)/(a2*c2-b2*d2);
    Hhat := (a2^2+b2^2-c2^2-d2^2)/(a2*b2-c2*d2);

    Khat := _computeKummer([Ehat, Fhat, Ghat, Hhat]);
    if timing then 
        return Khat, Cputime(t);
    else 
        return Khat;
    end if;
end function;
