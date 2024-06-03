/*
This file is used in Method 2 for scaling (in scalings.m) to compute the partition of monomials of degree N into those corresponding to the various coordinates of the isogeny (as described in Section 5.2.1 of the paper.

We stress that, for a fixed N, this can be done as a precomputation and so we don't count this in the total cost of the algorithm.

*/



function find_partition(N)
    /*
        INPUT: the degree N of the monomials in the equations for the (N,N)-isogeny

        OUTPUT: the partition of these monomials into four sets, where the i-th set gives the monomials in the i-th coordinate of the isogeny.    
    */

    num_mons := #MonomialsOfDegree(PolynomialRing(Rationals(), 4), N);
    pp := PolynomialRing(Rationals(), 4*num_mons);
    poly<X,Y,Z,T> := PolynomialRing(pp, 4);

    coord_transforms := [[X,Y,Z,T],[X,Y,-Z,-T],[X,-Y,Z,-T],[X,-Y,-Z,T]];

    mons := MonomialsOfDegree(poly, N);
    assert num_mons eq #mons;

    phi1 := poly!&+[pp.i*mons[i] : i in [1..num_mons]];
    phi2 := poly!&+[pp.(num_mons+i)*mons[i] : i in [1..num_mons]];
    phi3 := poly!&+[pp.(2*num_mons+i)*mons[i] : i in [1..num_mons]];
    phi4 := poly!&+[pp.(3*num_mons+i)*mons[i] : i in [1..num_mons]];

    phi_transforms := [[phi1,phi2,phi3,phi4],[phi1,phi2,-phi3,-phi4],[phi1,-phi2,phi3,-phi4],[phi1,-phi2,-phi3,phi4]];

    phi := phi_transforms[1];
    rels := [];

    for p in [1..4] do
        for i in [2..4] do
            F := Evaluate(phi[p], coord_transforms[i]) - phi_transforms[i][p];
            R := [MonomialCoefficient(F, m) : m in mons];
            for r in R do
                if r ne 0 and r notin rels then 
                    Append(~rels, r);
                end if;
            end for;
        end for;
    end for;

    // Using these relations to reduce the number of coefficients
    cs := [pp.i : i in [1..(4*num_mons)]];

    M:=[];
    for r in rels do
        Append(~M, [MonomialCoefficient(r, c) : c in cs]);
    end for;
    M:=Matrix(M);
    kermat := Basis(Nullspace(Transpose(M)));
    new_cs := [pp!0 : i in [1..(num_mons*4)]];
    for k in kermat do 
        i:=0;
        tmp := 0;
        while tmp eq 0 do
            i +:= 1;
            tmp := k[i];
        end while;
        for j in [1..(4*num_mons)] do
            if k[j] ne 0 then 
                new_cs[j] := cs[i];
            end if;
        end for;
    end for;

    phi1 := poly!&+[new_cs[i]*mons[i] : i in [1..num_mons]];
    phi2 := poly!&+[new_cs[num_mons+i]*mons[i] : i in [1..num_mons]];
    phi3 := poly!&+[new_cs[num_mons*2+i]*mons[i] : i in [1..num_mons]];
    phi4 := poly!&+[new_cs[num_mons*3+i]*mons[i] : i in [1..num_mons]];

    mons1 := Monomials(phi1);
    mons2 := Monomials(phi2);
    mons3 := Monomials(phi3);
    mons4 := Monomials(phi4);

    return [mons1, mons2, mons3, mons4];
end function; 


