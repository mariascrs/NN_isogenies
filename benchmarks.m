/*
This script contains the functions necessary to run the benchmarks used to create Figure 1, Figure 2, Table 1 and Table 2.
*/

load "functions.m";
load "isoNN.m";

function benchmark_fig1(num_primes, num_samples)
    // Experiments for Figure 1
    N := 7;

    // Generate the primes 
    primes := [];
    l2 := 10;
    "Generating the primes...";
    while #primes lt Min(num_primes, 20) do
    #primes;
        condition := false;
        while not condition do
            e := Random([2..l2]);
            p:=2^4*3^e*N-1;
            if IsPrime(p) and p notin primes then 
                condition := true;
            end if;
        end while;
        Append(~primes,p);
        l2 +:= 100;
    end while;
    primes := Sort(primes);
    logp := [Ceiling(Log(2,p)) : p in primes];
    "Generated ", #primes, "primes.\n";
    logp;

    // First with method 2
    "Getting results for Method 2.";
    av_times_m2 := [];
    method := 2;
    for p in primes do 
        "Working with: ", p; 
        Fp:=GF(p);
        Fp2<i>:=ExtensionField<Fp,x|x^2+1>;


        _<x>:=PolynomialRing(Fp2);

        C:=HyperellipticCurve(x^6+1);
        C:=random_walk(C);
        K,rosen:=KummerFromCurve(C);
        C:=HyperellipticCurve(x*(x-1)*(x-rosen[1])*(x-rosen[2])*(x-rosen[3]));
        J:=Jacobian(C);

        poly<X,Y,Z,T>:=PolynomialRing(Fp2,4);
        P := [X,Y,Z,T];

        times := [];
        for i in [1..num_samples] do
            repeat
                repeat
                    R:=((p+1) div N)*Random(J);
                until N*R eq J!0 and R ne J!0;
                repeat
                    S:=((p+1) div N)*Random(J);
                until N*S eq J!0 and S notin [J!0,R,-R];
            until WeilPairing(R,S,N) eq 1;

            R:=JtoK(R,rosen,K);
            S:=JtoK(S,rosen,K);

            "Sample: ", i;
            phi, ts := GetIsogeny(P, R, S, K, N, method: timing := true);
            Append(~times, ts[3][1]);
        end for;

        
        T := &+times/num_samples;
        Append(~av_times_m2, T);

        "Average:", T, "seconds.\n";
    end for;

    // Then  with method 3
    "Getting results for Method 3.";
    av_times_m3 := [];
    method := 3;
    for p in primes do 
        "Working with: ", p; 
        Fp:=GF(p);
        Fp2<i>:=ExtensionField<Fp,x|x^2+1>;


        _<x>:=PolynomialRing(Fp2);

        C:=HyperellipticCurve(x^6+1);
        C:=random_walk(C);
        K,rosen:=KummerFromCurve(C);
        C:=HyperellipticCurve(x*(x-1)*(x-rosen[1])*(x-rosen[2])*(x-rosen[3]));
        J:=Jacobian(C);

        poly<X,Y,Z,T>:=PolynomialRing(Fp2,4);
        P := [X,Y,Z,T];

        times := [];
        for i in [1..num_samples] do
            repeat
                repeat
                    R:=((p+1) div N)*Random(J);
                until N*R eq J!0 and R ne J!0;
                repeat
                    S:=((p+1) div N)*Random(J);
                until N*S eq J!0 and S notin [J!0,R,-R];
            until WeilPairing(R,S,N) eq 1;

            R:=JtoK(R,rosen,K);
            S:=JtoK(S,rosen,K);

            "Sample: ", i;
            phi, ts := GetIsogeny(P, R, S, K, N, method: timing := true);
            Append(~times, ts[3][1]);
        end for;

        
        T := &+times/num_samples;
        Append(~av_times_m3, T);

        "Average:", T, "seconds.\n";
        
    end for;

    return primes, logp, av_times_m2, av_times_m3;

end function;

function benchmark_fig2(logp, delta, num_samples : Ns := [5, 7, 11, 13, 17, 19] , primes := [])
    // Experiments for Figure 2

    // Generate the primes 
    
    if #primes ne #Ns then 
        "Generating the prime for each N of logarithm around", logp, "with an error on each side of", delta;
        primes := [];
        l2 := 10;

        for N in Ns do
            condition := false;
            while not condition do
                e := Random([2..logp]);
                p:=2^4*3^e*N-1;
                if IsPrime(p) and Ceiling(Log(2,p)) in [(logp-delta)..(logp+delta)] then 
                    condition := true;
                end if;
            end while;
            Append(~primes,p);
        end for;
        "Generated the primes.\n";
    else
        "Using pre-generated primes.\n";
    end if;

    logps := [Ceiling(Log(2,p)) : p in primes];

    // First with method 2
    "Getting results for Method 2.";
    av_times_m2 := [];
    method := 2;
    for j in [1..#primes] do 
        p := primes[j];
        N := Ns[j];

        "Working with N =", N; 

        Fp:=GF(p);
        Fp2<i>:=ExtensionField<Fp,x|x^2+1>;


        _<x>:=PolynomialRing(Fp2);

        C:=HyperellipticCurve(x^6+1);
        C:=random_walk(C);
        K,rosen:=KummerFromCurve(C);
        C:=HyperellipticCurve(x*(x-1)*(x-rosen[1])*(x-rosen[2])*(x-rosen[3]));
        J:=Jacobian(C);

        poly<X,Y,Z,T>:=PolynomialRing(Fp2,4);
        P := [X,Y,Z,T];

        times := [];
        for i in [1..num_samples] do
            repeat
                repeat
                    R:=((p+1) div N)*Random(J);
                until N*R eq J!0 and R ne J!0;
                repeat
                    S:=((p+1) div N)*Random(J);
                until N*S eq J!0 and S notin [J!0,R,-R];
            until WeilPairing(R,S,N) eq 1;

            R:=JtoK(R,rosen,K);
            S:=JtoK(S,rosen,K);

            "Sample:", i;
            phi, ts := GetIsogeny(P, R, S, K, N, method: timing := true);
            Append(~times, ts[3][1]);
        end for;

        T := &+times/num_samples;
        Append(~av_times_m2, T);

        "Average:", T, "seconds.\n";
    end for;

    // Then  with method 3
    "Getting results for Method 3.";
    av_times_m3 := [];
    method := 3;
    
    // Starting from N := 7 as we don't want to reproduce the N = 5 timings (we only use one scaling there)
    if Ns[1] eq 5 then 
        min:= 2;
    else
        min:= 1;
    end if;

    for j in [min..#primes] do 
        p := primes[j];
        N := Ns[j];

        "Working with N =", N; 

        Fp:=GF(p);
        Fp2<i>:=ExtensionField<Fp,x|x^2+1>;


        _<x>:=PolynomialRing(Fp2);

        C:=HyperellipticCurve(x^6+1);
        C:=random_walk(C);
        K,rosen:=KummerFromCurve(C);
        C:=HyperellipticCurve(x*(x-1)*(x-rosen[1])*(x-rosen[2])*(x-rosen[3]));
        J:=Jacobian(C);

        poly<X,Y,Z,T>:=PolynomialRing(Fp2,4);
        P := [X,Y,Z,T];

        times := [];
        for i in [1..num_samples] do
            repeat
                repeat
                    R:=((p+1) div N)*Random(J);
                until N*R eq J!0 and R ne J!0;
                repeat
                    S:=((p+1) div N)*Random(J);
                until N*S eq J!0 and S notin [J!0,R,-R];
            until WeilPairing(R,S,N) eq 1;

            R:=JtoK(R,rosen,K);
            S:=JtoK(S,rosen,K);

            "Sample:", i;
            phi, ts := GetIsogeny(P, R, S, K, N, method : timing := true);
            Append(~times, ts[3][1]);
        end for;

        T := &+times/num_samples;
        Append(~av_times_m3, T);

        "Average:", T, "seconds.\n";
    end for;

    return primes, logps, av_times_m2, av_times_m3;

end function;

function benchmark_table1(logp, delta, num_samples : Ns := [5, 7, 11, 13, 17, 19] , primes := [])
    // Experiments for Figure 2

    // Generate the primes 
    if #primes ne #Ns then 
        "Generating the prime for each N of logarithm around", logp, "with an error on each side of", delta;
        primes := [];
        l2 := 10;
        

        for N in Ns do
            condition := false;
            while not condition do
                e := Random([2..logp]);
                p:=2^4*3^e*N-1;
                if IsPrime(p) and Ceiling(Log(2,p)) in [(logp-delta)..(logp+delta)] then 
                    condition := true;
                end if;
            end while;
            Append(~primes,p);
        end for;
        "Generated the primes.\n";
    else
        "Using pre-generated primes.\n";
    end if;

    logps := [Ceiling(Log(2,p)) : p in primes];

    // First with method 2
    "Getting results for Method 2.";
    av_times_m2 := [];
    method := 2;
    for j in [1..#primes] do 
        p := primes[j];
        N := Ns[j];

        "Working with N =", N; 

        Fp:=GF(p);
        Fp2<i>:=ExtensionField<Fp,x|x^2+1>;


        _<x>:=PolynomialRing(Fp2);

        C:=HyperellipticCurve(x^6+1);
        C:=random_walk(C);
        K,rosen:=KummerFromCurve(C);
        C:=HyperellipticCurve(x*(x-1)*(x-rosen[1])*(x-rosen[2])*(x-rosen[3]));
        J:=Jacobian(C);

        poly<X,Y,Z,T>:=PolynomialRing(Fp2,4);
        P := [X,Y,Z,T];

        times := [];
        for i in [1..num_samples] do
            repeat
                repeat
                    R:=((p+1) div N)*Random(J);
                until N*R eq J!0 and R ne J!0;
                repeat
                    S:=((p+1) div N)*Random(J);
                until N*S eq J!0 and S notin [J!0,R,-R];
            until WeilPairing(R,S,N) eq 1;

            R:=JtoK(R,rosen,K);
            S:=JtoK(S,rosen,K);

            "Sample:", i;
            phi, ts := GetIsogeny(P, R, S, K, N, method: timing := true);
            Append(~times, ts);
        end for;

        av_times := [];
        for i in [1..4] do 
            T := &+[ts[i][1] : ts in times]/num_samples;
            Append(~av_times, T);
        end for;
        "Average:", av_times, ".\n";
        Append(~av_times_m2, av_times);
    end for;

    // Then with method 3
    "Getting results for Method 3.";
    av_times_m3 := [];
    method := 3;
    
    // Starting from N := 7 as we don't want to reproduce the N = 5 timings (we only use one scaling there)
    if Ns[1] eq 5 then 
        min:= 2;
    else
        min:= 1;
    end if;

    for j in [min..#primes] do 
        p := primes[j];
        N := Ns[j];

        "Working with N =", N; 

        Fp:=GF(p);
        Fp2<i>:=ExtensionField<Fp,x|x^2+1>;


        _<x>:=PolynomialRing(Fp2);

        C:=HyperellipticCurve(x^6+1);
        C:=random_walk(C);
        K,rosen:=KummerFromCurve(C);
        C:=HyperellipticCurve(x*(x-1)*(x-rosen[1])*(x-rosen[2])*(x-rosen[3]));
        J:=Jacobian(C);

        poly<X,Y,Z,T>:=PolynomialRing(Fp2,4);
        P := [X,Y,Z,T];

        times := [];
        for m in [1..num_samples] do
            repeat
                repeat
                    R:=((p+1) div N)*Random(J);
                until N*R eq J!0 and R ne J!0;
                repeat
                    S:=((p+1) div N)*Random(J);
                until N*S eq J!0 and S notin [J!0,R,-R];
            until WeilPairing(R,S,N) eq 1;

            R:=JtoK(R,rosen,K);
            S:=JtoK(S,rosen,K);

            "Sample:", m;
            phi, ts := GetIsogeny(P, R, S, K, N, method : timing := true);
            Append(~times, ts);

            av_times := [];
            for i in [1..4] do 
                T := &+[ts[i][1] : ts in times]/m;
                Append(~av_times, T);
            end for;
            "Average:", av_times, ".\n";
        end for;

        av_times := [];
        for i in [1..4] do 
            T := &+[ts[i][1] : ts in times]/num_samples;
            Append(~av_times, T);
        end for;
        "Average:", av_times, ".\n";
        Append(~av_times_m3, av_times);
    end for;

    return primes, logps, av_times_m2, av_times_m3;
end function;

function benchmark_getimage(logp, delta, num_samples : Ns := [7,11,13,17,19], primes := [])
    // Experiments for Table 2
    
    // Generate the primes 
    if #primes ne #Ns then 
        "Generating the prime for each N of logarithm around", logp, "with an error on each side of", delta;
        primes := [];
        l2 := 10;

        for N in Ns do
            condition := false;
            while not condition do
                e := Random([2..logp]);
                p:=2^4*3^e*N-1;
                if IsPrime(p) and Ceiling(Log(2,p)) in [(logp-delta)..(logp+delta)] then 
                    condition := true;
                end if;
            end while;
            Append(~primes,p);
            #primes;
        end for;
        "Generated the primes.\n";
    else
        "Using pre-generated primes.\n";
    end if;

    logps := [Ceiling(Log(2,p)) : p in primes];

    "Getting timings for GetImage.";
    av_times:= [];
    method := 0;

    for j in [1..#primes] do 
        p := primes[j];
        N := Ns[j];

        "Working with N =", N; 
        "With prime p s.t. Log(2,p) = ", Ceiling(Log(2,p));

        Fp:=GF(p);
        Fp2<i>:=ExtensionField<Fp,x|x^2+1>;


        _<x>:=PolynomialRing(Fp2);

        C:=HyperellipticCurve(x^6+1);
        C:=random_walk(C);
        K,rosen:=KummerFromCurve(C);
        C:=HyperellipticCurve(x*(x-1)*(x-rosen[1])*(x-rosen[2])*(x-rosen[3]));
        J:=Jacobian(C);

        poly<X,Y,Z,T>:=PolynomialRing(Fp2,4);
        P := [X,Y,Z,T];

        times := [];
        errors := 0; // counting number of errors we run into 
        for i in [1..num_samples] do
            repeat
                repeat
                    R:=((p+1) div N)*Random(J);
                until N*R eq J!0 and R ne J!0;
                repeat
                    S:=((p+1) div N)*Random(J);
                until N*S eq J!0 and S notin [J!0,R,-R];
            until WeilPairing(R,S,N) eq 1;

            R:=JtoK(R,rosen,K);
            S:=JtoK(S,rosen,K);

            "Sample:", i;
            phi, _ := GetIsogeny(P, R, S, K, N, method : timing := true);

            // Checking if GetImage will work (i.e., not landing on a product of ECs)
            image_thetas := Evaluate(phi, K[2]);
            try 
                image_kummer := KummerFromThetas(image_thetas);
                if &*image_kummer[2] ne 0 then 
                    _, t := GetImage(phi, K : timing:=true);
                    Append(~times, t);
                else 
                    print "Just ran into an errors: the image is likely a product of elliptic curves. Let's skip this sample.";
                    errors +:= 1;
                end if;
            catch e
                print "Just ran into an errors: the image is likely a product of elliptic curves. Let's skip this sample.";
                errors +:= 1;
            end try;

        
            "Average:", &+times/(i-errors), ".\n";
        end for;

        av_time := &+times/(num_samples-errors);
        errors;
        "Average:", av_time, ".\n";
        Append(~av_times, av_time);
    end for;

    return logps, av_times;

end function;


// Pregenerated primes for Figure 2, Table 1 and Table 2 - to reproduce experiments in paper
primes_logp100 := [2^4*3^58*5*29-1, 2^4*3^62*7-1, 2^4*3^55*11-1, 2^4*3^57*13-1, 2^4*3^54*17-1, 2^4*3^61*19-1];


// // Run to get Figure 1
// benchmark_fig1(17, 50);

// // Run to get Figure 2
// benchmark_fig2(100, 50, 50 : Ns := [5, 7, 11, 13, 17, 19] , primes := primes_logp100);

// // Run to get Table 1
// benchmark_table1(100, 50, 50 : Ns := [5, 7, 11, 13, 17, 19] , primes := primes_logp100);

// // Run to get Table 2
// benchmark_getimage(100, 50, 50 : Ns := [7, 11, 13, 17, 19], primes := primes_logp100);
