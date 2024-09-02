/*
    In this file we give an example of how to run the algorithms in this repository
    for computing (N,N)-isogenies between Fast Kummer surfaces

    For simplicity, we choose to work with superspecial Kummer surfaces defined over Fp2
            (indeed, in this way we can force there to be rational N-torsion)
*/

// Choose which N you would like to run
N:=5;

// Set up a prime 
// Here we ensure that we have rational N-torsion by setting p s.t. N | p+1, 
//      and we choose a large-ish prime (80 <= Log_2(p) <= 120)

p := 4;
repeat
    repeat 
        e := Random([20..100]);
        p:=2^4*3^e*N-1;
    until IsPrime(p);
until (Floor(Log(2,p)) in [80..120]);

load "functions.m";

load "isoNN.m";

Fp:=GF(p);
Fp2<i>:=ExtensionField<Fp,x|x^2+1>;
_<x>:=PolynomialRing(Fp2);


// Set up a Fast Kummer surface defined over Fp2
C:=HyperellipticCurve(x^6+1);
C:=random_walk(C);
K,rosen:=KummerFromCurve(C);
C:=HyperellipticCurve(x*(x-1)*(x-rosen[1])*(x-rosen[2])*(x-rosen[3]));
J:=Jacobian(C);


poly<X,Y,Z,T>:=PolynomialRing(Fp2,4);
P := [X,Y,Z,T];

// Generate the N-torsion points that will generate the kernel of the (N,N)-isogeny
"Generating Torsion Points...";
repeat
    repeat 
        R:=((p+1) div N)*Random(J);
    until N*R eq J!0 and R ne J!0;
    repeat
        S:=((p+1) div N)*Random(J);
    until N*S eq J!0 and S notin [J!0,R,-R];
until WeilPairing(R,S,N) eq 1;

// Map these points to the Fast Kummer surface
//     (as N is odd, we can use out JtoK function in 'functions.m')
R:=JtoK(R,rosen,K);
S:=JtoK(S,rosen,K);

"Generated the N torsion points where N is", N, "\n";


// Choose the method for the final scaling
//   - method := 0, this means that no final scaling is applied, and is used when running GetImage
//   - method := 2, method 2 is Scaling_GE as described in the paper 
//   - method := 3, method 3 is Scaling_sqrt as described in the paper 
//   - if N = 5, then we always run Scaling_5 (if method != 0)

method := 2;
if N ne 5 then 
    "Let's run GetIsogeny with method #", method, "!\n";
else
    "Let's run GetIsogeny with usual method for N = 5!\n";
end if;

// Running GetIsogeny
phi, times:= GetIsogeny(P,R, S, K, N, method : timing := true);

"Computed the isogeny in total time ", times[4][1], "seconds.\n";

// Checking the isogeny has the correct image
"Checking the isogeny...";
for i in [1..50] do 
    image_thetas := Evaluate(phi, K[2]);
    image_kummer := KummerFromThetas(image_thetas);

    image_P := Evaluate(phi, RandomKummerPoint(K));

    assert OnKummer(image_P, image_kummer);
end for;

"Isogeny has correct image (with overwhelming probability).\n";

// Running GetImage
"Now let's run GetImage!";
method := 0;
phi, times:= GetIsogeny(P,R, S, K, N, method : timing := true);
t1 := times[4][1];
Kim, t2 := GetImage(phi, K : timing:=true);
"Computed the image in total time ", t1+t2, "seconds.\n";

// Checking the image is correct 
"Checking the image...";
assert Kim eq Parent(Kim)!_computeKummer(image_kummer[1]);
"Image is correct!";

