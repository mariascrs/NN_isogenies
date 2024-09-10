interface(prettyprint, prettyprint=false):
with(linalg):
# We start with y^2 = x*(x-1)*(x-l)*(x-m)*(x-n) and map Kummer
# surface embedding (k1,k2,k3,k4) to (l1(D),l2(D),l3(D),l4(D))
# which simultaneously diagonalises addition by E1 and E2,
# where E1 = {infty,(0,0)} and E2 = {(1,0),(l,0)}.
# For D in Jac(C), let l(D) denote (l1(D),l2(D),l3(D),l4(D)).
#    Here:
# l :=
# (a^2-b^2+c^2-d^2)*(a^2-b^2-c^2+d^2)/(a^2+b^2+c^2+d^2)/(a^2+b^2-c^2-d^2):
# m := (a*c+b*d)*(a^2-b^2+c^2-d^2)/(a*c-b*d)/(a^2+b^2+c^2+d^2):
# n := (a*c+b*d)*(a^2-b^2-c^2+d^2)/(a*c-b*d)/(a^2+b^2-c^2-d^2):
#
# NOTE: What we call E1 in the article is E3 in this file: {(m,0),(n,0)},
#       what we call E2 in the article is E1 in this file: {infty,(0,0)},
#       what we call E3 in the article is E2 in this file: {(1,0),(l,0)}.
#
#   First consider a general curve of genus 2: 
# C: Y^2 = f6*X^6 + f5*X^5 + f4*X^4 + f3*X^3 + f2*X^2 + f1*X + f0.
# For any divisor class represented by {(x,y),(u,v)}
# (shorthand for the divisor class of (x,y) + (u,v) - infty+ - infty-),
# let (k1,k2,k3,k4) be the embedding of the Kummer surface (described,
# for example, in Chapter 3 of Prolegomena to a Middlebrow Arithmetic
# of Curves of Genus 2, by Cassels and Flynn) given by:
# k1 := 1:  k2 := x+u:  k3 := x*u:  k4 := (f0xu-2*y*v)/((x-u)^2): 
# ## where f0xu is:
# f0xu := 2*f0+f1*(x+u)+2*f2*(x*u)+f3*(x+u)*(x*u)
#              +2*f4*(x*u)^2+f5*(x+u)*(x*u)^2+2*f6*(x*u)^3:

# The following is our genus 2 curve y^2 = ccurve.
# We first defined the constants A,B,C,D,E,F,G,H in terms of
# our parameters a,b,c,d. 
unprotect(D):
A := (a^2+b^2+c^2+d^2):
B := (a^2+b^2-c^2-d^2):
C := (a^2-b^2+c^2-d^2):
D := (a^2-b^2-c^2+d^2):
# Note that there seems to be an error in the defn of E in
# Gaudry's article (where (a,b,c,d) is not even on given variety),
# so we are using the defn in the Magma files, which is correct.
E := a*b*c*d*( A*B*C*D/((a^2*d^2-b^2*c^2)
                 *(a^2*c^2-b^2*d^2)*(a^2*b^2-c^2*d^2)) ):
F := (a^4 - b^4 - c^4 + d^4)/(a^2*d^2 - b^2*c^2):
G := (a^4 - b^4 + c^4 - d^4)/(a^2*c^2 - b^2*d^2):
H := (a^4 + b^4 - c^4 - d^4)/(a^2*b^2 - c^2*d^2):
l := C*D/(A*B):
m := (a*c+b*d)*C/((a*c-b*d)*A):
n := (a*c+b*d)*D/((a*c-b*d)*B):

ccurve := x*(x-1)*(x-l)*(x-m)*(x-n):
DQ1 := x:
DQ2 := (x-1)*(x-l):
DQ3 := (x-m)*(x-n):
ccurve := DQ1*DQ2*DQ3:
subf := [ f0 = coeff(ccurve,x,0), f1 = coeff(ccurve,x,1),
         f2 = coeff(ccurve,x,2), f3 = coeff(ccurve,x,3),
         f4 = coeff(ccurve,x,4), f5 = coeff(ccurve,x,5),
         f6 = coeff(ccurve,x,6) ]:

# We now find a linear change of basis on the
# Kummer surface which simultaneously diagonalises addition by E1 and E2,
# where E1 = {infty, (0,0)} and E2 = {(1,0),(l,0)}.
#   We first look at addition by the 2-torsion elements: 
# Let E0 denote the identity element. 
# Recall that the Kummer embedding is: 
# k1 := 1:  k2 := x+u:  k3 := x*u:  k4 := (f0xu-2*y*v)/((x-u)^2):  
# ## where f0xu is: 
# f0xu := 2*f0+f1*(x+u)+2*f2*(x*u)+f3*(x+u)*(x*u) 
# +2*f4*(x*u)^2+f5*(x+u)*(x*u)^2+2*f6*(x*u)^3
# This is the same as: 
kk4 := ( 2*f0 + f1*k2/k1 + 2*f2*(k3/k1) + f3*(k2/k1)*(k3/k1) 
        + 2*f4*(k3/k1)^2 + f5*(k2/k1)*(k3/k1)^2 + 2*f6*(k3/k1)^3
         - 2*yv )/( (k2/k1)^2 - 4*(k3/k1) ):
# The Kummer equation is: 
kummeqn := 
 (k2^2-4*k1*k3)*k4^2
+ (-4*k1^3*f0-2*k1^2*k2*f1-4*k1^2*k3*f2-2*k1*k2*k3*f3-4*k1*k3^2*f4
  -2*k2*k3^2*f5-4*k3^3*f6)*k4
+ (-4*k1^4*f0*f2+k1^4*f1^2-4*k1^3*k2*f0*f3-2*k1^3*k3*f1*f3
   -4*k1^2*k2^2*f0*f4+4*k1^2*k2*k3*f0*f5-4*k1^2*k2*k3*f1*f4
   -4*k1^2*k3^2*f0*f6+2*k1^2*k3^2*f1*f5-4*k1^2*k3^2*f2*f4
   +k1^2*k3^2*f3^2-4*k1*k2^3*f0*f5+8*k1*k2^2*k3*f0*f6
   -4*k1*k2^2*k3*f1*f5+4*k1*k2*k3^2*f1*f6-4*k1*k2*k3^2*f2*f5
   -2*k1*k3^3*f3*f5-4*k2^4*f0*f6-4*k2^3*k3*f1*f6-4*k2^2*k3^2*f2*f6
   -4*k2*k3^3*f3*f6-4*k3^4*f4*f6+k3^4*f5^2):
# The following give the k-coordinates of E0, the identity element,
k1E0 := 0:
k2E0 := 0: 
k3E0 := 0: 
k4E0 := 1:
kE0 := [[k1E0],[k2E0],[k3E0],[k4E0]]:
# We define E1 to be the point of order 2 corresponding to the quadratic
# DQ1 (that is: {(rtDQ1a,0),(rtDQ1b,0)}, where rtDQ1a,rtDQ1b are the roots
# of DQ1). The following is (k1(E1),k2(E1),k3(E1),k4(E1)).
k1E1 := coeff(DQ1,x,2):  
k2E1 := -coeff(DQ1,x,1):  
k3E1 := coeff(DQ1,x,0):
if k1E1 = 0 then k4E1 := subs(op(subf),f5*k3E1^2/k2E1) else
k4E1 := factor(k1E1*subs(op(subf), k1=k1E1, k2=k2E1, k3=k3E1, yv = 0, kk4)) fi:
check := factor(subs(op(subf), k1=k1E1, k2=k2E1, k3=k3E1, k4=k4E1, kummeqn));
kE1 := [[k1E1],[k2E1],[k3E1],[k4E1]]:
# E2 corresponds to DQ2 in the same way. 
k1E2 := coeff(DQ2,x,2):  
k2E2 := -coeff(DQ2,x,1):  
k3E2 := coeff(DQ2,x,0):
k4E2 := factor(k1E2*subs(op(subf), k1=k1E2, k2=k2E2, k3=k3E2, yv = 0, kk4)):
check := factor(subs(op(subf), k1=k1E2, k2=k2E2, k3=k3E2, k4=k4E2, kummeqn));
kE2 := [[k1E2],[k2E2],[k3E2],[k4E2]]:
# E3 corresponds to DQ3 in the same way. 
k1E3 := coeff(DQ3,x,2):  
k2E3 := -coeff(DQ3,x,1):  
k3E3 := coeff(DQ3,x,0):
k4E3 := factor(k1E3*subs(op(subf), k1=k1E3, k2=k2E3, k3=k3E3, yv = 0, kk4)):
check := factor(subs(op(subf), k1=k1E3, k2=k2E3, k3=k3E3, k4=k4E3, kummeqn));
kE3 := [[k1E3],[k2E3],[k3E3],[k4E3]]:
# These are the roots of the sextic: 
# roots of DQ1:
rtDQ1a := infty:
rtDQ1b := 0:
# roots of DQ2:
rtDQ2a := 1:
rtDQ2b := l:
# roots of DQ3:
rtDQ3a := m:
rtDQ3b := n:
# Using this notation:
# E1 = {(rtDQ1a,0),(rtDQ1b,0)},
# E2 = {(rtDQ2a,0),(rtDQ2b,0)},
# E3 = {(rtDQ3a,0),(rtDQ3b,0)}.
#
# Now find the Kummer coordinates of the other 12 points of order 2. 
# DQ4 := (x-rtDQ1a)*(x-rtDQ2a), but rtDQ1a is infty, so in fact:
DQ4 := (x-rtDQ2a):
k1E4 := coeff(DQ4,x,2):  
k2E4 := -coeff(DQ4,x,1):  
k3E4 := coeff(DQ4,x,0):
if k1E4 = 0 then k4E4 := subs(op(subf),f5*k3E4^2/k2E4) else
k4E4 := factor(k1E4*subs(op(subf), k1=k1E4, k2=k2E4, k3=k3E4, yv = 0, kk4)) fi:
# The following checks that (k1E4,k2E4,k3E4,k4E4) satisfy the equation
# of the Kummer surface.
check := 
simplify(factor(subs(op(subf), k1=k1E4, k2=k2E4, k3=k3E4, k4=k4E4, kummeqn)));
kE4 := [[k1E4],[k2E4],[k3E4],[k4E4]]:
# DQ5 := (x-rtDQ1a)*(x-rtDQ2b), but rtDQ1a = infty, so in fact:
DQ5 := (x-rtDQ2b):
k1E5 := coeff(DQ5,x,2):  
k2E5 := -coeff(DQ5,x,1):  
k3E5 := coeff(DQ5,x,0):
if k1E5 = 0 then k4E5 := subs(op(subf),f5*k3E5^2/k2E5) else
k4E5 := factor(k1E5*subs(op(subf), k1=k1E5, k2=k2E5, k3=k3E5, yv = 0, kk4)) fi:
check := 
simplify(factor(subs(op(subf), k1=k1E5, k2=k2E5, k3=k3E5, k4=k4E5, kummeqn)));
kE5 := [[k1E5],[k2E5],[k3E5],[k4E5]]:
# DQ6 := (x-rtDQ1a)*(x-rtDQ3a), but rtDQ1a = infty, so in fact:
DQ6 := (x-rtDQ3a):
k1E6 := coeff(DQ6,x,2):  
k2E6 := -coeff(DQ6,x,1):  
k3E6 := coeff(DQ6,x,0):
if k1E6 = 0 then k4E6 := subs(op(subf),f5*k3E6^2/k2E6) else
k4E6 := factor(k1E6*subs(op(subf), k1=k1E6, k2=k2E6, k3=k3E6, yv = 0, kk4)) fi:
check := 
simplify(factor(subs(op(subf),k1=k1E6,k2=k2E6,k3=k3E6,k4=k4E6,kummeqn)));
kE6 := [[k1E6],[k2E6],[k3E6],[k4E6]]:
# DQ7 := (x-rtDQ1a)*(x-rtDQ3b), but rtDQ1a = infty, so in fact:
DQ7 := (x-rtDQ3b):
k1E7 := coeff(DQ7,x,2):
k2E7 := -coeff(DQ7,x,1):
k3E7 := coeff(DQ7,x,0):
if k1E7 = 0 then k4E7 := subs(op(subf),f5*k3E7^2/k2E7) else
k4E7 := factor(k1E7*subs(op(subf), k1=k1E7, k2=k2E7, k3=k3E7, yv = 0, kk4)) fi:
check :=
simplify(factor(subs(op(subf), k1=k1E7,k2=k2E7,k3=k3E7,k4=k4E7, kummeqn)));
kE7 := [[k1E7],[k2E7],[k3E7],[k4E7]]:
DQ8 := (x-rtDQ1b)*(x-rtDQ2a):
k1E8 := coeff(DQ8,x,2):
k2E8 := -coeff(DQ8,x,1):
k3E8 := coeff(DQ8,x,0):
k4E8 := factor(k1E8*subs(op(subf), k1=k1E8, k2=k2E8, k3=k3E8, yv = 0, kk4)):
check :=
simplify(factor(subs(op(subf), k1=k1E8, k2=k2E8, k3=k3E8, k4=k4E8, kummeqn)));
kE8 := [[k1E8],[k2E8],[k3E8],[k4E8]]:
DQ9 := (x-rtDQ1b)*(x-rtDQ2b):
k1E9 := coeff(DQ9,x,2):
k2E9 := -coeff(DQ9,x,1):
k3E9 := coeff(DQ9,x,0):
k4E9 := factor(k1E9*subs(op(subf), k1=k1E9, k2=k2E9, k3=k3E9, yv = 0, kk4)):
check :=
simplify(factor(subs(op(subf), k1=k1E9, k2=k2E9, k3=k3E9, k4=k4E9, kummeqn)));
kE9 := [[k1E9],[k2E9],[k3E9],[k4E9]]:
DQ10 := (x-rtDQ1b)*(x-rtDQ3a):
k1E10 := coeff(DQ10,x,2):
k2E10 := -coeff(DQ10,x,1):
k3E10 := coeff(DQ10,x,0):
k4E10 := factor(k1E10*subs(op(subf), k1=k1E10, k2=k2E10, k3=k3E10,yv = 0,kk4)):
check :=
simplify(factor(subs(op(subf),k1=k1E10,k2=k2E10,k3=k3E10,k4=k4E10,kummeqn)));
kE10 := [[k1E10],[k2E10],[k3E10],[k4E10]]:
DQ11 := (x-rtDQ1b)*(x-rtDQ3b):
k1E11 := coeff(DQ11,x,2):
k2E11 := -coeff(DQ11,x,1):
k3E11 := coeff(DQ11,x,0):
k4E11 := factor(k1E11*subs(op(subf), k1=k1E11, k2=k2E11, k3=k3E11,yv = 0,kk4)):
check :=
simplify(factor(subs(op(subf),k1=k1E11,k2=k2E11,k3=k3E11,k4=k4E11,kummeqn)));
kE11 := [[k1E11],[k2E11],[k3E11],[k4E11]]:
DQ12 := (x-rtDQ2a)*(x-rtDQ3a):
k1E12 := coeff(DQ12,x,2):
k2E12 := -coeff(DQ12,x,1):
k3E12 := coeff(DQ12,x,0):
k4E12 := factor(k1E12*subs(op(subf), k1=k1E12, k2=k2E12, k3=k3E12,yv = 0,kk4)):
check :=
simplify(factor(subs(op(subf),k1=k1E12,k2=k2E12,k3=k3E12,k4=k4E12,kummeqn)));
kE12 := [[k1E12],[k2E12],[k3E12],[k4E12]]:
DQ13 := (x-rtDQ2a)*(x-rtDQ3b):
k1E13 := coeff(DQ13,x,2):
k2E13 := -coeff(DQ13,x,1):
k3E13 := coeff(DQ13,x,0):
k4E13 := factor(k1E13*subs(op(subf), k1=k1E13, k2=k2E13, k3=k3E13,yv = 0,kk4)):
check :=
simplify(factor(subs(op(subf),k1=k1E13,k2=k2E13,k3=k3E13,k4=k4E13,kummeqn)));
kE13 := [[k1E13],[k2E13],[k3E13],[k4E13]]:
DQ14 := (x-rtDQ2b)*(x-rtDQ3a):
k1E14 := coeff(DQ14,x,2):
k2E14 := -coeff(DQ14,x,1):
k3E14 := coeff(DQ14,x,0):
k4E14 := factor(k1E14*subs(op(subf), k1=k1E14, k2=k2E14, k3=k3E14,yv = 0,kk4)):
check :=
simplify(factor(subs(op(subf),k1=k1E14,k2=k2E14,k3=k3E14,k4=k4E14,kummeqn)));
kE14 := [[k1E14],[k2E14],[k3E14],[k4E14]]:
DQ15 := (x-rtDQ2b)*(x-rtDQ3b):
k1E15 := coeff(DQ15,x,2):
k2E15 := -coeff(DQ15,x,1):
k3E15 := coeff(DQ15,x,0):
k4E15 := factor(k1E15*subs(op(subf), k1=k1E15, k2=k2E15, k3=k3E15,yv = 0,kk4)):
check :=
simplify(factor(subs(op(subf),k1=k1E15,k2=k2E15,k3=k3E15,k4=k4E15,kummeqn)));
kE15 := [[k1E15],[k2E15],[k3E15],[k4E15]]:
# We now set up the linear maps given by addition by the points of order 2.
# The following is the known formula (see pp.21,22 of Cassels and Flynn) 
# for y^2 = (g2*X^2 + g1*X + g0)*(h4*X^4 + h3*X^3 + h2*X^2 + h1*X + h0)
# which gives the linear map for addition by the order-2 point  
# corresponding to g2*X^2 + g1*X + g0. 
Wg := matrix(4,4):
Wg[1,1] := g2^2*h0+g0*g2*h2-g0^2*h4:
Wg[1,2] := g0*g2*h3-g0*g1*h4:
Wg[1,3] := g1*g2*h3-g1^2*h4+2*g0*g2*h4:
Wg[1,4] := g2:
Wg[2,1] := -g0*g2*h1-g0*g1*h2+g0^2*h3:
Wg[2,2] := g2^2*h0-g0*g2*h2+g0^2*h4:
Wg[2,3] := g2^2*h1-g1*g2*h2-g0*g2*h3:
Wg[2,4] := -g1:
Wg[3,1] := -g1^2*h0+2*g0*g2*h0+g0*g1*h1:
Wg[3,2] := -g1*g2*h0+g0*g2*h1:
Wg[3,3] := -g2^2*h0+g0*g2*h2+g0^2*h4:
Wg[3,4] := g0:
Wg[4,1] := -g1*g2^2*h0*h1+g1^2*g2*h0*h2+g0*g2^2*h1^2-4*g0*g2^2*h0*h2-
g0*g1*g2*h1*h2+g0*g1*g2*h0*h3-g0^2*g2*h1*h3:
Wg[4,2] := g1^2*g2*h0*h3-g1^3*h0*h4-2*g0*g2^2*h0*h3-g0*g1*g2*h1*h3+
4*g0*g1*g2*h0*h4+g0*g1^2*h1*h4-2*g0^2*g2*h1*h4:
Wg[4,3] := -g0*g2^2*h1*h3-g0*g1*g2*h2*h3+g0*g1*g2*h1*h4+g0*g1^2*h2*h4+
g0^2*g2*h3^2-4*g0^2*g2*h2*h4-g0^2*g1*h3*h4:
Wg[4,4] := -g2^2*h0-g0*g2*h2-g0^2*h4:
# 
# We now specialise this for each point of order 2. 
# The following array gives all possible quadratic factors of ccurve
# DQ1,..,DQ15, which correponding to the points of order 2: E1,..,E15. 
DQarray := [DQ1,DQ2,DQ3,DQ4,DQ5,DQ6,DQ7,DQ8,DQ9,DQ10,DQ11,DQ12,DQ13,DQ14,DQ15]:
# The following is the array of k-coordinates of the 2-torsion elements. 
EKarray := [ [[k1E0],[k2E0],[k3E0],[k4E0]],
             [[k1E1],[k2E1],[k3E1],[k4E1]],
             [[k1E2],[k2E2],[k3E2],[k4E2]],
             [[k1E3],[k2E3],[k3E3],[k4E3]],
             [[k1E4],[k2E4],[k3E4],[k4E4]],
             [[k1E5],[k2E5],[k3E5],[k4E5]],
             [[k1E6],[k2E6],[k3E6],[k4E6]],
             [[k1E7],[k2E7],[k3E7],[k4E7]],
             [[k1E8],[k2E8],[k3E8],[k4E8]],
             [[k1E9],[k2E9],[k3E9],[k4E9]],
             [[k1E10],[k2E10],[k3E10],[k4E10]],
             [[k1E11],[k2E11],[k3E11],[k4E11]],
             [[k1E12],[k2E12],[k3E12],[k4E12]],
             [[k1E13],[k2E13],[k3E13],[k4E13]],
             [[k1E14],[k2E14],[k3E14],[k4E14]],
             [[k1E15],[k2E15],[k3E15],[k4E15]] ]:
WE := [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]:
# The following uses the above Wg matrix, which described in general
# addition by a point of order 2, and we apply this for each
# of E1,..,E15, giving: WE, which is an array of matrices, such
# that WE[i] is the matrix representing addition by Ei. 
for i from 1 to 15 do
 gtemp := DQarray[i]: htemp := simplify(ccurve/gtemp):
 WEitemp := [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]:
 for j from 1 to 4 do for k from 1 to 4 do
  WEitemp[j][k] := 
    subs(g0=coeff(gtemp,x,0),g1=coeff(gtemp,x,1),g2=coeff(gtemp,x,2),
         h0=coeff(htemp,x,0),h1=coeff(htemp,x,1),h2=coeff(htemp,x,2),
         h3=coeff(htemp,x,3),h4=coeff(htemp,x,4), Wg[j,k]):
  WE[i] := [ WEitemp ]:
od od od:
# These following give the matrices of the linear maps for addition 
# by the 2-torsion points: E0,...,E15. 
WE0 := matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]):
WE1 := matrix(WE[1][1]):
WE2 := matrix(WE[2][1]):
WE3 := matrix(WE[3][1]):
WE4 := matrix(WE[4][1]):
WE5 := matrix(WE[5][1]):
WE6 := matrix(WE[6][1]):
WE7 := matrix(WE[7][1]):
WE8 := matrix(WE[8][1]):
WE9 := matrix(WE[9][1]):
WE10 := matrix(WE[10][1]):
WE11 := matrix(WE[11][1]):
WE12 := matrix(WE[12][1]):
WE13 := matrix(WE[13][1]):
WE14 := matrix(WE[14][1]):
WE15 := matrix(WE[15][1]):
WEarray := [WE0,WE1,WE2,WE3,WE4,WE5,WE6,WE7,
            WE8,WE9,WE10,WE11,WE12,WE13,WE14,WE15]:

# We now simultaneously diagonalise addition by E1,E2 (which will
# also diagonalise addition by E3, since E3 = E1 + E2).
# First find the eigenvectors of the matrices representing addition
# by E1,E2,E3.
evWE1 := eigenvects(WE1):
evWE2 := eigenvects(WE2):
evWE3 := eigenvects(WE3):
# This gives that WE1 has eigenvalue:
evalue1WE1 := (a*c+b*d)*C*D/((a*c-b*d)*A*B):
# with 2-dimensional eigenspace spanned by:
ev1aWE1 := matrix( [[-(a*c-b*d)*A*B/((a*c+b*d)*C*D)], [0], [1], [0]] ):
# The following double checks this:
check := factor(multiply(WE1,ev1aWE1)[1,1] - evalue1WE1*ev1aWE1[1,1]);
ev1bWE1 := matrix( [[0], [-(a*c-b*d)*A*B/((a*c+b*d)*C*D)], [0], [1]] ):
# The following double checks this:
check := factor(multiply(WE1,ev1bWE1)[1,1] - evalue1WE1*ev1bWE1[1,1]);
# and eigenvalue:
evalue2WE1 := -(a*c+b*d)*C*D/((a*c-b*d)*A*B):
# with 2-dimensional eigenspace spanned by:
ev2aWE1 := matrix( [[(a*c-b*d)*A*B/((a*c+b*d)*C*D)], [0], [1], [0]] ):
# The following double checks this:
check := factor(multiply(WE1,ev2aWE1)[1,1] - evalue2WE1*ev2aWE1[1,1]);
ev2bWE1 := matrix( [[0], [(a*c-b*d)*A*B/((a*c+b*d)*C*D)], [0], [1]] ):
# The following double checks this:
check := factor(multiply(WE1,ev2bWE1)[1,1] - evalue2WE1*ev2bWE1[1,1]);
#
# The above eigenvector command also gave that WE2 has eigenvalue:
evalue1WE2 := 4*(a*d-b*c)*(a*d+b*c)*(a*b-c*d)*(a*b+c*d)*C*D/(A*B*(a*c-b*d))^2:
# with 2-dimensional eigenspace spanned by:
ev1aWE2 := matrix( [[1], 
 [(a^5*d+2*a^2*b^3*c+a*b^4*d-a*c^4*d-a*d^5-2*b*c^3*d^2)/(a*d*A*B)],
                    [0], 
 [(a^5*d-2*a^2*b^3*c+a*b^4*d-a*c^4*d-a*d^5+2*b*c^3*d^2)*C*D*(a*c+b*d)^2
     /(a*d*A^2*B^2*(a*c-b*d)^2)]] ):
# The following double checks this:
check := factor(multiply(WE2,ev1aWE2)[1,1] - evalue1WE2*ev1aWE2[1,1]);

ev1bWE2 := matrix( [[0], 
 [(a^5*d-2*a^2*b^3*c+a*b^4*d-a*c^4*d-a*d^5+2*b*c^3*d^2)/(a*d*C*D)], 
                    [1], 
 [(a^5*d+2*a^2*b^3*c+a*b^4*d-a*c^4*d-a*d^5-2*b*c^3*d^2)/(a*d*A*B)]] ): 
# The following double checks this:
check := factor(multiply(WE2,ev1bWE2)[1,1] - evalue1WE2*ev1bWE2[1,1]);

# and eigenvalue:
evalue2WE2 := -4*(a*d-b*c)*(a*d+b*c)*(a*b-c*d)*(a*b+c*d)*C*D/(A*B*(a*c-b*d))^2:
# with 2-dimensional eigenspace spanned by:
ev2aWE2 := matrix( [[1], 
   [(a^4*b*c+2*a^3*b^2*d-2*a*c^2*d^3+b^5*c-b*c^5-b*c*d^4)/(b*c*A*B)], 
                    [0], 
  [(a^4*b*c-2*a^3*b^2*d+2*a*c^2*d^3+b^5*c-b*c^5-b*c*d^4)*C*D*(a*c+b*d)^2
   /(b*c*(a*c-b*d)^2*A^2*B^2)]] ):
# The following double checks this:
check := factor(multiply(WE2,ev2aWE2)[1,1] - evalue2WE2*ev2aWE2[1,1]);
ev2bWE2 := matrix( [[0], 
 [(a^4*b*c-2*a^3*b^2*d+2*a*c^2*d^3+b^5*c-b*c^5-b*c*d^4)/(b*c*C*D)], 
                    [1], 
 [(a^4*b*c+2*a^3*b^2*d-2*a*c^2*d^3+b^5*c-b*c^5-b*c*d^4)/(b*c*A*B)]] ):
# The following double checks this:
check := factor(multiply(WE2,ev2bWE2)[1,1] - evalue2WE2*ev2bWE2[1,1]);
# We now want the four simultaneous eigenvectors for WE1,WE2.
# Each 2-dim espace of WE1 and 2-dim e-space of WE2 will have a common vector.
# First find the vector which has evalue1WE1 for WE1 and evalue1WE2 for WE2. 
solve( { 
    rr*ev1aWE1[1,1] + ss*ev1bWE1[1,1] = rrr*ev1aWE2[1,1] + sss*ev1bWE2[1,1],
    rr*ev1aWE1[2,1] + ss*ev1bWE1[2,1] = rrr*ev1aWE2[2,1] + sss*ev1bWE2[2,1],
    rr*ev1aWE1[3,1] + ss*ev1bWE1[3,1] = rrr*ev1aWE2[3,1] + sss*ev1bWE2[3,1],
    rr*ev1aWE1[4,1] + ss*ev1bWE1[4,1] = rrr*ev1aWE2[4,1] + sss*ev1bWE2[4,1] },
         {rr,ss,rrr,sss} ):
# This gave: {rr = -1/2*d*(a*c-b*d)*A*B*ss
#                   /b/(a^4*d^2-2*a^2*b^2*c^2+b^4*d^2+c^4*d^2-d^6), 
# rrr = 1/2*A^2*B^2*(a*c-b*d)^2*d*ss
# /b/(a^4*d^2-2*a^2*b^2*c^2+b^4*d^2+c^4*d^2-d^6)/(a*c+b*d)/C/D, 
# ss = ss, 
# sss = -1/2*d*(a*c-b*d)*A*B*ss/b/(a^4*d^2-2*a^2*b^2*c^2+b^4*d^2+c^4*d^2-d^6)}
common4 := matadd( scalarmul(ev1aWE1,-1/2*d*(a*c-b*d)*A*B), 
          scalarmul(ev1bWE1,b*(a^4*d^2-2*a^2*b^2*c^2+b^4*d^2+c^4*d^2-d^6)) ):
common4[1,1] := factor(common4[1,1]):
common4[2,1] := factor(common4[2,1]):
common4[3,1] := factor(common4[3,1]):
common4[4,1] := factor(common4[4,1]):
# The above gives the following result:
common4 := matrix([[1/2*d*A^2*B^2*(a*c-b*d)^2/((a*c+b*d)*C*D)],
[-b*(a^4*d^2-2*a^2*b^2*c^2+b^4*d^2+c^4*d^2-d^6)*(a*c-b*d)*A*B/((a*c+b*d)*C*D)],
 [-1/2*d*(a*c-b*d)*A*B],
 [b*(a^4*d^2-2*a^2*b^2*c^2+b^4*d^2+c^4*d^2-d^6)]]):
# Multiply common4 through by 4*(a*c+b*d)*C*D/(A^2*B^2*(a*c-b*d)^2)
for i from 1 to 4 do common4[i,1] := 
  factor(common4[i,1]*4*(a*c+b*d)*C*D/(A^2*B^2*(a*c-b*d)^2)) od:
# The above gives the following result:
common4 := matrix([[2*d], 
  [-4*b*(a^4*d^2-2*a^2*b^2*c^2+b^4*d^2+c^4*d^2-d^6)/((a*c-b*d)*A*B)], 
  [-2*d*(a*c+b*d)*C*D/((a*c-b*d)*A*B)], 
  [4*b*(a^4*d^2-2*a^2*b^2*c^2+b^4*d^2+c^4*d^2-d^6)*(a*c+b*d)*C*D
     /((a*c-b*d)^2*A^2*B^2)]]):
# The following checks that common4 is indeed vector which 
# has evalue1WE1 for WE1 and evalue1WE2 for WE2. 
check := factor(multiply(WE1,common4)[1,1] - evalue1WE1*common4[1,1]);
check := factor(multiply(WE1,common4)[2,1] - evalue1WE1*common4[2,1]);
check := factor(multiply(WE1,common4)[3,1] - evalue1WE1*common4[3,1]);
check := factor(multiply(WE1,common4)[4,1] - evalue1WE1*common4[4,1]);
check := factor(multiply(WE2,common4)[1,1] - evalue1WE2*common4[1,1]);
check := factor(multiply(WE2,common4)[2,1] - evalue1WE2*common4[2,1]);
check := factor(multiply(WE2,common4)[3,1] - evalue1WE2*common4[3,1]);
check := factor(multiply(WE2,common4)[4,1] - evalue1WE2*common4[4,1]);

# Now find the vector which has evalue1WE1 for WE1 and evalue2WE2 for WE2. 
solve( { 
  rr*ev1aWE1[1,1] + ss*ev1bWE1[1,1] = rrr*ev2aWE2[1,1] + sss*ev2bWE2[1,1],
  rr*ev1aWE1[2,1] + ss*ev1bWE1[2,1] = rrr*ev2aWE2[2,1] + sss*ev2bWE2[2,1],
  rr*ev1aWE1[3,1] + ss*ev1bWE1[3,1] = rrr*ev2aWE2[3,1] + sss*ev2bWE2[3,1],
  rr*ev1aWE1[4,1] + ss*ev1bWE1[4,1] = rrr*ev2aWE2[4,1] + sss*ev2bWE2[4,1] },
         {rr,ss,rrr,sss} ):
# This gives
# rr = 1/2*b*(a*c-b*d)*A*B*ss/d/(a^4*b^2-2*a^2*c^2*d^2-b^6+b^2*c^4+b^2*d^4), 
# rrr = -1/2*A^2*B^2*(a*c-b*d)^2*b*ss
#       /d/(a^4*b^2-2*a^2*c^2*d^2-b^6+b^2*c^4+b^2*d^4)/(a*c+b*d)/C/D, 
# ss = ss, 
# sss = 1/2*b*(a*c-b*d)*A*B*ss/d/(a^4*b^2-2*a^2*c^2*d^2-b^6+b^2*c^4+b^2*d^4)
common2 := matadd( scalarmul(ev1aWE1,1/2*b*(a*c-b*d)*A*B),
          scalarmul(ev1bWE1,d*(a^4*b^2-2*a^2*c^2*d^2-b^6+b^2*c^4+b^2*d^4)) ):
# The above gives the following result:
common2 := matrix([[-1/2*b*A^2*B^2*(a*c-b*d)^2/((a*c+b*d)*C*D)],
[-d*(a^4*b^2-2*a^2*c^2*d^2-b^6+b^2*c^4+b^2*d^4)*(a*c-b*d)*A*B/((a*c+b*d)*C*D)],
  [1/2*b*(a*c-b*d)*A*B],
  [d*(a^4*b^2-2*a^2*c^2*d^2-b^6+b^2*c^4+b^2*d^4)]]):
# Multiply common2 through by 4*(a*c+b*d)*C*D/(A^2*B^2*(a*c-b*d)^2).
for i from 1 to 4 do common2[i,1] := 
  factor(common2[i,1]*4*(a*c+b*d)*C*D/(A^2*B^2*(a*c-b*d)^2)) od:
# The above gives the following result:
common2 := matrix([[-2*b], 
 [-4*d*(a^4*b^2-2*a^2*c^2*d^2-b^6+b^2*c^4+b^2*d^4)/((a*c-b*d)*A*B)], 
 [2*b*(a*c+b*d)*C*D/((a*c-b*d)*A*B)], 
 [4*d*(a^4*b^2-2*a^2*c^2*d^2-b^6+b^2*c^4+b^2*d^4)*(a*c+b*d)*C*D
   /((a*c-b*d)^2*A^2*B^2)]]):
# The following checks that common2 is the vector which 
# has evalue1WE1 for WE1 and evalue2WE2 for WE2. 
check := factor(multiply(WE1,common2)[1,1] - evalue1WE1*common2[1,1]);
check := factor(multiply(WE1,common2)[2,1] - evalue1WE1*common2[2,1]);
check := factor(multiply(WE1,common2)[3,1] - evalue1WE1*common2[3,1]);
check := factor(multiply(WE1,common2)[4,1] - evalue1WE1*common2[4,1]);
check := factor(multiply(WE2,common2)[1,1] - evalue2WE2*common2[1,1]);
check := factor(multiply(WE2,common2)[2,1] - evalue2WE2*common2[2,1]);
check := factor(multiply(WE2,common2)[3,1] - evalue2WE2*common2[3,1]);
check := factor(multiply(WE2,common2)[4,1] - evalue2WE2*common2[4,1]);

# Now find the vector which has evalue2WE1 for WE1 and evalue1WE2 for WE2.
solve( { 
  rr*ev2aWE1[1,1] + ss*ev2bWE1[1,1] = rrr*ev1aWE2[1,1] + sss*ev1bWE2[1,1],
  rr*ev2aWE1[2,1] + ss*ev2bWE1[2,1] = rrr*ev1aWE2[2,1] + sss*ev1bWE2[2,1],
  rr*ev2aWE1[3,1] + ss*ev2bWE1[3,1] = rrr*ev1aWE2[3,1] + sss*ev1bWE2[3,1],
  rr*ev2aWE1[4,1] + ss*ev2bWE1[4,1] = rrr*ev1aWE2[4,1] + sss*ev1bWE2[4,1] },
         {rr,ss,rrr,sss} ):
# This gives
# {rr = 1/2*a*(a*c-b*d)*A*B*ss/(c*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2)),
# rrr = 1/2*A^2*B^2*(a*c-b*d)^2*a*ss
#       /(c*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2)*(a*c+b*d)*C*D),
# ss = ss,
# sss = 1/2*a*(a*c-b*d)*A*B*ss/(c*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2))}
common1 := 
  matadd(scalarmul(ev2aWE1,
     1/2*a*(a*c-b*d)*A*B/(c*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2))),
           scalarmul(ev2bWE1,1)):
# The above gives the following result:
common1 := matrix( 
 [[1/2*a*(a*c-b*d)^2*A^2*B^2
   /(c*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2)*(a*c+b*d)*C*D)],
  [(a*c-b*d)*A*B/((a*c+b*d)*C*D)], 
  [1/2*a*(a*c-b*d)*A*B/(c*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2))],
  [1]] ):
# Multiply common2 through by 
# 4*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2)*c*(a*c+b*d)*C*D
#     /((a*c-b*d)^2*A^2*B^2)
for i from 1 to 4 do common1[i,1] := 
  factor(common1[i,1]
       *4*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2)*c*(a*c+b*d)*C*D
       /((a*c-b*d)^2*A^2*B^2)) od:
# The above gives the following result:
common1 := matrix([[2*a], 
 [4*c*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2)/((a*c-b*d)*A*B)], 
 [2*a*(a*c+b*d)*C*D/((a*c-b*d)*A*B)],
 [4*c*(a*c+b*d)*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2)*C*D
   /((a*c-b*d)^2*A^2*B^2)]]):
# The following checks that common1 is the vector which 
# has evalue2WE1 for WE1 and evalue1WE2 for WE2.
check := factor(multiply(WE1,common1)[1,1] - evalue2WE1*common1[1,1]);
check := factor(multiply(WE1,common1)[2,1] - evalue2WE1*common1[2,1]);
check := factor(multiply(WE1,common1)[3,1] - evalue2WE1*common1[3,1]);
check := factor(multiply(WE1,common1)[4,1] - evalue2WE1*common1[4,1]);
check := factor(multiply(WE2,common1)[1,1] - evalue1WE2*common1[1,1]);
check := factor(multiply(WE2,common1)[2,1] - evalue1WE2*common1[2,1]);
check := factor(multiply(WE2,common1)[3,1] - evalue1WE2*common1[3,1]);
check := factor(multiply(WE2,common1)[4,1] - evalue1WE2*common1[4,1]);
# Now find the vector which has evalue2WE1 for WE1 and evalue2WE2 for WE2.
solve( { 
  rr*ev2aWE1[1,1] + ss*ev2bWE1[1,1] = rrr*ev2aWE2[1,1] + sss*ev2bWE2[1,1],
  rr*ev2aWE1[2,1] + ss*ev2bWE1[2,1] = rrr*ev2aWE2[2,1] + sss*ev2bWE2[2,1],
  rr*ev2aWE1[3,1] + ss*ev2bWE1[3,1] = rrr*ev2aWE2[3,1] + sss*ev2bWE2[3,1],
  rr*ev2aWE1[4,1] + ss*ev2bWE1[4,1] = rrr*ev2aWE2[4,1] + sss*ev2bWE2[4,1] },
         {rr,ss,rrr,sss} ):
# This gives:
# {rr = 1/2*(a*c-b*d)*A*B*c*ss/(a*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4)),
# rrr = 1/2*A^2*B^2*(a*c-b*d)^2*c*ss
#   /(a*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4)*(a*c+b*d)*C*D),
# ss = ss, 
# sss = 1/2*(a*c-b*d)*A*B*c*ss/(a*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4))}
common3 := 
  matadd(scalarmul(ev2aWE1,
    1/2*(a*c-b*d)*A*B*c/(a*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4))),
           scalarmul(ev2bWE1,1)):
# The above gives the following result:
common3 := matrix( 
 [[1/2*c*(a*c-b*d)^2*A^2*B^2
     /(a*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4)*(a*c+b*d)*C*D)], 
 [(a*c-b*d)*A*B/((a*c+b*d)*(a^2-b^2+c^2-d^2)*(a^2-b^2-c^2+d^2))], 
 [1/2*c*(a*c-b*d)*A*B/(a*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4))],
 [1]] ):
# Multiply common3 through by 
# -4*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4)*a*(a*c+b*d)*C*D
#  /((a*c-b*d)^2*A^2*B^2)
for i from 1 to 4 do common3[i,1] := 
  factor(common3[i,1]
     *(-4*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4)*a*(a*c+b*d)*C*D
      /((a*c-b*d)^2*A^2*B^2))) od:
# The above gives the following result:
common3 := matrix([[-2*c], 
 [-4*a*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4)/((a*c-b*d)*A*B)],
 [-2*c*(a*c+b*d)*C*D/((a*c-b*d)*A*B)], 
 [-4*a*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4)*(a*c+b*d)*C*D
   /((a*c-b*d)^2*A^2*B^2)]]):
# The following checks that common1 is indeed the vector which 
# has evalue2WE1 for WE1 and evalue2WE2 for WE2.
check := factor(multiply(WE1,common3)[1,1] - evalue2WE1*common3[1,1]);
check := factor(multiply(WE1,common3)[2,1] - evalue2WE1*common3[2,1]);
check := factor(multiply(WE1,common3)[3,1] - evalue2WE1*common3[3,1]);
check := factor(multiply(WE1,common3)[4,1] - evalue2WE1*common3[4,1]);
check := factor(multiply(WE2,common3)[1,1] - evalue2WE2*common3[1,1]);
check := factor(multiply(WE2,common3)[2,1] - evalue2WE2*common3[2,1]);
check := factor(multiply(WE2,common3)[3,1] - evalue2WE2*common3[3,1]);
check := factor(multiply(WE2,common3)[4,1] - evalue2WE2*common3[4,1]);
# Now put the columns of common1,common2,common3,common4 into Pdiag,
# which gives the change of basis between our old Kummer coordiates
# and our new ones. Note that there are a number of possibilities,
# given that there are 16 different linear maps from K^gen to Gaudry
# (since one can compose with any of the linear maps given by
# addition of a point of order 2), but the following order of columns
# and scalar multiples of columns are chosen so that the identity
# maps to (X,Y,Z,T) = (a,b,c,d).
Pdiag := matrix( 
  [ [common3[1,1], common4[1,1], common1[1,1], common2[1,1]],
    [common3[2,1], common4[2,1], common1[2,1], common2[2,1]],
    [common3[3,1], common4[3,1], common1[3,1], common2[3,1]],
    [common3[4,1], common4[4,1], common1[4,1], common2[4,1]] ]):
# Let us now find what addition by our points of order 2 look
# like after this change of basis. We let WE0diag, WE1diag,...
# denote the matrices representing addition by E0,E1,...
# on the Kummer surface with respect to our new coordinates.
# Of course, addition by the identity is still given by the
# identity matrix:
WE0diag := matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]):
# The following is addition by E1:
WE1diag := multiply( inverse(Pdiag), multiply(WE1,Pdiag) ): 
for i from 1 to 4 do for j from 1 to 4 do
  WE1diag[i,j] := -factor(WE1diag[i,j]/evalue1WE1) od od:
# This gives:
# WE1diag := 
#   matrix([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]):
# which, reassuringly, is indeed diagonalised!
WE2diag := multiply( inverse(Pdiag), multiply(WE2,Pdiag) ): 
for i from 1 to 4 do for j from 1 to 4 do
  WE2diag[i,j] := -factor(WE2diag[i,j]/evalue1WE2) od od:
# This gives:
# WE2diag :=
#   matrix([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]]):
WE3diag := multiply( inverse(Pdiag), multiply(WE3,Pdiag) ): 
for i from 1 to 4 do for j from 1 to 4 do
  WE3diag[i,j] := factor(WE3diag[i,j]) od od:
const := WE3diag[1,1]: 
for i from 1 to 4 do for j from 1 to 4 do
  WE3diag[i,j] := factor(WE3diag[i,j]/const) od od:
# This gives:
# WE3diag :=
#   matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]):
# 
# Of course, addition by E4,..,E15 are not diagonalised, but we
# see below that they are greatly simplified:
WE4diag := simplify(multiply( inverse(Pdiag), multiply(WE4,Pdiag) )):
const := WE4diag[1,2]:
for i from 1 to 4 do for j from 1 to 4 do
  WE4diag[i,j] := simplify(WE4diag[i,j]/const) od od:
WE4diag := matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]):
WE5diag := simplify(multiply( inverse(Pdiag), multiply(WE5,Pdiag) )):
const := WE5diag[1,2]:
for i from 1 to 4 do for j from 1 to 4 do
  WE5diag[i,j] := simplify(WE5diag[i,j]/const) od od:
WE5diag := matrix([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]]):
WE6diag := simplify(multiply( inverse(Pdiag), multiply(WE6,Pdiag) )):
const := WE6diag[1,4]:
for i from 1 to 4 do for j from 1 to 4 do
  WE6diag[i,j] := simplify(WE6diag[i,j]/const) od od:
WE6diag := matrix([[0, 0, 0, 1], [0, 0, -1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]]):
WE7diag := simplify(multiply( inverse(Pdiag), multiply(WE7,Pdiag) )):
const := WE7diag[1,4]:
for i from 1 to 4 do for j from 1 to 4 do
  WE7diag[i,j] := simplify(WE7diag[i,j]/const) od od:
WE7diag := matrix([[0, 0, 0, 1], [0, 0, -1, 0], [0, -1, 0, 0], [1, 0, 0, 0]]):
WE8diag := simplify(multiply( inverse(Pdiag), multiply(WE8,Pdiag) )):
const := WE8diag[1,2]:
for i from 1 to 4 do for j from 1 to 4 do
  WE8diag[i,j] := simplify(WE8diag[i,j]/const) od od:
WE8diag := matrix([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]]):
WE9diag := simplify(multiply( inverse(Pdiag), multiply(WE9,Pdiag) )):
const := WE9diag[1,2]:
for i from 1 to 4 do for j from 1 to 4 do
  WE9diag[i,j] := simplify(WE9diag[i,j]/const) od od:
WE9diag := matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]):
WE10diag := simplify(multiply( inverse(Pdiag), multiply(WE10,Pdiag) )):
const := WE10diag[1,4]:
for i from 1 to 4 do for j from 1 to 4 do
  WE10diag[i,j] := simplify(WE10diag[i,j]/const) od od:
WE10diag := matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]):
WE11diag := simplify(multiply( inverse(Pdiag), multiply(WE11,Pdiag) )):
const := WE11diag[1,4]:
for i from 1 to 4 do for j from 1 to 4 do
  WE11diag[i,j] := simplify(WE11diag[i,j]/const) od od:
WE11diag := matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]]):
WE12diag := simplify(multiply( inverse(Pdiag), multiply(WE12,Pdiag) )):
const := WE12diag[1,3]:
for i from 1 to 4 do for j from 1 to 4 do
  WE12diag[i,j] := simplify(WE12diag[i,j]/const) od od:
WE12diag := matrix([[0, 0, 1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, -1, 0, 0]]):
WE13diag := simplify(multiply( inverse(Pdiag), multiply(WE13,Pdiag) )):
const := WE13diag[1,3]:
for i from 1 to 4 do for j from 1 to 4 do
  WE13diag[i,j] := simplify(WE13diag[i,j]/const) od od:
WE13diag := matrix([[0, 0, 1, 0], [0, 0, 0, -1], [-1, 0, 0, 0], [0, 1, 0, 0]]):
WE14diag := simplify(multiply( inverse(Pdiag), multiply(WE14,Pdiag) )):
const := WE14diag[1,3]:
for i from 1 to 4 do for j from 1 to 4 do
  WE14diag[i,j] := simplify(WE14diag[i,j]/const) od od:
WE14diag := matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]]):
WE15diag := simplify(multiply( inverse(Pdiag), multiply(WE15,Pdiag) )):
const := WE15diag[1,3]:
for i from 1 to 4 do for j from 1 to 4 do
  WE15diag[i,j] := simplify(WE15diag[i,j]/const) od od:
WE15diag := matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]):
# The following is the array of matrix which give addition by our 2-torsion
# points E0,..,E15, with respect to our new Kummer basis:
WEdiagarray := 
  [WE0diag,WE1diag,WE2diag,WE3diag,WE4diag,WE5diag,WE6diag,WE7diag,
   WE8diag,WE9diag,WE10diag,WE11diag,WE12diag,WE13diag,WE14diag,WE15diag]:
# 
# Let us denote the new coordinates as (l1,l2,l3,l4).
# We wish to find the l-coordinates of E0,..,E15.
# We shall first find inverse(Pdiag) and call it Pdiaginverse.
# We shall then look at:
# multiply(WEidiag, kEi) to get lEi.
#
# The following is the inverse of Pdiag.
Pdiaginverse := inverse(Pdiag):
for i from 1 to 4 do for j from 1 to 4 do
  Pdiaginverse[i,j] := 
  factor(inverse(Pdiag)[i,j]
  *16*(a^4*b^2*d^2-a^2*b^4*c^2-a^2*c^2*d^4+b^2*c^4*d^2)):
od od:
# We now find the l-coorinates of E0,..,E15.
multiply( Pdiaginverse, matrix(kE0) ):
lE0 := [ [factor(%[1,1]*(a*c+b*d)*C*D/A^2/B^2/(a*c-b*d)^2)], 
         [factor(%[2,1]*(a*c+b*d)*C*D/A^2/B^2/(a*c-b*d)^2)], 
         [factor(%[3,1]*(a*c+b*d)*C*D/A^2/B^2/(a*c-b*d)^2)], 
         [factor(%[4,1]*(a*c+b*d)*C*D/A^2/B^2/(a*c-b*d)^2)] ]:
l1E0 := lE0[1][1]: l2E0 := lE0[2][1]: l3E0 := lE0[3][1]: l4E0 := lE0[4][1]:
multiply( WE1diag, matrix(lE0) ):
lE1 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E1 := lE1[1][1]: l2E1 := lE1[2][1]: l3E1 := lE1[3][1]: l4E1 := lE1[4][1]:
multiply( WE2diag, matrix(lE0) ):
lE2 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E2 := lE2[1][1]: l2E2 := lE2[2][1]: l3E2 := lE2[3][1]: l4E2 := lE2[4][1]:
multiply( WE3diag, matrix(lE0) ):
lE3 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E3 := lE3[1][1]: l2E3 := lE3[2][1]: l3E3 := lE3[3][1]: l4E3 := lE3[4][1]:
multiply( WE4diag, matrix(lE0) ):
lE4 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E4 := lE4[1][1]: l2E4 := lE4[2][1]: l3E4 := lE4[3][1]: l4E4 := lE4[4][1]:
multiply( WE5diag, matrix(lE0) ):
lE5 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E5 := lE5[1][1]: l2E5 := lE5[2][1]: l3E5 := lE5[3][1]: l4E5 := lE5[4][1]:
multiply( WE6diag, matrix(lE0) ):
lE6 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E6 := lE6[1][1]: l2E6 := lE6[2][1]: l3E6 := lE6[3][1]: l4E6 := lE6[4][1]:
multiply( WE7diag, matrix(lE0) ):
lE7 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E7 := lE7[1][1]: l2E7 := lE7[2][1]: l3E7 := lE7[3][1]: l4E7 := lE7[4][1]:
multiply( WE8diag, matrix(lE0) ):
lE8 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E8 := lE8[1][1]: l2E8 := lE8[2][1]: l3E8 := lE8[3][1]: l4E8 := lE8[4][1]:
multiply( WE9diag, matrix(lE0) ):
lE9 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E9 := lE9[1][1]: l2E9 := lE9[2][1]: l3E9 := lE9[3][1]: l4E9 := lE9[4][1]:
multiply( WE10diag, matrix(lE0) ):
lE10 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E10:=lE10[1][1]: l2E10:=lE10[2][1]: l3E10:=lE10[3][1]: l4E10:=lE10[4][1]:
multiply( WE11diag, matrix(lE0) ):
lE11 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E11:=lE11[1][1]: l2E11:=lE11[2][1]: l3E11:=lE11[3][1]: l4E11:=lE11[4][1]:
multiply( WE12diag, matrix(lE0) ):
lE12 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E12:=lE12[1][1]: l2E12:=lE12[2][1]: l3E12:=lE12[3][1]: l4E12:=lE12[4][1]:
multiply( WE13diag, matrix(lE0) ):
lE13 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E13:=lE13[1][1]: l2E13:=lE13[2][1]: l3E13:=lE13[3][1]: l4E13:=lE13[4][1]:
multiply( WE14diag, matrix(lE0) ):
lE14 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E14:=lE14[1][1]: l2E14:=lE14[2][1]: l3E14:=lE14[3][1]: l4E14:=lE14[4][1]:
multiply( WE15diag, matrix(lE0) ):
lE15 := [ [%[1,1]], [%[2,1]], [%[3,1]], [%[4,1]] ]:
l1E15:=lE15[1][1]: l2E15:=lE15[2][1]: l3E15:=lE15[3][1]: l4E15:=lE15[4][1]:
# The following is the array of l-coodinates of all 2-torsion points:
ELarray := [ [[l1E0],[l2E0],[l3E0],[l4E0]],
             [[l1E1],[l2E1],[l3E1],[l4E1]],
             [[l1E2],[l2E2],[l3E2],[l4E2]],
             [[l1E3],[l2E3],[l3E3],[l4E3]],
             [[l1E4],[l2E4],[l3E4],[l4E4]],
             [[l1E5],[l2E5],[l3E5],[l4E5]],
             [[l1E6],[l2E6],[l3E6],[l4E6]],
             [[l1E7],[l2E7],[l3E7],[l4E7]],
             [[l1E8],[l2E8],[l3E8],[l4E8]],
             [[l1E9],[l2E9],[l3E9],[l4E9]],
             [[l1E10],[l2E10],[l3E10],[l4E10]],
             [[l1E11],[l2E11],[l3E11],[l4E11]],
             [[l1E12],[l2E12],[l3E12],[l4E12]],
             [[l1E13],[l2E13],[l3E13],[l4E13]],
             [[l1E14],[l2E14],[l3E14],[l4E14]],
             [[l1E15],[l2E15],[l3E15],[l4E15]] ]:

# The following provides a double check on correctness: 
checkmtxk := matrix(4,16):
for i from 1 to 4 do for j from 1 to 16 do
  checkmtxk[i,j] := EKarray[j][i][1]
od od:
checkmtxl := matrix(4,16):
for i from 1 to 4 do for j from 1 to 16 do
  checkmtxl[i,j] := ELarray[j][i][1]
od od:
checkmtxlvar := matrix(4,16):
for i from 1 to 4 do for j from 1 to 16 do
  checkmtxlvar[i,j] := multiply( WEdiagarray[j], matrix(lE0) )[i,1]
od od:

# The entries in the following check matrix are all zero
for i from 1 to 4 do for j from 1 to 16 do
   print( simplify( checkmtxlvar[i,j] - checkmtxl[i,j] )):
od od: 

# The following is the defining equation K^fast given by Gaudry.
Gaudry :=
 (X^4 + Y^4 + Z^4 + T^4) - F*(X^2*T^2+Y^2*Z^2)
    - G*(X^2*Z^2+Y^2*T^2) - H*(X^2*Y^2+Z^2*T^2) + 2*E*X*Y*Z*T:
Gaudryvar := subs(X=l1, Y=l2, Z=l3, T=l4, Gaudry):
# Let's also see what the Kummer equation becomes in l-coordinates.
kummereqnl := subs( op(subf),
  k1 = multiply( Pdiag, matrix([[l1],[l2],[l3],[l4]]) )[1,1],
  k2 = multiply( Pdiag, matrix([[l1],[l2],[l3],[l4]]) )[2,1],
  k3 = multiply( Pdiag, matrix([[l1],[l2],[l3],[l4]]) )[3,1],
  k4 = multiply( Pdiag, matrix([[l1],[l2],[l3],[l4]]) )[4,1],
 kummeqn):
kummereqnl := simplify(factor(kummereqnl)):
# The following checks that the transformed kummer is the
# same as Gaudry.
check := factor(Gaudryvar 
         + 1/4096*B^6*A^6*(a*c-b*d)^5/D^2/(a*c+b*d)^3/(a*d+b*c)^3
            /(a*b-c*d)^3/(a*b+c*d)^3/(a*d-b*c)^3/C^2*kummereqnl);

# So our linear map from Kummer(D) to Gaudry is Pdiag,
# which confirms the claim in the article.

# The following are the Gaudry biquadratic forms:
B11gv := (TP^2+XP^2+YP^2+ZP^2)*(TQ^2+XQ^2+YQ^2+ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2+(1/
2)*c^2+(1/2)*d^2))+(-TP^2+XP^2+YP^2-ZP^2)*(-TQ^2+XQ^2+YQ^2-ZQ^2)/(4*((1/2)*a^2+(
1/2)*b^2-(1/2)*c^2-(1/2)*d^2))+(-TP^2+XP^2-YP^2+ZP^2)*(-TQ^2+XQ^2-YQ^2+ZQ^2)/(4*
((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2))+(TP^2+XP^2-YP^2-ZP^2)*(TQ^2+XQ^2-YQ^2
-ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)):

B22gv := (TP^2+XP^2+YP^2+ZP^2)*(TQ^2+XQ^2+YQ^2+ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2+(1/
2)*c^2+(1/2)*d^2))+(-TP^2+XP^2+YP^2-ZP^2)*(-TQ^2+XQ^2+YQ^2-ZQ^2)/(4*((1/2)*a^2+(
1/2)*b^2-(1/2)*c^2-(1/2)*d^2))-(-TP^2+XP^2-YP^2+ZP^2)*(-TQ^2+XQ^2-YQ^2+ZQ^2)/(4*
((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2))-(TP^2+XP^2-YP^2-ZP^2)*(TQ^2+XQ^2-YQ^2
-ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)):

B33gv := (TP^2+XP^2+YP^2+ZP^2)*(TQ^2+XQ^2+YQ^2+ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2+(1/
2)*c^2+(1/2)*d^2))-(-TP^2+XP^2+YP^2-ZP^2)*(-TQ^2+XQ^2+YQ^2-ZQ^2)/(4*((1/2)*a^2+(
1/2)*b^2-(1/2)*c^2-(1/2)*d^2))+(-TP^2+XP^2-YP^2+ZP^2)*(-TQ^2+XQ^2-YQ^2+ZQ^2)/(4*
((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2))-(TP^2+XP^2-YP^2-ZP^2)*(TQ^2+XQ^2-YQ^2
-ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)):

B44gv := (TP^2+XP^2+YP^2+ZP^2)*(TQ^2+XQ^2+YQ^2+ZQ^2)/(4*((1/2)*a^2+(1/2)*b^2+(1/
2)*c^2+(1/2)*d^2))-(-TP^2+XP^2+YP^2-ZP^2)*(-TQ^2+XQ^2+YQ^2-ZQ^2)/(4*((1/2)*a^2+(
1/2)*b^2-(1/2)*c^2-(1/2)*d^2))-(-TP^2+XP^2-YP^2+ZP^2)*(-TQ^2+XQ^2-YQ^2+ZQ^2)/(4*
((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2))+(TP^2+XP^2-YP^2-ZP^2)*(TQ^2+XQ^2-YQ^2
-ZQ^2)/(4*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)):

B12gv := (2*(a*b*(TP*TQ*ZP*ZQ+XP*XQ*YP*YQ)-c*d*(TP*XQ*YQ*ZP+TQ*XP*YP*ZQ)))/(((1/
2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)-
((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d
^2)):

B13gv := (2*(a*c*(TP*TQ*YP*YQ+XP*XQ*ZP*ZQ)-b*d*(TP*XQ*YP*ZQ+TQ*XP*YQ*ZP)))/(((1/
2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2)-
((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d
^2)):

B14gv := (2*(a*d*(TP*TQ*XP*XQ+YP*YQ*ZP*ZQ)-b*c*(TP*XP*YQ*ZQ+TQ*XQ*YP*ZP)))/(((1/
2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)-
((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d
^2)):

B23gv := (2*(b*c*(TP*TQ*XP*XQ+YP*YQ*ZP*ZQ)-a*d*(TP*XP*YQ*ZQ+TQ*XQ*YP*ZP)))/(((1/
2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2)-
((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d
^2)):

B24gv := (2*(b*d*(TP*TQ*YP*YQ+XP*XQ*ZP*ZQ)-a*c*(TP*XQ*YP*ZQ+TQ*XP*YQ*ZP)))/(((1/
2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)-
((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d
^2)):

B34gv := (2*(c*d*(TP*TQ*ZP*ZQ+XP*XQ*YP*YQ)-a*b*(TP*XQ*YQ*ZP+TQ*XP*YP*ZQ)))/(((1/
2)*a^2-(1/2)*b^2+(1/2)*c^2-(1/2)*d^2)*((1/2)*a^2-(1/2)*b^2-(1/2)*c^2+(1/2)*d^2)-
((1/2)*a^2+(1/2)*b^2+(1/2)*c^2+(1/2)*d^2)*((1/2)*a^2+(1/2)*b^2-(1/2)*c^2-(1/2)*d
^2)):

Bgv := matrix( [
        [B11gv,B12gv,B13gv,B14gv],
        [B12gv,B22gv,B23gv,B24gv],
        [B13gv,B23gv,B33gv,B34gv],
        [B14gv,B24gv,B34gv,B44gv]
       ]):

# The following is a more succinct way of writing them. 
# The could be even more succinctly written if we were to include
# a notation for things like TP^2+XP^2+YP^2+ZP^2, etc.
# These might be more elegant with some permutation of a,b,c,d,
# which could be labelled e1,e2,e3,e4 in some order.
# Note that the following g1,g2,g3,g4 are the same as A,B,C,D.
g1 := a^2+b^2+c^2+d^2:
g2 := a^2+b^2-c^2-d^2:
g3 := a^2-b^2+c^2-d^2:
g4 := a^2-b^2-c^2+d^2:
B11g := (XP^2+YP^2+ZP^2+TP^2)*(XQ^2+YQ^2+ZQ^2+TQ^2)/(4*g1)
       +(XP^2+YP^2-ZP^2-TP^2)*(XQ^2+YQ^2-ZQ^2-TQ^2)/(4*g2)
       +(XP^2-YP^2+ZP^2-TP^2)*(XQ^2-YQ^2+ZQ^2-TQ^2)/(4*g3)
       +(XP^2-YP^2-ZP^2+TP^2)*(XQ^2-YQ^2-ZQ^2+TQ^2)/(4*g4):

B22g := (XP^2+YP^2+ZP^2+TP^2)*(XQ^2+YQ^2+ZQ^2+TQ^2)/(4*g1)
       +(XP^2+YP^2-ZP^2-TP^2)*(XQ^2+YQ^2-ZQ^2-TQ^2)/(4*g2)
       -(XP^2-YP^2+ZP^2-TP^2)*(XQ^2-YQ^2+ZQ^2-TQ^2)/(4*g3)
       -(XP^2-YP^2-ZP^2+TP^2)*(XQ^2-YQ^2-ZQ^2+TQ^2)/(4*g4):

B33g := (XP^2+YP^2+ZP^2+TP^2)*(XQ^2+YQ^2+ZQ^2+TQ^2)/(4*g1)
       -(XP^2+YP^2-ZP^2-TP^2)*(XQ^2+YQ^2-ZQ^2-TQ^2)/(4*g2)
       +(XP^2-YP^2+ZP^2-TP^2)*(XQ^2-YQ^2+ZQ^2-TQ^2)/(4*g3)
       -(XP^2-YP^2-ZP^2+TP^2)*(XQ^2-YQ^2-ZQ^2+TQ^2)/(4*g4):

B44g := (XP^2+YP^2+ZP^2+TP^2)*(XQ^2+YQ^2+ZQ^2+TQ^2)/(4*g1)
       -(XP^2+YP^2-ZP^2-TP^2)*(XQ^2+YQ^2-ZQ^2-TQ^2)/(4*g2)
       -(XP^2-YP^2+ZP^2-TP^2)*(XQ^2-YQ^2+ZQ^2-TQ^2)/(4*g3)
       +(XP^2-YP^2-ZP^2+TP^2)*(XQ^2-YQ^2-ZQ^2+TQ^2)/(4*g4):

B12g := 4*(a*b*(XP*YP*XQ*YQ+ZP*TP*ZQ*TQ)-c*d*(XP*YP*ZQ*TQ+ZP*TP*XQ*YQ))
          /(g1*g2-g3*g4):

B13g := 4*(a*c*(XP*ZP*XQ*ZQ+YP*TP*YQ*TQ)-b*d*(XP*ZP*YQ*TQ+YP*TP*XQ*ZQ))
          /(g1*g3-g2*g4):

B14g := 4*(a*d*(XP*TP*XQ*TQ+YP*ZP*YQ*ZQ)-b*c*(XP*TP*YQ*ZQ+YP*ZP*XQ*TQ))
          /(g1*g4-g2*g3):

B23g := 4*(b*c*(XP*TP*XQ*TQ+YP*ZP*YQ*ZQ)-a*d*(XP*TP*YQ*ZQ+YP*ZP*XQ*TQ))
          /(g2*g3-g1*g4):

B24g := 4*(b*d*(XP*ZP*XQ*ZQ+YP*TP*YQ*TQ)-a*c*(XP*ZP*YQ*TQ+YP*TP*XQ*ZQ))
          /(g2*g4-g1*g3):

B34g := 4*(c*d*(XP*YP*XQ*YQ+ZP*TP*ZQ*TQ)-a*b*(ZP*TP*XQ*YQ+XP*YP*ZQ*TQ))
          /(g3*g4-g1*g2):

Bg := matrix( [
        [B11g,B12g,B13g,B14g],
        [B12g,B22g,B23g,B24g],
        [B13g,B23g,B33g,B34g],
        [B14g,B24g,B34g,B44g]
       ]):

check := factor(B11gv - 2*B11g);
check := factor(B22gv - 2*B22g);
check := factor(B33gv - 2*B33g);
check := factor(B44gv - 2*B44g);
check := factor(B12gv - 2*B12g);
check := factor(B13gv - 2*B13g);
check := factor(B14gv - 2*B14g);
check := factor(B23gv - 2*B23g);
check := factor(B24gv - 2*B24g);
check := factor(B34gv - 2*B34g);

# We should be above to develop division polynomials from these.
# Natural to define:
# phi_X^(0) = a, phi_Y^(0) = b, phi_Z^(0) = c, phi_T^(0) = d,
# phi_X^(1) = X, phi_Y^(1) = Y, phi_Z^(1) = Z, phi_T^(1) = T,
# and then inductively define:
# phi_X^(2N) = B11( (phi_X^(N),phi_Y^(N),phi_Z^(N),phi_T^(N)),
#                      (phi_X^(N),phi_Y^(N),phi_Z^(N),phi_T^(N)) )/a
# phi_Y^(2N) = B22( (phi_X^(N),phi_Y^(N),phi_Z^(N),phi_T^(N)),
#                      (phi_X^(N),phi_Y^(N),phi_Z^(N),phi_T^(N)) )/b
# phi_Z^(2N) = B33( (phi_X^(N),phi_Y^(N),phi_Z^(N),phi_T^(N)),
#                      (phi_X^(N),phi_Y^(N),phi_Z^(N),phi_T^(N)) )/c
# phi_T^(2N) = B44( (phi_X^(N),phi_Y^(N),phi_Z^(N),phi_T^(N)),
#                      (phi_X^(N),phi_Y^(N),phi_Z^(N),phi_T^(N)) )/d
# and:
# phi_X^(2N+1) = B11( (phi_X^(N+1),phi_Y^(N+1),phi_Z^(N+1),phi_T^(N+1)),
#                      (phi_X^(N),phi_Y^(N),phi_Z^(N),phi_T^(N)) )/X
# phi_Y^(2N+1) = B22( (phi_X^(N+1),phi_Y^(N+1),phi_Z^(N+1),phi_T^(N+1)),
#                      (phi_X^(N),phi_Y^(N),phi_Z^(N),phi_T^(N)) )/Y
# phi_Z^(2N+1) = B33( (phi_X^(N+1),phi_Y^(N+1),phi_Z^(N+1),phi_T^(N+1)),
#                      (phi_X^(N),phi_Y^(N),phi_Z^(N),phi_T^(N)) )/Z
# phi_T^(2N+1) = B44( (phi_X^(N+1),phi_Y^(N+1),phi_Z^(N+1),phi_T^(N+1)),
#                      (phi_X^(N),phi_Y^(N),phi_Z^(N),phi_T^(N)) )/T
# However, when we try this for phi_X^(3),phi_Y^(3),phi_Z^(3),phi_T^(3),
# we don't get that these are polynomials, so maybe these need to
# be adjusted. When they are adjusted, it seems likely that
# for N odd, the effect of addition by the points of order 2
# can easily be shown to transfer over transparently (not merely
# projectively); for example: X <--> Y, Z <--> T will
# have the effect: phi_X^(N) <--> phi_Y^(N), phi_Z^(N) <--> phi_T^(N)
# as affine polynomails, and the terms involved in these
# partition as described for the (N,N)-isogeny.
X2P := subs(XP=X,YP=Y,ZP=Z,TP=T, XQ=X,YQ=Y,ZQ=Z,TQ=T, B11g)/a:
Y2P := subs(XP=X,YP=Y,ZP=Z,TP=T, XQ=X,YQ=Y,ZQ=Z,TQ=T, B22g)/b:
Z2P := subs(XP=X,YP=Y,ZP=Z,TP=T, XQ=X,YQ=Y,ZQ=Z,TQ=T, B33g)/c:
T2P := subs(XP=X,YP=Y,ZP=Z,TP=T, XQ=X,YQ=Y,ZQ=Z,TQ=T, B44g)/d:
X3P := subs(XP=X2P,YP=Y2P,ZP=Z2P,TP=T2P, XQ=X,YQ=Y,ZQ=Z,TQ=T, B11g)/X:
Y3P := subs(XP=X2P,YP=Y2P,ZP=Z2P,TP=T2P, XQ=X,YQ=Y,ZQ=Z,TQ=T, B22g)/Y:
Z3P := subs(XP=X2P,YP=Y2P,ZP=Z2P,TP=T2P, XQ=X,YQ=Y,ZQ=Z,TQ=T, B33g)/Z:
T3P := subs(XP=X2P,YP=Y2P,ZP=Z2P,TP=T2P, XQ=X,YQ=Y,ZQ=Z,TQ=T, B44g)/T:
# However, X3P,Y3P,Z3P,T3P are not polynomials, so need to work
# modulo the Gaudry equation!
Xa := X*X3P:
Xb := factor(subs(X=0,Xa)):
Gaudrymult := factor(Xb/subs(X=0,Gaudry)):
X3Pfinal := factor(Xa - Gaudrymult*Gaudry)/X:
Ya := Y*Y3P:
Yb := factor(subs(Y=0,Ya)):
Gaudrymult := factor(Yb/subs(Y=0,Gaudry)):
Y3Pfinal := factor(Ya - Gaudrymult*Gaudry)/Y:
Za := Z*Z3P:
Zb := factor(subs(Z=0,Za)):
Gaudrymult := factor(Zb/subs(Z=0,Gaudry)):
Z3Pfinal := factor(Za - Gaudrymult*Gaudry)/Z:
Ta := T*T3P:
Tb := factor(subs(T=0,Ta)):
Gaudrymult := factor(Tb/subs(T=0,Gaudry)):
T3Pfinal := factor(Ta - Gaudrymult*Gaudry)/T:
# We see the X3P has only the degree 9 monomials which are invariant
# under addition by the following pts of order 2, with maps
#   (X,Y,Z,T) |--> (X,Y,Z,T).
#   (X,Y,Z,T) |--> (X,Y,-Z,-T).
#   (X,Y,Z,T) |--> (X,-Y,Z,-T).
#   (X,Y,Z,T) |--> (X,-Y,-Z,T).
# Similarly Y3P has only degree 9 monomials which are left
# unchanged by the first and second and negated by the third and fourth.
# Similarly Z3P has only degree 9 monomials which are left
# unchanged by the first and third and negated by the second and fourth.
# Similarly T3P has only degree 9 monomials which are left
# unchanged by the first and fourth and negated by the second and third.
#   This is as we would expect, since mult by 3 takes any given
# point of order 2 to itself.
#   Also, the effect of adding the other points of order 2 gives
# the same result transparently on X3P, Y3P, Z3P, T3P as
# affine polynomials (not merely projectively). For example:
#   X <--> Y, Z <--> T literally swaps X3P <--> Y3P, Z3P <--> T3P.
subs(X=x,Y=y,Z=z,T=t,Y3Pfinal):
check := factor(subs(x=Y,y=X,z=T,t=Z,%) - X3Pfinal);
#
#
# The following are the same biquadratic forms, but with:
#     (XP,YP,ZP,TP) replaced by (v1,v2,v3,v4)
# and (XQ,YQ,ZQ,TQ) replaced by (w1,w2,w3,w4).
# and g1,g2,g3,g4 replaced by A,B,C,D.
B11h := (v1^2+v2^2+v3^2+v4^2)*(w1^2+w2^2+w3^2+w4^2)/(4*g1)
       +(v1^2+v2^2-v3^2-v4^2)*(w1^2+w2^2-w3^2-w4^2)/(4*g2)
       +(v1^2-v2^2+v3^2-v4^2)*(w1^2-w2^2+w3^2-w4^2)/(4*g3)
       +(v1^2-v2^2-v3^2+v4^2)*(w1^2-w2^2-w3^2+w4^2)/(4*g4):

B22h := (v1^2+v2^2+v3^2+v4^2)*(w1^2+w2^2+w3^2+w4^2)/(4*g1)
       +(v1^2+v2^2-v3^2-v4^2)*(w1^2+w2^2-w3^2-w4^2)/(4*g2)
       -(v1^2-v2^2+v3^2-v4^2)*(w1^2-w2^2+w3^2-w4^2)/(4*g3)
       -(v1^2-v2^2-v3^2+v4^2)*(w1^2-w2^2-w3^2+w4^2)/(4*g4):

B33h := (v1^2+v2^2+v3^2+v4^2)*(w1^2+w2^2+w3^2+w4^2)/(4*g1)
       -(v1^2+v2^2-v3^2-v4^2)*(w1^2+w2^2-w3^2-w4^2)/(4*g2)
       +(v1^2-v2^2+v3^2-v4^2)*(w1^2-w2^2+w3^2-w4^2)/(4*g3)
       -(v1^2-v2^2-v3^2+v4^2)*(w1^2-w2^2-w3^2+w4^2)/(4*g4):

B44h := (v1^2+v2^2+v3^2+v4^2)*(w1^2+w2^2+w3^2+w4^2)/(4*g1)
       -(v1^2+v2^2-v3^2-v4^2)*(w1^2+w2^2-w3^2-w4^2)/(4*g2)
       -(v1^2-v2^2+v3^2-v4^2)*(w1^2-w2^2+w3^2-w4^2)/(4*g3)
       +(v1^2-v2^2-v3^2+v4^2)*(w1^2-w2^2-w3^2+w4^2)/(4*g4):

B12h := 4*(a*b*(v1*v2*w1*w2+v3*v4*w3*w4)-c*d*(v1*v2*w3*w4+v3*v4*w1*w2))
          /(g1*g2-g3*g4):

B13h := 4*(a*c*(v1*v3*w1*w3+v2*v4*w2*w4)-b*d*(v1*v3*w2*w4+v2*v4*w1*w3))
          /(g1*g3-g2*g4):

B14h := 4*(a*d*(v1*v4*w1*w4+v2*v3*w2*w3)-b*c*(v1*v4*w2*w3+v2*v3*w1*w4))
          /(g1*g4-g2*g3):

B23h := 4*(b*c*(v1*v4*w1*w4+v2*v3*w2*w3)-a*d*(v1*v4*w2*w3+v2*v3*w1*w4))
          /(g2*g3-g1*g4):

B24h := 4*(b*d*(v1*v3*w1*w3+v2*v4*w2*w4)-a*c*(v1*v3*w2*w4+v2*v4*w1*w3))
          /(g2*g4-g1*g3):

B34h := 4*(c*d*(v1*v2*w1*w2+v3*v4*w3*w4)-a*b*(v3*v4*w1*w2+v1*v2*w3*w4))
          /(g3*g4-g1*g2):

Bh := matrix( [
        [B11h,B12h,B13h,B14h],
        [B12h,B22h,B23h,B24h],
        [B13h,B23h,B33h,B34h],
        [B14h,B24h,B34h,B44h]
       ]):

check := B11h - subs(XP=v1,YP=v2,ZP=v3,TP=v4,XQ=w1,YQ=w2,ZQ=w3,TQ=w4,B11g);
check := B22h - subs(XP=v1,YP=v2,ZP=v3,TP=v4,XQ=w1,YQ=w2,ZQ=w3,TQ=w4,B22g);
check := B33h - subs(XP=v1,YP=v2,ZP=v3,TP=v4,XQ=w1,YQ=w2,ZQ=w3,TQ=w4,B33g);
check := B44h - subs(XP=v1,YP=v2,ZP=v3,TP=v4,XQ=w1,YQ=w2,ZQ=w3,TQ=w4,B44g);
check := B12h - subs(XP=v1,YP=v2,ZP=v3,TP=v4,XQ=w1,YQ=w2,ZQ=w3,TQ=w4,B12g);
check := B13h - subs(XP=v1,YP=v2,ZP=v3,TP=v4,XQ=w1,YQ=w2,ZQ=w3,TQ=w4,B13g);
check := B14h - subs(XP=v1,YP=v2,ZP=v3,TP=v4,XQ=w1,YQ=w2,ZQ=w3,TQ=w4,B14g);
check := B23h - subs(XP=v1,YP=v2,ZP=v3,TP=v4,XQ=w1,YQ=w2,ZQ=w3,TQ=w4,B23g);
check := B24h - subs(XP=v1,YP=v2,ZP=v3,TP=v4,XQ=w1,YQ=w2,ZQ=w3,TQ=w4,B24g);
check := B34h - subs(XP=v1,YP=v2,ZP=v3,TP=v4,XQ=w1,YQ=w2,ZQ=w3,TQ=w4,B34g);

for i from 1 to 4 do for j from 1 to 4 do
Pdiag[i,j] := -Pdiag[i,j]*A^2*B^2*(a*c-b*d)^2/2 od od:

# As this point, we give a summary of the l-basis information.
#   Recall that the k-coordinates k1,k2,k3,k4 were our original
# standard Kummer coordinates, and we have changed bases
# to l-coordinates l1,l2,l3,l4, where:
#  ( l1 )                  ( k1 )
#  ( l2 )  =  Pdiaginverse ( k2 )
#  ( l3 )                  ( k3 )
#  ( l4 )                  ( k4 )
# where:
# Pdiaginverse := matrix(
# [[-2*c*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2),
# a*(a*c-b*d)*A*B,
#-2*(a*c-b*d)*A*B*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2)*c/((a*c+b*d)*C*D),
# a*A^2*B^2*(a*c-b*d)^2/((a*c+b*d)*C*D)],
# [2*d*(a^4*b^2-2*a^2*c^2*d^2-b^6+b^2*c^4+b^2*d^4),
# -b*(a*c-b*d)*A*B,
#-2*d*A*B*(a^4*b^2-2*a^2*c^2*d^2-b^6+b^2*c^4+b^2*d^4)*(a*c-b*d)/((a*c+b*d)*C*D),
# b*(a*c-b*d)^2*A^2*B^2/((a*c+b*d)*C*D)],
# [-2*a*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4), 
# c*(a*c-b*d)*A*B, 
#-2*a*(a*c-b*d)*A*B*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4)/((a*c+b*d)*C*D),
# c*A^2*B^2*(a*c-b*d)^2/((a*c+b*d)*C*D)], 
# [-2*b*(a^4*d^2-2*a^2*b^2*c^2+b^4*d^2+c^4*d^2-d^6), 
# -d*(a*c-b*d)*A*B, 
# 2*b*(a^4*d^2-2*a^2*b^2*c^2+b^4*d^2+c^4*d^2-d^6)*(a*c-b*d)*A*B/((a*c+b*d)*C*D),
# d*(a*c-b*d)^2*A^2*B^2/((a*c+b*d)*C*D)]]):
#
# Pdiag := matrix(
# [[c*(a*c-b*d)^2*A^2*B^2,
# -d*(a*c-b*d)^2*A^2*B^2,
# -a*(a*c-b*d)^2*A^2*B^2,
# b*(a*c-b*d)^2*A^2*B^2],
# [2*a*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4)*(a*c-b*d)*A*B,
#  2*b*(a^4*d^2-2*a^2*b^2*c^2+b^4*d^2+c^4*d^2-d^6)*(a*c-b*d)*A*B,
#  -2*c*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2)*(a*c-b*d)*A*B,
#  2*d*(a^4*b^2-2*a^2*c^2*d^2-b^6+b^2*c^4+b^2*d^4)*(a*c-b*d)*A*B],
# [c*(a*c+b*d)*(a*c-b*d)*A*B*C*D,
#  d*(a*c+b*d)*(a*c-b*d)*A*B*C*D,
#  -a*(a*c+b*d)*(a*c-b*d)*A*B*C*D,
#  -b*(a*c+b*d)*(a*c-b*d)*A*B*C*D],
# [2*a*(a^4*c^2-2*a^2*b^2*d^2+b^4*c^2-c^6+c^2*d^4)*(a*c+b*d)*C*D,
#  -2*b*(a^4*d^2-2*a^2*b^2*c^2+b^4*d^2+c^4*d^2-d^6)*(a*c+b*d)*C*D,
#  -2*c*(a^6-a^2*b^4-a^2*c^4-a^2*d^4+2*b^2*c^2*d^2)*(a*c+b*d)*C*D,
#  -2*d*(a^4*b^2-2*a^2*c^2*d^2-b^6+b^2*c^4+b^2*d^4)*(a*c+b*d)*C*D]]):
#
# In these l-coordinates, the 2-torsion points are:
# lE0,lE1,lE2,lE3,lE4,lE5,lE6,lE7,lE8,lE9,lE10,lE11,lE12,lE13,lE14,
# all in the array: ELarray.
# Note that in particular (these are column vectors),
# the above commands gave:
# lE0 := [[a], [b], [c], [d]]:
# lE1 := [[a], [-b], [c], [-d]]:
# lE2 := [[a], [-b], [-c], [d]]:
# lE3 := [[a], [b], [-c], [-d]]:

# Addition by these 2-torsion points (in l-coordinates) are:
# WE0diag, WE1diag, WE2diag, WE3diag, WE4diag, WE5diag, WE6diag, WE7diag, 
# WE8diag, WE9diag, WE10diag, WE11diag, WE12diag, WE13diag, WE14diag, WE15diag
# all in the array: WEdiagarray.
# Note that in particular:
# WE0diag := matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]):
# WE1diag := matrix([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]):
# WE2diag := matrix([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]]):
# WE3diag := matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]):


# We now look at linear maps : Gaudry --> some Gaudry.
# First note that addition by any point of order 2 gives
# a linear map from Gaudry to itself.
# Also, if there exists i in K then (X,Y,i*Z,i*T) gives
# a map from Gaudry_{a,b,c,d} to Gaudry_{a,b,c,-d}, where
# the equation has E negated.
# If we start with Gaudry_{a,b,c,d} and, for example, map:
# (X,Y,Z,T) |--> (X+Y+Z+T, X+Y-Z-T, X-Y+Z-T, X-Y-Z+T)
# then this maps Gaudry_{a,b,c,d} to 
# Gaudry_{a+b+c+d,a+b-c-d,a-b+c-d,a-b-c+d}.
subs(X=XX,Y=YY,Z=ZZ,T=TT, Gaudry):
subs(XX = X+Y+Z+T, YY=X+Y-Z-T, ZZ=X-Y+Z-T, TT=X-Y-Z+T, %):
makemonic := factor(coeff(%,X,4)):
Gaudrynew1 := %%/makemonic:
(a+b+c+d)^4 + (a+b-c-d)^4 - (a-b+c-d)^4 - (a-b-c+d)^4:
# The above reveals an interpretation of: a^3*b+3*a^2*c*d+...
# and might give an even cleaner way of writing some structures.
subs(a=aa,b=bb,c=cc,d=dd,E):
Enew := factor(subs(aa=a+b+c+d, bb=a+b-c-d, cc=a-b+c-d, dd=a-b-c+d, %)):
subs(a=aa,b=bb,c=cc,d=dd,F):
Fnew := factor(subs(aa=a+b+c+d, bb=a+b-c-d, cc=a-b+c-d, dd=a-b-c+d, %)):
subs(a=aa,b=bb,c=cc,d=dd,G):
Gnew := factor(subs(aa=a+b+c+d, bb=a+b-c-d, cc=a-b+c-d, dd=a-b-c+d, %)):
subs(a=aa,b=bb,c=cc,d=dd,H):
Hnew := factor(subs(aa=a+b+c+d, bb=a+b-c-d, cc=a-b+c-d, dd=a-b-c+d, %)):
Gaudrynew2 :=
 (X^4 + Y^4 + Z^4 + T^4) - Fnew*(X^2*T^2+Y^2*Z^2)
    - Gnew*(X^2*Z^2+Y^2*T^2) - Hnew*(X^2*Y^2+Z^2*T^2) + 2*Enew*X*Y*Z*T:
check := factor(Gaudrynew1 - Gaudrynew2);
# The above shows that (X,Y,Z,T) |-> (X+Y+Z+T, X+Y-Z-T, X-Y+Z-T, X-Y-Z+T)
# is a map : Gaudry_{a,b,c,d} --> Gaudry_{a+b+c+d,a+b-c-d,a-b+c-d,a-b-c+d}
# as claimed in the article.