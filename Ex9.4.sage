################################################################################

# Chabauty example
R.<x> = QQ[]
f = (x^2 - x - 1) * (x^2 - 2) * (x^2 + x - 1) # f mod5 = (x + 2)^2 * (x + 3)^2 * (x^2 + 3)
X = HyperellipticCurve(f) # LMFDB label 3200.f.819200.1
genus2reduction(0,f)
### Reduction data about this proper smooth genus 2 curve:
		y^2 = x^6 - 5*x^4 + 7*x^2 - 2
	A Minimal Equation (away from 2):
		y^2 = x^6 - 5*x^4 + 7*x^2 - 2
	Minimal Discriminant (away from 2):  800
	Conductor (away from 2): 25
	Local Data:
		p=2
		(potential) stable reduction:  (V), j1+j2=0, j1*j2=0
		p=5
		(potential) stable reduction:  (III)
		reduction at p: [I{1-1-0}] page 179, (1), f=2

################################################################################

k = GF(5)
fk = f.change_ring(k)
valuations = []
T.<t> = PowerSeriesRing(k, 't')
x = t + 1
y = fk(x)^(1/2) # for the point (1,1)
omega = x*x.derivative()/(2*y)
val = omega.valuation()
valuations.append(val)
x = t + 1
y = -fk(x)^(1/2) # for the point (1,-1)
omega = x*x.derivative()/(2*y)
val = omega.valuation()
valuations.append(val)
x = t + 4
y = fk(x)^(1/2) # for the point (4,1)
omega = x*x.derivative()/(2*y)
val = omega.valuation()
valuations.append(val)
x = t + 4
y = -fk(x)^(1/2) # for the point (4,-1)
omega = x*x.derivative()/(2*y)
val = omega.valuation()
valuations.append(val)
x = 1/t
y = fk(x)^(1/2) # for the point infty^+
omega = x*x.derivative()/(2*y)
val = omega.valuation()
valuations.append(val)
x = 1/t
y = -fk(x)^(1/2) # for the point infty^-
omega = x*x.derivative()/(2*y)
val = omega.valuation()
valuations.append(val)
valuations
### [0, 0, 0, 0, 0, 0]

################################################################################

R.<x> = QQ[]
f = x^2-2
fder = f.derivative()
K = Qp(5,60)
L.<a> = K.ext(x^2-2) # unramified extension
def Log(z): return L(z).log(p_branch=0, change_frac=True) # Iwasawa branch
N = 80 # truncation level
T.<t> = PowerSeriesRing(L, 't', default_prec = 4*N)
l = (1-5/4*t)^(-1/2)
UU = 0
for i in range(N):
    UU += l.list()[i]*((t^2)^(N-1-i))
U = UU/t^(2*N-1) # be careful with t

################################################################################

# first pole reduction: poles with x-coordinate 1/2
x = t+1/2 # local coordinates at x
y = f(x).sqrt() # local coordiantes at y
w =[U*U(t+1)*x^i*x.derivative()/(2*y) for i in range(2)]
c1 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(2)] # coefficients of 1/(x-1/2)*dx/2y
w = [w[i] - c1[i]*(1/t*x.derivative()/(2*y)) for i in range(2)]
ED = [0,0] + [(t*fder(x)-2*i*f(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,2*N-1)] # exact differentials
FF = [[] for i in range(2)]
for i in range(2):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
R.<x,y> = QQ[]
F1 = [0 for i in range(2)]
for i in range(2):
    for j in range(2*N-2):
        F1[i] += FF[i][j]*y/(x^(2*N-2-j)) # be careful with the x-coordinate: x--> x-1/2

################################################################################

# test
x = t+1/2
y = f(x).sqrt()
[(U*U(t+1)*x^i*x.derivative()/(2*y) - F1[i](t,y).derivative() - c1[i]*(1/t*x.derivative()/(2*y))).valuation() for i in range(2)]
### [0, 0]

################################################################################

# second pole reduction: poles with x-coordinate -1/2
x = t-1/2 # local coordinates at x
y = f(x).sqrt() # local coordiantes at y
w = [U(t-1)*U*x^i*x.derivative()/(2*y) - F1[i](t-1,y).derivative() - c1[i]*(1/(t-1)*x.derivative()/(2*y)) for i in range(2)]
c2 = [w[i].residue()/(1/t*x.derivative()/(2*y)).residue() for i in range(2)] # coefficients of 1/(x+1/2)*dx/2y
w = [w[i] - c2[i]*(1/t*x.derivative()/(2*y)) for i in range(2)]
ED = [0,0] + [(t*fder(x)-2*i*f(x))/(t^(i+1))*x.derivative()/(2*y) for i in range(1,2*N-1)] # exact differentials
FF = [[] for i in range(2)]
for i in range(2):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
R.<x,y> = QQ[]
F2 = [0 for i in range(2)]
for i in range(2):
    for j in range(2*N-2):
        F2[i] += FF[i][j]*y/(x^(2*N-2-j)) # be careful with the x-coordinate: x--> x+1/2

################################################################################

# test
x = t-1/2
y = f(x).sqrt()
[U(t-1)*U*x^i*x.derivative()/(2*y) == F1[i](t-1,y).derivative() + c1[i]*(1/(t-1)*x.derivative()/(2*y)) + F2[i](t,y).derivative() + c2[i]*(1/t*x.derivative()/(2*y)) for i in range(2)]
### [True, True]

################################################################################

def U(z):
    A = z[0] - 1/2
    AA = K(1-(5/4)/(A^2)).sqrt()
    B = z[0] + 1/2
    BB = K(1-(5/4)/(B^2)).sqrt()
    return (z[0],z[1]/(A*AA*B*BB))
def Int(i,z1,z2): # integral of omega_i from z1 to z2
    exact_part1 = F1[i](U(z2)[0]-1/2,U(z2)[1]) - F1[i](U(z1)[0]-1/2,U(z1)[1])
    exact_part2 = F2[i](U(z2)[0]+1/2,U(z2)[1]) - F2[i](U(z1)[0]+1/2,U(z1)[1])
    S = L(-7).sqrt()
    t1 = U(z1)[0] + U(z1)[1]
    t2 = U(z2)[0] + U(z2)[1]
    third_kind_part1 = c1[i]*1/S*(Log((t2-1/2-S/2)/(t2-1/2+S/2)) - Log((t1-1/2-S/2)/(t1-1/2+S/2)))
    third_kind_part2 = c2[i]*1/S*(Log((t2+1/2-S/2)/(t2+1/2+S/2)) - Log((t1+1/2-S/2)/(t1+1/2+S/2)))
    return exact_part1 + exact_part2 + third_kind_part1 + third_kind_part2

################################################################################

# Some tests:
P,Q,R,S = X(-1,-1),X(1,1),X(-1,1),X(1,-1) # All finite rational points on X
Int(0,P,Q) + Int(0,Q,S) == Int(0,P,S)
Int(1,P,S) + Int(1,S,Q) == Int(1,P,Q)
Int(1,P,Q) + Int(1,Q,S) + Int(1,S,R) == Int(1,P,R)
Int(1,P,Q) + Int(1,Q,S) + Int(1,S,R) + Int(1,R,P) == 0
Int(0,R,Q) + Int(0,Q,S) == Int(0,R,P) + Int(0,P,S)
### True
### True
### True
### True
### True

################################################################################

zero = [Int(0,P,Q), Int(0,P,R), Int(0,P,S), Int(0,Q,R), Int(0,Q,S), Int(0,R,S)]
one = [Int(1,P,Q), Int(1,P,R), Int(1,P,S), Int(1,Q,R), Int(1,Q,S), Int(1,R,S)]

################################################################################

zero[4]
### 2*5 + 5^4 + 3*5^6 + 2*5^7 + 2*5^8 + 4*5^9 + 5^10 + 4*5^11 + 2*5^12 + 4*5^13 + 2*5^14 + 3*5^15 + 5^16 + 3*5^17 + 4*5^18 + 3*5^19 + 3*5^20 + 4*5^21 + 3*5^22 + 3*5^23 + 3*5^24 + 5^25 + 5^26 + 3*5^27 + 3*5^30 + 5^31 + 2*5^32 + 4*5^35 + 4*5^36 + 4*5^37 + 4*5^38 + 5^39 + 2*5^40 + 3*5^41 + 4*5^43 + 4*5^44 + 5^45 + 5^47 + 3*5^48 + 2*5^49 + 2*5^50 + 4*5^51 + 3*5^52 + 3*5^55 + 3*5^56 + 2*5^58 + 5^59 + O(5^60)

################################################################################
