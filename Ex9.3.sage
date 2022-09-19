################################################################################

# An example discussed in OSU
R.<x> = QQ['x']
X = HyperellipticCurve(x*(x-13)*(x-169)*(x-1)*(x-14)*(x-27)*(x-4))
K = Qp(13,25)
XK = X.change_ring(K)
L.<a> = K.extension(x^4-13)
XL = X.change_ring(L)
(T0,T13,T169,T1,T14,T27,T4) = (XK(0,0),XK(13,0),XK(169,0),XK(1,0),XK(14,0),XK(27,0),XK(4,0)) # finite W's points
(P1,P2,P3,P4) = (XK.lift_x(2),XK.lift_x(20/7),XK.lift_x(-13/12),XK.lift_x(169/14)) # ref pts in vertices
(R1,R2,(R3,R4)) = (XL.lift_x(a^2+1),XL.lift_x(a^2),XL.lift_x(13*a^2,all=True)) # ref pts in edges
def Log(z): return L(z).log(p_branch=0, change_frac=True) # Iwasawa branch
N = 40 # truncation level

################################################################################

# Integration on \pi⁻¹(U1)
R.<x> = QQ[]
H = HyperellipticCurve(x*(x-1)*(x-4))
g = H.hyperelliptic_polynomials()[0]
gder = g.derivative()
(HK,HL) = (H.change_ring(K),H.change_ring(L))
T.<u> = PowerSeriesRing(K, 'u', default_prec = 2*N)
l1 = ((1-13*u)*(1-169*u))^(-1/2)
T.<v> = PowerSeriesRing(K, 'v', default_prec = 2*N)
l2 = ((1-13*v)*(1-26*v))^(-1/2)
T.<u,v> = PowerSeriesRing(K, 'u,v', default_prec = 2*N)
l = T(l1)*T(l2)
l = l.truncate(N)

################################################################################

# first pole reduction
S = HK(0,0)
x,y = HK.local_coord(S,2*N+3)
ED = [0] + [(x*gder(x)-2*i*g(x))/(x^(i+1))*x.derivative()/(2*y) for i in range(1,N+1)]
w = [l(1/x,1/(x-1))*x^(i-1)/(x-1)*x.derivative()/(2*y) for i in range(3)]
FF = [[] for i in range(3)]
for i in range(3):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val/2].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val/2].list()[0])*ED[-val/2]
R.<x,y> = QQ[]
F1 = [0 for i in range(3)]
for i in range(3):
    for j in range(N-i):
        F1[i] += FF[i][j]*y/(x^(N-i-j)) # first exact part

################################################################################

# test
x,y = HK.local_coord(S,2*N+3)
[(l(1/x,1/(x-1))*x^(i-1)/(x-1)*x.derivative()/(2*y) - F1[i](x,y).derivative()).valuation() for i in range(3)]
### [0,0,0]

################################################################################

# second pole reduction
S = HK(1,0)
x,y = HK.local_coord(S,2*N+3)
ED = [0] + [((x-1)*gder(x)-2*i*g(x))/((x-1)^(i+1))*x.derivative()/(2*y) for i in range(1,N+1)]
w = [l(1/x,1/(x-1))*x^(i-1)/(x-1)*x.derivative()/(2*y) - F1[i](x,y).derivative() for i in range(3)]
FF = [[] for i in range(3)]
for i in range(3):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val/2].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val/2].list()[0])*ED[-val/2]
R.<x,y> = QQ[]
F2 = [0 for i in range(3)]
for i in range(3):
    for j in range(N):
        F2[i] += FF[i][j]*y/(x^(N-j)) # second exact part # BE CAREFUL WITH T --> X-1

################################################################################

# test
x,y = HK.local_coord(S,2*N+3)
[(l(1/x,1/(x-1))*x^(i-1)/(x-1)*x.derivative()/(2*y) - F1[i](x,y).derivative() - F2[i](x-1,y).derivative()).valuation() for i in range(3)]
### [0,0,0]

################################################################################

# de Rham coefficients
x,y = HK.local_coordinates_at_infinity(2)
omega = [x^i*x.derivative()/(2*y) for i in range(2)]
w = [l(1/x,1/(x-1))*x^(i-1)/(x-1)*x.derivative()/(2*y) - F1[i](x,y).derivative() - F2[i](x-1,y).derivative() for i in range(3)]
C = [[] for i in range(3)]
for i in range(3):
    while w[i] !=0:
        val = w[i].valuation()
        C[i] = [w[i].list()[0]/omega[-val/2].list()[0]] + C[i]
        w[i] = w[i] - (w[i].list()[0]/omega[-val/2].list()[0])*omega[-val/2]

################################################################################

# test
[l.truncate(N)(1/x,1/(x-1))*x^(i-1)/(x-1)*x.derivative()/(2*y) - F1[i](x,y).derivative() - F2[i](x-1,y).derivative() - C[i][0]*omega[0] - C[i][1]*omega[1] for i in range(3)]
### [O(t^4), O(t^4), O(t^4)]

################################################################################

# result: i-th (truncated) form = dF1[i] + dF2[i] + C[i][0]*omega0 + C[i][1]*omega1

################################################################################

def H(z):
    if z.scheme().base_ring() == K:
        return HK(z[0],z[1]/(z[0]*(z[0]-1)*((1-13/z[0])*(1-169/z[0])*(1-13/(z[0]-1))*(1-26/(z[0]-1))).square_root()))
    else:
        return HL(z[0],z[1]/(z[0]*(z[0]-1)*((1-13/z[0])*(1-169/z[0])*(1-13/(z[0]-1))*(1-26/(z[0]-1))).square_root()))
(HT4,HP1,HR1,HR2) = (H(T4),H(P1),H(R1),H(R2))
def EP(i,z1,z2): return F1[i](z2[0],z2[1]) + F2[i](z2[0]-1,z2[1]) - F1[i](z1[0],z1[1]) - F2[i](z1[0]-1,z1[1]) # exact part

################################################################################

T = HK(1,0)
x,y = HK.local_coord(T, 3*N)
HR1_to_T_dR = [-(((C[i][0] + C[i][1]*x)*x.derivative()/(2*y)).integral()).truncate()(HR1[1]) for i in range(3)]
x,y = HK.monsky_washnitzer_gens()
w = HK.invariant_differential()
T_to_HP1_dR = [HK.coleman_integral(C[i][0]*w+C[i][1]*x*w,T,HP1) for i in range(3)]
R1_to_P1 = [EP(i,HR1,HP1) + HR1_to_T_dR[i] + T_to_HP1_dR[i] for i in range(3)]
T = HK(0,0)
x,y = HK.local_coord(T,3*N)
T_to_HR2_dR = [(((C[i][0] + C[i][1]*x)*x.derivative()/(2*y)).integral()).truncate()(HR2[1]) for i in range(3)]
x,y = HK.monsky_washnitzer_gens()
w = HK.invariant_differential()
HP1_to_T_dR = [HK.coleman_integral(C[i][0]*w+C[i][1]*x*w,HP1,T) for i in range(3)]
P1_to_R2 = [EP(i,HP1,HR2) + HP1_to_T_dR[i] + T_to_HR2_dR[i] for i in range(3)]
Q = XK.lift_x(7)
Q_to_P1 = [EP(i,H(Q),HP1) + HK.coleman_integral(C[i][0]*w+C[i][1]*x*w,H(Q),HP1) for i in range(3)]

################################################################################

# Integration on \pi⁻¹(U2)
R.<x> = QQ['x']
H = HyperellipticCurve(x*(x-1)*(x-2))
g = H.hyperelliptic_polynomials()[0]
gder = g.derivative()
(HK,HL) = (H.change_ring(K),H.change_ring(L))
T.<t> = PowerSeriesRing(K, 't', default_prec = 2*N)
l = (((13*t+1)*(13*t-12)*(13*t-168)*(13*t-3)).square_root())^-1
l = l.truncate(N)
U = [(13*t+1)^i*l for i in range(3)]

################################################################################

# pole reduction at infinity
ED = [0,0] + [t^i*gder(t) + 2*i*t^(i-1)*g(t) for i in range(N+1)]
w = [U[i] for i in range(3)]
FF = [[] for i in range(3)]
for i in range(3):
    while w[i].degree() > 1:
        deg = w[i].degree()
        FF[i].append(w[i].list()[deg]/ED[deg].list()[deg])
        w[i] = w[i] - (w[i].list()[deg]/ED[deg].list()[deg])*ED[deg]
R.<x,y> = QQ[]
F = [0 for i in range(3)]
for i in range(3):
    for j in range(N-2+i):
        F[i] += FF[i][j]*(x^(N-3+i-j)*y) # exact part
f = [0 for i in range(3)]
for i in range(3):
    for j in range(N-2+i):
        f[i] += FF[i][j]*ED[N-1+i-j]

################################################################################

# test
[(U[i] - f[i]).degree() for i in range(3)]
### [1, 1, 1]

################################################################################

# de Rham coefficients
x,y = HK.local_coordinates_at_infinity(3*N)
w = [U[i](x)*x.derivative()/(2*y) - F[i](x,y).derivative() for i in range(3)]
omega = [x^i*x.derivative()/(2*y) for i in range(2)]
C = [[] for i in range(3)]
for i in range(3):
    while w[i] !=0:
        val = w[i].valuation()
        C[i] = [w[i].list()[0]/omega[-val/2].list()[0]] + C[i]
        w[i] = w[i] - (w[i].list()[0]/omega[-val/2].list()[0])*omega[-val/2]

################################################################################

# test
[U[i] == f[i] + C[i][0] + C[i][1]*t for i in range(3)]
### [True, True, True]

################################################################################

# result: i-th (truncated) form = dF[i] + C[i][0]*omega0 + C[i][1]*omega1

################################################################################

def H(z):
    if z[1]==0:
        return HK((z[0]-1)/13,0)
    else:
        return HL((z[0]-1)/13,z[1]/(13*a^2*(z[0]*(z[0]-13)*(z[0]-169)*(z[0]-4)).square_root()))
(HT1,HT14,HT27,HP2,HR1) = (H(T1),H(T14),H(T27),H(P2),H(R1))
def EP(i,z1,z2): return 1/a^2*(F[i](z2[0],z2[1]) - F[i](z1[0],z1[1])) # exact part

################################################################################

T = HK(2,0)
x,y = HK.local_coord(T,3*N)
HP2_to_T_dR = [-1/a^2*(((C[i][0] + C[i][1]*x)*x.derivative()/(2*y)).integral()).truncate()(HP2[1]) for i in range(3)]
HR1rational = HK.lift_x(1/169) # a Qp rational point in disc of HR1
x,y = HK.monsky_washnitzer_gens()
w = HK.invariant_differential()
T_to_HR1rational_dR = [1/a^2*(HK.coleman_integral(C[i][0]*w+C[i][1]*x*w,T,HR1rational)) for i in range(3)]
tiny_integrals_dR = HL.tiny_integrals_on_basis(HR1rational,HR1)
HR1rational_to_HR1_dR = [1/a^2*(C[i][0]*tiny_integrals_dR[0] + C[i][1]*tiny_integrals_dR[1]) for i in range(3)]
P2_to_R1 = [EP(i,HP2,HR1) + HP2_to_T_dR[i] + T_to_HR1rational_dR[i] + HR1rational_to_HR1_dR[i] for i in range(3)]
T1_to_P2 = [EP(i,T,HP2) - HP2_to_T_dR[i] for i in range(3)] # T1_to_P2 = HT1_to_HP2 = HT1_to_T + T_to_HP2 = T_to_HP2

################################################################################

# Integration on \pi⁻¹(U3)
R.<x> = QQ['x']
g = x-13
gder = g.derivative()
T.<t> = PowerSeriesRing(K, 't', default_prec = 2*N)
l = [(1-169*t)^(-1/2), (((t-1)*(t-14)*(t-27)*(t-4)).square_root())^-1]
UU=0
for i in range(N):
    UU += l[0].list()[i]*t^(N-1-i)
U = [UU/t^(N-i)*l[1].truncate(N) for i in range(3)]

################################################################################

# first pole reduction
T.<t> = PowerSeriesRing(L, 't', default_prec = 2*N)
(xS,yS) = (t,(t-13).sqrt()) # local interpolation around x=0
wS = [U[i]*xS.derivative()/(2*yS) for i in range(3)]
ED = [0,0] + [(t*gder(xS)-2*i*g(xS))/(t^(i+1))*xS.derivative()/(2*yS) for i in range(1,N)]
c = [wS[i].residue()/(1/t*xS.derivative()/(2*yS)).residue() for i in range(3)]
w = [wS[i] - c[i]*(1/t*xS.derivative()/(2*yS)) for i in range(3)]
FF = [[] for i in range(3)]
for i in range(3):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] -(w[i].list()[0]/ED[-val].list()[0])*ED[-val]
R.<x,y> = QQ[]
F1 = [0 for i in range(3)]
for i in range(3):
    for j in range(N-1-i):
        F1[i] += FF[i][j]*y/(x^(N-1-i-j)) # first exact part
EDx = [0,0] + [(t*gder(t)-2*i*g(t))/t^(i+1) for i in range(1,N)]
f1 = [0 for i in range(3)]
for i in range(3):
    for j in range(N-1-i):
        f1[i] += FF[i][j]*EDx[N-i-j]

################################################################################

# test
[(U[i] - f1[i] - c[i]*t^-1).valuation() for i in range(3)]
### [0, 0, 0]

################################################################################

# pole reduction at infinity
ED = [t^i*gder(t) + 2*i*t^(i-1)*g(t) for i in range(N+1)]
w = [U[i] - f1[i] - c[i]*t^-1 for i in range(3)]
FF = [[] for i in range(3)]
for i in range(3):
    while w[i] != 0:
        deg = w[i].degree()
        FF[i].append(w[i].list()[deg]/ED[deg].list()[deg])
        w[i] = w[i] - (w[i].list()[deg]/ED[deg].list()[deg])*ED[deg]
R.<x,y> = QQ[]
F2 = [0 for i in range(3)]
for i in range(3):
    for j in range(N-1+i):
        F2[i] += FF[i][j]*(x^(N-2+i-j)*y) # second exact part
f2 = [0 for i in range(3)]
for i in range(3):
    for j in range(N-1+i):
        f2[i] += FF[i][j]*ED[N-2+i-j]

################################################################################

# test
[U[i] == c[i]*t^-1 + f1[i] + f2[i] for i in range(3)]
### [True, True, True]

################################################################################

# result: i-th (truncated) form = dF[i] + c[i]*1/x*dx/2y

################################################################################

def EP(i,z1,z2): return F1[i](z2[0],z2[1]) + F2[i](z2[0],z2[1]) - F1[i](z1[0],z1[1]) - F2[i](z1[0],z1[1]) # exact part
def H(z): return (z[0],z[1]/(z[0]*(1-169/z[0]).square_root()*((z[0]-1)*(z[0]-14)*(z[0]-27)*(z[0]-4)).square_root()))
I = K(-1).square_root()
def Int(i,z1,z2): return EP(i,H(z1),H(z2)) + c[i]/(2*I*a^2)*(Log((H(z2)[1]-I*a^2)/(H(z2)[1]+I*a^2))-Log((H(z1)[1]-I*a^2)/(H(z1)[1]+I*a^2)))
R2_to_P3 = [Int(i,R2,P3) for i in range(3)]
P3_to_R3 = [Int(i,P3,R3) for i in range(3)]
R4_to_P3 = [Int(i,R4,P3)  for i in range(3)]
T13_to_P3 = [Int(i,T13,P3)  for i in range(3)]

################################################################################

# Integration on \pi⁻¹(U4)
R.<x> = QQ['x']
g = x*(x-1)
gder = g.derivative()
T.<t> = PowerSeriesRing(L, 't', default_prec = 2*N)
l = (((169*t-13)*(169*t-1)*(169*t-14)*(169*t-27)*(169*t-4)).square_root())^-1
U = [(169*t)^i*l.truncate(N) for i in range(3)]

################################################################################

# pole reduction at infinity
ED = [0] + [t^i*gder(t) + 2*i*t^(i-1)*g(t) for i in range(N+1)]
w = [U[i] for i in range(3)]
FF = [[] for i in range(3)]
for i in range(3):
    while w[i].degree() > 0:
        deg = w[i].degree()
        FF[i].append(w[i].list()[deg]/ED[deg].list()[deg])
        w[i] = w[i] - (w[i].list()[deg]/ED[deg].list()[deg])*ED[deg]
R.<x,y> = QQ[]
F = [0 for i in range(3)]
for i in range(3):
    for j in range(N-1+i):
        F[i] += FF[i][j]*(x^(N-2+i-j)*y) # exact part
f = [0 for i in range(3)]
for i in range(3):
    for j in range(N-1+i):
        f[i] += FF[i][j]*ED[N-1+i-j]

################################################################################

# test
[(U[i] - f[i]).degree() for i in range(3)]
### [0, 0, 0]

################################################################################

# test
c = [(U[i] - f[i])[0] for i in range(3)]
[U[i] == f[i] + c[i] for i in range(3)]
### [True, True, True]

################################################################################

# result: i-th (truncated) form = dF[i] + c[i]*dx/2y

################################################################################

def EP(i,z1,z2): return F[i](z2[0],z2[1]) - F[i](z1[0],z1[1])
def H(z): return (z[0]/169,z[1]/(169*L((z[0]-13)*(z[0]-1)*(z[0]-14)*(z[0]-27)*(z[0]-4)).square_root()))
def Int(i,z1,z2): return EP(i,H(z1),H(z2)) + c[i]/2*(Log(H(z2)[0]+H(z2)[1]-1/2)-Log(H(z1)[0]+H(z1)[1]-1/2))
R3_to_P4 = [Int(i,R3,P4) for i in range(3)]
P4_to_R4 = [Int(i,P4,R4) for i in range(3)]
S = XK.lift_x(13*169)
P4_to_S = [Int(i,P4,S) for i in range(3)]

################################################################################

BC_along_gamma = [Q_to_P1[i] + P1_to_R2[i] + R2_to_P3[i] + P3_to_R3[i] + R3_to_P4[i] + P4_to_S[i] for i in range(3)]
period_integral = [P3_to_R3[i] + R3_to_P4[i] + P4_to_R4[i] + R4_to_P3[i] for i in range(3)]
CT = [period_integral[i]/2 for i in range(3)]
AB_Q_to_S = [BC_along_gamma[i] - CT[i] for i in range(3)]

################################################################################

AB_Q_to_S[1] + AB_Q_to_S[2]
### 11 + 2*a^4 + 6*a^8 + 8*a^12 + 9*a^16 + 9*a^20 + 4*a^24 + 4*a^28 + 4*a^32 + a^36 + 8*a^40 + 12*a^44 + 11*a^48 + 5*a^56 + 6*a^60 + 11*a^64 + 2*a^76 + 2*a^80 + 8*a^81 + 12*a^83 + 3*a^84 + 3*a^87 + 2*a^88 + 2*a^89 + 2*a^91 + 8*a^92 + 10*a^95 + 5*a^96 + 12*a^99 + O(a^100)

################################################################################
