################################################################################

# Hyperelliptic curve y²=x(x-1)(x-2)(x-3)(x-7)
R.<x> = QQ[]
X = HyperellipticCurve(x*(x-1)*(x-2)*(x-3)*(x-7))
K = Qp(7,30)
XK = X.change_ring(K)
L.<a> = K.extension(x^2-7)
XL = X.change_ring(L)
(T0,T1,T2,T3,T7) = (XK(0,0),XK(1,0),XK(2,0),XK(3,0),XK(7,0)) # finite Weierstrass points
(P1,P2) = (XK.lift_x(-1), XK.lift_x(14)) # reference pts in vertices
(R1,R2) = XL.lift_x(a,all=True) # reference pts in edges
def Log(z): return L(z).log(p_branch=0, change_frac=True) # Iwasawa branch

################################################################################

# Integration on \pi⁻¹(U1)
R.<x> = QQ[]
H = HyperellipticCurve((x-1)*(x-2)*(x-3))
g = H.hyperelliptic_polynomials()[0]
gder = g.derivative()
(HK,HL) = (H.change_ring(K),H.change_ring(L))
N1 = 40 # truncation level
T.<t> = PowerSeriesRing(K, 't', default_prec = 2*N1)
l = (1-7*t)^(-1/2)
UU=0
for i in range(N1):
    UU += l.list()[i]*x^(N1-1-i) # polynomial on numerator
U = [UU/(x^(N1-i)) for i in range(2)]
(S,iS) = HK.lift_x(0,all=True) # poles of the form
xS,yS = HK.local_coord(S,2*N1)
omega = [(xS)^i*xS.derivative()/(2*yS) for i in range(2)] # basis for first DeRham cohomology (around S)
wS =[U[i](xS)*xS.derivative()/(2*yS) for i in range(2)]
ED = [(xS*gder(xS)-2*i*g(xS))/(xS^(i+1))*xS.derivative()/(2*yS) for i in range(1,N1)]
d = [0,0] + ED
c = [wS[i].residue()/(1/xS*xS.derivative()/(2*yS)).residue() for i in range(2)] # coefficients of 1/x*dx/2y
w = [wS[i] - c[i]*(1/xS*xS.derivative()/(2*yS)) for i in range(2)]
FF = [[] for i in range(2)]
for i in range(2):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/d[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/d[-val].list()[0])*d[-val]
R.<x,y> = QQ[]
F = [0 for i in range(2)]
for i in range(2):
    for j in range(N1-1-i):
        F[i] += FF[i][j]*y/(x^(N1-1-i-j)) # exact parts
C = [[] for i in range(2)]
w = [wS[i] - F[i](xS,yS).derivative() - c[i]*(1/xS*xS.derivative()/(2*yS)) for i in range(2)]
for i in range(2):
    while w[i] !=0:
        val = w[i].valuation()
        C[i].append(w[i].list()[0]/omega[val].list()[val])
        w[i] = w[i] - (w[i].list()[0]/omega[val].list()[val])*omega[val] # deRham coefficients
def H(z):
    if z.scheme().base_ring() == K:
        return HK(z[0],z[1]/(z[0]*(1-7/z[0]).square_root()))
    else:
        return HL(z[0],z[1]/(z[0]*(1-7/z[0]).square_root()))
(HP1,HP2,HR1,HR2) = (H(P1),H(P2),H(R1),H(R2))
(HT1,HT2,HT3) = (H(T1),H(T2),H(T3))
(AP1,AP2) = HK.lift_x(7,all=True) # auxiliary points
x1,y1,z1 = HL.local_analytic_interpolation(AP1,HR1,prec = 2*N1)
x2,y2,z2 = HL.local_analytic_interpolation(HR2,AP2,prec = 2*N1)
def EP(i,z1,z2):
    return F[i](z2[0],z2[1]) - F[i](z1[0],z1[1])
exact_part = [EP(i,AP1,HR1) for i in range(2)]
third_kind_part = [c[i]/(2*S[1])*(((g(0)-g(x1))/(y1*x1*(S[1]+y1))*x1.derivative()).integral().truncate()(1) + Log(HR1[0]/AP1[0])) for i in range(2)]
deRham_part = [(C[i][0]*x1.derivative()/(2*y1)+C[i][1]*x1*x1.derivative()/(2*y1)).integral().truncate()(1) for i in range(2)]
AP1_to_HR1 = [exact_part[i] + third_kind_part[i] + deRham_part[i] for i in range(2)]
exact_part = [EP(i,HP1,AP1) for i in range(2)]
(div1,div2) = ([(1,S),(-1,iS)],[(1,AP1),(-1,HP1)])
HK.init_height(div1,div2)
third_kind_part = [c[i]/(2*S[1])*(HK.omega_integral(div1,div2)) for i in range(2)]
x,y = HK.monsky_washnitzer_gens()
w = HK.invariant_differential()
deRham_part = [HK.coleman_integral(C[i][0]*w+C[i][1]*x*w,HP1,AP1) for i in range(2)]
HP1_to_AP1 = [exact_part[i] + third_kind_part[i] + deRham_part[i] for i in range(2)]
P1_to_R1 = [HP1_to_AP1[i] + AP1_to_HR1[i] for i in range(2)]
exact_part = [EP(i,HR2,AP2) for i in range(2)]
third_kind_part = [-c[i]/(2*S[1])*(((g(0)-g(x2))/(y2*x2*(iS[1]+y2))*x2.derivative()).integral().truncate()(1) + Log(AP2[0]/HR2[0])) for i in range(2)]
deRham_part = [(C[i][0]*x2.derivative()/(2*y2)+C[i][1]*x2*x2.derivative()/(2*y2)).integral().truncate()(1) for i in range(2)]
HR2_to_AP2 = [exact_part[i] + third_kind_part[i] + deRham_part[i] for i in range(2)]
exact_part = [EP(i,AP2,HP1) for i in range(2)]
(div1,div2) = ([(1,iS),(-1,S)],[(1,AP2),(-1,HP1)])
HK.init_height(div1,div2)
third_kind_part = [c[i]/(2*S[1])*(HK.omega_integral(div1,div2)) for i in range(2)]
x,y = HK.monsky_washnitzer_gens()
w = HK.invariant_differential()
deRham_part = [HK.coleman_integral(C[i][0]*w+C[i][1]*x*w,AP2,HP1) for i in range(2)]
AP2_to_HP1 = [exact_part[i] + third_kind_part[i] + deRham_part[i] for i in range(2)]
R2_to_P1 = [HR2_to_AP2[i] + AP2_to_HP1[i] for i in range(2)]

################################################################################

# integration in \pi⁻¹(U2)
N2 = 40
T.<t> = PowerSeriesRing(K, 't', default_prec = 2*N2)
w = 0
for i in range(N2):
    w += binomial(-1/2,i)*(((t+7/2)^2/(2*t)-1)*((t+7/2)^2/(2*t)-2)*((t+7/2)^2/(2*t)-3)-1)^i
omega = [((t+7/2)^2/(2*t))^i*w*1/(2*t) for i in range(2)]
res = [omega[i].residue() for i in range(2)]
int = [(omega[i] - res[i]*t^-1).integral() for i in range(2)]
def t(z): return z[0] + z[1]/(L((z[0]-1)*(z[0]-2)*(z[0]-3)).sqrt()) - 7/2
def Int(i,z1,z2): return int[i](t(z2))-int[i](t(z1)) + res[i]*(Log(t(z2))-Log(t(z1)))
R1_to_P2 = [Int(i,R1,P2) for i in range(2)]
P2_to_R2 = [Int(i,P2,R2) for i in range(2)]

################################################################################

period_integral = [P1_to_R1[i] + R1_to_P2[i] + P2_to_R2[i] + R2_to_P1[i] for i in range(2)]
correction_term = [period_integral[i]/2 for i in range(2)]

################################################################################

iHP1 = HK(HP1[0],-HP1[1])
exact_part = [EP(i,iHP1,HP1) for i in range(2)]
(div1,div2) = ([(1,S),(-1,iS)],[(1,HP1),(-1,iHP1)])
HK.init_height(div1,div2)
third_kind_part = [c[i]/(2*S[1])*(HK.omega_integral(div1,div2)) for i in range(2)]
x,y = HK.monsky_washnitzer_gens()
w = HK.invariant_differential()
deRham_part = [HK.coleman_integral(C[i][0]*w+C[i][1]*x*w,iHP1,HP1) for i in range(2)]
T1_to_P1 = [(exact_part[i] + third_kind_part[i] + dirham_part[i])/2 for i in range(2)]
P2_to_T0 = [Int(i,P2,T0) for i in range(2)]
BC_along_gamma = [T1_to_P1[i] + P1_to_R1[i] + R1_to_P2[i] + P2_to_T0[i] for i in range(2)] # gamma is the path
AB_T1_to_T0 = [BC_along_gamma[i] - correction_term[i] for i in range(2)]

################################################################################

(AB_T1_to_T0[0],AB_T1_to_T0[1])
### (O(a^30), O(a^32))

################################################################################
