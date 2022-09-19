################################################################################

# Elliptic Curve with LMFDB label 272.b2
R.<x> = QQ[]
X = EllipticCurve([0, 0, 0, -91, 330]) # eq'n: y²=(x-6)(x-5)(x+11), M-W group: Z x Z/2Z x Z/2Z
K = Qp(17,12)
XK = X.change_ring(K)
(P,T1,T2,T3) = (XK(-3,24),XK(6,0),XK(5,0),XK(-11,0)) # gen's of M-W group
L.<a> = K.extension(x^2-17)
XL = X.change_ring(L)
T.<t> = PowerSeriesRing(K, 't')
(P1,P2) = (XK.lift_x(1), XK.lift_x(-28)) # reference pts in vertices
(R1,R2) = XL.lift_x(a+6,all=True) # reference pts in edges
def Log(z): return L(z).log(p_branch=0, change_frac=True) # Iwasawa branch

################################################################################

# integration in \pi⁻¹(U1)
N1 = 25
w1 = 0
for i in range(N1):
    w1 += binomial(-1/2,i)*((t-17)^2/(4*t))^i
w1 = 1/(2*t)*w1
r1 = w1.residue()
int1 = (w1 - r1*t^-1).integral()
def t1(z): return 17*(z[1]-(z[0]-6)*L(1+17/(z[0]-6)).sqrt())/(z[1]+(z[0]-6)*L(1+17/(z[0]-6)).sqrt())
def Int1(z1,z2): return int1(t1(z2))-int1(t1(z1)) + r1*(Log(t1(z2))-Log(t1(z1)))

################################################################################

# integration in \pi⁻¹(U2)
N2 = 25
w2 = 0
for i in range(N2):
    w2 += binomial(-1/2,i)*((t-17/2)^2/(2*t))^i
w2 = 1/(2*t)*w2
r2 = w2.residue()
int2 = (w2 - r2*t^-1).integral()
def t2(z): return z[0]-6 + z[1]/(L(z[0]-5).sqrt()) + 17/2
def Int2(z1,z2): return int2(t2(z2))-int2(t2(z1)) + r2*(Log(t2(z2))-Log(t2(z1)))

################################################################################

period_integral = Int1(P1,R1) + Int2(R1,P2) + Int2(P2,R2) + Int1(R2,P1)
CT = (period_integral)/2 # correction term

################################################################################

(Q,S) = (XK(7,6),XK(23,102))
BC_along_gamma = Int1(Q,R1) + Int2(R1,S) # gamma is the path
AB_Q_to_S = BC_along_gamma - CT; AB_Q_to_S
### 12*a^2 + 8*a^4 + 15*a^6 + 9*a^8 + 16*a^10 + 8*a^12 + 15*a^14 + 16*a^16 + 10*a^18 + O(a^20)

################################################################################

# consistincy check
formal_log = XK.formal_group().log()
p = 17
c = X.tamagawa_number(p)
mS = c*(p-1)*S
nQ = (p-1)*Q
logS = 1/(c*(p-1))*formal_log.truncate()(-mS[0]/mS[1])
logQ = 1/(p-1)*formal_log.truncate()(-nQ[0]/nQ[1])
logS - logQ
### 12*17 + 8*17^2 + 15*17^3 + 9*17^4 + 16*17^5 + 8*17^6 + 15*17^7 + 16*17^8 + 10*17^9 + 12*17^10 + O(17^11)

################################################################################

# test on points whose difference is torsion
LLL = []
for i in range(100):
    LLL.append(Int1(T2+i*P,R1) + Int2(R1,T3+i*P) - CT)
    LLL.append(Int1(T2+i*P,R1) + Int2(R1,T1+i*P) - CT)

LLL
### [O(a^20), O(a^20), O(a^20), O(a^20), O(a^20), O(a^20), O(a^20), O(a^20), O(a^20), O(a^20), O(a^20), O(a^20), O(a^20), O(a^20), O(a^20), O(a^20), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^10), O(a^10), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14), O(a^14)]

################################################################################
