import scipy.integrate as integrate
import numpy as np
np.set_printoptions(precision=16)
sqrt = lambda z: (z)**(1/2);

RELAMBDA = 0.8458;
IMLAMBDA = -0.045185;
KAPPA = -2;

al =-1.11384351264428;
alp =-1.11516534249106;
be =-0.477251571512785;
bet =-0.47598538151937;

w1 = lambda z: sqrt((-1)*sqrt(z**2/3 +KAPPA/3 - RELAMBDA/3*z**(-1))*(3*z**2 +KAPPA- (z**2/3 +KAPPA/3 - RELAMBDA/3*z**(-1)))- IMLAMBDA)*(2*z/3 + RELAMBDA/3*z**(-2))/sqrt(z**2/3 +KAPPA/3  - RELAMBDA/3*z**(-1))/2;
S1 = integrate.quad(w1, alp, al);
w2 = lambda z: -((z**2/3 +KAPPA/3  - RELAMBDA/3*z**(-1))**(1/2)*(3*z**2 +KAPPA- (z**2/3 +KAPPA/3 - RELAMBDA/3*z**(-1)))- IMLAMBDA)**(1/2)*(2*z/3 + RELAMBDA/3*z**(-2))/(z**2/3 +KAPPA/3  - RELAMBDA/3*z**(-1))**(1/2)/2;
S2 = integrate.quad(w2, be, bet);
w3 = lambda z: -sqrt( sqrt((z**3 +KAPPA*z -RELAMBDA)**2 + IMLAMBDA**2)+IMLAMBDA)/sqrt(2);
S3 = integrate.quad(w3, -1, al);
w4 = lambda z: sqrt( sqrt((z**3 +KAPPA*z -RELAMBDA)**2 + IMLAMBDA**2)+IMLAMBDA)/sqrt(2);
S4 = integrate.quad(w4, al, be);

xi = S1[0] - S3[0];
eta = S1[0] + S4[0]+ S2[0];
