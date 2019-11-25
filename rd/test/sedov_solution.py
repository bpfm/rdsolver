import numpy as np
import matplotlib.pyplot as plt

t = 0.5
n = 200

energy = 1.0                  # 0.47
rho1 = 10.0

lam = 1.1527

gamma = 5.0 / 3.0
nu = 3.0

gammap1 = gamma + 1
gammam1 = gamma - 1
gammam2 = gamma - 2

nup2 = nu + 2

alpha2 = -gammam1 / (2 * gammam1 + nu)

alpha1 = nup2 * gamma / (2 + nu * gammam1) * (2 * nu * (-gammam2) / (gamma * nup2**2) - alpha2)

alpha3 = nu / (2 * gammam1 + nu)

alpha4 = alpha1 * nup2 / (-gammam2)

alpha5 = 2.0 / gammam2



rho2 = gammap1 / gammam1 * rho1
r2 = lam * (energy / rho1)**(1.0 / nup2) * t**(2.0 / nup2)


V0 = 2 / (nup2 * gamma)
V1 = 4 / (nup2 * gammap1)

V = np.arange(n)

V = V / (float(n) - 1) * (V1 - V0) + V0

r = r2 * (nup2 * gammap1 / 4 * V)**(-2 / nup2) * (gammap1 / gammam1 *(nup2 * gamma / 2 * V - 1))**(-alpha2) * (nup2 * gammap1 / (nup2 * gammap1 - 2 * (2 + nu * gammam1)) * (1 - (2 + nu * gammam1) / 2 * V))**(-alpha1)

rho = rho2 * (gammap1 / gammam1 * (nup2 * gamma / 2 * V - 1))**alpha3 * (gammap1 / gammam1 * (1 - nup2 / 2 * V))**alpha5 * (nup2 * gammap1 / (nup2 * gammap1 - 2 * (2 + nu * gammam1)) * (1 - (2 + nu * gammam1) / 2 * V))**alpha4

plt.scatter(r,rho)
plt.show()