import numpy as np
from matplotlib import pyplot as pp

# mcp-pmt parameters
Vs = 2400
C = 2.86e-13
Lmcp = 1.5e-3
G = 0.78
k = 1.16e4

Qs = C*Vs
L = k*Lmcp


# time and space discretization
sizes = (50, 50)
times = np.linspace(0, 5e-8, sizes[1])
xcoord = np.linspace(0, L, sizes[0])

i0 = np.repeat(1.75e-5, sizes[1]) # constant input signal
psi0 = np.zeros(sizes[0])   # unsaturated channel

# Q0 calculation
Q0 = np.zeros(sizes[1])
for i in range(sizes[1]):
    Q0[i] = np.trapz(i0[:i], times[:i])

# g0 calculation
g0 = np.zeros(sizes)
for x in range(sizes[0]):
    for t in range(sizes[1]):
        g0[x, t] = np.exp(G * xcoord[x])


gs = [g0]

QQss = []
QwQss = []
psis = []

for i in np.concatenate((np.linspace(0.1, 1, 10), np.ones(30))):
    g = gs[-1]

    # calculate Q/Qs
    QQs = np.zeros(sizes)
    for x in range(sizes[0]):
        for t in range(sizes[1]):
            QQs[x, t] = np.trapz(i0[:t]*g[x, :t], times[:t])/Qs

    QQss.append(QQs)

    # calculate Qw/Qs
    QwQs = np.zeros(sizes[1])
    for t in range(sizes[1]):
        QwQs[t] = np.trapz(QQs[:, t], xcoord)/L

    QwQss.append(QwQs)

    # calculate psi
    psi = np.zeros(sizes)
    for x in range(sizes[0]):
        for t in range(sizes[1]):
            psi[x, t] = psi0[x] + (QwQs[t] - QQs[x, t])

    psis.append(psi)

    # calculate new g
    for x in range(sizes[0]):
        for t in range(sizes[1]):
            g[x, t] = np.exp(G*xcoord[x] + np.trapz(np.log(1 + psi[:x, t]), xcoord[:x]))

    gs.append(g)


print(gs[-1])