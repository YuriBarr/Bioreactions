# O intuito deste codigo é simular a reação:
# D + R <-> DR <-> DR*
# Onde D representa uma droga, R um receptor, DR o complexo droga-receptor e DR* o complexo ativo.
# No codigo abaixo utilizamos k1 é a taxa na qual o DR é formado e k2 é a taxa que DR é disassociado, enquanto que β e α são as taxas te ativação e desativação.

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parâmetros
k1 = 0.1   # taxa na qual o DR é formado
k2 = 1.1    # taxa de disassociação
beta = 1.0   # taxa de ativação
alfa = 0.5 # taxa de desativação

# Sistema de EDOs
def lotka_volterra(t, z):
    D, R, DR, DR_a = z
    dDdt = k2*DR - k1*D*R
    dRdt = k2*DR - k1*D*R
    dDRdt = k1*D*R + alfa*DR_a - k2*DR - beta*DR
    dDR_adt = beta*DR - alfa*DR_a
    return [dDdt, dRdt, dDRdt, dDR_adt]

# Condições iniciais
d0 = 20
r0 = 20
dr0 = 0
dra0 = 0
z0 = [d0, r0, dr0, dra0]

# Intervalo de tempo
t_span = (0, 30)
t_eval = np.linspace(0, 30, 1000)

# Resolver
sol = solve_ivp(lotka_volterra, t_span, z0, t_eval=t_eval)

# Plot temporal
plt.figure()
plt.plot(sol.t, sol.y[0], label='Droga')
plt.plot(sol.t, sol.y[1], label='Receptor livre')
plt.plot(sol.t, sol.y[2], label='droga-receptor')
plt.plot(sol.t, sol.y[3], label='ativação')
plt.xlabel('Tempo')
plt.ylabel('concentração')
plt.legend()
plt.show()

# Espaço de fase
plt.figure()
plt.plot(sol.y[2], sol.y[3])
plt.xlabel('Complexo droga-receptor')
plt.ylabel('Efeito')
plt.title('Espaço de Fase')
plt.show()
