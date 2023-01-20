# modulos
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.integrate import solve_ivp
from cmath import rect

from matplotlib import pyplot as plt

# parametros
import parametros as p
sqrt2 = np.sqrt(2)
fl = "{:15.5e}"

def main():
    # matriz no vacuo
    H0re, H0im = calculaMatriz(p.num_nu)

    # densidade eletronica
    r, elecdens = np.loadtxt(p.elecdens_file, comments='#', unpack=True)
    elecdens = sqrt2 * p.G_F * p.N_A * elecdens * (p.inv_cm_to_eV)**3
    D = UnivariateSpline(r, elecdens, k = 3, s = 0)
    # distribuicao de energia
    E, p_E = np.loadtxt("./8b-energy.txt", comments='#', unpack=True)
    E = p.MeV * E   # MeV = 10**6 eV
    p_E = p_E / sum(p_E)    # normalizando p_E
    r_prod, p_prod = np.loadtxt("./8b-distr.txt", comments='#', unpack=True)
    p_prod = p_prod / sum(p_prod)   # normalizando p_prod

    # Gerando as distribuicoes
    #r_ini = r[0]; r_fin = r[-1]
    for t_ini in np.linspace(0.9, 1, num = 200):
        for _ in range(16): # sorteamos 16 energias para cada dado de raio
            # sorteando energias de acordo com p_E
            energ = np.random.choice(E, p=p_E)

            # EDO
            y0 = np.array([p.re1, p.re2, p.re3, p.im1, p.im2, p.im3])
            sol = solve_ivp(func, (t_ini, p.t_fin), y0,
                            method=p.metodo, args = (H0re, H0im, D, energ),   # Runge-Kutta ordem 8
                            atol=p.eps_abs, rtol=p.eps_rel)
            # t = sol.t
            y = np.transpose(sol.y)

            # printando a probabilidade de sobrevivencia ao sair do Sol
            print((4*fl).format(t_ini, energ, sobrev(y[-1]), norm(y[-1])))
            # o output sao tres colunas da forma
            #                   t_ini  energ  p_sobrev       norma ?= 1

            ### printar todos os steps
            # d  = "{:" + str(int(np.log10(length)) + 1) + "d}" # numero de iteracoes
            # length = len(t)
            # for k in range(length):
            #     print((d + 8*fl).format(k, t[k], y[k][0], y[k][1], y[k][2],
            #                             y[k][3], y[k][4], y[k][5], norm(y[k])))

def teste2gen():
    # matriz no vacuo
    H0re, H0im = calculaMatriz(p.num_nu)

    dm2 = p.dm2_21
    th = graustorad(p.theta12)

    # densidade eletronica
    r, elecdens = np.loadtxt(p.elecdens_test, comments='#', unpack=True)
    elecdens = sqrt2 * p.G_F * p.N_A * elecdens * (p.inv_cm_to_eV)**3
    D = UnivariateSpline(r, elecdens, k = 3, s = 0)
    elecdens_media = elecdens[0]

    # distribuicao de energia
    E, p_E = np.loadtxt("./8b-energy.txt", comments='#', unpack=True)
    energ = media(E)   # MeV = 10**6 eV


    cteD = 2 * energ * elecdens_media
    sin2THm = (dm2**2 * np.sin(2*th))**2 / ( (cteD - dm2 * np.cos(2*th))**2 + (dm2 * np.sin(2*th))**2 )
    Delta = np.sqrt( (cteD - dm2 * np.cos(2*th))**2 + (dm2 * np.sin(2*th))**2 )
    P_exact = lambda t : 1 - sin2THm * np.sin(Delta * t / (4*energ))**2

    # EDO
    y0 = np.array([p.re1, p.re2, p.re3, p.im1, p.im2, p.im3])
    sol = solve_ivp(func, (0, p.t_fin), y0,
                    method=p.metodo, args = (H0re, H0im, D, energ),   # Runge-Kutta ordem 8
                    atol=p.eps_abs, rtol=p.eps_rel)
    R = sol.t
    y = np.transpose(sol.y)
    sobrev_array = np.array([sobrev(y[k]) for k in range(len(y))])

    #for k in range(30):
    #    print(R[k], sobrev_array[k], P_exact(R[k]))

    plt.plot(R, P_exact(R), label=r'$P_{exact}$')
    plt.plot(R, sobrev_array, label=r'$P_{num}$')
    plt.xlabel(r'$t$', fontsize=20)
    plt.ylabel(r'$P_e(t)$', fontsize=20)
    #plt.ylim(0, 1)
    plt.legend(fontsize=14)
    plt.title(r'2 $\nu$ - Interpol vs Analitico')
    plt.savefig("2genteste.png", dpi=300, format='png', bbox_inches="tight")

# https://stackoverflow.com/questions/4265988/generate-random-numbers-with-a-given-numerical-distribution

def norm(vec):
    sum = 0.0
    for ele in vec:
        sum += ele*ele
    return sum


def sobrev(state):  # state = (re1, re2, re3, im1, im2, im3)
    return state[0]*state[0] + state[3]*state[3]


def calculaMatriz(num_nu):
    th12 = graustorad(p.theta12)

    if num_nu == 2:
        th13, th23, d_cp = 0.0, 0.0, 0.0
    else:
        th13, th23, d_cp = graustorad(p.theta13), graustorad(p.theta23), graustorad(p.deltacp)

    s12, c12 = np.sin(th12), np.cos(th12)
    s13, c13 = np.sin(th13), np.cos(th13)
    s23, c23 = np.sin(th23), np.cos(th23)

    # Matriz mudança de base, da base de massas para base de sabores:

    U = np.array([[c12*c23, s12*c13, rect(s13, -d_cp)],
                  [-s12*c23-rect(c12*s23*s13, d_cp), c12*c23-rect(s12*s23*s13, d_cp),  s23* c13],
                  [s12*s23-rect(c12*c23*s13, d_cp), -c12*s23-rect(s12*c23*s13,d_cp), c23*c13]])

    M = np.zeros((3,3), dtype=complex)
    # Aqui podemos ter que mudar a energ
    # Ordenamento (depois):
    d1, d2, d3 = -p.dm2_21, 0.0, p.dm2_32
    M[0][0] = d1
    M[1][1] = d2
    M[2][2] = d3

    H0 = np.matmul(np.matmul(U, M), np.conj(np.transpose(U)))
    H0re = np.real(H0)
    H0im = np.imag(H0)
    return H0re, H0im


def func(t, y, H0re, H0im, D, energ):
    c =  1.0 / (2.0 * energ)
    f = np.array([c * (  H0im[0][0]*y[0] + H0im[0][1]*y[1] + H0im[0][2]*y[2]  +  H0re[0][0]*y[3] + H0re[0][1]*y[4] + H0re[0][2]*y[5]) + D(t)*y[3],
                  c * (  H0im[1][0]*y[0] + H0im[1][1]*y[1] + H0im[1][2]*y[2]  +  H0re[1][0]*y[3] + H0re[1][1]*y[4] + H0re[1][2]*y[5]),
                  c * (  H0im[2][0]*y[0] + H0im[2][1]*y[1] + H0im[2][2]*y[2]  +  H0re[2][0]*y[3] + H0re[2][1]*y[4] + H0re[2][2]*y[5]),
                  c * (- H0re[0][0]*y[0] - H0re[0][1]*y[1] - H0re[0][2]*y[2]  +  H0im[0][0]*y[3] + H0im[0][1]*y[4] + H0im[0][2]*y[5]) - D(t)*y[0],
                  c * (- H0re[1][0]*y[0] - H0re[1][1]*y[1] - H0re[1][2]*y[2]  +  H0im[1][0]*y[3] + H0im[1][1]*y[4] + H0im[1][2]*y[5]),
                  c * (- H0re[2][0]*y[0] - H0re[2][1]*y[1] - H0re[2][2]*y[2]  +  H0im[2][0]*y[3] + H0im[2][1]*y[4] + H0im[2][2]*y[5])])
    return f


def graustorad(a):
    return (a * np.pi) / 180


def media(L):
    return sum(L)/len(L)


if __name__ == "__main__":
    if p.teste_nu2:
        teste2gen()
    else:
        main()
