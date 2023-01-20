# modulos
import numpy as np
#from scipy.interpolate import UnivariateSpline
from scipy.integrate import solve_ivp
from cmath import rect

from matplotlib import pyplot as plt
from matplotlib import rc
plt.style.use('bmh')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# parametros
import parametros as p
fl = "{:15.5e}"

def main(part, rang):
    # matriz no vacuo
    H0re_2E, H0im_2E = calculaMatriz(p.num_nu)  # Ainda sem dividir por 2 * E

    # densidade eletronica
    #r, log10elecdens = np.loadtxt(p.elecdens_file, comments='#', unpack=True)
    #D = interpolate_data(r, log10elecdens, p.R_sun)

    #energ = 6.44e+6 # energia mais provavel em eV
    #print((1 / (2*energ)) * H0re_2E)
    #print("D(1.5985) = ", D(1.5985))
    #quit()

    # distribuicao do ponto de producao
    #r_prod, p_prod = np.loadtxt("./8b-distr.txt", comments='#', unpack=True)
    #p_prod = p_prod / sum(p_prod)   # normalizando p_prod

    #energias = np.linspace(0.2e+6, 15e+6, 30)
    #lista_sobrev = []
    for t_ini in [0]:
        for _ in range(1):
            #energ = 6.44e+6 # energia mais provavel em eV
            H0re = 2 * p.factor * H0re_2E   # fator 2 que salva a vida
            H0im = 2 * p.factor * H0im_2E

            # EDO
            y0 = np.array([p.re1, p.re2, p.re3, p.im1, p.im2, p.im3])
            sol = solve_ivp(func, (t_ini, p.t_fin), y0,
                            method=p.metodo, args = (H0re, H0im),
                            atol=p.eps_abs, rtol=p.eps_rel)#, jac=jacob)
            t = sol.t
            y = np.transpose(sol.y)

            #'''
            # printando a probabilidade de sobrevivencia ao sair do Sol
            #print((4*fl).format(t_ini, energ, sobrev(y[-1]), norm(y[-1])))
            #lista_sobrev.append(sobrev(y[-1]))
            # o output sao tres colunas da forma
            #                   t_ini  energ  p_sobrev       norma ?= 1
            #'''

            '''
            # printar todos os steps
            length = len(t)
            d  = "{:" + str(int(np.log10(length)) + 1) + "d}"   # numero de iteracoes
            for k in range(length):
                print((d + 8*fl).format(k, t[k], y[k][0], y[k][1], y[k][2],
                                        y[k][3], y[k][4], y[k][5], norm(y[k])))
            '''

            #'''
            # fazendo o plot da probabilidade de P_e(t)
            sobrev_array = np.array([sobrev(y[k]) for k in range(len(y))])
            mu_array = np.array([mu(y[k]) for k in range(len(y))])
            tau_array = np.array([tau(y[k]) for k in range(len(y))])
            plt.plot(t , sobrev_array, label=p.label_part(part, 'e'))
            plt.plot(t , mu_array, label=p.label_part(part, '\\mu'))
            plt.plot(t , tau_array, label=p.label_part(part, '\\tau'))
            plt.xlabel(r'$L/E$ (km/GeV)', fontsize=20)
            plt.ylabel(r'Probabilidade', fontsize=20)
            plt.ylim(0, 1.05)
            plt.legend(fontsize=16)
            plt.title(p.titulo(part))
            plt.savefig(str(p.num_nu)+"nu-"+p.metodo+"-"+part+rang+".png", dpi=300, format='png', bbox_inches="tight")
            #'''

    '''
    # fazendo o plot da sobrevivencia em funcao da energia
    array_sobrev = np.array(lista_sobrev)
    plt.plot(energias / 1e+6, array_sobrev, label=r'$P_{ee}$')
    plt.xlabel(r'$E$ (MeV)', fontsize=20)
    plt.ylabel(r'$P_{ee}$', fontsize=20)
    plt.ylim(0, 1.05)
    plt.legend(fontsize=14)
    plt.title(r'Probabilidade de sobrevivência de $\nu_e$ para ' + str(p.num_nu) + r' gerações')
    plt.savefig("sobrev-"+str(p.num_nu)+"nu.png", dpi=300, format='png', bbox_inches="tight")
    '''


#def interpolate_data(r, log10elecdens, R_sun):
#    u = p.solar * r
#    cte = R_sun * p.sqrt2 * p.G_FxN_A * (p.gamma)**3
#    v = cte * ( 10 ** log10elecdens )
#    D = UnivariateSpline(u, v, k = 3, s = 0)
#    return D

def norm(vec):
    return sum(vec * vec)

def sobrev(state):  # state = (re1, re2, re3, im1, im2, im3)
    return state[0]*state[0] + state[3]*state[3]

def mu(state):  # state = (re1, re2, re3, im1, im2, im3)
    return state[1]*state[1] + state[4]*state[4]

def tau(state):  # state = (re1, re2, re3, im1, im2, im3)
    return state[2]*state[2] + state[5]*state[5]


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
                  [-s12*c23-rect(c12*s23*s13, d_cp), c12*c23-rect(s12*s23*s13, d_cp), s23*c13],
                  [s12*s23-rect(c12*c23*s13, d_cp), -c12*s23-rect(s12*c23*s13,d_cp), c23*c13]])

    M = np.zeros((3,3), dtype=complex)
    # Aqui podemos ter que mudar a energ
    # Ordenamento (depois):
    d1, d2, d3 = -p.dm2_21, 0, p.dm2_32
    M[0][0] = d1
    M[1][1] = d2
    M[2][2] = d3

    H0 = np.matmul(np.matmul(U, M), np.conj(np.transpose(U)))
    H0re = np.real(H0)
    H0im = np.imag(H0)
    return H0re, H0im


def func(t, y, H0re, H0im):
    t = t
    f = np.array([  H0im[0][0]*y[0] + H0im[0][1]*y[1] + H0im[0][2]*y[2]  +  H0re[0][0]*y[3] + H0re[0][1]*y[4] + H0re[0][2]*y[5],#+ D(t)*y[3],
                    H0im[1][0]*y[0] + H0im[1][1]*y[1] + H0im[1][2]*y[2]  +  H0re[1][0]*y[3] + H0re[1][1]*y[4] + H0re[1][2]*y[5],
                    H0im[2][0]*y[0] + H0im[2][1]*y[1] + H0im[2][2]*y[2]  +  H0re[2][0]*y[3] + H0re[2][1]*y[4] + H0re[2][2]*y[5],
                  - H0re[0][0]*y[0] - H0re[0][1]*y[1] - H0re[0][2]*y[2]  +  H0im[0][0]*y[3] + H0im[0][1]*y[4] + H0im[0][2]*y[5],#- D(t)*y[0],
                  - H0re[1][0]*y[0] - H0re[1][1]*y[1] - H0re[1][2]*y[2]  +  H0im[1][0]*y[3] + H0im[1][1]*y[4] + H0im[1][2]*y[5],
                  - H0re[2][0]*y[0] - H0re[2][1]*y[1] - H0re[2][2]*y[2]  +  H0im[2][0]*y[3] + H0im[2][1]*y[4] + H0im[2][2]*y[5]])
    return f


#def jacob(t, y, H0re, H0im, D):
#    j = np.array([[ H0im[0][0]       ,  H0im[0][1],  H0im[0][2], H0re[0][0]+D(t), H0re[0][1], H0re[0][2]],
#                  [ H0im[1][0]       ,  H0im[1][1],  H0im[1][2], H0re[1][0]     , H0re[1][1], H0re[1][2]],
#                  [ H0im[2][0]       ,  H0im[2][1],  H0im[2][2], H0re[2][0]     , H0re[2][1], H0re[2][2]],
#                  [-H0re[0][0] - D(t), -H0re[0][1], -H0re[0][2], H0im[0][0]     , H0im[0][1], H0im[0][2]],
#                  [-H0re[1][0]       , -H0re[1][1], -H0re[1][2], H0im[1][0]     , H0im[1][1], H0im[1][2]],
#                  [-H0re[2][0]       , -H0re[2][1], -H0re[2][2], H0im[2][0]     , H0im[2][1], H0im[2][2]]])
#    return j


def graustorad(a):
    return (a * np.pi) / 180


if __name__ == "__main__":
    if p.re1 == 1.0 or p.im1 == 1.0:
        part = r'e'
    elif p.re2 == 1.0 or p.im2 == 1.0:
        part = r'\mu'
    else:
        part = r'\tau'

    if p.t_fin - p.t_ini < 10000:
        rang = '_short'
    else:
        rang = ''

    main(part, rang)
