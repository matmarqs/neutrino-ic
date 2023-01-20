# modulos
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import odeint
from cmath import rect

# parametros
import parametros as p
sqrt2 = np.sqrt(2)
NUM_IT = round((p.t_fin - p.t_ini) / p.passo)   # numero de iteracoes na EDO

def main():
    # matriz no vacuo
    H0re, H0im = calculaMatriz()

    # interpolacao
    r, elecdens = np.loadtxt("./elecdens.txt", comments='#', unpack=True)
    N_e = interp1d(r, elecdens, kind='cubic')
    # N = len(r)

    # EDO
    t  = np.array([p.t_ini + k * p.passo for k in range(NUM_IT + 1)])
    y0 = np.array([p.re1, p.re2, p.re3, p.im1, p.im2, p.im3])
    y  = odeint(func, y0, t,   # talvez trocar pro ivp
                 atol=p.eps_abs, rtol=p.eps_rel,
                 args = (H0re, H0im, N_e))

    # printando
    d  = "{:" + str(int(np.log10(NUM_IT)) + 1) + "d}"
    fl = "{:15.5e}"
    for k in range(NUM_IT + 1):
        print((d + 8*fl).format(k, t[k], y[k][0], y[k][1], y[k][2],
                                    y[k][3], y[k][4], y[k][5], norm(y[k])))
'''
rtol, atol
The input parameters rtol and atol determine the error control performed by the solver.
The solver will control the vector, e, of estimated local errors in y, according to an inequality
of the form max-norm of (e / ewt) <= 1, where ewt is a vector of positive error weights computed as
ewt = rtol * abs(y) + atol. rtol and atol can be either vectors the same length as y or scalars.
Defaults to 1.49012e-8.
'''


def norm(vec):
    sum = 0.0
    for ele in vec:
        sum += ele*ele
    return sum


def calculaMatriz():
    th12 = graustorad(p.theta12)

    if p.num_nu == 2:
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
    M[0][0] = d1 / (2.0 * p.energ)
    M[1][1] = d2 / (2.0 * p.energ)
    M[2][2] = d3 / (2.0 * p.energ)

    H0 = np.matmul(np.matmul(U, M), np.conj(np.transpose(U)))
    H0re = np.real(H0)
    H0im = np.imag(H0)
    return H0re, H0im


def D(t, N_e):
    return sqrt2 * p.G_F * N_e(t)


def func(y, t, H0re, H0im, N_e):
    f = np.array([  H0im[0][0]*y[0] + H0im[0][1]*y[1] + H0im[0][2]*y[2]  +  H0re[0][0]*y[3] + H0re[0][1]*y[4] + H0re[0][2]*y[5] + D(t, N_e)*y[3],
                    H0im[1][0]*y[0] + H0im[1][1]*y[1] + H0im[1][2]*y[2]  +  H0re[1][0]*y[3] + H0re[1][1]*y[4] + H0re[1][2]*y[5],
                    H0im[2][0]*y[0] + H0im[2][1]*y[1] + H0im[2][2]*y[2]  +  H0re[2][0]*y[3] + H0re[2][1]*y[4] + H0re[2][2]*y[5],
                  - H0re[0][0]*y[0] - H0re[0][1]*y[1] - H0re[0][2]*y[2]  +  H0im[0][0]*y[3] + H0im[0][1]*y[4] + H0im[0][2]*y[5] - D(t, N_e)*y[0],
                  - H0re[1][0]*y[0] - H0re[1][1]*y[1] - H0re[1][2]*y[2]  +  H0im[1][0]*y[3] + H0im[1][1]*y[4] + H0im[1][2]*y[5],
                  - H0re[2][0]*y[0] - H0re[2][1]*y[1] - H0re[2][2]*y[2]  +  H0im[2][0]*y[3] + H0im[2][1]*y[4] + H0im[2][2]*y[5]])
    return f


def graustorad(a):
    return (a * np.pi) / 180


if __name__ == "__main__":
    main()
