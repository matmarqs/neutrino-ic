def teste2gen():
    # matriz no vacuo
    H0re, H0im = calculaMatriz(p.num_nu)

    dm2 = p.dm2_21
    th = graustorad(p.theta12)

    # densidade eletronica
    r, elecdens = np.loadtxt(p.elecdens_test, comments='#', unpack=True)
    r, elecdens, D = interpolate_data(r, elecdens)

    # distribuicao de energia
#   E, p_E = np.loadtxt("./8b-energy.txt", comments='#', unpack=True)
#   energ = p.MeV * media(E)   # MeV = 10**6 eV
    energ = 6.44e+6 # unidade unit_E

    cteD = 2 * energ * elecdens_media
    sin2THm = (dm2**2 * np.sin(2*th))**2 / ( (cteD - dm2 * np.cos(2*th))**2 + (dm2 * np.sin(2*th))**2 )
    Delta = np.sqrt( (cteD - dm2 * np.cos(2*th))**2 + (dm2 * np.sin(2*th))**2 )
    P_exact = lambda t : 1 - sin2THm * np.sin(Delta * t / (4*energ))**2

    # EDO
    y0 = np.array([p.re1, p.re2, p.re3, p.im1, p.im2, p.im3])
    sol = solve_ivp(func, (0, 0.2), y0,
                    method=p.metodo, args = (H0re, H0im, D, energ),   # Runge-Kutta ordem 8
                    atol=p.eps_abs, rtol=p.eps_rel)
    R = sol.t
    y = np.transpose(sol.y)
    sobrev_array = np.array([sobrev(y[k]) for k in range(len(y))])
    print(elecdens_media)

    plt.plot(R, P_exact(R), label=r'$P_{exact}$')
    plt.plot(R, sobrev_array, label=r'$P_{num}$')
    plt.xlabel(r'$t$', fontsize=20)
    plt.ylabel(r'$P_e(t)$', fontsize=20)
    #plt.ylim(0, 1)
    plt.legend(fontsize=14)
    plt.title(r'2 $\nu$ - Interpol vs Analitico')
    plt.savefig("2genteste.png", dpi=300, format='png', bbox_inches="tight")

# https://stackoverflow.com/questions/4265988/generate-random-numbers-with-a-given-numerical-distribution

