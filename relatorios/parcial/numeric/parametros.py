# condicao inicial
re1 = 1.0
re2 = 0.0
re3 = 0.0
im1 = 0.0
im2 = 0.0
im3 = 0.0

t_fin = 35000
factor = 1.2669326791371    # 1.27 da conversão de unidades

# constantes (todas listadas aqui estao potencias de eV)
sqrt2 = 1.4142135623730951
#solar = 1e+3    # plotaremos para r in [0, solar * R_sun]
#gamma = 1.23984198e-04  # 1/cm = gamma * eV
## Temos que as constantes G_F e N_A so aparecem juntas da forma G_F * N_A
## G_F = 1.1663787e-23   # eV^{-2}
## N_A = 6.02214076e+23  # Avogadro
#G_FxN_A = 7.0240967108658126    # 1.1663787 * 6.02214076
#R_sun = (1/solar) * (1/gamma) * 6.957e+10   # R_sun/1000 em eV^{-1}; 6.957e+10 = raio solar em cm

# condicoes
#elecdens_file = "elecdens.txt"
num_nu = 3  # quantos neutrinos, 2 ou 3

# simulacao
metodo = 'DOP853'   # DOP853, Radau ou BDF (para o IVP)
t_ini = 0

# dados do NuFit
dm2_21 = 7.42e-5    # eV^2
dm2_32 = 2.515e-3   # eV^2
# todos os angulos em graus
theta12 = 33.44
theta23 = 49.2
theta13 = 8.57
deltacp = 194.0

# parametros da EDO
eps_abs = 1e-13  # na ordem de 1e-7 eh minimo para nao dar ruim
eps_rel = 1e-13  # na ordem de 1e-7 eh minimo para nao dar ruim

def sci(x):
    return "{:.8e}".format(x)

def label_part(part, final):
    return r'$P(\nu_{' + part + r'} \to \nu_{' + final + r'})$'

def titulo(part):
    return r'Probabilidade de $\nu_' + part + r'$ para ' + str(num_nu) + r' gerações'
