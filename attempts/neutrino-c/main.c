#include "main.h"
#include "param.h"

/* DEFINITIONS */
void multiplyConst(double *Ne, int N);

/* check if norm = 1, only for debugging */
double norm(double Psi[DIM]) {
    double calc = 0;
    for (int i = 0; i < DIM; i++)
        calc += Psi[i] * Psi[i];
    return sqrt(calc);
}

/* check if norm = 1, only for debugging */
double sobrev(double Psi[DIM]) {
    return Psi[0]*Psi[0] + Psi[3]*Psi[3];
}


int main(int argc, char *argv[]) {
    FILE *elecdens = fopen("elecdens.txt", "r");
    FILE *energy = fopen("8b-energy.txt", "r");
    FILE *distr = fopen("8b-distr.txt", "r");

                            /*************************/
                            /***   READING DATA    ***/
                            /*************************/

    double *x, *Ne;
    int N = readalloc(elecdens, &x, &Ne, 4800); /* we have 4726 points of data */
    multiplyConst(Ne, N);
    double *E, *p_E;
    double *r0, *p_r0;
    /*int num_E = */readalloc(energy, &E, &p_E, 900);
    /*int num_r = */readalloc(distr, &r0, &p_r0, 1300);
    fclose(elecdens); fclose(energy); fclose(distr);

                            /*************************/
                            /***   INITIALIZING    ***/
                            /*************************/

    /* generating all parameters and allocating memory */
    par *params = genAllocParams(x, Ne, N);

    /* sistema de EDO */            /* jac */
    gsl_odeiv2_system ode_sys = {func, NULL, DIM, params};
    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&ode_sys, METHOD, PASSO, EPS_ABS, EPS_REL);

    /* declaring */
    int status = GSL_SUCCESS;
    double ti, t, t_inic;
    char *format = "%15.5e%15.5e%15.5e%15.5e\n";


                            /*************************/
                            /***        EDO        ***/
                            /*************************/

    int r_ini = 1200, E_ini = 400;
    for (int i_r = r_ini; i_r < r_ini + 6; i_r++) {
        t_inic = r0[i_r];
        for (int i_E = E_ini; i_E < E_ini + 6; i_E++) {
            t = t_inic; params->energ = E[i_E];
            double Psi[DIM] = { RE1, RE2, RE3, IM1, IM2, IM3 }; /* initial condition */
            long num_it = lround((T_FINAL - t_inic) / PASSO);
            //int it_width = (int) log10(num_it) + 1;
            for (long i = 0; i <= num_it; i++) {
                ti = i * PASSO + t_inic;
                if ((status = gsl_odeiv2_driver_apply(driver, &t, ti, Psi)) != GSL_SUCCESS) {
                    printf ("Error, return value = %d\n", status);
                    break;
                }
                //printf("%*d%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e%15.5e\n", it_width,
                //        i,t,Psi[0],Psi[1],Psi[2],Psi[3],Psi[4],Psi[5],norm(Psi));
            }
            printf(format, t_inic, params->energ, sobrev(Psi), norm(Psi));
        }
    }



                            /*************************/
                            /* FREEING THE RESOURCES */
                            /*************************/

    /* params */
    freeParams(params);
    /* EDO */
    gsl_odeiv2_driver_free(driver);
    /* data */
    free(x); free(Ne);
    free(E); free(p_E);
    free(r0); free(p_r0);

    if (status != GSL_SUCCESS)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;

}

void multiplyConst(double *Ne, int N) {
    for (int i = 0; i < N; i++)
        Ne[i] *= M_SQRT2 * G_F * N_A;
}

/* gera os parametros que precisamos e retorna o pointer dele */
par *genAllocParams(double *x, double *Ne, int N) {
    par *params = malloc(sizeof(par));
    params->acc = gsl_interp_accel_alloc();
    params->spline = gsl_spline_alloc(gsl_interp_steffen, N);
    gsl_spline_init(params->spline, x, Ne, N);
    params->energ = 1.0;    /* so inicializando, a energia ainda sera iterada */

    /* defining the parameters */
    double th12 = DEG_TO_RAD(THETA12),
           th23, th13, d_CP;

    /* checking number of neutrinos and setting parameters accordingly */
    if (NUM_NU == 2)
        th23 = th13 = d_CP = 0.0;
    else
        th23 = DEG_TO_RAD(THETA23), th13 = DEG_TO_RAD(THETA13), d_CP = DEG_TO_RAD(DELTACP);

    /* calculating mixing matrix */
    double s12 = sin(th12), c12 = cos(th12),
           s23 = sin(th23), c23 = cos(th23),
           s13 = sin(th13), c13 = cos(th13);
    gsl_complex
    U11=REAL(c12 * c13),                        U12=REAL(s12 * c13),                        U13=POLAR(s13,-d_CP),
    U21=ADD(POLAR(-c12*s23*s13,d_CP),-s12*c23), U22=ADD(POLAR(-s12*s23*s13,d_CP),c12*c23),  U23=REAL(s23*c13),
    U31=ADD(POLAR(-c12*c23*s13,d_CP),s12*s23),  U32=ADD(POLAR(-s12*c23*s13,d_CP),-c12*s23), U33=REAL(c23*c13);

    gsl_matrix_complex *U = gsl_matrix_complex_alloc(GERAC, GERAC);
    MAT(U, 0, 0, U11);  MAT(U, 0, 1, U12);  MAT(U, 0, 2, U13);
    MAT(U, 1, 0, U21);  MAT(U, 1, 1, U22);  MAT(U, 1, 2, U23);
    MAT(U, 2, 0, U31);  MAT(U, 2, 1, U32);  MAT(U, 2, 2, U33);

    gsl_matrix_complex *M = gsl_matrix_complex_alloc(GERAC, GERAC);
    gsl_matrix_complex_set_zero(M);
    double d1 = -DM2_21, d2 = 0.0, d3 = DM2_32; /* de acordo com o paper da Gonzalez */
    MAT(M, 0, 0, REAL(d1));
                            MAT(M, 1, 1, REAL(d2));
                                                    MAT(M, 2, 2, REAL(d3));

    params->h0 = gsl_matrix_complex_alloc(GERAC, GERAC);
    action(U, M, params->h0);   /* h0 = U . M . U^dagger */
    /* h0 = H0_re + i H0_im */
    params->H0_re = gsl_matrix_alloc(GERAC, GERAC); params->H0_im = gsl_matrix_alloc(GERAC, GERAC);

    /* colocando as partes reais e imaginarias */
    for (int i = 0; i < GERAC; i++) {
        for (int j = 0; j < GERAC; j++) {
            gsl_matrix_set(params->H0_re, i, j, GSL_REAL(gsl_matrix_complex_get(params->h0, i, j)));
            gsl_matrix_set(params->H0_im, i, j, GSL_IMAG(gsl_matrix_complex_get(params->h0, i, j)));
        }
    }

    /* liberando a memoria das matrizes auxiliares */
    gsl_matrix_complex_free(M);
    gsl_matrix_complex_free(U);

    return params;
}

void freeParams(par *params) {
    /* matrices */
    gsl_matrix_free(params->H0_re); gsl_matrix_free(params->H0_im);
    gsl_matrix_complex_free(params->h0);
    /* interpol */
    gsl_spline_free(params->spline); gsl_interp_accel_free(params->acc);
    free(params);
}

/* dPsi/dt = f(Psi, t) */
int func(double t, const double Psi[], double f[], void *params) {
    par *param = (par *) params;
    gsl_matrix *H0_re = param->H0_re;
    gsl_matrix *H0_im = param->H0_im;
    double c = 1.0 / (2.0 * param->energ);
                                /* notice: the DENSITY TERM always corresponds to R11 = H0re(0, 0) */
/*                  J11                 J12                 J13                   R11                 R12                 R13               DENSITY TERM  */
f[0] = c * (  H0im(0, 0)*Psi[0] + H0im(0, 1)*Psi[1] + H0im(0, 2)*Psi[2]  +  H0re(0, 0)*Psi[3] + H0re(0, 1)*Psi[4] + H0re(0, 2)*Psi[5]) + D(t, param)*Psi[3];
/*                  J21                 J22                 J23                   R21                 R22                 R23       */
f[1] = c * (  H0im(1, 0)*Psi[0] + H0im(1, 1)*Psi[1] + H0im(1, 2)*Psi[2]  +  H0re(1, 0)*Psi[3] + H0re(1, 1)*Psi[4] + H0re(1, 2)*Psi[5]);
/*                  J31                 J32                 J33                   R31                 R32                 R33       */
f[2] = c * (  H0im(2, 0)*Psi[0] + H0im(2, 1)*Psi[1] + H0im(2, 2)*Psi[2]  +  H0re(2, 0)*Psi[3] + H0re(2, 1)*Psi[4] + H0re(2, 2)*Psi[5]);
/*                 -R11                -R12                -R13                   J11                 J12                 J13               DENSITY TERM  */
f[3] = c * (- H0re(0, 0)*Psi[0] - H0re(0, 1)*Psi[1] - H0re(0, 2)*Psi[2]  +  H0im(0, 0)*Psi[3] + H0im(0, 1)*Psi[4] + H0im(0, 2)*Psi[5]) - D(t, param)*Psi[0];
/*                 -R21                -R22                -R23                   J21                 J22                 J23       */
f[4] = c * (- H0re(1, 0)*Psi[0] - H0re(1, 1)*Psi[1] - H0re(1, 2)*Psi[2]  +  H0im(1, 0)*Psi[3] + H0im(1, 1)*Psi[4] + H0im(1, 2)*Psi[5]);
/*                 -R31                -R32                -R33                   J31                 J32                 J33       */
f[5] = c * (- H0re(2, 0)*Psi[0] - H0re(2, 1)*Psi[1] - H0re(2, 2)*Psi[2]  +  H0im(2, 0)*Psi[3] + H0im(2, 1)*Psi[4] + H0im(2, 2)*Psi[5]);
                                /*******************************************************************/
    return GSL_SUCCESS;
}

/* calculates D = sqrt(2) . G_F . N_e, where N_e = exp(log N_e) is interpolated */
double D(double t, void *params) {  /* void pointer because this defines a gsl_function */
    par *param = (par *) params;    /* for numeric differentiation, we use gsl_function */
    return gsl_spline_eval(param->spline, t, param->acc);
}

/* evaluates C = A B A^dagger */
void action(const gsl_matrix_complex *A, const gsl_matrix_complex *B, gsl_matrix_complex *C) {
    gsl_matrix_complex *A_mult_B = gsl_matrix_complex_alloc(A->size1, B->size2);
    /* zgemm does: C = alpha op(A) op(B) + beta C */
                                                  /*  alpha                          beta  */
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,   GSL_COMPLEX_ONE, A, B,        GSL_COMPLEX_ZERO, A_mult_B); /* A_mult_B = A . B */
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, A_mult_B, A, GSL_COMPLEX_ZERO, C);        /* C = A.B.A^dagger */

    gsl_matrix_complex_free(A_mult_B);  /* freeing the memory of the auxiliary variable */
}

/* reads data from standard input, storages them in x_ptr and y_ptr and returns N */
int readalloc(FILE *stream, double **x_ptr, double **y_ptr, int chunk) {
    char *line = NULL;
    size_t len = 0;
    int i = 0, N = 0;
    double *x = (double *) malloc(chunk * sizeof(double));
    double *y = (double *) malloc(chunk * sizeof(double));
    /* getting input */
    while (getline(&line, &len, stream) != -1) {
        if (line[0] == '#')
            continue;
        else {
            sscanf(line, "%lf%lf", &x[N], &y[N]);
            i++; N++;
            if (i > chunk-1) {
                x = (double *) realloc(x, (N+chunk)*sizeof(double));
                y = (double *) realloc(y, (N+chunk)*sizeof(double));
                i = 0;
            }
        }
    }
    /* resizing the arrays to correct size */
    *x_ptr = (double *) realloc(x, N * sizeof(double));
    *y_ptr = (double *) realloc(y, N * sizeof(double));
    /* freeing resources */
    free(line);
    return N;
}
