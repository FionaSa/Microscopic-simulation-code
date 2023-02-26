#pragma once

int N_particules_total = 0;
double epislon_etoile = 0.2;
double r_etoile = 3;
int N_sym = 27;
int R_cut = 10;
double L = 30;
int m = 18;
double CONSTANTE_R = 0.00199;
double Ndl;
double T0 = 300;
double k = 1.380649E-23;
int m_step = 5;

struct Particle
{
    double x, y, z;
};

struct Forces
{
    double fx, fy, fz;
};

struct Cinetique
{
    double px, py, pz;
};

struct Cinetique *cin;

extern inline double calcul_energie(double r_frac6, double r_frac3);
double verlet(struct Particle *particle, int print);
void berendsen(struct Particle *particle, int print);
void non_periodique(struct Particle *particle, int print);
void periodique(struct Particle *particle, int print);

struct Forces *forces_np;

struct Forces *forces_p;