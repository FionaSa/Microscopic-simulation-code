#pragma once

int N_particules_total = 0;
#define epislon_etoile 0.2
#define r_etoile 3
#define N_sym 27
#define R_cut 10
#define L 30
#define m 18
#define CONSTANTE_R 0.00199
double Ndl;
#define T0 300
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