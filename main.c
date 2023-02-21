#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

inline double calcul_energie(double r_frac6, double r_frac3)
{
    return epislon_etoile * (r_frac6 - 2 * r_frac3);
}

double non_periodique(struct Particle *particle)
{
    double energie = 0;

    forces_np = malloc(sizeof(struct Forces) * N_particules_total);

    double sum = 0;

    // init forces
    for (int i = 0; i < N_particules_total; i++)
    {
        forces_np[i].fx = 0;
        forces_np[i].fy = 0;
        forces_np[i].fz = 0;
    }

    for (int i = 0; i < N_particules_total; i++)
    {
        for (int j = i + 1; j < N_particules_total; j++)

        {

            double x = particle[i].x - particle[j].x;

            double y = particle[i].y - particle[j].y;

            double z = particle[i].z - particle[j].z;

            double rij2 = (x) * (x) + (y) * (y) + (z) * (z);

            double r_frac = (r_etoile * r_etoile) / rij2;

            // Calcul énergie

            double r_frac3 = r_frac * r_frac * r_frac;

            double r_frac6 = r_frac3 * r_frac3;

            energie += calcul_energie(r_frac6, r_frac3);

            // Calcul forces

            double r_frac4 = r_frac * r_frac * r_frac * r_frac;

            double r_frac7 = r_frac4 * r_frac * r_frac * r_frac;

            double calc = -48 * epislon_etoile * (r_frac7 - r_frac4);

            forces_np[i].fx += calc * (x);
            forces_np[i].fy += calc * (y);
            forces_np[i].fz += calc * (z);

            forces_np[j].fx += calc * -(x);
            forces_np[j].fy += calc * -(y);
            forces_np[j].fz += calc * -(z);
        }
    }

    energie *= 4;

    //! Mettre dans la boucle principale
    for (int i = 0; i < N_particules_total; i++)
    {
        sum += forces_np[i].fx + forces_np[i].fy + forces_np[i].fz;
    }
    printf("Somme des forces = %e %lf\n", sum, sum);

    return energie;
}

double calcul_energie_periodique(struct Particle *particle)
{
    double Sym[] = {0, 0, 0, 0, 0, L, 0, L, 0, 0, L, L, L, 0, 0, L, 0, L, L, L, 0, L, L, L, L, L, -L, L, -L, L, L, -L, -L, -L, L, L, -L, L, -L, -L, -L, L, -L, -L, -L, 0, 0, -L, 0, -L, 0, 0, -L, -L, -L, 0, 0, -L, 0, -L, -L, -L, 0, L, 0, -L, L, -L, 0, -L, 0, L, 0, L, -L, 0, -L, L, -L, L, 0};

    double energie = 0;

    double forces = 0;

    for (int p = 0; p < N_sym; p = p + 3)
    {

        for (int i = 0; i < N_particules_total; i++)
        {
            for (int j = i + 1; j < N_particules_total; j++)
            {

                double x = particle[j].x + Sym[p];
                double y = particle[j].y + Sym[p + 1];
                double z = particle[j].z + Sym[p + 2];

                double rij = (particle[i].x - x) * (particle[i].x - x) + (particle[i].y - y) * (particle[i].y - y) + (particle[i].z - z) * (particle[i].z - z);

                // energie += epislon_etoile * (1.0 / rij12 - 2 / rij6);

                if (rij < (R_cut * R_cut))
                {
                    double r_frac = (r_etoile * r_etoile) / rij;

                    double r_frac3 = r_frac * r_frac * r_frac;

                    double r_frac6 = r_frac3 * r_frac3;

                    energie += (r_frac6 - 2 * r_frac3);

                    double r_frac4 = r_frac * r_frac * r_frac * r_frac;

                    double r_frac7 = r_frac4 * r_frac * r_frac * r_frac;

                    double calc = -48 * epislon_etoile * (r_frac7 - r_frac4);

                    forces += calc * (x);
                    forces += calc * (y);
                    forces += calc * (z);

                    forces += calc * -(x);
                    forces += calc * -(y);
                    forces += calc * -(z);
                }
            }
        }
    }

    energie *= (4 * epislon_etoile) / 2;

    printf("Energie périodique = %lf Forces périodiques = %lf %e\n", energie, forces, forces);

    return energie;
}

double fonction_signe(double num)
{
    if (num <= 0)
        return 1;
    else
        return -1;
}

double calcul_energie_cinetique(struct Particle *particle)
{
    srand(time(NULL));
    Ndl = 3 * (N_particules_total - 3);

    double energie = 0;

    cin = malloc(sizeof(cin) * N_particules_total + 1);

    double Px = 0, Py = 0, Pz = 0;

    double vi = sqrt(((Ndl * 2) / (Ndl * m)) * CONSTANTE_R * T0);

    double kb = 0.5 * Ndl * m * pow(vi, 2);

    printf("vi = %e, k = %e kb = %e vrai = %e\n", vi, k, kb, CONSTANTE_R * T0);

    for (int i = 0; i < N_particules_total; i++)
    {
        cin[i].px = fonction_signe(0.5 - ((double)rand() / (double)(RAND_MAX))) * ((double)rand() / (double)(RAND_MAX)) * vi;
        Px += cin[i].px;
        cin[i].py = fonction_signe(0.5 - ((double)rand() / (double)(RAND_MAX))) * ((double)rand() / (double)(RAND_MAX)) * vi;
        Py += cin[i].py;

        cin[i].pz = fonction_signe(0.5 - ((double)rand() / (double)(RAND_MAX))) * ((double)rand() / (double)(RAND_MAX)) * vi;
        Pz += cin[i].pz;

        printf("Temp = %lf %lf %lf %lf  \n", cin[i].px, cin[i].py, cin[i].pz, vi);
    }

    double temperature = 0;

    double a = (1 / (2 * 0.0001 * 4.186));

    double energie_cin = 1;

    Px /= N_particules_total;

    Py /= N_particules_total;

    Pz /= N_particules_total;

    double energie0 = Ndl * CONSTANTE_R * T0;

    for (int i = 0; i < N_particules_total; i++)
    {
        double x = cin[i].px;

        double y = cin[i].py;

        double z = cin[i].pz;

        //  printf("x = %lf , y = %lf , z = %lf\n", x, y, z);

        energie += ((x * x) + (y * y) + (z * z)) / m;
    }

    energie *= a;

    double RAPPORT = sqrt(energie0 / energie);

    energie = 0;

    for (int i = 0; i < N_particules_total; i++)
    {
        double x = cin[i].px;

        double y = cin[i].py;

        double z = cin[i].pz;
        // energie_cin = a * energie;

        energie += ((x * x) + (y * y) + (z * z)) / m;
        temperature = (1 / (Ndl * CONSTANTE_R)) * energie_cin;
        energie_cin = a * energie;

        // 3
        cin[i + 1].px *= RAPPORT;
        cin[i + 1].py *= RAPPORT;
        cin[i + 1].pz *= RAPPORT;

        // 4
        cin[i + 1].px -= Px;
        cin[i + 1].py -= Py;
        cin[i + 1].pz -= Pz;
    }
    // printf("Temp = %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf\n", cin[i].px, cin[i].py, cin[i].pz, temperature, energie_cin, Px, Py, Pz, energie, a);

    /*cin[i + 1].px *= sqrt(temperature / T0);
    cin[i + 1].py *= sqrt(temperature / T0);
    cin[i + 1].pz *= sqrt(temperature / T0);*/

    printf("ENERGIE CINETIQUE\n\tEnergie: %lf\n\tTempérature: %lf\n", energie_cin, temperature);

    return energie_cin;
}

double calcul_energie_cinetique_berendsen(struct Particle *particle)
{
    srand(time(NULL));
    Ndl = 3 * (N_particules_total - 3);

    double energie = 0;

    double temperature = 0;

    double energie_cin = 0;

    // (1.0 / (2.0 * 0.0001 * 4.186))

    double b = 1194.4577162;

    for (int i = 0; i < N_particules_total; i++)
    {

        double x = cin[i].px;
        double y = cin[i].py;
        double z = cin[i].pz;

        energie += ((x * x) + (y * y) + (z * z)) / m;

        energie_cin = b * energie;
        temperature = (1 / (Ndl * CONSTANTE_R)) * energie_cin;

        // 6
        double a = 0.01 * ((T0 / temperature) - 1);

        cin[i + 1].px += a * cin[i + 1].px;
        cin[i + 1].py += a * cin[i + 1].py;
        cin[i + 1].pz += a * cin[i + 1].pz;
    }

    printf("BERENDSEN\n\tEnergie cinétique: %lf\n", energie_cin);
    printf("\tTempérature: %lf\n", temperature);

    return energie;
}

struct Particle *read_file(struct Particle *particle)
{
    FILE *in = fopen("./data/particule.xyz", "r");
    if (in == NULL)
        return perror("Failed: "), -1;

    double x, y, z;

    if (fscanf(in, " %d %d", &x, &y) == EOF)
        return -1;

    int i = 0;

    while (fscanf(in, " %d %lf %lf %lf", &x, &x, &y, &z) != EOF)
    {
        (particle[i]).x = x;
        (particle[i]).y = y;
        (particle[i]).z = z;
        // printf("%lf %lf %lf\n", (particle[i]).x, (particle[i]).y, (particle[i]).z);
        i++;
        particle = realloc(particle, (i + 1) * (sizeof(struct Particle)));
    }

    N_particules_total = i;

    printf("Particules totales = %d\n", N_particules_total);

    return particle;
}

int main(int argc, char const *argv[])
{
    struct Particle *particle = malloc(sizeof(*particle));

    if ((particle = read_file(particle)) == -1)
        return printf("Erreur lecture des coordonnées\n"), -1;

    printf("%lf\n", particle[0].x);

    printf("Energie = %lf\n", non_periodique(particle));
    printf("Energie périodique = %lf\n", calcul_energie_periodique(particle));

    double energie_cin = calcul_energie_cinetique(particle);

    double energie_b = calcul_energie_cinetique_berendsen(particle);

    return 0;
}
