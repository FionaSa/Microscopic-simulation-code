#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <errno.h>

void sauvegarde_deb(FILE *fp, int xdim, int ydim, int zdim)
{

    fprintf(fp, "CRYST1 %d %d %d 90.00 90.00 90.00 P 1\n", xdim, ydim, zdim);
}
void sauvegarde_deb2(FILE *fp, int iteration)
{
    fprintf(fp, "MODEL %d\n", iteration);
}

void sauvegarde_mid(FILE *fp, double x, double y, double z, char *k, int j)
{
    fprintf(fp, "ATOM %s C 0 %lf %lf %lf MRES\n", k, x, y, z);
}

void sauvegarde_fin(FILE *fp)
{
    fprintf(fp, "TER\nENDMDL\n");
}

void simulation_p(int iterations, struct Particle *particles)
{
    FILE *fp = fopen("output_p.pdb", "w");

    if (!fp)
        perror("Erreur alloc file");

    sauvegarde_deb(fp, L, L, L);

    double frac = 1 / m;

    double vi = sqrt(((Ndl * 2) / (Ndl * m)) * CONSTANTE_R * T0);

    for (int i = 1; i <= iterations; i++)
    {
        sauvegarde_deb2(fp, i);
        for (int j = 0; j < N_particules_total; j++)
        {
            // a_i = 1/m_i * ∑_(j≠i) f_ij

            double accx = frac * forces_p[j].fx;
            double accy = frac * forces_p[j].fy;
            double accz = frac * forces_p[j].fz;

            char *label = malloc((sizeof(char)) * 10);

            sprintf(label, "%dX%d", i, j);

            // r(t+Δt) = r(t) + v(t)Δt + 0.5a(t)Δt^2
            // Δt = i
            double x = particles[j].x + vi * (i) + 0.5 * accx * (i * i);
            double y = particles[j].y + vi * (i) + 0.5 * accy * (i * i);
            double z = particles[j].z + vi * (i) + 0.5 * accz * (i * i);

            sauvegarde_mid(fp, x, y, z, label, j);

            free(label);
        }
        sauvegarde_fin(fp);
    }

    fclose(fp);
}

void simulation_np(int iterations, struct Particle *particles)
{
    FILE *fp = fopen("output_np.pdb", "w");

    if (!fp)
        perror("Erreur alloc file");

    sauvegarde_deb(fp, 0, 0, 0);

    double frac = 1 / m;

    double vi = sqrt(((Ndl * 2) / (Ndl * m)) * CONSTANTE_R * T0);

    for (int i = 1; i <= iterations; i++)
    {
        sauvegarde_deb2(fp, i);
        for (int j = 0; j < N_particules_total; j++)
        {
            // a_i = 1/m_i * ∑_(j≠i) f_ij

            double accx = frac * forces_np[j].fx;
            double accy = frac * forces_np[j].fy;
            double accz = frac * forces_np[j].fz;

            char *label = malloc((sizeof(char)) * 10);

            sprintf(label, "%dX%d", i, j);

            // r(t+Δt) = r(t) + v(t)Δt + 0.5a(t)Δt^2
            // Δt = i
            double x = particles[j].x + vi * (i) + 0.5 * accx * (i * i);
            double y = particles[j].y + vi * (i) + 0.5 * accy * (i * i);
            double z = particles[j].z + vi * (i) + 0.5 * accz * (i * i);

            sauvegarde_mid(fp, x, y, z, label, j);

            free(label);
        }
        sauvegarde_fin(fp);
    }

    fclose(fp);
}

inline double calcul_energie(double r_frac6, double r_frac3)
{
    return epislon_etoile * (r_frac6 - 2 * r_frac3);
}

double non_periodique(struct Particle *particle)
{
    double energie = 0;

    forces_np = malloc(sizeof(struct Forces) * N_particules_total);
    forces_p = malloc(sizeof(struct Forces) * N_particules_total);

    double sum = 0;

    // init forces
    for (int i = 0; i < N_particules_total; i++)
    {
        forces_np[i].fx = 0;
        forces_p[i].fx = 0;
        forces_np[i].fy = 0;
        forces_p[i].fy = 0;
        forces_np[i].fz = 0;
        forces_p[i].fz = 0;
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

            double tmp = calc * (x);

            forces_np[i].fx += tmp;

            sum += tmp;

            tmp = calc * (y);

            forces_np[i].fy += tmp;

            sum += tmp;

            tmp = calc * (z);

            forces_np[i].fz += tmp;

            sum += tmp;

            tmp = calc * -(x);

            forces_np[j].fx += tmp;

            sum += tmp;

            tmp = calc * -(y);

            forces_np[j].fy += tmp;

            sum += tmp;

            tmp = calc * -(z);

            forces_np[j].fz += tmp;

            sum += tmp;
        }
    }

    energie *= 4;

    printf("\033[32mNON PERIODIQUE\n\tEnergie: \033[37m%lf\033[38;5;54m J\n\t\033[32mSomme des forces = \033[37m%lf\033[38;5;54m N , \033[37m%e\033[38;5;54m N\n", energie, sum, sum);

    return energie;
}

double periodique(struct Particle *particle)
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

                if (rij < (R_cut * R_cut))
                {
                    double r_frac = (r_etoile * r_etoile) / rij;

                    double r_frac3 = r_frac * r_frac * r_frac;

                    double r_frac6 = r_frac3 * r_frac3;

                    energie += (r_frac6 - 2 * r_frac3);

                    double r_frac4 = r_frac * r_frac * r_frac * r_frac;

                    double r_frac7 = r_frac4 * r_frac * r_frac * r_frac;

                    double calc = -48 * epislon_etoile * (r_frac7 - r_frac4);

                    double tmp;

                    tmp = calc * (x);

                    forces += tmp;
                    forces_p[i].fx += tmp;

                    tmp = calc * (y);

                    forces += tmp;
                    forces_p[i].fy += tmp;

                    tmp = calc * (z);

                    forces += tmp;
                    forces_p[i].fz += tmp;

                    tmp = calc * -(x);

                    forces += tmp;
                    forces_p[j].fx += tmp;

                    tmp = calc * -(y);

                    forces += tmp;
                    forces_p[j].fy += tmp;

                    tmp = calc * -(z);

                    forces += tmp;
                    forces_p[j].fz += tmp;
                }
            }
        }
    }

    energie *= (4 * epislon_etoile) / 2;

    printf("\033[38;5;40mPERIODIQUE\n\tEnergie = \033[37m%lf\033[38;5;54m J\n\t\033[38;5;40mSomme des forces = \033[37m%lf\033[38;5;54m N, \033[37m%e\033[38;5;54m N\n", energie, forces, forces);

    return energie;
}

double fonction_signe(double num)
{
    if (num <= 0.5)
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

    // printf("vi = %e, k = %e kb = %e vrai = %e\n", vi, k, kb, CONSTANTE_R * T0);

    for (int i = 0; i < N_particules_total; i++)
    {
        cin[i].px = fonction_signe(((double)rand() / (double)(RAND_MAX))) * ((double)rand() / (double)(RAND_MAX)) * vi;
        Px += cin[i].px;

        cin[i].py = fonction_signe(((double)rand() / (double)(RAND_MAX))) * ((double)rand() / (double)(RAND_MAX)) * vi;
        Py += cin[i].py;

        cin[i].pz = fonction_signe(((double)rand() / (double)(RAND_MAX))) * ((double)rand() / (double)(RAND_MAX)) * vi;
        Pz += cin[i].pz;
        double x = cin[i].px;

        double y = cin[i].py;

        double z = cin[i].pz;

        //  printf("x = %lf , y = %lf , z = %lf\n", x, y, z);

        energie += ((x * x) + (y * y) + (z * z)) / m;
        //  printf("Temp = %lf %lf %lf %lf  \n", cin[i].px, cin[i].py, cin[i].pz, vi);
    }

    double temperature = 0;

    double tmp = 1194.4577162;

    double energie_cin = 1;

    Px /= N_particules_total;

    Py /= N_particules_total;

    Pz /= N_particules_total;

    double energie0 = Ndl * CONSTANTE_R * T0;

    energie *= tmp;

    double RAPPORT = sqrt(energie0 / energie);

    // printf("%lf\n", RAPPORT);

    energie = 0;

    for (int i = 0; i < N_particules_total; i++)
    {
        double x = cin[i].px;

        double y = cin[i].py;

        double z = cin[i].pz;
        // energie_cin = a * energie;

        energie += ((x * x) + (y * y) + (z * z)) / m;
        temperature = (1 / (Ndl * CONSTANTE_R)) * energie_cin;
        energie_cin = tmp * energie;

        // 4
        cin[i + 1].px -= Px;
        cin[i + 1].py -= Py;
        cin[i + 1].pz -= Pz;

        // 3
        cin[i + 1].px *= RAPPORT;
        cin[i + 1].py *= RAPPORT;
        cin[i + 1].pz *= RAPPORT;

        // printf("Temp = %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf\n", cin[i].px, cin[i].py, cin[i].pz, temperature, energie_cin, Px, Py, Pz, energie, tmp);
    }

    /*cin[i + 1].px *= sqrt(temperature / T0);
    cin[i + 1].py *= sqrt(temperature / T0);
    cin[i + 1].pz *= sqrt(temperature / T0);*/

    printf("\033[036mENERGIE CINETIQUE\n\033[036m\tEnergie: \033[37m%lf\033[38;5;54m J\n\t\033[036mTempérature: \033[37m%lf\033[38;5;54m K\n", energie_cin, temperature);

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

    double tmp = 1194.4577162;

    for (int i = 0; i < N_particules_total; i++)
    {

        double x = cin[i].px;
        double y = cin[i].py;
        double z = cin[i].pz;

        energie += ((x * x) + (y * y) + (z * z)) / m;

        energie_cin = tmp * energie;
        temperature = (1 / (Ndl * CONSTANTE_R)) * energie_cin;

        // 6
        double a = 0.01 * ((T0 / temperature) - 1);

        cin[i + 1].px += a * cin[i + 1].px;
        cin[i + 1].py += a * cin[i + 1].py;
        cin[i + 1].pz += a * cin[i + 1].pz;
    }

    printf("\033[34mBERENDSEN\n\tEnergie cinétique: \033[37m%lf\033[38;5;54m J\n", energie_cin);
    printf("\t\033[34mTempérature:\033[37m %lf\033[38;5;54m K\n", temperature);

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
        i++;
        particle = realloc(particle, (i + 1) * (sizeof(struct Particle)));
    }

    N_particules_total = i;

    printf("Particules totales = %d\n", N_particules_total);

    fclose(in);

    return particle;
}

int main(int argc, char const *argv[])
{
    struct Particle *particle = malloc(sizeof(*particle));

    if ((particle = read_file(particle)) == -1)
        return printf("Erreur lecture des coordonnées\n"), -1;

    int iterations = 1000;

    if (argc > 1)
        iterations = atoi(argv[1]);

    // printf("%lf\n", particle[0].x);

    non_periodique(particle);
    periodique(particle);

    calcul_energie_cinetique(particle);

    calcul_energie_cinetique_berendsen(particle);

    simulation_np(iterations, particle);
    simulation_p(iterations, particle);

    return 0;
}
