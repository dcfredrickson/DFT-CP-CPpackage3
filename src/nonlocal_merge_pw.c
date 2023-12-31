/*     nonlocal_byatom, part of the Fredrickson Group Chemical Pressure Package  

                  Copyright (C) 2012, by Daniel C. Fredrickson

                    Last modified:  Mar. 28, 2012

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>

#define CUT3D_COMMAND "cut3d"
#define PI 3.14159265
#define TWO_PI 6.28318531
#define ONE_TWO_PI 0.159154943
#define PI_5_4 4.182513398
#define MAX_YAEHMOPFILES 10
#define LONGEST_FILENAME 100
#define DOSPOINTS_MAX 30000
#define NGX_MAX 155
#define NGY_MAX 155
#define NGZ_MAX 305
#define NATOMS_MAX 300

#define TBANDS_MAX 3000
#define BANDS_MAX 1000
#define NORBS_MAX 3000
#define MAX_SYM 1000
#define KPOINTS_MAX 1000
#define WAVES_MAX 301000
#define NTYPES_MAX 10

struct psp_file {
    double rloc[NTYPES_MAX];
    double rrs[NTYPES_MAX];
    double rrp[NTYPES_MAX];
    double rrd[NTYPES_MAX];
    double rrf[NTYPES_MAX];
    double cc1[NTYPES_MAX];
    double cc2[NTYPES_MAX];
    double cc3[NTYPES_MAX];
    double cc4[NTYPES_MAX];
    double h[3][3][4][NTYPES_MAX];
    double h11s[NTYPES_MAX];
    double h22s[NTYPES_MAX];
    double h33s[NTYPES_MAX];
    double h11p[NTYPES_MAX];
    double h22p[NTYPES_MAX];
    double h33p[NTYPES_MAX];
    double h11d[NTYPES_MAX];
    double h22d[NTYPES_MAX];
    double h33d[NTYPES_MAX];
    double h11f[NTYPES_MAX];
    double h22f[NTYPES_MAX];
    double h33f[NTYPES_MAX];
    double k11p[NTYPES_MAX];
    double k22p[NTYPES_MAX];
    double k33p[NTYPES_MAX];
    double k11d[NTYPES_MAX];
    double k22d[NTYPES_MAX];
    double k33d[NTYPES_MAX];
    double k11f[NTYPES_MAX];
    double k22f[NTYPES_MAX];
    double k33f[NTYPES_MAX];
    double k12p[NTYPES_MAX];
    double k23p[NTYPES_MAX];
    double k13p[NTYPES_MAX];
    double k12d[NTYPES_MAX];
    double k23d[NTYPES_MAX];
    double k13d[NTYPES_MAX];
    double k12f[NTYPES_MAX];
    double k23f[NTYPES_MAX];
    double k13f[NTYPES_MAX];
}
psp;

struct wfk_file {
    int nbands;
    int nkpts;
    int npw[2][KPOINTS_MAX];
    int nspinor[2][KPOINTS_MAX];
    int nband[2][KPOINTS_MAX];
    int kg[WAVES_MAX][3];
    double eigen[BANDS_MAX];
    double occ[BANDS_MAX];
    double cg[WAVES_MAX][BANDS_MAX][2];
}
wfk1;

double nonlocalE_byatom[3][NATOMS_MAX][6];

struct XSFfile {
    char systemname[100];
    double cellvolume;
    double scale_xyz;
    double cella_x, cella_y, cella_z;
    double cellb_x, cellb_y, cellb_z;
    double cellc_x, cellc_y, cellc_z;
    int atomicno[NATOMS_MAX];
    double Xcart[NATOMS_MAX];
    double Ycart[NATOMS_MAX];
    double Zcart[NATOMS_MAX];
    int NIONS;
    int NGX, NGY, NGZ;
    double ** * grid;
    double VoxelV;
    double NELECT;
    double epsatm[NATOMS_MAX];
    int epsatm_map[NATOMS_MAX];
    int epsatm_map_type[NATOMS_MAX];
    double nelectrons;
}
den0, nlden1, nlden2;

double cgt[2][WAVES_MAX][2];
double ga[WAVES_MAX], gb[WAVES_MAX], gc[WAVES_MAX];
double gx[WAVES_MAX], gy[WAVES_MAX], gz[WAVES_MAX];
double costheta[WAVES_MAX];
double theta[WAVES_MAX];
double phi[WAVES_MAX];
double g[WAVES_MAX];
double Plm[WAVES_MAX][4][4];
/* p1=projp(l,i,g[pw1],pspin,atom_type,cellvolume); */
double p[WAVES_MAX][4][3][NTYPES_MAX];
int typat[NATOMS_MAX];

FILE * inputfile;
FILE * outputfile;
char filename[LONGEST_FILENAME];
int nkpt;

int AllocDbl(struct XSFfile * gridin) {
    /* called by: */
    /* calls: none */
    int gridx;
    int gridy;
    int gridz;
    int jx, jy, jz;
    gridx = gridin -> NGX;
    gridy = gridin -> NGY;
    gridz = gridin -> NGZ;
    gridin -> grid = (double ** * ) malloc(gridx * sizeof(double ** ));
    for (jx = 0; jx < gridx; jx++) {
        gridin -> grid[jx] = (double ** ) malloc(gridy * sizeof(double * ));
        for (jy = 0; jy < gridy; jy++) {
            gridin -> grid[jx][jy] = (double * ) malloc(gridz * sizeof(double));
        }
    }
    return 0;
}

void finish_line(FILE * f3) {
    int cont = 0;
    char check;
    while (cont == 0) {
        check = getc(f3);
        if ((check == 10) || (check == EOF)) cont = 1;
    }
}

void outputXSF(struct XSFfile * XSFIN, struct XSFfile * XSFOUT, FILE * f2) {
    int jx, jy, jz;
    int j;
    int line_counter;
    fprintf(f2, " DIM-GROUP\n");
    fprintf(f2, " 3  1\n");
    fprintf(f2, " PRIMVEC\n");
    fprintf(f2, "%17.10lf%19.10lf%19.10lf\n", XSFIN -> cella_x, XSFIN -> cella_y, XSFIN -> cella_z);
    fprintf(f2, "%17.10lf%19.10lf%19.10lf\n", XSFIN -> cellb_x, XSFIN -> cellb_y, XSFIN -> cellb_z);
    fprintf(f2, "%17.10lf%19.10lf%19.10lf\n", XSFIN -> cellc_x, XSFIN -> cellc_y, XSFIN -> cellc_z);
    fprintf(f2, " PRIMCOORD\n");
    fprintf(f2, "%12d%3d\n", XSFIN -> NIONS, 1);
    for (j = 0; j < XSFIN -> NIONS; j++) {
        fprintf(f2, "%9d%20.10lf%20.10lf%20.10lf\n", XSFIN -> atomicno[j], XSFIN -> Xcart[j], XSFIN -> Ycart[j], XSFIN -> Zcart[j]);
    }
    fprintf(f2, " ATOMS\n");
    for (j = 0; j < XSFIN -> NIONS; j++) {
        fprintf(f2, "%9d%20.10lf%20.10lf%20.10lf\n", XSFIN -> atomicno[j], XSFIN -> Xcart[j], XSFIN -> Ycart[j], XSFIN -> Zcart[j]);
    }
    fprintf(f2, " BEGIN_BLOCK_DATAGRID3D\n");
    fprintf(f2, " datagrids\n");
    fprintf(f2, " DATAGRID_3D_DENSITY\n");
    fprintf(f2, "%12d%12d%12d\n", XSFIN -> NGX, XSFIN -> NGY, XSFIN -> NGZ);
    fprintf(f2, " 0.0 0.0 0.0\n");
    fprintf(f2, "%17.10lf%19.10lf%19.10lf\n", XSFIN -> cella_x, XSFIN -> cella_y, XSFIN -> cella_z);
    fprintf(f2, "%17.10lf%19.10lf%19.10lf\n", XSFIN -> cellb_x, XSFIN -> cellb_y, XSFIN -> cellb_z);
    fprintf(f2, "%17.10lf%19.10lf%19.10lf\n", XSFIN -> cellc_x, XSFIN -> cellc_y, XSFIN -> cellc_z);
    line_counter = 0;
    for (jz = 0; jz < XSFIN -> NGZ; jz++) {
        for (jy = 0; jy < XSFIN -> NGY; jy++) {
            for (jx = 0; jx < XSFIN -> NGX; jx++) {
                line_counter++;
                fprintf(f2, "%20.10lf", XSFOUT -> grid[jx][jy][jz]);
                if (line_counter == 6) {
                    fprintf(f2, "\n");
                    line_counter = 0;
                }
            }
        }
    }
    fprintf(f2, " END_DATAGRID_3D\n");
    fprintf(f2, " END_BLOCK_DATAGRID3D\n");

}

void xyz2sph(double X, double Y, double Z, double * r, double * theta, double * phi) {
    /* Z = r cos(theta), X = r sin(theta) cos(phi) Y = r sin(theta)sin(phi) */
    * r = pow((pow(X, 2.0) + pow(Y, 2.0) + pow(Z, 2.0)), 0.5);
    if ( * r > 0.0) {
        * theta = acos(Z / * r);
        if (fabs(sin( * theta)) > 0.0) {
            * phi = acos(X / ( * r * sin( * theta)));
            if (X / ( * r * sin( * theta)) > 1.0) {
                * phi = acos(1.0);
            }
            if (X / ( * r * sin( * theta)) < -1.0) {
                * phi = acos(-1.0);
            }
            if (Y < 0.0) * phi = - * phi;
        } else * phi = 0;
    } else {
        * phi = 0.0;
        * theta = 0.0;
    }
}

double projp(int l, int i, double g1, struct psp_file * pspin, int typat, double cellvolume) {
    //                                                  p1=projp(l,i,g1,pspin,atom_type,cellvolume);
    double p = 0.0;
    typat--;
    if ((l == 0) && (i == 0)) {
        p = 4.0 * pow(2.0 * pow(pspin -> rrs[typat], 3.0), 0.5) * PI_5_4 / (pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrs[typat], 2.0)));
        //	       printf("p = %lf   rrs = %lf g1 = %lf atom_type = %d,  cell_volume = %lf\n",p, pspin->rrs[0],g1,typat,cellvolume);
        //               exit(0);
    }
    if ((l == 0) && (i == 1)) {
        p = 8.0 * pow(2.0 * pow(pspin -> rrs[typat], 3.0) / 15.0, 0.5) * PI_5_4 * (3.0 - pow(g1 * pspin -> rrs[typat], 2.0)) / (pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrs[typat], 2.0)));
    }
    if ((l == 0) && (i == 2)) {
        p = 16.0 * pow(2.0 * pow(pspin -> rrs[typat], 3.0) / 105.0, 0.5) * PI_5_4 * (15.0 - 10.0 * pow(g1 * pspin -> rrs[typat], 2.0) + pow(g1 * pspin -> rrs[typat], 4.0)) / (3.0 * pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrs[typat], 2.0)));
    }
    if ((l == 1) && (i == 0)) {
        p = 8.0 * pow(pow(pspin -> rrp[typat], 5.0) / 3.0, 0.5) * PI_5_4 * (g1) / (pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrp[typat], 2.0)));
    }
    if ((l == 1) && (i == 1)) {
        p = 16.0 * pow(pow(pspin -> rrp[typat], 5.0) / 105.0, 0.5) * PI_5_4 * (g1) * (5.0 - pow(g1 * pspin -> rrp[typat], 2.0)) / (pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrp[typat], 2.0)));
    }
    if ((l == 1) && (i == 2)) {
        p = 32.0 * pow(pow(pspin -> rrp[typat], 5.0) / 1155.0, 0.5) * PI_5_4 * (g1) * (35.0 - 14.0 * pow(g1 * pspin -> rrp[typat], 2.0) + pow(g1 * pspin -> rrp[typat], 4.0)) / (3.0 * pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrp[typat], 2.0)));
    }
    if ((l == 2) && (i == 0)) {
        p = 8.0 * pow(2 * pow(pspin -> rrd[typat], 7.0) / 15.0, 0.5) * PI_5_4 * (g1 * g1) / (pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrd[typat], 2.0)));
    }
    if ((l == 2) && (i == 1)) {
        p = 16.0 * pow(2 * pow(pspin -> rrd[typat], 7.0) / 105.0, 0.5) * PI_5_4 * (g1 * g1) * (7.0 - pow(g1 * pspin -> rrd[typat], 2.0)) / (3.0 * pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrd[typat], 2.0)));
    }
    if ((l == 3) && (i == 0)) {
        p = 16.0 * pow(pow(pspin -> rrf[typat], 9.0) / 105.0, 0.5) * PI_5_4 * (g1 * g1 * g1) / (pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrf[typat], 2.0)));
    }
    //        if(p<0.0) printf("YIKES!");
    return (p);
}

double projp_r(int l, int i, double r, struct psp_file * pspin, int typat) {
    double p = 0.0;
    double r_l;
    if (l == 0) r_l = pspin -> rrs[typat];
    if (l == 1) r_l = pspin -> rrp[typat];
    if (l == 2) r_l = pspin -> rrd[typat];
    if (l == 3) r_l = pspin -> rrf[typat];
    if (r_l > 0.0) {
        p = sqrt(2.0) * pow(r, 1.0 * l + 2.0 * (1.0 * i + 1.0) - 2.0) * exp(-pow(r, 2.0) / (2.0 * pow(r_l, 2.0))) / (pow(r_l, 1.0 * l + 2.0 * (1.0 * i + 1.0) - 0.5) * sqrt(gsl_sf_gamma(l + 2.0 * (1.0 * i + 1.0) - 0.5)));
    }
    return (p);
}

void read_calc_den(char filename[200], struct XSFfile * XSFIN) {
    FILE * f2;
    int k;
    char codvsn[110];
    char title[132];
    char test;
    int headform;
    int fform;
    int bandtot;
    int date;
    int intxc;
    int ixc;
    int natom;
    int atomno, pw1, pw2;
    int ngfftx;
    int ngffty;
    int ngfftz;
    int npsp;
    int j = 0;
    int i, m, l;
    int jx, jy, jz;
    int stop = 0;
    int nspden;
    int nsppol;
    int nsym;
    int ntypat;
    int occopt;
    int pertcase;
    int usepaw;
    double ecut;
    double ecutdg;
    double ecutsm;
    double ecut_eff;
    double qptnx;
    double qptny;
    double qptnz;
    double rprimd_ax;
    double rprimd_ay;
    double rprimd_az;
    double rprimd_bx;
    double rprimd_by;
    double rprimd_bz;
    double rprimd_cx;
    double rprimd_cy;
    double rprimd_cz;
    double ax_star;
    double ay_star;
    double az_star;
    double bx_star;
    double by_star;
    double bz_star;
    double cx_star;
    double cy_star;
    double cz_star;
    double ox, oy, oz;
    double theta1, phi1, g1;
    double theta2, phi2, g2;
    double stmbias;
    double tphysel;
    double tsmear;
    double znuclpsp;
    double zionpsp;
    int pspso;
    int pspdat;
    int pspcod;
    int pspxc;
    int type;
    int lmn_size;
    int atom_type;
    int usewvl;
    int istwfk[KPOINTS_MAX];
    int istwfkv;
    int nband[BANDS_MAX];
    int nbandv;
    int npwarr[KPOINTS_MAX];
    int npwarrv;
    int so_psp[NATOMS_MAX];
    int symafm[MAX_SYM];
    int symrel[3][3][MAX_SYM];
    double kpt[3][KPOINTS_MAX];
    double occ;
    double tnons[3][MAX_SYM];
    double znucltypat[NATOMS_MAX];
    double wtk[KPOINTS_MAX];
    double wtkv;
    double residm;
    double x, y, z;
    double etotal, fermie;
    double Enonlocal_temp;
    double ga1, gb1, gc1;
    double gx1, gy1, gz1;
    double ga2, gb2, gc2;
    double gx2, gy2, gz2;
    double xred[3][NATOMS_MAX];
    double kptv;
    double php;
    int npw;
    int pw;
    int nspinor;
    int nband_temp;
    int kptno;
    int kx, ky, kz;
    int band;
    double normalization;
    double eigen, occ_temp, cg;
    double cellvolume, p1, p2;
    double atom_phase0;

    f2 = fopen(filename, "rb+");
    if (f2 == NULL) {
        printf("%s not found.\n", filename);
        exit(0);
    }
    /*   READ HEADER OF WFK FILE    */
    j = 0;
    fread( & j, sizeof(int), 1, f2);
    j = fread(codvsn, sizeof(char), 6, f2);
    fread( & headform, sizeof(int), 1, f2);
    fread( & fform, sizeof(int), 1, f2);
    //    printf("%s %d %d\n",codvsn,headform,fform);
    fread( & j, sizeof(int), 1, f2);
    fread( & j, sizeof(int), 1, f2);
    fread( & bandtot, sizeof(int), 1, f2);
    /*    int date;*/
    fread( & date, sizeof(int), 1, f2);
    //    printf("%d %d\n",bandtot,date);
    /*    int intxc;  */
    fread( & intxc, sizeof(int), 1, f2);
    /*    int ixc;*/
    fread( & ixc, sizeof(int), 1, f2);
    /*    int natom;*/
    fread( & natom, sizeof(int), 1, f2);
    XSFIN -> NIONS = natom;
    /*    int ngfftx;*/
    fread( & ngfftx, sizeof(int), 1, f2);
    /*    int ngffty;*/
    fread( & ngffty, sizeof(int), 1, f2);
    /*    int ngfftz;*/
    fread( & ngfftz, sizeof(int), 1, f2);
    //    printf("ngfft = %d x %d x %d\n",ngfftx,ngffty,ngfftz);
    XSFIN -> NGX = ngfftx + 1;
    XSFIN -> NGY = ngffty + 1;
    XSFIN -> NGZ = ngfftz + 1;
    AllocDbl(XSFIN);
    /*    int nkpt;*/
    fread( & nkpt, sizeof(int), 1, f2);
    /*    int nspden;*/
    fread( & nspden, sizeof(int), 1, f2);
    /*    int nspinor;*/
    fread( & nspinor, sizeof(int), 1, f2);
    /*    int nsppol;*/
    fread( & nsppol, sizeof(int), 1, f2);
    /*    int nsym;*/
    fread( & nsym, sizeof(int), 1, f2);
    /*    int npsp;*/
    fread( & npsp, sizeof(int), 1, f2);
    /*    int ntypat;*/
    fread( & ntypat, sizeof(int), 1, f2);
    /*    int occopt;*/
    fread( & occopt, sizeof(int), 1, f2);
    /*    int pertcase;*/
    fread( & pertcase, sizeof(int), 1, f2);
    /*    int usepaw;*/
    fread( & usepaw, sizeof(int), 1, f2);
    /*    double ecut;*/
    fread( & ecut, sizeof(double), 1, f2);
    /*     double ecutdg;*/
    fread( & ecutdg, sizeof(double), 1, f2);
    /*    double ecutsm;*/
    fread( & ecutsm, sizeof(double), 1, f2);
    /*    double ecut_eff;*/
    fread( & ecut_eff, sizeof(double), 1, f2);
    /*    double qptnx;*/
    fread( & qptnx, sizeof(double), 1, f2);
    /*    double qptny;*/
    fread( & qptny, sizeof(double), 1, f2);
    /*    double qptnz;*/
    fread( & qptnz, sizeof(double), 1, f2);
    /*    double rprimd_ax;*/
    fread( & rprimd_ax, sizeof(double), 1, f2);
    XSFIN -> cella_x = rprimd_ax * 0.52917720859;
    /*    double rprimd_ay;*/
    fread( & rprimd_ay, sizeof(double), 1, f2);
    XSFIN -> cella_y = rprimd_ay * 0.52917720859;
    /*    double rprimd_az;*/
    fread( & rprimd_az, sizeof(double), 1, f2);
    XSFIN -> cella_z = rprimd_az * 0.52917720859;
    /*    double rprimd_bx;*/
    fread( & rprimd_bx, sizeof(double), 1, f2);
    XSFIN -> cellb_x = rprimd_bx * 0.52917720859;
    /*    double rprimd_by;*/
    fread( & rprimd_by, sizeof(double), 1, f2);
    XSFIN -> cellb_y = rprimd_by * 0.52917720859;
    /*    double rprimd_bz;*/
    fread( & rprimd_bz, sizeof(double), 1, f2);
    XSFIN -> cellb_z = rprimd_bz * 0.52917720859;
    /*    double rprimd_cx;*/
    fread( & rprimd_cx, sizeof(double), 1, f2);
    XSFIN -> cellc_x = rprimd_cx * 0.52917720859;
    /*    double rprimd_cy;*/
    fread( & rprimd_cy, sizeof(double), 1, f2);
    XSFIN -> cellc_y = rprimd_cy * 0.52917720859;
    /*    double rprimd_cz;*/
    fread( & rprimd_cz, sizeof(double), 1, f2);
    XSFIN -> cellc_z = rprimd_cz * 0.52917720859;
    /*    double stmbias;*/
    cellvolume = (rprimd_ax * (rprimd_by * rprimd_cz - rprimd_bz * rprimd_cy) - rprimd_ay * (rprimd_bx * rprimd_cz - rprimd_bz * rprimd_cx) + rprimd_az * (rprimd_bx * rprimd_cy - rprimd_by * rprimd_cx));
    //    printf(" cella = %lf %lf %lf \n",rprimd_ax,rprimd_ay,rprimd_az);
    //    printf(" cellb = %lf %lf %lf \n",rprimd_bx,rprimd_by,rprimd_bz);
    //    printf(" cellc = %lf %lf %lf \n\n",rprimd_cx,rprimd_cy,rprimd_cz);
    //    printf(" cell volume = %lf\n\n",cellvolume);
    XSFIN -> cellvolume = cellvolume;
    XSFIN -> VoxelV = XSFIN -> cellvolume / ((XSFIN -> NGX - 1) * (XSFIN -> NGY - 1) * (XSFIN -> NGZ - 1));
    ax_star = 2 * PI * (rprimd_by * rprimd_cz - rprimd_cy * rprimd_bz) / cellvolume;
    ay_star = -2 * PI * (rprimd_bx * rprimd_cz - rprimd_cx * rprimd_bz) / cellvolume;
    az_star = 2 * PI * (rprimd_bx * rprimd_cy - rprimd_cx * rprimd_by) / cellvolume;
    bx_star = 2 * PI * (rprimd_cy * rprimd_az - rprimd_ay * rprimd_cz) / cellvolume;
    by_star = -2 * PI * (rprimd_cx * rprimd_az - rprimd_ax * rprimd_cz) / cellvolume;
    bz_star = 2 * PI * (rprimd_cx * rprimd_ay - rprimd_ax * rprimd_cy) / cellvolume;
    cx_star = 2 * PI * (rprimd_ay * rprimd_bz - rprimd_by * rprimd_az) / cellvolume;
    cy_star = -2 * PI * (rprimd_ax * rprimd_bz - rprimd_bx * rprimd_az) / cellvolume;
    cz_star = 2 * PI * (rprimd_ax * rprimd_by - rprimd_bx * rprimd_ay) / cellvolume;
    //    printf(" cella* = %lf %lf %lf \n",ax_star,ay_star,az_star);
    //    printf(" cellb* = %lf %lf %lf \n",bx_star,by_star,bz_star);
    //    printf(" cellc* = %lf %lf %lf \n\n",cx_star,cy_star,cz_star);
    //    printf(" a.a*=%lf, a.b*=%lf, a.c*=%lf\n",rprimd_ax*ax_star+rprimd_ay*ay_star+rprimd_az*az_star,rprimd_ax*bx_star+rprimd_ay*by_star+rprimd_az*bz_star,rprimd_ax*cx_star+rprimd_ay*cy_star+rprimd_az*cz_star);
    //    printf(" b.a*=%lf, b.b*=%lf, b.c*=%lf\n",rprimd_bx*ax_star+rprimd_by*ay_star+rprimd_bz*az_star,rprimd_bx*bx_star+rprimd_by*by_star+rprimd_bz*bz_star,rprimd_bx*cx_star+rprimd_by*cy_star+rprimd_bz*cz_star);
    //    printf(" c.a*=%lf, c.b*=%lf, c.c*=%lf\n\n",rprimd_cx*ax_star+rprimd_cy*ay_star+rprimd_cz*az_star,rprimd_cx*bx_star+rprimd_cy*by_star+rprimd_cz*bz_star,rprimd_cx*cx_star+rprimd_cy*cy_star+rprimd_cz*cz_star);
    fread( & stmbias, sizeof(double), 1, f2);
    /*    double tphysel;*/
    fread( & tphysel, sizeof(double), 1, f2);
    /*    double tsmear;*/
    fread( & tsmear, sizeof(double), 1, f2);
    /*    int usewvl;*/
    fread( & usewvl, sizeof(int), 1, f2);
    //    printf("natoms = %d   ecut = %lf  tsmear = %lf  occopt = %d \n",natom,ecut,tsmear, occopt);   
    fread( & j, sizeof(int), 1, f2);
    fread( & j, sizeof(int), 1, f2);
    /*    int istwfk [KPOINTS_MAX]; */
    //    printf("nkpt = %d    nsppol = %d   npsp = %d  ntypat = %d \n",nkpt,nsppol,npsp,ntypat);
    for (j = 0; j < nkpt; j++) {
        fread( & istwfkv, sizeof(int), 1, f2);
    }
    /*    int nband [BANDS_MAX]; */
    for (j = 0; j < (nkpt * nsppol); j++) {
        fread( & nbandv, sizeof(int), 1, f2);
    }
    /*    int npwarr [KPOINTS_MAX]; */
    for (j = 0; j < (nkpt); j++) {
        fread( & npwarrv, sizeof(int), 1, f2);
    }
    /*    int so_psp [NATOMS_MAX]; */
    for (j = 0; j < (npsp); j++) {
        fread( & so_psp[j], sizeof(int), 1, f2);
    }
    /*    int symafm [MAX_SYM]; */
    for (j = 0; j < (nsym); j++) {
        fread( & symafm[j], sizeof(int), 1, f2);
    }
    /*    int symrel [3][3][MAX_SYM]; */
    for (j = 0; j < (nsym); j++) {
        fread( & symrel[0][0][j], sizeof(int), 1, f2);
        fread( & symrel[1][0][j], sizeof(int), 1, f2);
        fread( & symrel[2][0][j], sizeof(int), 1, f2);
        fread( & symrel[0][1][j], sizeof(int), 1, f2);
        fread( & symrel[1][1][j], sizeof(int), 1, f2);
        fread( & symrel[2][1][j], sizeof(int), 1, f2);
        fread( & symrel[0][2][j], sizeof(int), 1, f2);
        fread( & symrel[1][2][j], sizeof(int), 1, f2);
        fread( & symrel[2][2][j], sizeof(int), 1, f2);
    }
    /*    int typat [NATOMS_MAX]; */
    for (j = 0; j < (natom); j++) {
        fread( & typat[j], sizeof(int), 1, f2);
        //	 printf("typeat = %d\n",typat[j]);         
    }
    /*    double kpt [3][KPOINTS_MAX]; */
    for (j = 0; j < (nkpt); j++) {
        fread( & kptv, sizeof(double), 1, f2);
        fread( & kptv, sizeof(double), 1, f2);
        fread( & kptv, sizeof(double), 1, f2);
        //         printf("kpoint %d:  %lf %lf %lf \n",j+1,kpt[0][j],kpt[1][j],kpt[2][j]);
    }
    /*    double occ(TBANDS_MAX); */
    for (j = 0; j < (bandtot); j++) {
        fread( & occ, sizeof(double), 1, f2);
    }
    /*    double tnons [3][MAX_SYM]; */
    for (j = 0; j < (nsym); j++) {
        fread( & tnons[0][j], sizeof(double), 1, f2);
        fread( & tnons[1][j], sizeof(double), 1, f2);
        fread( & tnons[2][j], sizeof(double), 1, f2);
    }

    /*    double znucltypat[NATOMS_MAX];  */
    for (j = 0; j < (ntypat); j++) {
        fread( & znucltypat[j], sizeof(double), 1, f2);
    }

    for (j = 0; j < natom; j++) {
        type = typat[j] - 1;
        //          printf("atom %d type=%d\n",j, type);
        XSFIN -> atomicno[j] = (int) znucltypat[type];
    }
    /*    double wtk[KPOINTS_MAX]; */
    for (j = 0; j < (nkpt); j++) {
        fread( & wtkv, sizeof(double), 1, f2);
    }
    fread( & j, sizeof(int), 1, f2);
    for (k = 0; k < (npsp); k++) {
        fread( & j, sizeof(int), 1, f2);
        fread(title, sizeof(char), 132, f2);
        //          printf("     %s\n",title);
        fread( & znuclpsp, sizeof(double), 1, f2);
        //          printf("     %lf ",znuclpsp);
        fread( & zionpsp, sizeof(double), 1, f2);
        //          printf("%lf ",zionpsp);
        fread( & pspso, sizeof(int), 1, f2);
        fread( & pspdat, sizeof(int), 1, f2);
        fread( & pspcod, sizeof(int), 1, f2);
        fread( & pspxc, sizeof(int), 1, f2);
        fread( & lmn_size, sizeof(int), 1, f2);
        //          printf("%d %d %d %d %d \n",pspso,pspdat,pspcod,pspxc,lmn_size);
        fread( & j, sizeof(int), 1, f2);
    }
    if (usepaw == 0) {
        fread( & j, sizeof(int), 1, f2);
        fread( & residm, sizeof(double), 1, f2);
        //        printf("     residm = %lf\n",residm);
        //          printf("Enter origin x y z: ");
        for (k = 0; k < natom; k++) {
            fread( & x, sizeof(double), 1, f2);
            fread( & y, sizeof(double), 1, f2);
            fread( & z, sizeof(double), 1, f2);
            //              printf("     Atom %d:  %lf %lf %lf \n",k+1,x,y,z);
            XSFIN -> Xcart[k] = (x) * XSFIN -> cella_x + (y) * XSFIN -> cellb_x + (z) * XSFIN -> cellc_x;
            XSFIN -> Ycart[k] = (x) * XSFIN -> cella_y + (y) * XSFIN -> cellb_y + (z) * XSFIN -> cellc_y;
            XSFIN -> Zcart[k] = (x) * XSFIN -> cella_z + (y) * XSFIN -> cellb_z + (z) * XSFIN -> cellc_z;
        }
        fread( & etotal, sizeof(double), 1, f2);
        fread( & fermie, sizeof(double), 1, f2);
        //          printf("     Etotal = %lf    FermiE = %lf\n\n",etotal,fermie);
        fread( & j, sizeof(int), 1, f2);
    } else {
        printf("Yikes!  usepaw!=0.  We haven't written code for this case yet. \n");
    }
    /*   HEADER FINISHED -  READ DENSITY DATA   */
    fread( & j, sizeof(int), 1, f2);
    for (jz = 0; jz < XSFIN -> NGZ; jz++) {
        for (jy = 0; jy < XSFIN -> NGY; jy++) {
            for (jx = 0; jx < XSFIN -> NGX; jx++) {

                if ((jx < XSFIN -> NGX - 1) && (jy < XSFIN -> NGY - 1) && (jz < XSFIN -> NGZ - 1)) {
                    fread( & eigen, sizeof(double), 1, f2);
                    XSFIN -> grid[jx][jy][jz] = eigen;
                }
            }
        }
    }
    fread( & j, sizeof(int), 1, f2);
    fclose(f2);
    for (jz = 0; jz < XSFIN -> NGZ; jz++) {
        for (jy = 0; jy < XSFIN -> NGY; jy++) {
            for (jx = 0; jx < XSFIN -> NGX; jx++) {
                if (jx == XSFIN -> NGX - 1) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[0][jy][jz];
                if (jy == XSFIN -> NGY - 1) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[jx][0][jz];
                if (jz == XSFIN -> NGZ - 1) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[jx][jy][0];
                if ((jx == XSFIN -> NGX - 1) && (jy == XSFIN -> NGY - 1)) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[0][0][jz];
                if ((jx == XSFIN -> NGX - 1) && (jz == XSFIN -> NGZ - 1)) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[0][jy][0];
                if ((jy == XSFIN -> NGY - 1) && (jz == XSFIN -> NGZ - 1)) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[jx][0][0];
                if ((jx == XSFIN -> NGX - 1) && (jy == XSFIN -> NGY - 1) && (jz == XSFIN -> NGZ - 1)) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[0][0][0];
            }
        }
    }
}

void NL2XSF(struct XSFfile * XSFIN, struct XSFfile * XSFOUT1) {
    int jx, jy, jz, jx1, jy1, jz1, ha, hb, hc;
    int j, i1, i2;
    int l;
    int line_counter;
    int nvoxels_atm;
    int check;
    int stop;
    int stop2;
    int step;
    double r_max;
    double elementname[4];
    char outfilename[1000];
    char str[1000];
    int NGX, NGY, NGZ, k;
    double cell_volume, dist, dist2;
    double rcut, rcut2;
    double xf, yf, zf;
    double voxel_centerx;
    double voxel_centery;
    double voxel_centerz;
    double deltar;
    double Wnl_atom[NATOMS_MAX][4];
    double nelect_core;
    int n_epsatm_values = 0;
    int step1, step2;
    double averageE1_vs_r[10001];
    double averageE2_vs_r[10001];
    double sigmaE1_vs_r[10001];
    double sigmaE2_vs_r[10001];
    double NL1, NL2, NELECT;
    int averageP_vs_r_npoints[10001];
    int kptno;
    int do_map = 0;
    int atomic_no;
    double Z;
    double tolerance = 0.5 * pow(XSFIN -> VoxelV * pow(0.52917720859, 3.0), 0.33333333333333); //0.005-->0.2-->0.4-->0.5
    double Z_ion_temp;
    double E_core_up, E_core_down;
    int Zuse[100];
    double r_by_Z[100];
    FILE * f2;
    FILE * f4;
    char profile_file[200];
    double cella, cellb, cellc;
    NGX = XSFIN -> NGX;
    NGY = XSFIN -> NGY;
    NGZ = XSFIN -> NGZ;
    r_max = 4.0;
    cella = pow(XSFIN -> cella_x * XSFIN -> cella_x + XSFIN -> cella_y * XSFIN -> cella_y + XSFIN -> cella_z * XSFIN -> cella_z, 0.5) / (NGX - 1);
    cellb = pow(XSFIN -> cellb_x * XSFIN -> cellb_x + XSFIN -> cellb_y * XSFIN -> cellb_y + XSFIN -> cellb_z * XSFIN -> cellb_z, 0.5) / (NGY - 1);
    cellc = pow(XSFIN -> cellc_x * XSFIN -> cellc_x + XSFIN -> cellc_y * XSFIN -> cellc_y + XSFIN -> cellc_z * XSFIN -> cellc_z, 0.5) / (NGZ - 1);
    stop = 0;
    XSFOUT1 -> NGX = NGX;
    XSFOUT1 -> NGY = NGY;
    XSFOUT1 -> NGZ = NGZ;
    AllocDbl(XSFOUT1);
    for (jz = 0; jz < NGZ; jz++) {
        for (jy = 0; jy < NGY; jy++) {
            for (jx = 0; jx < NGX; jx++) {
                XSFOUT1 -> grid[jx][jy][jz] = 0.0;
            }
        }
    }
    XSFOUT1 -> VoxelV = XSFIN -> VoxelV;
    for (j = 0; j < XSFIN -> NIONS; j++) {
        Wnl_atom[j][0] = 0.0;
        Wnl_atom[j][1] = 0.0;
        Wnl_atom[j][2] = 0.0;
        Wnl_atom[j][3] = 0.0;
        for (jz = -NGZ / 2; jz < XSFIN -> NGZ - 1 + NGZ / 2; jz++) {
            zf = (jz * 1.000) / ((NGZ - 1) * 1.000);
            // printf("Here: %lf\n",zf); 
            for (jy = -NGY / 2; jy < XSFIN -> NGY - 1 + NGY / 2; jy++) {
                yf = (jy * 1.000) / ((NGY - 1) * 1.000);
                for (jx = -NGX / 2; jx < XSFIN -> NGX - 1 + NGX / 2; jx++) {
                    xf = (jx * 1.000) / ((NGX - 1) * 1.000);
                    voxel_centerx = (xf) * XSFIN -> cella_x + (yf) * XSFIN -> cellb_x + (zf) * XSFIN -> cellc_x;
                    voxel_centery = (xf) * XSFIN -> cella_y + (yf) * XSFIN -> cellb_y + (zf) * XSFIN -> cellc_y;
                    voxel_centerz = (xf) * XSFIN -> cella_z + (yf) * XSFIN -> cellb_z + (zf) * XSFIN -> cellc_z;
                    dist2 = (voxel_centerx - XSFIN -> Xcart[j]) * (voxel_centerx - XSFIN -> Xcart[j]) + (voxel_centery - XSFIN -> Ycart[j]) * (voxel_centery - XSFIN -> Ycart[j]) + (voxel_centerz - XSFIN -> Zcart[j]) * (voxel_centerz - XSFIN -> Zcart[j]);
                    if (dist2 < 1.05 * r_max * r_max) {
                        jx1 = (jx + NGX - 1) % (NGX - 1);
                        jy1 = (jy + NGY - 1) % (NGY - 1);
                        jz1 = (jz + NGZ - 1) % (NGZ - 1);
                        dist = pow(dist2, 0.5) / 0.52917720859;
                        for (i1 = 0; i1 < 3; i1++) {
                            for (i2 = 0; i2 < 3; i2++) {
                                for (l = 0; l < 4; l++) {
                                    Wnl_atom[j][l] += XSFIN -> grid[jx1][jy1][jz1] * psp.h[i1][i2][l][typat[j] - 1] * projp_r(l, i1, dist, & psp, typat[j] - 1) * projp_r(l, i2, dist, & psp, typat[j] - 1) * XSFIN -> VoxelV;
                                }
                            }
                        }
                    }
                } /*  end loop on jx */
            } /*  end loop on jy */
        } /*  end loop on jz */
    } /*  end loop on j */
    for (j = 0; j < XSFIN -> NIONS; j++) {
        printf("Weights:  (s) %lf (p) %lf (d) %lf (f) %lf\n", Wnl_atom[j][0], Wnl_atom[j][1], Wnl_atom[j][2], Wnl_atom[j][3]);
        for (jz = -NGZ / 2; jz < XSFIN -> NGZ - 1 + NGZ / 2; jz++) {
            zf = (jz * 1.000) / ((NGZ - 1) * 1.000);
            for (jy = -NGY / 2; jy < XSFIN -> NGY - 1 + NGY / 2; jy++) {
                yf = (jy * 1.000) / ((NGY - 1) * 1.000);
                for (jx = -NGX / 2; jx < XSFIN -> NGX - 1 + NGX / 2; jx++) {
                    xf = (jx * 1.000) / ((NGX - 1) * 1.000);
                    voxel_centerx = (xf) * XSFIN -> cella_x + (yf) * XSFIN -> cellb_x + (zf) * XSFIN -> cellc_x;
                    voxel_centery = (xf) * XSFIN -> cella_y + (yf) * XSFIN -> cellb_y + (zf) * XSFIN -> cellc_y;
                    voxel_centerz = (xf) * XSFIN -> cella_z + (yf) * XSFIN -> cellb_z + (zf) * XSFIN -> cellc_z;
                    dist2 = (voxel_centerx - XSFIN -> Xcart[j]) * (voxel_centerx - XSFIN -> Xcart[j]) + (voxel_centery - XSFIN -> Ycart[j]) * (voxel_centery - XSFIN -> Ycart[j]) + (voxel_centerz - XSFIN -> Zcart[j]) * (voxel_centerz - XSFIN -> Zcart[j]);
                    if (dist2 < 1.05 * r_max * r_max) {
                        jx1 = (jx + NGX - 1) % (NGX - 1);
                        jy1 = (jy + NGY - 1) % (NGY - 1);
                        jz1 = (jz + NGZ - 1) % (NGZ - 1);
                        dist = pow(dist2, 0.5) / 0.52917720859;
                        if (nonlocalE_byatom[0][j][0] != 0.0) {
                            for (i1 = 0; i1 < 3; i1++) {
                                for (i2 = 0; i2 < 3; i2++) {
                                    for (l = 0; l < 4; l++) {
                                        if (nonlocalE_byatom[0][j][l] != 0.0) {
                                            XSFOUT1 -> grid[jx1][jy1][jz1] += nonlocalE_byatom[0][j][l] * XSFIN -> grid[jx1][jy1][jz1] * psp.h[i1][i2][l][typat[j] - 1] * projp_r(l, i1, dist, & psp, typat[j] - 1) * projp_r(l, i2, dist, & psp, typat[j] - 1) / Wnl_atom[j][l];
                                        }
                                    }
                                }
                            }
                        }
                    }
                } /*  end loop on jx */
            } /*  end loop on jy */
        } /*  end loop on jz */
    } /*  end loop on j */
    NELECT = 0.0;
    NL1 = 0.0;
    for (jz = 0; jz < NGZ - 1; jz++) {
        for (jy = 0; jy < NGY - 1; jy++) {
            for (jx = 0; jx < NGX - 1; jx++) {
                NELECT += XSFIN -> grid[jx][jy][jz] * XSFIN -> VoxelV;
                NL1 += XSFOUT1 -> grid[jx][jy][jz] * XSFIN -> VoxelV;
            }
        }
    }
    printf("Total electron count in DS2_DEN file:  %lf\n", NELECT);
    printf("Total nonlocal energy for NL.xsf:  %lf\n", NL1);
    for (jy = 0; jy < XSFIN -> NGY - 1; jy++) {
        for (jx = 0; jx < XSFIN -> NGX - 1; jx++) {
            XSFOUT1 -> grid[jx][jy][NGZ - 1] = XSFOUT1 -> grid[jx][jy][0];
        }
    }
    for (jz = 0; jz < XSFIN -> NGZ - 1; jz++) {
        for (jx = 0; jx < XSFIN -> NGX - 1; jx++) {
            XSFOUT1 -> grid[jx][NGY - 1][jz] = XSFOUT1 -> grid[jx][0][jz];
        }
    }
    for (jz = 0; jz < XSFIN -> NGZ - 1; jz++) {
        for (jy = 0; jy < XSFIN -> NGY - 1; jy++) {
            XSFOUT1 -> grid[NGX - 1][jy][jz] = XSFOUT1 -> grid[0][jy][jz];
        }
    }
    for (jz = 0; jz < XSFIN -> NGZ - 1; jz++) {
        XSFOUT1 -> grid[NGX - 1][NGY - 1][jz] = XSFOUT1 -> grid[0][0][jz];
    }
    for (jy = 0; jy < XSFIN -> NGY - 1; jy++) {
        XSFOUT1 -> grid[NGX - 1][jy][NGZ - 1] = XSFOUT1 -> grid[0][jy][0];
    }
    for (jx = 0; jx < XSFIN -> NGX - 1; jx++) {
        XSFOUT1 -> grid[jx][NGY - 1][NGZ - 1] = XSFOUT1 -> grid[jx][0][0];
    }
    XSFOUT1 -> grid[NGX - 1][NGY - 1][NGZ - 1] = XSFOUT1 -> grid[0][0][0];
}

void read_psp_data(char filename[200], struct psp_file * pspin) {
    FILE * f2;
    int jx, jy, jz, ha, hb, hc;
    int j, l;
    int line_counter;
    int nvoxels_atm;
    int check;
    int stop = 0;
    int stop2;
    char str[1000];
    int NGX, NGY, NGZ, k, ixc = -10;
    double cell_volume, dist;
    double rcut;
    double xf, yf, zf;
    double voxel_centerx;
    double voxel_centery;
    double voxel_centerz;
    int n_epsatm_values = 0;
    double Z;
    int typat;
    double Z_ion_temp;
    double Z_ion[NATOMS_MAX];
    int Z_values[NATOMS_MAX];
    double epsatm_values[NATOMS_MAX];
    f2 = fopen(filename, "r");
    typat = 0;
    ixc = 1;
    while (stop == 0) {
        check = fscanf(f2, "%s", str);
        if (check == EOF) {
            stop = 1;
        }
        if (strcmp(str, "ixc") == 0) {
            fscanf(f2, "%d", & ixc);
            printf("XC correlation functional = %d\n", ixc);
        }
        if (strcmp(str, "ETOT") == 0) stop = 1;
        if ((strcmp(str, "pspini:") == 0) && (ixc == 1)) {
            finish_line(f2);
            finish_line(f2);
            finish_line(f2);
            fscanf(f2, "%s", str);
            fscanf(f2, "%lf %lf", & Z, & Z_ion[typat]);
            Z_values[typat] = Z;
            stop2 = 0;
            while (stop2 == 0) {
                check = fscanf(f2, "%s", str);

                if (strcmp(str, "rloc=") == 0) {
                    fscanf(f2, "%lf", & pspin -> rloc[typat]);
                }
                if (strcmp(str, "cc1") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc1[typat]);
                }
                if (strcmp(str, "cc2") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc2[typat]);
                }
                if (strcmp(str, "cc3") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc3[typat]);
                }
                if (strcmp(str, "cc4") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc4[typat]);
                }
                if (strcmp(str, "rrs") == 0) {
                    fscanf(f2, " = %lf", & pspin -> rrs[typat]);
                }
                if (strcmp(str, "h11s=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[0][0][0][typat]);
                }
                if (strcmp(str, "h22s=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[1][1][0][typat]);
                }
                if (strcmp(str, "h33s=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[2][2][0][typat]);
                }
                if (strcmp(str, "rrp") == 0) {
                    fscanf(f2, " = %lf", & pspin -> rrp[typat]);
                }
                if (strcmp(str, "h11p=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[0][0][1][typat]);
                }
                if (strcmp(str, "h22p=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[1][1][1][typat]);
                }
                if (strcmp(str, "h33p=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[2][2][1][typat]);
                }
                if (strcmp(str, "k11p=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k11p[typat]);
                }
                if (strcmp(str, "k22p=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k22p[typat]);
                }
                if (strcmp(str, "k33p=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k33p[typat]);
                }
                if (strcmp(str, "rrd") == 0) {
                    fscanf(f2, " = %lf", & pspin -> rrd[typat]);
                }
                if (strcmp(str, "h11d=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[0][0][2][typat]);
                }
                if (strcmp(str, "h22d=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[1][1][2][typat]);
                }
                if (strcmp(str, "h33d=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[2][2][2][typat]);
                }
                if (strcmp(str, "k11d=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k11d[typat]);
                }
                if (strcmp(str, "k22d=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k22d[typat]);
                }
                if (strcmp(str, "k33d=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k33d[typat]);
                }
                if (strcmp(str, "rrf") == 0) {
                    fscanf(f2, " = %lf", & pspin -> rrf[typat]);
                }
                if (strcmp(str, "h11f=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[0][0][3][typat]);
                }
                if (strcmp(str, "h22f=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[1][1][3][typat]);
                }
                if (strcmp(str, "h33f=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[2][2][3][typat]);
                }
                if (strcmp(str, "k11f=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k11f[typat]);
                }
                if (strcmp(str, "k22f=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k22f[typat]);
                }
                if (strcmp(str, "k33f=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k33f[typat]);
                }

                if (strcmp(str, "COMMENT") == 0) {
                    stop2 = 1;
                }
            }
            printf("PSEUDOPOTENTIAL %d\n", typat + 1);
            printf("   %lf  %lf \n", Z, Z_ion[typat]);
            printf("   rloc= %lf \n", pspin -> rloc[typat]);
            printf("   cc1 = %lf; cc2 = %lf; cc3 = %lf; cc4 = %lf \n", pspin -> cc1[typat], pspin -> cc2[typat], pspin -> cc3[typat], pspin -> cc4[typat]);
            printf("   rrs = %lf; h11s= %lf; h22s= %lf; h33s= %lf \n", pspin -> rrs[typat], pspin -> h[0][0][0][typat], pspin -> h[1][1][0][typat], pspin -> h[2][2][0][typat]);
            printf("   rrp = %lf; h11p= %lf; h22p= %lf; h33p= %lf \n", pspin -> rrp[typat], pspin -> h[0][0][1][typat], pspin -> h[1][1][1][typat], pspin -> h[2][2][1][typat]);
            printf("   rrp = %lf; k11p= %lf; k22p= %lf; k33p= %lf \n", pspin -> rrp[typat], pspin -> k11p[typat], pspin -> k22p[typat], pspin -> k33p[typat]);
            printf("   rrd = %lf; h11d= %lf; h22d= %lf; h33d= %lf \n", pspin -> rrd[typat], pspin -> h[0][0][2][typat], pspin -> h[1][1][2][typat], pspin -> h[2][2][2][typat]);
            printf("   rrd = %lf; k11d= %lf; k22d= %lf; k33d= %lf \n", pspin -> rrd[typat], pspin -> k11d[typat], pspin -> k22d[typat], pspin -> k33d[typat]);
            printf("   rrf = %lf; h11f= %lf; h22f= %lf; h33f= %lf \n", pspin -> rrf[typat], pspin -> h[0][0][3][typat], pspin -> h[1][1][3][typat], pspin -> h[2][2][3][typat]);
            printf("   rrf = %lf; k11f= %lf; k22f= %lf; k33f= %lf \n", pspin -> rrf[typat], pspin -> k11f[typat], pspin -> k22f[typat], pspin -> k33f[typat]);
            printf("\n");
            pspin -> h[0][1][0][typat] = -0.5 * pow(3.0 / 5.0, 0.5) * pspin -> h[1][1][0][typat];
            pspin -> h[1][0][0][typat] = pspin -> h[0][1][0][typat];
            pspin -> h[0][2][0][typat] = 0.5 * pow(5.0 / 21.0, 0.5) * pspin -> h[2][2][0][typat];
            pspin -> h[2][0][0][typat] = pspin -> h[0][2][0][typat];
            pspin -> h[1][2][0][typat] = -0.5 * pow(100.0 / 63.0, 0.5) * pspin -> h[2][2][0][typat];
            pspin -> h[2][1][0][typat] = pspin -> h[1][2][0][typat];
            pspin -> h[0][1][1][typat] = -0.5 * pow(5.0 / 7.0, 0.5) * pspin -> h[1][1][1][typat];
            pspin -> h[1][0][1][typat] = pspin -> h[0][1][1][typat];
            pspin -> h[0][2][1][typat] = (1.0 / 6.0) * pow(35.0 / 11.0, 0.5) * pspin -> h[2][2][1][typat];
            pspin -> h[2][0][1][typat] = pspin -> h[0][2][1][typat];
            pspin -> h[1][2][1][typat] = -(1.0 / 6.0) * 14 * pow(1.0 / 11.0, 0.5) * pspin -> h[2][2][1][typat];
            pspin -> h[2][1][1][typat] = pspin -> h[1][2][1][typat];
            pspin -> h[0][1][2][typat] = -(0.5) * pow(7.0 / 9.0, 0.5) * pspin -> h[1][1][2][typat];
            pspin -> h[1][0][2][typat] = pspin -> h[0][1][2][typat];
            pspin -> h[0][2][2][typat] = (0.5) * pow(63.0 / 143.0, 0.5) * pspin -> h[2][2][2][typat];
            pspin -> h[2][0][2][typat] = pspin -> h[0][2][2][typat];
            pspin -> h[1][2][2][typat] = -(0.5) * 18.0 * pow(1.0 / 143.0, 0.5) * pspin -> h[2][2][2][typat];
            pspin -> h[2][1][2][typat] = pspin -> h[1][2][2][typat];
            typat++;
        }
        if ((strcmp(str, "pspini:") == 0) && (ixc == 11)) {
            finish_line(f2);
            finish_line(f2);
            finish_line(f2);
            fscanf(f2, "%s", str);
            fscanf(f2, "%lf %lf", & Z, & Z_ion[typat]);
            Z_values[typat] = Z;
            stop2 = 0;
            l = -1;
            while (stop2 == 0) {
                check = fscanf(f2, "%s", str);
                if (strcmp(str, "angular") == 0) l++;
                if (strcmp(str, "rloc=") == 0) {
                    fscanf(f2, "%lf", & pspin -> rloc[typat]);
                }
                if (strcmp(str, "cc(1:1)=") == 0) {
                    fscanf(f2, " %lf", & pspin -> cc1[typat]);
                }
                if (strcmp(str, "cc2") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc2[typat]);
                }
                if (strcmp(str, "cc3") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc3[typat]);
                }
                if (strcmp(str, "cc4") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc4[typat]);
                }
                if (strcmp(str, "r(l)") == 0) {
                    if (l == 0) fscanf(f2, " = %lf", & pspin -> rrs[typat]);
                    if (l == 1) fscanf(f2, " = %lf", & pspin -> rrp[typat]);
                    if (l == 2) fscanf(f2, " = %lf", & pspin -> rrd[typat]);
                    if (l == 3) fscanf(f2, " = %lf", & pspin -> rrf[typat]);
                }
                if (strcmp(str, "h11,") == 0) {
                    fscanf(f2, "%s", str);
                    fscanf(f2, "%s", str);
                    if (l == 0) fscanf(f2, " = %lf %lf %lf", & pspin -> h[0][0][0][typat], & pspin -> h[0][1][0][typat], & pspin -> h[0][2][0][typat]);
                    if (l == 1) fscanf(f2, " = %lf %lf %lf", & pspin -> h[0][0][1][typat], & pspin -> h[0][1][1][typat], & pspin -> h[0][2][1][typat]);
                    if (l == 2) fscanf(f2, " = %lf %lf %lf", & pspin -> h[0][0][2][typat], & pspin -> h[0][1][2][typat], & pspin -> h[0][2][2][typat]);
                    if (l == 3) fscanf(f2, " = %lf %lf %lf", & pspin -> h[0][0][3][typat], & pspin -> h[0][1][3][typat], & pspin -> h[0][2][3][typat]);
                }
                if (strcmp(str, "h22,") == 0) {
                    fscanf(f2, "%s", str);
                    if (l == 0) fscanf(f2, " = %lf %lf", & pspin -> h[1][1][0][typat], & pspin -> h[1][2][0][typat]);
                    if (l == 1) fscanf(f2, " = %lf %lf", & pspin -> h[1][1][1][typat], & pspin -> h[1][2][1][typat]);
                    if (l == 2) fscanf(f2, " = %lf %lf", & pspin -> h[1][1][2][typat], & pspin -> h[1][2][2][typat]);
                    if (l == 3) fscanf(f2, " = %lf %lf", & pspin -> h[1][1][3][typat], & pspin -> h[1][2][3][typat]);
                }
                if (strcmp(str, "h33") == 0) {
                    if (l == 0) fscanf(f2, " = %lf", & pspin -> h[2][2][0][typat]);
                    if (l == 1) fscanf(f2, " = %lf", & pspin -> h[2][2][1][typat]);
                    if (l == 2) fscanf(f2, " = %lf", & pspin -> h[2][2][2][typat]);
                    if (l == 3) fscanf(f2, " = %lf", & pspin -> h[2][2][3][typat]);
                }
                if (strcmp(str, "k11,") == 0) {
                    fscanf(f2, "%s", str);
                    fscanf(f2, "%s", str);
                    //                    if(l==0) fscanf(f2,"= %lf %lf %lf",&pspin->k11s[typat],&pspin->k12s[typat],&pspin->k13s[typat]);
                    if (l == 1) fscanf(f2, " = %lf %lf %lf", & pspin -> k11p[typat], & pspin -> k12p[typat], & pspin -> k13p[typat]);
                    if (l == 2) fscanf(f2, " = %lf %lf %lf", & pspin -> k11d[typat], & pspin -> k12d[typat], & pspin -> k13d[typat]);
                    if (l == 3) fscanf(f2, " = %lf %lf %lf", & pspin -> k11f[typat], & pspin -> k12f[typat], & pspin -> k13f[typat]);
                }
                if (strcmp(str, "k22,") == 0) {
                    fscanf(f2, "%s", str);
                    //                    if(l==0) fscanf(f2,"= %lf %lf",&pspin->k22s[typat],&pspin->k23s[typat]);
                    if (l == 1) fscanf(f2, " = %lf %lf", & pspin -> k22p[typat], & pspin -> k23p[typat]);
                    if (l == 2) fscanf(f2, " = %lf %lf", & pspin -> k22d[typat], & pspin -> k23d[typat]);
                    if (l == 3) fscanf(f2, " = %lf %lf", & pspin -> k22f[typat], & pspin -> k23f[typat]);
                }
                if (strcmp(str, "k33") == 0) {
                    fscanf(f2, "%s", str);
                    //                   if(l==0) fscanf(f2,"= %lf",&pspin->k33s[typat]);
                    if (l == 1) fscanf(f2, " = %lf", & pspin -> k33p[typat]);
                    if (l == 2) fscanf(f2, " = %lf", & pspin -> k33d[typat]);
                    if (l == 3) fscanf(f2, " = %lf", & pspin -> k33f[typat]);
                }
                if (strcmp(str, "COMMENT") == 0) {
                    stop2 = 1;
                }
            }
            printf("PSEUDOPOTENTIAL %d\n", typat + 1);
            printf("   %lf  %lf \n", Z, Z_ion[typat]);
            printf("   rloc= %lf \n", pspin -> rloc[typat]);
            printf("   cc1 = %lf; cc2 = %lf; cc3 = %lf; cc4 = %lf \n", pspin -> cc1[typat], pspin -> cc2[typat], pspin -> cc3[typat], pspin -> cc4[typat]);
            printf("   rrs = %lf; h11s= %lf; h22s= %lf; h33s= %lf \n", pspin -> rrs[typat], pspin -> h[0][0][0][typat], pspin -> h[1][1][0][typat], pspin -> h[2][2][0][typat]);
            printf("   rrp = %lf; h11p= %lf; h22p= %lf; h33p= %lf \n", pspin -> rrp[typat], pspin -> h[0][0][1][typat], pspin -> h[1][1][1][typat], pspin -> h[2][2][1][typat]);
            printf("   rrp = %lf; k11p= %lf; k22p= %lf; k33p= %lf \n", pspin -> rrp[typat], pspin -> k11p[typat], pspin -> k22p[typat], pspin -> k33p[typat]);
            printf("   rrd = %lf; h11d= %lf; h22d= %lf; h33d= %lf \n", pspin -> rrd[typat], pspin -> h[0][0][2][typat], pspin -> h[1][1][2][typat], pspin -> h[2][2][2][typat]);
            printf("   rrd = %lf; k11d= %lf; k22d= %lf; k33d= %lf \n", pspin -> rrd[typat], pspin -> k11d[typat], pspin -> k22d[typat], pspin -> k33d[typat]);
            printf("   rrf = %lf; h11f= %lf; h22f= %lf; h33f= %lf \n", pspin -> rrf[typat], pspin -> h[0][0][3][typat], pspin -> h[1][1][3][typat], pspin -> h[2][2][3][typat]);
            printf("   rrf = %lf; k11f= %lf; k22f= %lf; k33f= %lf \n", pspin -> rrf[typat], pspin -> k11f[typat], pspin -> k22f[typat], pspin -> k33f[typat]);
            printf("\n");
            //               pspin->h[0][1][0][typat]=-0.5*pow(3.0/5.0,0.5)*pspin->h[1][1][0][typat];
            pspin -> h[1][0][0][typat] = pspin -> h[0][1][0][typat];
            //               pspin->h[0][2][0][typat]=0.5*pow(5.0/21.0,0.5)*pspin->h[2][2][0][typat];
            pspin -> h[2][0][0][typat] = pspin -> h[0][2][0][typat];
            //               pspin->h[1][2][0][typat]=-0.5*pow(100.0/63.0,0.5)*pspin->h[2][2][0][typat];
            pspin -> h[2][1][0][typat] = pspin -> h[1][2][0][typat];
            //               pspin->h[0][1][1][typat]=-0.5*pow(5.0/7.0,0.5)*pspin->h[1][1][1][typat];
            pspin -> h[1][0][1][typat] = pspin -> h[0][1][1][typat];
            //              pspin->h[0][2][1][typat]=(1.0/6.0)*pow(35.0/11.0,0.5)*pspin->h[2][2][1][typat];
            pspin -> h[2][0][1][typat] = pspin -> h[0][2][1][typat];
            //               pspin->h[1][2][1][typat]=-(1.0/6.0)*14*pow(1.0/11.0,0.5)*pspin->h[2][2][1][typat];
            pspin -> h[2][1][1][typat] = pspin -> h[1][2][1][typat];
            //               pspin->h[0][1][2][typat]=-(0.5)*pow(7.0/9.0,0.5)*pspin->h[1][1][2][typat];
            pspin -> h[1][0][2][typat] = pspin -> h[0][1][2][typat];
            //               pspin->h[0][2][2][typat]=(0.5)*pow(63.0/143.0,0.5)*pspin->h[2][2][2][typat];
            pspin -> h[2][0][2][typat] = pspin -> h[0][2][2][typat];
            //               pspin->h[1][2][2][typat]=-(0.5)*18.0*pow(1.0/143.0,0.5)*pspin->h[2][2][2][typat];
            pspin -> h[2][1][2][typat] = pspin -> h[1][2][2][typat];
            printf(" h11s h12s h13s = %lf %lf %lf \n", pspin -> h[0][0][0][typat], pspin -> h[0][1][0][typat], pspin -> h[0][2][0][typat]);
            printf("      h22s h23s =     %lf %lf \n", pspin -> h[1][1][0][typat], pspin -> h[1][2][0][typat]);
            printf("           h33s =     %lf %lf \n", pspin -> h[2][2][0][typat]);
            printf(" h11p h12p h13p = %lf %lf %lf \n", pspin -> h[0][0][1][typat], pspin -> h[0][1][1][typat], pspin -> h[0][2][1][typat]);
            printf("      h22p h23p =     %lf %lf \n", pspin -> h[1][1][1][typat], pspin -> h[1][2][1][typat]);
            printf("           h33p =     %lf %lf \n", pspin -> h[2][2][1][typat]);
            printf(" h11d h12d h13d = %lf %lf %lf \n", pspin -> h[0][0][2][typat], pspin -> h[0][1][2][typat], pspin -> h[0][2][2][typat]);
            printf("      h22d h23d =     %lf %lf \n", pspin -> h[1][1][2][typat], pspin -> h[1][2][2][typat]);
            printf("           h33d =     %lf %lf \n", pspin -> h[2][2][2][typat]);
            printf(" h11f h12f h13f = %lf %lf %lf \n", pspin -> h[0][0][3][typat], pspin -> h[0][1][3][typat], pspin -> h[0][2][3][typat]);
            printf("      h22f h23f =     %lf %lf \n", pspin -> h[1][1][3][typat], pspin -> h[1][2][3][typat]);
            printf("           h33f =     %lf %lf \n", pspin -> h[2][2][3][typat]);
            typat++;
        }
    }
    fclose(f2);
}

int main(int argc, char * argv[]) {
    int j;
    int nparams;
    FILE * f2;
    FILE * f4;
    double dV;
    int choice;
    char tim0[100];
    char tim1[100];
    char tim2[100];
    char KDENfile1[100];
    char KDENfile2[100];
    char POTfile1[100];
    char POTfile2[100];
    char DENfile0[100];
    char DENfile1[100];
    char DENfile2[100];
    char VHXCfile1[100];
    char VHXCfile2[100];
    char VHAfile1[100];
    char VHAfile2[100];
    char systemcmd[200];
    char WFKfilename[100];
    char DENfilename[100];
    char WFKfilename1[100];
    char WFKfilename2[100];
    char NLfilename1[100];
    char NLfilename2[100];
    char outfilename[100];
    char outfilename2[100];
    char tmp1[100], tmp2[100];
    double res;
    double V1;
    double V2;
    int stop;
    int k_min, k_max, k;
    double Ealpha_map1, Ealpha_map2, Ealpha_nomap1, Ealpha_nomap2;
    double temp1, temp2, temp3, temp4, temp5;
    double P_entropy;
    double P_Ewald;
    double P_nonlocal;
    double P_Ealpha;
    double min_occup;
    int dtsetnum;
    if (argc > 4) {
        strcpy(outfilename, argv[1]);
        strcpy(WFKfilename, argv[2]);
        sscanf(argv[4], "%lf", & min_occup);
        sscanf(argv[3], "%d", & dtsetnum);
    } else {
        printf("Usage:  nonlocal14 <_out file> <_o base> <dtsetnum> <minimum band occupancy>\n");
        exit(0);
    }
    read_psp_data(outfilename, & psp);
    sprintf(WFKfilename1, "%s_o_DS%d_WFK", WFKfilename, dtsetnum);
    sprintf(NLfilename1, "%s_o_DS%d_NL.xsf", WFKfilename, dtsetnum);
    strcpy(DENfilename, WFKfilename);
    strcat(DENfilename, "_o_DS2_DEN");
    printf("_HERE!\n");
    read_calc_den(DENfilename, & den0);
    printf("HERE@!\n");
    for (j = 0; j < den0.NIONS; j++) {
        nonlocalE_byatom[0][j][5] = 0;
        nonlocalE_byatom[0][j][0] = 0;
        nonlocalE_byatom[0][j][1] = 0;
        nonlocalE_byatom[0][j][2] = 0;
        nonlocalE_byatom[0][j][3] = 0;
    }
    stop = 0;
    k = 0;
    while (stop == 0) {
        k++;
        sprintf(outfilename2, "%s_o_DS%d_NL-%d.log", WFKfilename, dtsetnum, k);
        f2 = fopen(outfilename2, "r");
        if (f2 == NULL) {
            printf("%s not found. \n", outfilename2);
            stop = 1;
        } else {
            finish_line(f2);
            printf("Reading nonlocal components from pw group %d...\n", k);
            for (j = 0; j < den0.NIONS; j++) {
                //           fscanf(f2,"%s %s %lf %lf %lf %lf %lf",tmp1,tmp2,&nonlocalE_byatom[0][j][5],&nonlocalE_byatom[0][j][0],&nonlocalE_byatom[0][j][1],&nonlocalE_byatom[0][j][2],&nonlocalE_byatom[0][j][3]);
                fscanf(f2, "%s %s %lf %lf %lf %lf %lf", tmp1, tmp2, & temp1, & temp2, & temp3, & temp4, & temp5);
                nonlocalE_byatom[0][j][5] += temp1;
                nonlocalE_byatom[0][j][0] += temp2;
                nonlocalE_byatom[0][j][1] += temp3;
                nonlocalE_byatom[0][j][2] += temp4;
                nonlocalE_byatom[0][j][3] += temp5;
                finish_line(f2);
                printf("    atom %d:  nonlocalE = %13.8lf = ", j + 1, nonlocalE_byatom[0][j][5]);
                printf(" (s) %13.8lf + ", nonlocalE_byatom[0][j][0]);
                printf(" (p) %13.8lf + ", nonlocalE_byatom[0][j][1]);
                printf(" (d) %13.8lf + ", nonlocalE_byatom[0][j][2]);
                printf(" (f) %13.8lf \n", nonlocalE_byatom[0][j][3]);
            }
            fclose(f2);
        }
    }
    temp1 = 0;
    for (j = 0; j < den0.NIONS; j++) {
        temp1 += nonlocalE_byatom[0][j][5];
    }
    printf("  Total NL energy:  %lf\n", temp1);
    printf("Proceeding to map generation...\n");
    NL2XSF( & den0, & nlden1);
    f4 = fopen(NLfilename1, "w");
    outputXSF( & den0, & nlden1, f4);
    fclose(f4);
}
