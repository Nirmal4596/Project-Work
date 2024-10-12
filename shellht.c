#include<stdio.h>
#include<conio.h>
#include<stdlib.h>
#include<math.h>

FILE *fin,*fout,*fst,*flm;
int nel, nsn, nmat, ox, oe, oz, npr = 2, *n, *ne, *lmc, *matid, nnl, nnt,dnnt,i,j,k,    /*npr: no. of material property per material id*/
    *nt, *dnt, ix, ie, iz, ii, jj, ndn /* no of degrees of freedom per node*/, it, jt, nct, nrt, nr, nc, iel, nng, p, q, r;
double ts,Tx,*dT,*ang,*r0,*y,**kel,*pm, N[9], dN[9][4],*NN,**nn,**dN1,**dN2,**dN3,**B,**R,**RT, **KR, K[4][4], th = 0,/*  change for 3d */
       det, **ke, **kh, *qs, *qq, *qgc, *f, c, *s, *T, xm, ym, zm, xdh, ydh, zdh, J[4][4], Jinv[4][4], *xt, *yt, *zt, *t, *xb, *yb, *zb, xi0, eta0,
       XIn[] = { 0, -1, 0, 1, 1, 1, 0, -1, -1 }, ETAn[] = { 0, -1, -1, -1, 0, 1, 1, 1, 0 }, vf, cnst, k1, k2,pf,pi,rad,teta,vin /* volume fraction*/;
char type;


#include "C:\\Users\\Nirmal\\Desktop\\shellht\\shellht\\arraypgm.c"
#include "C:\\Users\\Nirmal\\Desktop\\shellht\\shellht\\matrixpgm.c"
void main()
{
	/*void shape(void);*/
	/*void jacobian(void);*/
	int factorial(int);
	double hs, hg, *hf, *hb;
	int *hgc, *hsc;
	fin = fopen("inputs.txt", "r");
	fout = fopen("outputs", "w+");
	fst = fopen("sts", "w+");
	flm = fopen("lms", "w+");

	if (fin == 0)
	{
		printf("\ndat : File not found ! \n\a");
		exit(0);
	}


	fprintf(fout, "FE program : Degenerated Shell Element.");
	fprintf(fout, "\n---------------------------------------\n\n");
	printf("FE program : Degenerated Shell Element.");
	printf("\n---------------------------------------\n\n");

	fscanf(fin, "%c", &type);

	fscanf(fin, "%d %d %d %d %d %d %d", &nsn, &nel, &nmat, &ox, &oe, &oz,&ndn);
	/* nsn - total number of nodes.    */
	/* t   - Shell thickness.          */
	/* nel - total number of elements. */
	/* nmat no. of materials            */
	/* ox oe oz: Gauss integration order */
	

	fprintf(fout, "\nTotal number of nodes    : %d", nsn);
	fprintf(fout, "\nTotal number of elements : %d", nel);
	fprintf(fout, "\nTotal number of material sets: %d", nmat);
	fprintf(fout, "\nGauss Integration order : %d %d %d", ox, oe, oz);
	fprintf(fout, "\nHierarchial Temperature approximation order: %d", ndn - 1);
	/* fprintf(fout, "\nShell thickness          : %f\n", t);  */

	/* Memory allocation for nodal data:
	x, y, z co-ordinate;
	Shell thickness at each node;
	material data.                                    */

	n = (int *)calloc(nsn + 1, sizeof(int));     /* n for node number storage*/
	ne = (int *)calloc(nel + 1, sizeof(int));	 /* ne for element number storage*/
	lmc = (int*)calloc(1 + nel * 8, sizeof(int));
	r0 = (double*)calloc(nsn + 1, sizeof(double));
	ang = (double*)calloc(nsn + 1, sizeof(double));
	y = (double*)calloc(nsn + 1, sizeof(double));
	xt = (double *)calloc(nsn + 1, sizeof(double));
	xb = (double *)calloc(nsn + 1, sizeof(double));
	yt = (double *)calloc(nsn + 1, sizeof(double));
	yb = (double *)calloc(nsn + 1, sizeof(double));
	zt = (double *)calloc(nsn + 1, sizeof(double));
	zb = (double *)calloc(nsn + 1, sizeof(double));
	t = (double *)calloc(nsn + 1, sizeof(double));
	matid = (int *)calloc(nel + 1, sizeof(int));
	/*pm = (double *)calloc(nmat*npr, sizeof(double));*/

	if (type == 'p')
	{
		fprintf(fout, "\n\n\nFlat plate\n");
		printf("\n\n\nFlat plate\n");
		i = 1;
		/* Read nodal data*/
		do{
			fscanf(fin, "%d %lf %lf %lf %lf", (n + i), (xt + i), (yt + i), (zt + i), (t + i));
			*(xb + i) = *(xt + i);
			*(yb + i) = *(yt + i);
			*(zb + i) = *(zt + i) - *(t + i);
			i++;
		} while (i <= nsn);
	}

	pi = 4.0 * atan(1.0);
	rad = pi / 180;
	if (type == 's')
	{
		fprintf(fout, "\n\n\nShell\n");
		printf("\n\n\nShell\n");
		i = 1;
		do{
			fscanf(fin, "%d %lf %lf %lf %lf", (n + i), (r0 + i), (ang+i),(t + i), (y + i));
			i++;
		} while (i <= nsn);

		j = 1;
		do{
			fprintf(fout, "\n%d  %.4lf  %.4lf  %.4lf  %.4lf", *(n + j), *(r0 + j), *(ang + j), *(t + j), *(y + j));
			j++;
		} while (j <= nsn);

		i = 1;
		do{
			teta = *(ang + i)*rad;
			*(xb + i) = (*(r0 + i)- *(t+i)/2)*sin(teta);
			*(xt + i) = (*(r0 + i) + *(t + i) / 2)*sin(teta);
			*(zb + i) = (*(r0 + i)- *(t+i)/2)*cos(teta)/*-*(r0+i)*/;
			*(zt + i) = (*(r0 + i) + *(t + i)/2)*cos(teta)/* - *(r0 + i)*/;
			/*(zb + i) = -(*(r0 + i) - *(t + i) / 2)*cos(teta) + *(r0 + i);
			*(zt + i) = -(*(r0 + i) + *(t + i) / 2)*cos(teta) + *(r0 + i); */
			*(yt + i) = *(yb + i) = *(y + i);
			/*(zt + i) = *(zb + i) = *(y + i)*(-1);
			*(yb + i) = ((*(r0 + i) - (*(t + i) / 2))*cos(teta));
			*(yt + i) = ((*(r0 + i) + (*(t + i) / 2))*cos(teta));*/
			i++;
		} while (i <= nsn);

	}

	/* Read element connectivity & material ID */
	int lmcr[9];
	for (i = 1; i <= nel;)
	{
		fscanf(fin, "%d", &ne[i]);
		for (j = 1; j <= 8;)
		{
			fscanf(fin, "%d", &lmcr[j]);
			*(lmc + (i - 1) * 8 + j) = lmcr[j];   // lmc[(i-1)*8 + j] and *(lmc + (i - 1) * 8 + j) are same
			j++;
		}
		fscanf(fin, "%d", &matid[i]);
		i++;
	}

	/* Write nodal data in the file*/
	fprintf(fout, " \n\nNode no.  Co-ordinates Top: X  Y  Z  Thickness  Cordinates Bottom: X Y Z ");
	i = 1;
	do{
		fprintf(fout, "\n%d  \t %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf", *(n + i), *(xt + i), *(yt + i), *(zt + i), *(t + i), *(xb + i), *(yb + i), *(zb + i));
		i++;
	} while (i <= nsn);

	/* Print connectivity data in file*/
	fprintf(fout, "\n\nNode no. \t Element Connectivity    Material ID \n");
	for (i = 1; i <= nel;)
	{
		fprintf(fout, "\n %d\t", ne[i]);
		for (j = 1; j <= 8;)
		{
			fprintf(fout, " %d", *(lmc + (i - 1) * 8 + j));
			j++;
		}
		fprintf(fout, "\t   %d", matid[i]);
		i++;
	}

	/* Boundary data entry */

	/*nnt: no. of node in which Temperature is specified*/
	/*dnnt: no. of node in which Temperature derivative is specified*/
	/* nt: Node no. in which Temperature is specified*/
	/*dnt: Node no. in which Temperature derivative is specified*/
	/*T: Temperature */
	/*dT: Temperature derivative*/

	fscanf(fin, "%d %d", &nnt,&dnnt);
	/* Memory allocation to read temperature boundary condition*/
	nt = ialloc1(nnt + 1);
	T = dalloc1(nnt + 1);
	dnt = ialloc1(dnnt + 1);
	dT = dalloc1(dnnt + 1);
	int nno, cd,dcd;
	double temp, dtemp;
	i = 1;
	j = 1;
	k = 1;
	fprintf(fout, "\n\nNode no.  Condition  Temperature");
	do{
		fscanf(fin, "%d %d %lf %d %lf", &nno, &cd, &temp, &dcd, &dtemp);
		if (cd != 0)
		{
			nt[j] = nno;
			T[j] = temp;
			j++;
		}
		if (dcd != 0)
		{
			dnt[k] = nno;
			dT[k] = dtemp;
			k++;
		}
		fprintf(fout, "\n  %d    %d    %.2lf    %d    %.2lf", nno, cd, temp, dcd, dtemp);  /* prints to file, the temeperature boundary condition at all nodes*/
		i++;
	} while (i <= nsn);

	i = 1;
	fprintf(fout, "\n\n Temperature Boundary condition Specified");
	if (j > 1)
	{
		do
		{
			fprintf(fout, "\n %d \t %.2lf", nt[i], T[i]);
			i++;
		} while (i <= nnt);
	}

	i = 1;
	fprintf(fout, "\n\n Temperature Gradient Boundary condition Specified");
	if (k > 1)
	{
		do
		{
			fprintf(fout, "\n %d \t %.2lf", dnt[i], dT[i]);
			i++;
		} while (i <= dnnt);
	}
	/* Read flux boundary conditions */

	hgc = ialloc1(nel + 1);		/*hgc: Specifies if heat generation is present in an element:-  1: present; 0: No heat gen.*/
	hf = dalloc1(nel + 1);		/*hf: Heat flux at the surface*/
	hsc = ialloc1(nel + 1);		/*hsc: Specifies if surface convective is present in an element:-  1: present; 0: no convection*/
	/*hb = (double *)calloc(nel+1, sizeof(double));		hb: Boundary convective eo-efficient*/


	i = 1;
	do
	{
		fscanf(fin, "%d %d %lf %d", &ne[i], &hgc[i], &hf[i], &hsc[i]);
		i++;
	} while (i <= nel);

	/* Write boundary conditions to the file*/
	j = 1;
	fprintf(fout, "\n\nElement no.  hgc \t hf \t hsc");
	do
	{
		fprintf(fout, "\n %d   \t%d \t%.2lf \t%d ", ne[j], hgc[j], hf[j], hsc[j]);  /* Stored from 1 level*/
		j++;
	} while (j <= nel);

	/* Read boundary convection data*/

	/*                                          */

	/* Read material properties*/

	fscanf(fin, "%lf %lf %lf %lf %lf %lf", &k1, &k2, &vin, &ts, &hg, &hs);
	/*print material data*/
	fprintf(fout, "\n\nConduction co-efficient top: %.2lf \nConduction co-efficient bottom: %.2lf\nVolume fraction index: %.2lf", k1, k2, vin);
	fprintf(fout, "\nAmbient Temperature: %.2lf \nHeat generation: %.2lf \nSurface convective co-eff: %.2lf", ts, hg, hs);

	/* Bandwidth calculation */
	int nmin, nmax, nbw = 0;
	for (i = 1; i <= nel; i++)
	{
		nmin = *(lmc + (i - 1) * 8 + 1);
		nmax = *(lmc + (i - 1) * 8 + 1);
		for (j = 1; j <= 8; ++j)
		{
			if ((*(lmc + (i - 1) * 8 + j) < nmin))
				nmin = *(lmc + (i - 1) * 8 + j);

			if ((nmax < *(lmc + (i - 1) * 8 + j)))
				nmax = *(lmc + (i - 1) * 8 + j);
		}

		if (nbw < ((nmax*ndn) - (nmin*ndn - (ndn-1)) + 1))
			nbw = (nmax*ndn) - (nmin*ndn - (ndn-1)) + 1;
	}
	
	fprintf(fout, " \n\nBandwidth: %d", nbw);
	for (i = 1; i <= nel * 8; i++)
	{
		fprintf(flm, "\n %d \t %d", i, *(lmc + i));
	}

	/* Element stiffness matrix*/

	kel = dalloc2(8 * ndn + 1, 8 * ndn + 1);                    /* Stifness matrix per gauss point*/
	B = dalloc2(4, 1 + 8 * ndn);
	R = dalloc2(4, 1 + 8 * ndn);                               /* Strain displacement matrix*/
	RT = dalloc2(1 + 8 * ndn, 4);
	KR = dalloc2(4, 1 + 8 * ndn);                              /*    Store K*B    */
	/*initmat(K,2,2);*/      /*   material property matrix  */
	ke = dalloc2(1 + 8 * ndn, 1 + 8 * ndn);                    /* Element stiffness matrix */
	kh = dalloc2(1 + 8 * ndn, 1 + 8 * ndn);                  /* Element convective stiffness matrix*/
	qgc = dalloc1(1 + 8 * ndn);                                /* Element heat generation vector*/
	qs = dalloc1(1 + 8 * ndn);                                  /* Element surface convection heat vector*/
	qq = dalloc1(1 + 8 * ndn);                                  /* Element heat flux load vector*/
	f = dalloc1(nsn * ndn);                                 /* assembled global load vector*/
	/*c = hs;*/
	nn = dalloc2(9, 1 + ndn);
	NN = dalloc1(8 * ndn + 1);
	dN1 = dalloc2(9, 1 + ndn);
	dN2 = dalloc2(9, 1 + ndn);
	dN3 = dalloc2(9, 1 + ndn);


	/* Allocate memory for banded stiffness matrix*/
	s = dalloc1(nsn*nbw*ndn);   /* stored from level 0*/

	/* Initialize Gauss Integration points*/

	/*xi = dalloc2(7, 7);
	wxi = dalloc2(7, 7);
	eta = xi;
	weta = wxi;
	zeta = xi;
	wzeta = wxi;*/

	/*initmat(xi, 6, 6);
	initmat(wxi, 6, 6);*/
	double xi[5][5], eta[5][5], zeta[5][5], wxi[5][5], weta[5][5], wzeta[5][5];
	initmat(xi, 4, 4);  initmat(eta, 4, 4); initmat(zeta, 4, 4);  initmat(wxi, 4, 4);  initmat(weta, 4, 4); initmat(wzeta, 4, 4);
	zeta[1][1] = eta[1][1] = xi[1][1] = 0.0;
	zeta[2][1] = eta[2][1] = xi[2][1] = -0.5773502691896;  zeta[2][2] = eta[2][2] = xi[2][2] = 0.5773502691896;
	zeta[3][1] = eta[3][1] = xi[3][1] = -0.7745966692415;  zeta[3][2] = eta[3][2] = xi[3][2] = 0.0;
	zeta[3][3] = eta[3][3] = xi[3][3] = 0.7745966692415;
	zeta[4][1] = eta[4][1] = xi[4][1] = -0.339981; zeta[4][3] = eta[4][3] = xi[4][3] = 0.339981;
	zeta[4][2] = eta[4][2] = xi[4][2] = -0.861136; zeta[4][4] = eta[4][4] = xi[4][4] = 0.861136;

	wzeta[1][1] = weta[1][1] = wxi[1][1] = 2.0;
	wzeta[2][1] = weta[2][1] = wxi[2][1] = 1.0;  wzeta[2][2] = weta[2][2] = wxi[2][2] = 1.0;
	wzeta[3][1] = wzeta[3][3] = weta[3][1] = weta[3][3] = wxi[3][1] = wxi[3][3] = 0.55555555555556; wzeta[3][2] = weta[3][2] = wxi[3][2] = 0.8888888888889;
	wzeta[4][1] = weta[4][1] = wxi[4][1] = 0.652145; wzeta[4][3] = weta[4][3] = wxi[4][3] = 0.652145;
	wzeta[4][2] = weta[4][2] = wxi[4][2] = 0.347855; wzeta[4][4] = weta[4][4] = wxi[4][4] = 0.347855;
	

	int count = 0;
	for (iel = 1; iel <= nel; ++iel)
	{                                   /* Element loop starts*/

		initmat(kh, 8 * ndn, 8 * ndn);
		initmat(kel, 8 * ndn, 8 * ndn);
		initmat(ke, 8 * ndn, 8 * ndn);
		initmat(R, 3, 8 * ndn);
		initmat(B, 3, 8 * ndn);
		initmat(KR, 3, 8 * ndn);
		initmat(dN1, 8, ndn);
		initmat(dN2, 8, ndn);
		initmat(dN3, 8,ndn);
		initvec(qgc, 8 * ndn);
		initvec(qs, 8 * ndn);
		initvec(qq, 8 * ndn);
		initvec(NN, 8 * ndn);
		fprintf(fst, "\n\nElement No. %d\n", iel);
		for (ix = 1; ix <= ox; ++ix)      /* Gauss loop starts*/
		{
			for (ie = 1; ie <= oe; ++ie)
			{
				for (iz = 1; iz <= oz; ++iz)
				{

					/*shape();*/
					for (nnl = 1; nnl <= 8; ++nnl)                      /* Shape functions ordered in clockwise direction continuously*/
					{
						xi0 = XIn[nnl] * xi[ox][ix];
						eta0 = ETAn[nnl] * eta[oe][ie];
						if ((nnl == 1) || (nnl == 3) || (nnl == 5) || (nnl == 7))
						{
							N[nnl] = 0.25*(1 + xi0)*(1 + eta0)*(xi0 + eta0 - 1.0);
							dN[nnl][1] = 0.25*(1 + eta0)*(2.0*xi0 + eta0)*XIn[nnl];
							dN[nnl][2] = 0.25*(1 + xi0)*(xi0 + 2.0*eta0)*ETAn[nnl];
							dN[nnl][3] = 0.0;
						}
						else if ((nnl == 2) || (nnl == 6))
						{
							N[nnl] = 0.5*(1.0 - xi[ox][ix] * xi[ox][ix])*(1.0 + eta0);
							dN[nnl][1] = -xi[ox][ix] * (1.0 + eta0);
							dN[nnl][2] = 0.5*(1.0 - xi[ox][ix] * xi[ox][ix])*ETAn[nnl];
							dN[nnl][3] = 0.0;
						}
						else if ((nnl == 4) || (nnl == 8))
						{
							N[nnl] = 0.5*(1.0 - eta[oe][ie] * eta[oe][ie])*(1.0 + xi0);
							dN[nnl][1] = 0.5*(1.0 - eta[oe][ie] * eta[oe][ie])*XIn[nnl];
							dN[nnl][2] = -eta[oe][ie] * (1.0 + xi0);
							dN[nnl][3] = 0.0;
						}
					}

					fprintf(fst, "\n\n %d %d %d\n\n", ix, ie, iz);    /* Print shape function and derivatives*/
					for (nnl = 1; nnl <= 8; nnl++)
					{
						fprintf(fst, "  %lf ", N[nnl]);
					}
					fprintf(fst, "\n\tdN/dXi\tdN/dEta\n");
					for (nnl = 1; nnl <= 8; nnl++)
					{
						fprintf(fst, "\n\t%0.5lf \t%0.5lf", dN[nnl][1], dN[nnl][2]);

					}

					/*jacobian();*/
					initmat(J, 3, 3);
					fprintf(fst, "  \n  ");
					for (nnl = 1; nnl <= 8; ++nnl)
					{
						nng = *(lmc + (iel - 1) * 8 + nnl);
						xm = (xt[nng] + xb[nng]) / 2.0;
						ym = (yt[nng] + yb[nng]) / 2.0;
						zm = (zt[nng] + zb[nng]) / 2.0;
						xdh = (xt[nng] - xb[nng]) / 2.0;
						ydh = (yt[nng] - yb[nng]) / 2.0;
						zdh = (zt[nng] - zb[nng]) / 2.0;
						fprintf(fst, "  %d  %.3lf  %.3lf  %.3lf  %.3lf  %.3lf  %.3lf",nng, xm, ym, zm,xdh,ydh,zdh);

						/*J[1][1] += dN[nnl][1] * (xm + zeta[oz][iz] * xdh);
						J[1][2] += dN[nnl][1] * (ym + zeta[oz][iz] * ydh);
						J[1][3] += dN[nnl][1] * (zm + zeta[oz][iz] * zdh);
						J[2][1] += dN[nnl][2] * (xm + zeta[oz][iz] * xdh);
						J[2][2] += dN[nnl][2] * (ym + zeta[oz][iz] * ydh);
						J[2][3] += dN[nnl][2] * (zm + zeta[oz][iz] * zdh);
						J[3][1] += N[nnl] * xdh;
						J[3][2] += N[nnl] * ydh;
						J[3][3] += N[nnl] * zdh;*/

						/* Jacobian at mid surface */
						J[1][1] += dN[nnl][1] * (xm);
						J[1][2] += dN[nnl][1] * (ym);
						J[1][3] += dN[nnl][1] * (zm);
						J[2][1] += dN[nnl][2] * (xm);
						J[2][2] += dN[nnl][2] * (ym);
						J[2][3] += dN[nnl][2] * (zm);
						J[3][1] += N[nnl] * xdh;
						J[3][2] += N[nnl] * ydh;
						J[3][3] += N[nnl] * zdh;
					}

					det = J[1][1] * (J[2][2] * J[3][3] - J[3][2] * J[2][3])
						- J[1][2] * (J[2][1] * J[3][3] - J[3][1] * J[2][3])
						+ J[1][3] * (J[2][1] * J[3][2] - J[3][1] * J[2][2]);

					fprintf(fst, "\n");          /* print Jacobian*/
					for (i = 1; i <= 3; i++)
					{
						fprintf(fst, "\n");
						for (j = 1; j <= 3; j++)
						{
							fprintf(fst, "  %lf \t  ", J[i][j]);
						}
					}

					fprintf(fst, "\n\n%lf\n", det);           /* print determinant of  Jacobian*/

					

					if (det <= 0.0)
					{
						fprintf(fst, "\nElement no. %d", iel);
						fprintf(fst, "\nJacobian determinant is zero or negative !");
						fprintf(fst, "\n Jdet = %lf ", det);
						exit(0);
					}


					Jinv[1][1] = (J[2][2] * J[3][3] - J[3][2] * J[2][3]);
					Jinv[2][1] = -(J[2][1] * J[3][3] - J[3][1] * J[2][3]);
					Jinv[3][1] = (J[2][1] * J[3][2] - J[2][2] * J[3][1]);
					Jinv[1][2] = -(J[1][2] * J[3][3] - J[3][2] * J[1][3]);
					Jinv[2][2] = (J[1][1] * J[3][3] - J[3][1] * J[1][3]);
					Jinv[3][2] = -(J[1][1] * J[3][2] - J[3][1] * J[1][2]);
					Jinv[1][3] = (J[1][2] * J[2][3] - J[2][2] * J[1][3]);
					Jinv[2][3] = -(J[1][1] * J[2][3] - J[2][1] * J[1][3]);
					Jinv[3][3] = (J[1][1] * J[2][2] - J[2][1] * J[1][2]);
					for (i = 1; i <= 3; ++i)
					{
						for (j = 1; j <= 3; ++j)
						{
							Jinv[i][j] = Jinv[i][j] / det;
						}
					}

					
					for (i = 1; i <= 3; i++)       /* print inverse of Jacobian*/
					{
						fprintf(fst, "\n");
						for (j = 1; j <= 3; j++)
						{
							fprintf(fst, "  %lf \t", Jinv[i][j]);
						}
					}
					
					/* Hierarchical approximation function for temperature degree of freedom */
					for (nnl = 1; nnl <= 8; nnl++)
					{
						nng = *(lmc + 8 * (iel - 1) + nnl);

						for (i = 1; i <= ndn; i++)
						{
							xi0 = XIn[nnl] * xi[ox][ix];
							eta0 = ETAn[nnl] * eta[oe][ie];
							if ((nnl == 1) || (nnl == 3) || (nnl == 5) || (nnl == 7))
							{
								nn[nnl][i] = 0.25*(1 + xi0)*(1 + eta0)*(xi0 + eta0 - 1.0)*pow(zeta[oz][iz], i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);

							}
							else if ((nnl == 2) || (nnl == 6))
							{
								nn[nnl][i] = 0.5*(1.0 - xi[ox][ix] * xi[ox][ix])*(1.0 + eta0)*pow(zeta[oz][iz], i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);

							}
							else if ((nnl == 4) || (nnl == 8))
							{
								nn[nnl][i] = 0.5*(1.0 - eta[oe][ie] * eta[oe][ie])*(1.0 + xi0)*pow(zeta[oz][iz], i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);

							}
						}
					}

					/* Assembling Hierarchical Approximation function in 1D array */
					for (nnl = 1; nnl <= 8; nnl++)
					{
						for (i = 1; i <= ndn; i++)
						{
							NN[(nnl - 1)*ndn + i] = nn[nnl][i];
						}
					}

					/* derivative of hierarchical temperature approximation function*/
					for (nnl = 1; nnl <= 8; nnl++)
					{
						nng = *(lmc + 8 * (iel - 1) + nnl);

						for (i = 1; i <= ndn; i++)
						{
							xi0 = XIn[nnl] * xi[ox][ix];
							eta0 = ETAn[nnl] * eta[oe][ie];
							if ((nnl == 1) || (nnl == 3) || (nnl == 5) || (nnl == 7))
							{
								dN1[nnl][i] = 0.25*(1 + eta0)*(2.0*xi0 + eta0)*XIn[nnl] * pow(zeta[oz][iz], i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);
								dN2[nnl][i] = 0.25*(1 + xi0)*(xi0 + 2.0*eta0)*ETAn[nnl] * pow(zeta[oz][iz], i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);
								dN3[nnl][i] = 0.25*(1 + xi0)*(1 + eta0)*(xi0 + eta0 - 1.0)*pow(zeta[oz][iz], i - 2)*pow(t[nng] / 2, i - 1) / factorial(i - 2);
								dN3[nnl][1] = 0.0;
							}
							else if ((nnl == 2) || (nnl == 6))
							{
								dN1[nnl][i] = -xi[ox][ix] * (1.0 + eta0)*pow(zeta[oz][iz], i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);
								dN2[nnl][i] = 0.5*(1.0 - xi[ox][ix] * xi[ox][ix])*ETAn[nnl] * pow(zeta[oz][iz], i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);
								dN3[nnl][i] = 0.5*(1.0 - xi[ox][ix] * xi[ox][ix])*(1.0 + eta0)*pow(zeta[oz][iz], i - 2)*pow(t[nng] / 2, i - 1) / factorial(i - 2);
								dN3[nnl][1] = 0.0;
							}
							else if ((nnl == 4) || (nnl == 8))
							{
								dN1[nnl][i] = 0.5*(1.0 - eta[oe][ie] * eta[oe][ie])*XIn[nnl] * pow(zeta[oz][iz], i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);
								dN2[nnl][i] = -eta[oe][ie] * (1.0 + xi0)*pow(zeta[oz][iz], i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);
								dN3[nnl][i] = 0.5*(1.0 - eta[oe][ie] * eta[oe][ie])*(1.0 + xi0)*pow(zeta[oz][iz], i - 2)*pow(t[nng] / 2, i - 1) / factorial(i - 2);
								dN3[nnl][1] = 0.0;
							}
						}
					}

					fprintf(fst, "\n\n Hierarchical Function\n");
					for (i = 1; i <= ndn * 8; i++)
					{
						fprintf(fst, "  %.4lf   ", NN[i]);
					}


					/* B matrix*/

					for (nnl = 1; nnl <= 8; nnl++)
					{
						for (i = 1; i <= ndn; i++)
						{
							B[1][(nnl - 1)*ndn + i] = dN1[nnl][i];
							B[2][(nnl - 1)*ndn + i] = dN2[nnl][i];
							B[3][(nnl - 1)*ndn + i] = dN3[nnl][i];
						}
					}

					/* R matrix*/
					matmul(Jinv, B, R, 3, 3, 8*ndn);

					/*for (nnl = 1; nnl <= 8; nnl++)
					{
					for (i = 1; i <= ndn; i++)
					{
					R[1][(nnl - 1)*ndn + i] = Jinv[1][1] * B[1][(nnl - 1)*ndn + i] + Jinv[1][2] * B[2][(nnl - 1)*ndn + i] + Jinv[1][3] * B[3][(nnl - 1)*ndn + i];
					R[2][(nnl - 1)*ndn + i] = Jinv[2][1] * B[1][(nnl - 1)*ndn + i] + Jinv[2][2] * B[2][(nnl - 1)*ndn + i] + Jinv[2][3] * B[3][(nnl - 1)*ndn + i];
					R[3][(nnl - 1)*ndn + i] = Jinv[3][1] * B[1][(nnl - 1)*ndn + i] + Jinv[3][2] * B[2][(nnl - 1)*ndn + i] + Jinv[3][3] * B[3][(nnl - 1)*ndn + i];
					}
					}*/


					fprintf(fst, "\n\nB matrix");    /* Print B matrix*/
					for (j = 1; j <= 3; j++)
					{
						fprintf(fst, "\n");
						for (i = 1; i <= 8 * ndn; i++)
						{
							fprintf(fst, "  %5.5lf  ", B[j][i]);
						}
					}

					fprintf(fst, "\n\nR matrix");    /* Print R matrix*/
					for (j = 1; j <= 3; j++)
					{
						fprintf(fst, "\n");
						for (i = 1; i <= 8 * ndn; i++)
						{
							fprintf(fst, "  %5.5lf  ", R[j][i]);
						}
					}

					/* KR matrix*/
					cnst = (zeta[oz][iz] * 1.0 + 1.0) / 2.0;
					pf = pow(cnst, vin);
					vf = k2 + pf*(k1 - k2); /*vin- volume fraction index */
					K[1][1] = K[2][2] = K[3][3] = vf; K[1][2] = K[1][3] = K[2][1] = K[2][3] = K[3][1] = K[3][2] = 0;

					matmul(K, R, KR, 3, 3, 8*ndn);

					/*for (i = 1; i <= 8*ndn; i++)  /* [k]*[R]*/
					/*{
					KR[1][i] = K[1][1] * R[1][i] + K[1][2] * R[2][i] + K[1][3] * R[3][i];
					KR[2][i] = K[2][1] * R[1][i] + K[2][2] * R[2][i] + K[2][3] * R[3][i];
					KR[3][i] = K[3][1] * R[1][i] + K[3][2] * R[2][i] + K[3][3] * R[3][i];
					}*/

					/*for (j = 1; j <= 3; j++)
					{
						for (i = 1; i <= 8 * ndn; i++)  /* [k]*[R]*/
						/*{
							KR[j][i] = k2* R[j][i];
						}
					}*/

					fprintf(fst, "\n\n Volume fraction: %lf\t%lf", pf, vf);


					/* Transpose of R matrix*/
					for (i = 1; i <= 3; i++)
					{
						for (j = 1; j <= 8 * ndn; j++)
						{
							RT[j][i] = R[i][j];
						}
					}

					/* Elemental Conduction Matrix */
					tranmul(R, KR, kel, 3, 8*ndn, 8*ndn);



					/*for (i = 1; i <= 8*ndn; i++)
					{
					for (j = 1; j <= 8*ndn; j++)
					{

					kel[i][j] = R[1][i] * KR[1][j]+ R[2][i] * KR[2][j] + R[3][i] * KR[3][j];

					}
					}*/
					

					/* Element conduction stiffness matrix*/
					for (i = 1; i <= 8 * ndn; i++)
					for (j = 1; j <= 8 * ndn; j++)
					{
						ke[i][j] += kel[i][j] * det*wxi[ox][ix] * weta[oe][ie] * wzeta[oz][iz];
					}

					/* Equivalent thickness for each element*/
					th = 0;
					for (nnl = 1; nnl <= 8; nnl++)
					{
						nng = *(lmc + 8 * (iel - 1) + nnl);
						th += t[nng];
					}
					th = th / 8;

					/* Heat generation matrix*/

					if (hgc[iel] != 0)  /* Heat generation is present if hgc>0 */
						c = hg;
					else
						c = 0.0;
					for (i = 1; i <= 8 * ndn; i++)
					{
						qgc[i] += c*NN[i] * det*wxi[ox][ix] * weta[oe][ie] * wzeta[oz][iz];  /* Load vector due to heat generation*/
					}

					/* Convection load matrix*/
					for (nnl = 1; nnl <= 8; nnl++)
					{
						nng = *(lmc + 8 * (iel - 1) + nnl);

						for (i = 1; i <= ndn; i++)
						{
							xi0 = XIn[nnl] * xi[ox][ix];
							eta0 = ETAn[nnl] * eta[oe][ie];
							if ((nnl == 1) || (nnl == 3) || (nnl == 5) || (nnl == 7))
							{
								nn[nnl][i] = 0.25*(1 + xi0)*(1 + eta0)*(xi0 + eta0 - 1.0)*pow(1, i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);

							}
							else if ((nnl == 2) || (nnl == 6))
							{
								nn[nnl][i] = 0.5*(1.0 - xi[ox][ix] * xi[ox][ix])*(1.0 + eta0)*pow(1, i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);

							}
							else if ((nnl == 4) || (nnl == 8))
							{
								nn[nnl][i] = 0.5*(1.0 - eta[oe][ie] * eta[oe][ie])*(1.0 + xi0)*pow(1, i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);

							}
						}
					}

					for (nnl = 1; nnl <= 8; nnl++)
					{
						for (i = 1; i <= ndn; i++)
						{
							NN[(nnl - 1)*ndn + i] = nn[nnl][i];
						}
					}

					if (hsc[iel] != 0)  
						c = hs;
					else
						c = 0.0;
					for (i = 1; i <= 8 * ndn; i++)
					{
						qs[i] += c*ts*NN[i] * det*wxi[ox][ix] * weta[oe][ie];  /* Load vector due to convective heat transfer*/
					}

					/* Elemental Convective stiffness matrix*/
					if (hsc[iel] != 0)
						c = hs;
					else
						c = 0.0;
					for (i = 1; i <= 8 * ndn; i++)
					{
						for (j = 1; j <= 8 * ndn; j++)
						{

							kh[i][j] += c * det * NN[i] * NN[j] * wxi[ox][ix] * weta[oe][ie];
						}
					}

					/* External Heat flux*/
					for (nnl = 1; nnl <= 8; nnl++)
					{
						nng = *(lmc + 8 * (iel - 1) + nnl);

						for (i = 1; i <= ndn; i++)
						{
							xi0 = XIn[nnl] * xi[ox][ix];
							eta0 = ETAn[nnl] * eta[oe][ie];
							if ((nnl == 1) || (nnl == 3) || (nnl == 5) || (nnl == 7))
							{
								nn[nnl][i] = 0.25*(1 + xi0)*(1 + eta0)*(xi0 + eta0 - 1.0)*pow(-1, i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);

							}
							else if ((nnl == 2) || (nnl == 6))
							{
								nn[nnl][i] = 0.5*(1.0 - xi[ox][ix] * xi[ox][ix])*(1.0 + eta0)*pow(-1, i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);

							}
							else if ((nnl == 4) || (nnl == 8))
							{
								nn[nnl][i] = 0.5*(1.0 - eta[oe][ie] * eta[oe][ie])*(1.0 + xi0)*pow(-1, i - 1)*pow(t[nng] / 2, i - 1) / factorial(i - 1);

							}
						}
					}

					for (nnl = 1; nnl <= 8; nnl++)
					{
						for (i = 1; i <= ndn; i++)
						{
							NN[(nnl - 1)*ndn + i] = nn[nnl][i];
						}
					}

					for (i = 1; i <= 8 * ndn; i++)
					{
						qq[i] += hf[iel] * NN[i] * det*wxi[ox][ix] * weta[oe][ie];  /* Load vector due to external heat flux*/
					}


				}
			}
		}   /*Gauss integration loop ends*/

		fprintf(fst, "\n\n conduction matrix");
		for (i = 1; i <= 8 * ndn; i++)
		{
			fprintf(fst, "\n");
			for (j = 1; j <= 8 * ndn; j++)
			{
				fprintf(fst, "  %5.3lf  ", ke[i][j]);
			}
		}

		fprintf(fst, "\n\n Conv. matrix");
		for (i = 1; i <= 8 * ndn; i++)
		{
			fprintf(fst, "\n");
			for (j = 1; j <= 8 * ndn; j++)
			{
				fprintf(fst, "  %5.5lf  ", kh[i][j]);
			}
		}

		for (i = 1; i <= 8 * ndn; i++)
		{
			for (j = 1; j <= 8 * ndn; j++)
			{
				ke[i][j] += kh[i][j];
			}
		}

		fprintf(fst, "\n\n Element stiffness matrix");
		for (i = 1; i <= 8 * ndn; i++)
		{
			fprintf(fst, "\n");
			for (j = 1; j <= 8 * ndn; j++)
			{
				fprintf(fst, "  %5.2lf  ", ke[i][j]);
			}
		}
		fprintf(fst, "\n\n Element load vectors");
		fprintf(fst, "\n qgc\t    qs\t    qq\t    qe  ");
		for (i = 1; i <= 8 * ndn; i++)
		{
			fprintf(fst, "\n");
			fprintf(fst, "  %5.4lf %5.4lf %5.4lf  %5.4lf", qgc[i], qs[i], qq[i], qgc[i] + qs[i] - qq[i]);

		}
		
		/*if (iel >= 11)goto loop;
		/* Assembling global stiffness matrix*/
		/*for (ii = 0; ii < 8; ii++) {
			nrt = ndn*(*(lmc + 8 * (iel - 1) + ii + 1) - 1);
			for (it = 0; it < ndn; it++) {
				nr = nrt + it;
				i = ndn * ii + it;
				for (jj = 0; jj < 8; jj++) {
					nct = ndn * (*(lmc + 8 * (iel - 1) + jj + 1) - 1);
					for (jt = 0; jt < ndn; jt++) {
						j = ndn * jj + jt;
						nc = nct + jt - nr;
						if (nc >= 0)
							s[nbw*nr + nc] = s[nbw*nr + nc] + ke[i][j];
					}
				}
				f[nr] = f[nr] + qgc[i] + qs[i] + qq[i];
			}
		}*/


		/* New */
		for (ii = 1; ii <=8; ii++) {
			nrt = ndn*(*(lmc + 8 * (iel - 1) + ii ) - 1);
			for (it = 1; it <=ndn; it++) {
				nr = nrt + it;
				i = ndn * (ii-1) + it;
				for (jj = 1; jj <=8; jj++) {
					nct = ndn * (*(lmc + 8 * (iel - 1) + jj ) - 1);
					for (jt = 1; jt <= ndn; jt++) {
						j = ndn * (jj-1) + jt;
						nc = nct + jt - nr;
						if (nc >= 0)
							s[nbw*(nr-1) + nc] = s[nbw*(nr-1) + nc] + ke[i][j];
					}
				}
				f[nr-1] = f[nr-1] + qgc[i] + qs[i] - qq[i];
			}
		}


		/*for (ii = 0; ii < nen; ii++) {
			nrt = ndn * (noc[nen*n + ii] - 1);
			for (it = 0; it < ndn; it++) {
				nr = nrt + it;
				i = ndn * ii + it;
				for (jj = 0; jj < nen; jj++) {
					nct = ndn * (noc[nen*n + jj] - 1);
					for (jt = 0; jt < ndn; jt++) {
						j = ndn * jj + jt;
						nc = nct + jt - nr;
						if (nc >= 0)
							s[nbw*nr + nc] = s[nbw*nr + nc] + se[i][j];
					}
				}
				f[nr] = f[nr] + tl[i];
			}
		}*/

		count++;
		fprintf(fout, "\n%d\t%d", count,iel);
		
	/*loop:
		;*/
	}   /* Element loop ends*/

	fprintf(flm, "\n\n Global Stiffness Matrix\n");   /* Print Assembled global stiffness matrix*/
	for (i = 0; i < nsn*ndn; i++)
	{
		for (j = 0; j < nbw; j++)
		{
			fprintf(flm, " %.3lf ", s[nbw*i + j]);
		}
		fprintf(flm, "\n");
	}

	fprintf(flm, "\n\nGlobal load vector\n"); /*Global load vector*/
	for (i = 0; i < nsn*ndn; i++)
	{
		fprintf(flm, "%.5lf", f[i]);
		fprintf(flm, "\n");
	}

	/* Determining Temperature vector order*/
	/*int *lm,cn;
	lm = ialloc1(nsn + 1);
	lm[1] = *(lmc);
	for (iel = 0; iel < nel; iel++)
	{
	for (nnl = 0; nnl < 8; nnl++)
	{
	cn = *(lmc + 8 * iel + nnl);
	for (i = 1; i <= j; i++)
	{
	if (lm[i]!=cn)

	}
	}
	}*/

	/* Modifying Stiffness and  load for temperature boundary condition*/
	double pe = 0.0;
	for (k = 0; k < nsn*ndn; k++)
	{
		pe += s[nbw*(k)] / (nsn*ndn);
	}
	pe = pe*10000000;
	
		for (i = 1; i <= nnt; i++)
		{
			
				s[nbw*(nt[i] - 1)*ndn] += pe;
				f[(nt[i] - 1)*ndn] +=  pe*T[i];

			
		}
		

		for (i = 1; i <= dnnt; i++)
		{

			s[nbw*(dnt[i]*ndn-1)] += pe;
			f[(dnt[i]*ndn-1)] += pe*dT[i];


		}

		fprintf(flm, "\n\n Global Stiffness Matrix\n");   /* Print Assembled global stiffness matrix*/
		for (i = 0; i < nsn*ndn; i++)
		{
			for (j = 0; j < nbw; j++)
			{
				fprintf(flm, " %.3lf ", s[nbw*i + j]);
			}
			fprintf(flm, "\n");
		}

		fprintf(flm, "\n\nGlobal load vector\n"); /*Global load vector*/
		for (i = 0; i < nsn*ndn; i++)
		{
			fprintf(flm, "%.5lf", f[i]);
			fprintf(flm, "\n");
		}

		
		/* ----- band solver ----- */
		int n1, nk, i1, j1, kk, nq;
		double c1;
		nq = nsn*ndn;
		n1 = nq - 1;

		/* --- forward elimination --- */
		for (k = 1; k <= n1; k++) {
			nk = nq - k + 1;
			if (nk > nbw)
				nk = nbw;
			for (i = 2; i <= nk; i++) {
				c1 = s[nbw*(k - 1) + i - 1] / s[nbw*(k - 1)];
				i1 = k + i - 1;
				for (j = i; j <= nk; j++) {
					j1 = j - i + 1;
					s[nbw*(i1 - 1) + j1 - 1] = s[nbw*(i1 - 1) + j1 - 1] - c1 * s[nbw*(k - 1) + j - 1];
				}
				f[i1 - 1] = f[i1 - 1] - c1 * f[k - 1];
			}
		}

		/*s[nbw*(nq - 1)] = 0.01;*/
		/* --- back-substitution --- */
		f[nq - 1] = f[nq - 1] / s[nbw*(nq - 1)];
		for (kk = 1; kk <= n1; kk++) {
			k = nq - kk;
			c1 = 1 / s[nbw*(k - 1)];
			f[k - 1] = c1 * f[k - 1];
			nk = nq - k + 1;
			if (nk > nbw)
				nk = nbw;
			for (j = 2; j <= nk; j++) {
				f[k - 1] = f[k - 1] - c1 * s[nbw*(k - 1) + j - 1] * f[k + j - 2];
			}
		}

	
		fprintf(flm, "\n\n");
		for (i = 0; i < nsn*ndn; i++)
		{
			fprintf(flm, "\n%d\t%.2lf", i+1,f[i]);

		}

		/* Temperature at mid surface*/

		fprintf(flm, "\nTemperature at mid surface\n");
		for (i = 0; i < nsn; i++)
		{
			j = i*ndn;
			fprintf(flm, "\n%d\t%.2lf", i+1, f[j]);
		}

		/* Temperature at top*/
		fprintf(flm, "\nTemperature at Top surface\n");
		for (i = 0; i <nsn; i++)
		{
			j = i*ndn;
			k = j + 1;
			Tx = f[j] + 1 * (0.2)*f[k];
			fprintf(flm, "\n%d\t%.2lf", i + 1, Tx);

		}
		
		/* Temperature at Bottom*/
		fprintf(flm, "\nTemperature at Bottom surface\n");
		for (i = 0; i <nsn; i++)
		{
			j = i*ndn;
			k = j + 1;
			Tx = f[j] - 1 *(0.2)*f[k];
			fprintf(flm, "\n%d\t%.2lf", i + 1, Tx);

		}

		/* Temperature plot*/
		fprintf(flm, "\n\nTemperature plot- Mid");
		for (i = 0; i <nel; i++)
		{
			nng = *(lmc + 8 * (i) + 6);
			fprintf(flm, "\n%d\t%0.3lf", nng, f[(nng-1)*ndn]);
			
		}

		fprintf(flm, "\n\nTemperature plot- Top");
		for (i = 0; i <nel; i++)
		{
			nng = *(lmc + 8 * (i)+6);
			Tx = f[(nng-1)*ndn] + 1 * (0.2)*f[nng*ndn-1];
			fprintf(flm, "\n%d\t%0.3lf", nng, Tx);

		}

		fprintf(flm, "\n\nTemperature plot- bottom");
		for (i = 0; i <nel; i++)
		{
			nng = *(lmc + 8 * (i)+6);
			Tx = f[(nng - 1)*ndn] - 1 * (0.2)*f[nng*ndn - 1];
			fprintf(flm, "\n%d\t%0.3lf", nng, Tx);

		}


		fclose(fin);
		fclose(fout);
		fclose(fst);
		fclose(flm);
}

int factorial(int it)
{
	int fact=1, zz;
	if (it > 1) 
	{
		for(zz = 1; zz <= it; zz++)
			fact = fact*zz;
	}
	return(fact);

}