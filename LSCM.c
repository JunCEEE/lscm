// LSCM.c -- Large Single Crystal Maker
// gcc LSCM.c -std=c99 -o LSCM
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <math.h>


struct args_t
{
	char* output; // -o output filename
	char* input; // -r input filename
	int outMode; // output mode --cfg --custom
	double a;
	double b;
	double c;
	int mode;
	int is_s;
	double diameter;
	int nc[3];
}args={NULL,NULL,0,-1,-1,-1};

static const char *optString = "p:a:x:y:z:r:fbcndo:s:h?";


static const struct option longOpts[] = {
    { "custom", no_argument, NULL,0},
    { "cfg", no_argument, NULL,1},
    { "data", no_argument, NULL,5},
    { "pa", required_argument, NULL,2},
    { "pb", required_argument, NULL,3},
    { "pc", required_argument, NULL,4},
};

void get_curr_time()
{
	time_t curr_time;
	struct tm * timeinfo;
	time(&curr_time);
	timeinfo = localtime(&curr_time);
	printf("%s", asctime(timeinfo));
}

void writeAtoms(char* output, int *nc, double *lc, int una, double *x, double *y, double *z, int *symbol, int ne)
{
	unsigned long long na = (unsigned long long)  una*nc[0]*nc[1]*nc[2];
	FILE *fout = fopen(output,"w");
	double box[3];


	double center[3];
	if (args.outMode == 0)
	{
		printf("%llu atoms\n",na);
	//cout << nc[0] << " " << nc[1] << " " << nc[2] << " " << endl; 
	//cout << "Atom number = " << na << endl; 
		fprintf(fout,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS:\n");
		fprintf(fout,"%llu\n",na);
		fprintf(fout,"ITEM: BOX BOUNDS pp pp pp\n"); 

		double center[3];
		for (int i = 0; i < 3; i++)
		{
			fprintf(fout,"0 %f\n",(nc[i])*lc[i]);
			center[i] = nc[i]*lc[i]/2;
		}
		fprintf(fout,"ITEM: ATOMS type x y z\n");
	}
	else if (args.outMode == 1)
	{
		for (int i = 0; i < 3; i++)
		{
			box[i] = nc[i]*lc[i];
			center[i] = box[i]/2;
		}
		fprintf(fout,"Number of particles = %llu\n",na);
		fprintf(fout,"A = 1 Angstrom (basic length-scale)\n"); 
		fprintf(fout,"H0(1,1) = %f A\n",box[0]);
		fprintf(fout,"H0(1,2) = 0 A\n"); 
		fprintf(fout,"H0(1,3) = 0 A\n"); 
		fprintf(fout,"H0(2,1) = 0 A\n"); 
		fprintf(fout,"H0(2,2) = %f A\n",box[1]); 
		fprintf(fout,"H0(2,3) = 0 A\n"); 
		fprintf(fout,"H0(3,1) = 0 A\n"); 
		fprintf(fout,"H0(3,2) = 0 A\n"); 
		fprintf(fout,"H0(3,3) = %f A\n",box[2]); 
		fprintf(fout,".NO_VELOCITY.\nentry_count = 3\n");
		fprintf(fout,"63.55\nCu\n");
	}
	else if (args.outMode == 2)
	{
		printf("%llu atoms\n",na);
	//cout << nc[0] << " " << nc[1] << " " << nc[2] << " " << endl; 
	//cout << "Atom number = " << na << endl; 
		fprintf(fout,"# LAMMPS data file written by LSCM\n");
		fprintf(fout,"%llu atoms\n",na);
		fprintf(fout,"%i atom types\n",ne); 

		double center[3];
		for (int i = 0; i < 3; i++)
		{
			center[i] = nc[i]*lc[i]/2;
		}
		fprintf(fout,"0.0 %f xlo xhi\n",(nc[0])*lc[0]);
		fprintf(fout,"0.0 %f ylo yhi\n",(nc[1])*lc[1]);
		fprintf(fout,"0.0 %f zlo zhi\n",(nc[2])*lc[2]);

		fprintf(fout,"\n");

		fprintf(fout,"Atoms # atomic\n");

		fprintf(fout,"\n");
	}

	
	unsigned long long dna = 0;
	double xx, yy, zz;
	int n = 0;
	for (int i = 0; i < nc[0]; i++)
		for (int j = 0; j < nc[1]; j++)
			for (int k = 0; k < nc[2]; k++)
			{
				for (int l = 0; l < una; l++)
				{
					xx = (i + x[l])*lc[0];
					yy = (j + y[l])*lc[1];
					zz = (k + z[l])*lc[2];
					if (args.is_s && ( (xx-center[0])*(xx-center[0]) + (yy-center[1])*(yy-center[1]) + (zz-center[2])*(zz-center[2]) > args.diameter*args.diameter/4 ) ) 
					{
						dna++;
						continue;
					}
					else 
					{
						if (args.outMode == 0)
						{
							fprintf(fout,"%i %f %f %f\n",symbol[l],xx,yy,zz);
						}
						else if (args.outMode == 1)
						{
							fprintf(fout,"%f %f %f\n",xx/box[0],yy/box[1],zz/box[2]);
						}
						else if (args.outMode == 2)
						{
							n++;
							fprintf(fout,"%i %i %f %f %f\n",n,symbol[l],xx,yy,zz);
						}
					}
				}
			}
	
	fclose(fout);
	if (args.is_s)
	{
		printf("%llu atoms deleted, %llu atoms left.\n",dna,na-dna);
		char buffer[255];
		if (args.outMode == 0)
			sprintf(buffer,"sed -i '4c%llu' %s\n",na-dna,output);
		else if (args.outMode == 1)
			sprintf(buffer,"sed -i '2cNumber of particles = %llu' %s\n",na-dna,output);
		printf("%s",buffer);
		system(buffer);
	}
}

// void writeAtomsTri(char* output, int *nc, double *lc, int una, double *x, double *y, double *z, int *symbol)
// {
// 	unsigned long long na = (unsigned long long)  una*nc[0]*nc[1]*nc[2];
// 	FILE *fout = fopen(output,"w");
// 	printf("%llu atoms\n",na);
// 	//cout << nc[0] << " " << nc[1] << " " << nc[2] << " " << endl; 
// 	//cout << "Atom number = " << na << endl; 
// 	fprintf(fout,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS:\n");
// 	fprintf(fout,"%llu\n",na);
// 	fprintf(fout,"ITEM: BOX BOUNDS pp pp pp\n"); 
	
// 	double center[3];
// 	for (int i = 0; i < 3; i++)
// 	{
// 		fprintf(fout,"0 %f\n",(nc[i])*lc[i]);
// 		center[i] = nc[i]*lc[i]/2;
// 	}


// 	fprintf(fout,"ITEM: ATOMS type x y z\n");

	
// 	unsigned long long dna = 0;
// 	double xx, yy, zz;
// 	for (int i = 0; i < nc[0]; i++)
// 		for (int j = 0; j < nc[1]; j++)
// 			for (int k = 0; k < nc[2]; k++)
// 			{
// 				for (int l = 0; l < una; l++)
// 				{
// 					xx = (i + x[l])*lc[0];
// 					yy = (j + y[l])*lc[1];
// 					zz = (k + z[l])*lc[2];
// 					if (args.is_s && ( (xx-center[0])*(xx-center[0]) + (yy-center[1])*(yy-center[1]) + (zz-center[2])*(zz-center[2]) > args.diameter*args.diameter/4 ) ) 
// 					{
// 						dna++;
// 						continue;
// 					}
// 					else 
// 					{
// 						fprintf(fout,"%i %f %f %f\n",symbol[l],xx,yy,zz);
// 					}
// 				}
// 			}
	
// 	fclose(fout);
// 	if (args.is_s)
// 	{
// 		printf("%llu atoms deleted, %llu atoms left.\n",dna,na-dna);
// 		char buffer[255];
// 		sprintf(buffer,"sed -i '4c%llu' %s\n",na-dna,output);
// 		printf("%s",buffer);
// 		system(buffer);
// 	}
// }

void NaCl(char* output, int* nc, double a)
{
	printf("NaCl\n");
	int una = 8;
	int ne = 2;
	double x[8] = {0.0,0.5,0.5,0.0,0.5,0.0,0.0,0.5};
	double y[8] = {0.0,0.5,0.0,0.5,0.5,0.0,0.5,0.0};
	double z[8] = {0.0,0.0,0.5,0.5,0.5,0.5,0.0,0.0};
	int symbol[8] = {1,1,1,1,2,2,2,2}; 
	double lc[3] = {a,a,a};
	writeAtoms(output,nc,lc,una,x,y,z,symbol,ne);
}


void SP(char* output, int* nc, double a,double b, double c)
{
	printf("SP: Simple Cubic\n");
	int una = 1;
	int ne = 1;
	double x[1] = {0.0};
	double y[1] = {0.0};
	double z[1] = {0.0};
	int symbol[1] = {1}; 
	double lc[3] = {a,b,c};
	writeAtoms(output,nc,lc,una,x,y,z,symbol,ne);
}

void FCC(char* output, int* nc, double a)
{
	printf("FCC\n");
	int una = 4;
	int ne = 1;
	double x[4] = {0.0,0.5,0.5,0.0};
	double y[4] = {0.0,0.5,0.0,0.5};
	double z[4] = {0.0,0.0,0.5,0.5};
	int symbol[4] = {1,1,1,1};
	double lc[3] = {a,a,a};
	writeAtoms(output,nc,lc,una,x,y,z,symbol,ne);
}

void BCC(char* output, int* nc, double a)
{
	printf("BCC\n");
	int una = 2;
	int ne = 1;
	double x[2] = {0.0,0.5};
	double y[2] = {0.0,0.5};
	double z[2] = {0.0,0.5};
	int symbol[2] = {1,1};
	double lc[3] = {a,a,a};
	writeAtoms(output,nc,lc,una,x,y,z,symbol,ne);
}

void HCP(double c, char* output, int* nc, double a)
{
	printf("HCP\n");
	printf("a =  b = %f\n",a);
	printf("c = %f\n",c);
	double lc[3];
	lc[0] = a * sqrt(3);
	lc[1] = a;
	lc[2] = c;
	int una = 4;
	int ne = 1;
	double x[4] = {0.33333,0.16667,0.66667,0.833330};
	double y[4] = {0.0,0.5,0.0,0.5};
	double z[4] = {0.25,0.75,0.75,0.25};
	int symbol[4] = {1,1,1,1};
	writeAtoms(output,nc,lc,una,x,y,z,symbol,ne);
}

void Diamond(char* output, int* nc, double a)
{
	printf("Diamond\n");
	int una = 8;
	int ne = 1;
	double x[8] = {0.00,0.50,0.50,0.00,0.25,0.75,0.75,0.25};
	double y[8] = {0.00,0.50,0.00,0.50,0.25,0.75,0.25,0.75};
	double z[8] = {0.00,0.00,0.50,0.50,0.25,0.25,0.75,0.75};
	int symbol[8] = {1,1,1,1,1,1,1,1};
	double lc[3] = {a,a,a};
	writeAtoms(output,nc,lc,una,x,y,z,symbol,ne);
}

struct RR_data
{
	int una;
	int ne;
	double* x;
	double* y;
	double* z;
	int* symbol;
	double aa[3];
	double lo[3];
}xyz_data;

void read_data(char *input)
{
	char buffer[255];
	int n = 0;
	int na,ne;
	int r;
	double xlo,xhi,ylo,yhi,zlo,zhi;
	double aa[3];
    
    FILE * fin = fopen(input,"r");
    int skip_na = 0;
    int skip_aa = 0;
    // bool skip_box = 0;

	while (fgets(buffer,255,fin) != NULL )
	{
		if (!skip_na)
		{
			r = sscanf(buffer,"%i %s",&na,buffer);
			if (r == 2)
			{
				fgets(buffer,255,fin);
				sscanf(buffer,"%i %*s",&ne);
				skip_na = 1;
				continue;
			}
		}
		else if (!skip_aa)
		{
			r = sscanf(buffer,"%lf %lf %s %s",&xlo,&xhi,buffer,buffer);
			if (r == 4)
			{
				aa[0] = xhi-xlo;
				fgets(buffer,255,fin);
				sscanf(buffer,"%lf %lf",&ylo,&yhi);
				aa[1] = yhi-ylo;
				fgets(buffer,255,fin);
				sscanf(buffer,"%lf %lf",&zlo,&zhi);
				aa[2] = zhi-zlo;
				skip_aa = 1;
				continue;
			}
		}
	    else if (strncmp(buffer,"Atoms",5) == 0)
		{
			//fgets(buffer,255,fin);
			break;
		}
	}

	double *x = malloc(na * sizeof(double));
	double *y = malloc(na * sizeof(double));
	double *z = malloc(na * sizeof(double));
	int *symbol = malloc(na * sizeof(int));

	while (fgets(buffer,255,fin) != NULL )
	{
		//cout << buffer;
		if(buffer[0] == '\n' || buffer[0] == '\r' || buffer[0]== '\0') continue;  

		sscanf(buffer,"%*s %i""%lf" "%lf" "%lf",&symbol[n],&x[n],&y[n],&z[n]);
		if (symbol[n] > ne)
		{
			printf("Error: only %i element types defined, %i element types detected!\n",ne,symbol[n]);
			exit(0);
		}
		n++;
	}

	fclose(fin);

	if (n < na)
	{
		printf("Error: only %i atoms found. %i atoms declared in this file.\n",n,na);
		exit(0);
	}

	xyz_data.una = na;
	xyz_data.ne = ne;
	xyz_data.x = x;
	xyz_data.y = y;
	xyz_data.z = z;
	xyz_data.symbol = symbol;
	xyz_data.aa[0] = aa[0];
	xyz_data.aa[1] = aa[1];
	xyz_data.aa[2] = aa[2];
	xyz_data.lo[0] = xlo;
	xyz_data.lo[1] = ylo;
	xyz_data.lo[2] = zlo;
}



void RR(char* input, char* output, int* nc)
{
	printf("Read and replicate\n");
	read_data(input);

	unsigned long long na = (unsigned long long)  xyz_data.una*nc[0]*nc[1]*nc[2];
	FILE *fout = fopen(output,"w");
	printf("%llu atoms\n",na);
	//cout << nc[0] << " " << nc[1] << " " << nc[2] << " " << endl; 
	//cout << "Atom number = " << na << endl; 
	fprintf(fout,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS:\n");
	fprintf(fout,"%llu\n",na);
	fprintf(fout,"ITEM: BOX BOUNDS pp pp pp\n"); 
	
	double center[3];
	for (int i = 0; i < 3; i++)
	{
		fprintf(fout,"%f %f\n",xyz_data.lo[i],xyz_data.lo[i]+nc[i]*xyz_data.aa[i]);
	}
	fprintf(fout,"ITEM: ATOMS type x y z\n");

	
	unsigned long long dna = 0;
	double xx, yy, zz;
	for (int i = 0; i < nc[0]; i++)
		for (int j = 0; j < nc[1]; j++)
			for (int k = 0; k < nc[2]; k++)
			{
				for (int l = 0; l < xyz_data.una; l++)
				{
					xx = xyz_data.x[l] + i*xyz_data.aa[0];
					yy = xyz_data.y[l] + j*xyz_data.aa[1];
					zz = xyz_data.z[l] + k*xyz_data.aa[2];
					fprintf(fout,"%i %f %f %f\n",xyz_data.symbol[l],xx,yy,zz);
				}
			}
	
	fclose(fout);
}

void displayUsage(char* bin)
{
	puts("v1.2");
	printf("Usage: %s  [-p c]  [-f] [-b] [-d] [-n] [-a lattice_constant] [-r file.data] -x nx -y ny -z nz [--cfg] [--custom] [--data] -o output\n", bin);
	printf("Example: %s -f -a 3.615 -x 50 -y 50 -z 50  -o singleCu.custom\n",bin);
	printf("Example: %s -p 5.21033 -a 3.20927 -x 30 -y 50 -z 30 --cfg -o singleMg.cfg\n",bin);
	puts( "    -c, cP: simple cubic " );
	puts( "    -s, Spherical shape " );
	puts( "    -f, FCC " );
	puts( "    -b, BCC " );
	puts( "    -d, Diamond " );
	puts( "    -n, Nacl " );
	puts( "    -p, HCP " );
	puts( "    -r, Read and replicate LAMMPS DATA" );
	puts( "    -o, outputfile" );
	puts( "    -h, print help" );
	puts( "    --cfg, output in .cfg format" );
	puts( "    --custom, output in LAMMPS .custom format" );
	puts( "    --data, output in .data format" );
}



int main(int argc, char* argv[])
{

	int longIndex;
	int opt;
    while ((opt = getopt_long (argc, argv, optString, longOpts, &longIndex)) != -1)
    {
    	switch(opt)
    	{
    		case 'x': args.nc[0] = atof(optarg); break;
    		case 'y': args.nc[1] = atof(optarg); break;
    		case 'z': args.nc[2] = atof(optarg); break;
    		case 'a': args.a = atof(optarg); break;
    		case 'f': args.mode = 1; break;
    		case 'b': args.mode = 2; break;
    		case 'n': args.mode = 3; break;
    		case 'd': args.mode = 4; break;
    		case 'c': args.mode = 7; break;
    		case 'p': args.mode = 6;args.input = optarg; break;
    		case 'r': args.mode = 5;args.input = optarg; break;
    		case 'o': args.output = optarg;break;
    		case 's': args.is_s = 1; args.diameter = atof(optarg); break;
			case 'h': // just go to '?'
			case 0: args.outMode = 0; break; // custom
			case 1: args.outMode = 1; break; // cfg
			case 5: args.outMode = 2; break; // data
			case 2: args.a = atof(optarg); break; // Lattice parameter a
			case 3: args.b = atof(optarg); break; // Lattice parameter b
			case 4: args.c = atof(optarg); break; // Lattice parameter c
			case '?': 
			default:
    			{
	    			displayUsage(argv[0]); 
	    		}
				return(0);
    			break;
    	}
    }

    if (args.mode == 5 && args.nc[0] == 0 && args.nc[1] == 0 && args.nc[2] == 0)
    {
    	displayUsage(argv[0]);
    	return(0); 
    }

	get_curr_time();
	clock_t begin = clock();

	switch(args.mode)
	{
		case 1 : FCC(args.output,args.nc,args.a); break;
		case 2 : BCC(args.output,args.nc,args.a); break;
		case 3 : NaCl(args.output,args.nc,args.a); break;
		case 4 : Diamond(args.output,args.nc,args.a); break;
		case 7 : if(args.b == -1 || args.c == -1)  {SP(args.output,args.nc,args.a,args.a,args.a);} else {SP(args.output,args.nc,args.a,args.b,args.c);}      break;
		case 5 : RR(args.input,args.output,args.nc); break;
		case 6 : HCP(atof(args.input),args.output,args.nc,args.a); break;
		default:
			{
    			displayUsage(argv[0]); 
    		}
			return(0);
			break;
	}


	
    get_curr_time();
    clock_t end = clock();
    double elapsed_secs = (double) (end - begin) / CLOCKS_PER_SEC;
    printf("Time elapsed: %f s\n",elapsed_secs);
    return 0;
}
