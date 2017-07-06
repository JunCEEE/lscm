// LSCM.c -- Large Single Crystal Maker
// gcc LSCM.c -std=c99 -o LSCM
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>


struct args_t
{
	char* output; // -o output filename
	int mode;
	int is_s;
	double diameter;
}args={NULL,0,0,0};

static const char *optString = "fbndo:s:h?";

void get_curr_time()
{
	time_t curr_time;
	struct tm * timeinfo;
	time(&curr_time);
	timeinfo = localtime(&curr_time);
	printf("%s", asctime(timeinfo));
}

void writeAtoms(char* output, int *nc, double a, int una, double *x, double *y, double *z, int *symbol)
{
	unsigned long long na = (unsigned long long)  una*nc[0]*nc[1]*nc[2];
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
		fprintf(fout,"0 %f\n",(nc[i])*a);
		center[i] = nc[i]*a/2;
	}
	fprintf(fout,"ITEM: ATOMS type x y z\n");

	
	unsigned long long dna = 0;
	double xx, yy, zz;
	for (int i = 0; i < nc[0]; i++)
		for (int j = 0; j < nc[1]; j++)
			for (int k = 0; k < nc[2]; k++)
			{
				for (int l = 0; l < una; l++)
				{
					xx = (i + x[l])*a;
					yy = (j + y[l])*a;
					zz = (k + z[l])*a;
					if (args.is_s && ( (xx-center[0])*(xx-center[0]) + (yy-center[1])*(yy-center[1]) + (zz-center[2])*(zz-center[2]) > args.diameter*args.diameter/4 ) ) 
					{
						dna++;
						continue;
					}
					else 
					{
						fprintf(fout,"%i %f %f %f\n",symbol[l],xx,yy,zz);
					}
				}
			}
	
	fclose(fout);
	if (args.is_s)
	{
		printf("%llu atoms deleted, %llu atoms left.\n",dna,na-dna);
		char buffer[255];
		sprintf(buffer,"sed -i '4c%llu' %s\n",na-dna,output);
		printf("%s",buffer);
		system(buffer);
	}
}

void NaCl(char* output, int* nc, double a)
{
	printf("NaCl\n");
	int una = 8;
	double x[8] = {0.0,0.5,0.5,0.0,0.5,0.0,0.0,0.5};
	double y[8] = {0.0,0.5,0.0,0.5,0.5,0.0,0.5,0.0};
	double z[8] = {0.0,0.0,0.5,0.5,0.5,0.5,0.0,0.0};
	int symbol[8] = {1,1,1,1,2,2,2,2}; 
	writeAtoms(output,nc,a,una,x,y,z,symbol);
}



void FCC(char* output, int* nc, double a)
{
	printf("FCC\n");
	int una = 4;
	double x[4] = {0.0,0.5,0.5,0.0};
	double y[4] = {0.0,0.5,0.0,0.5};
	double z[4] = {0.0,0.0,0.5,0.5};
	int symbol[4] = {1,1,1,1};
	writeAtoms(output,nc,a,una,x,y,z,symbol);
}

void BCC(char* output, int* nc, double a)
{
	printf("BCC\n");
	int una = 2;
	double x[2] = {0.0,0.5};
	double y[2] = {0.0,0.5};
	double z[2] = {0.0,0.5};
	int symbol[2] = {1,1};
	writeAtoms(output,nc,a,una,x,y,z,symbol);
}

void Diamond(char* output, int* nc, double a)
{
	printf("Diamond\n");
	int una = 8;
	double x[8] = {0.00,0.50,0.50,0.00,0.25,0.75,0.75,0.25};
	double y[8] = {0.00,0.50,0.00,0.50,0.25,0.75,0.25,0.75};
	double z[8] = {0.00,0.00,0.50,0.50,0.25,0.25,0.75,0.75};
	int symbol[8] = {1,1,1,1,1,1,1,1};
	writeAtoms(output,nc,a,una,x,y,z,symbol);
}


void displayUsage(char* bin)
{
	printf("Usage: %s nx ny nz a [structure] -o output\n", bin);
	printf("Example: %s 50 50 50 3.615 -f -o singleCu.custom\n",bin);
	puts( "    -s, Spherical shape " );
	puts( "    -f, FCC " );
	puts( "    -b, BCC " );
	puts( "    -d, Diamond " );
	puts( "    -n, Nacl " );
	puts( "    -o, outputfile" );
	puts( "    -h, print help" );
}



int main(int argc, char* argv[])
{
	if (argc < 8)
	{
		displayUsage(argv[0]);
		return 0;
	}

	int nc[3] = {atoi(argv[1]),atoi(argv[2]),atoi(argv[3])};
	double a = atof(argv[4]);

	int opt;
    while ((opt = getopt (argc, argv, optString)) != -1)
    {
    	switch(opt)
    	{
    		case 'f': args.mode = 1; break;
    		case 'b': args.mode = 2; break;
    		case 'n': args.mode = 3; break;
    		case 'd': args.mode = 4; break;
    		case 'o': args.output = optarg; break;
    		case 's': args.is_s = 1; args.diameter = atof(optarg); break;
			case 'h': // just go to '?'
			case '?': 
			default:
    			{
	    			displayUsage(argv[0]); 
	    		}
				return(0);
    			break;
    	}
    }

	get_curr_time();
	clock_t begin = clock();

	switch(args.mode)
	{
		case 1 : FCC(args.output,nc,a); break;
		case 2 : BCC(args.output,nc,a); break;
		case 3 : NaCl(args.output,nc,a); break;
		case 4 : Diamond(args.output,nc,a); break;
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
