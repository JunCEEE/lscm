// LSCM.c -- Large Single Crystal Maker
// gcc LSCM.c -o LSCM
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


void get_curr_time()
{
	time_t curr_time;
	struct tm * timeinfo;
	time(&curr_time);
	timeinfo = localtime(&curr_time);
	printf("%s", asctime(timeinfo));
}

void NaCl(char* output, int* nc, double a)
{
	FILE *fout = fopen(output,"w");
	printf("NaCl\n");
	double x[8] = {0.0,0.5,0.5,0.0,0.5,0.0,0.0,0.5};
	double y[8] = {0.0,0.5,0.0,0.5,0.5,0.0,0.5,0.0};
	double z[8] = {0.0,0.0,0.5,0.5,0.5,0.5,0.0,0.0};
	int symbol[8] = {1,1,1,1,2,2,2,2}; 
	unsigned long long na = (unsigned long long)  8*nc[0]*nc[1]*nc[2];
	printf("%llu atoms\n",na);
	//cout << nc[0] << " " << nc[1] << " " << nc[2] << " " << endl; 
	//cout << "Atom number = " << na << endl; 
	fprintf(fout,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS:\n");
	fprintf(fout,"%llu\n",na);
	fprintf(fout,"ITEM: BOX BOUNDS pp pp pp\n"); 
	for (int i = 0; i < 3; i++)
		fprintf(fout,"0 %f\n",(nc[i])*a);
	fprintf(fout,"ITEM: ATOMS type x y z\n");
	
	double xx, yy, zz;
	for (int i = 0; i < nc[0]; i++)
		for (int j = 0; j < nc[1]; j++)
			for (int k = 0; k < nc[2]; k++)
			{
				for (int l = 0; l < 8; l++)
				{
					xx = (i + x[l])*a;
					yy = (j + y[l])*a;
					zz = (k + z[l])*a;
					fprintf(fout,"%i %f %f %f\n",symbol[l],xx,yy,zz);
				}
			}
}

void FCC(char* output, int* nc, double a)
{
	FILE *fout = fopen(output,"w");
	printf("FCC\n");
	double x[4] = {0.0,0.5,0.5,0.0};
	double y[4] = {0.0,0.5,0.0,0.5};
	double z[4] = {0.0,0.0,0.5,0.5};
	int symbol[4] = {1,1,1,1}; 
	unsigned long long na = (unsigned long long)  4*nc[0]*nc[1]*nc[2];
	printf("%llu atoms\n",na);
	//cout << nc[0] << " " << nc[1] << " " << nc[2] << " " << endl; 
	//cout << "Atom number = " << na << endl; 
	fprintf(fout,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS:\n");
	fprintf(fout,"%llu\n",na);
	fprintf(fout,"ITEM: BOX BOUNDS pp pp pp\n"); 
	for (int i = 0; i < 3; i++)
		fprintf(fout,"0 %f\n",(nc[i])*a);
	fprintf(fout,"ITEM: ATOMS type x y z\n");
	
	double xx, yy, zz;
	for (int i = 0; i < nc[0]; i++)
		for (int j = 0; j < nc[1]; j++)
			for (int k = 0; k < nc[2]; k++)
			{
				for (int l = 0; l < 4; l++)
				{
					xx = (i + x[l])*a;
					yy = (j + y[l])*a;
					zz = (k + z[l])*a;
					fprintf(fout,"%i %f %f %f\n",symbol[l],xx,yy,zz);
				}
			}
}

void BCC(char* output, int* nc, double a)
{
	FILE *fout = fopen(output,"w");
	printf("BCC\n");
	double x[2] = {0.0,0.5};
	double y[2] = {0.0,0.5};
	double z[2] = {0.0,0.5};
	int symbol[2] = {1,1}; 
	unsigned long long na = (unsigned long long)  2*nc[0]*nc[1]*nc[2];
	printf("%llu atoms\n",na);
	//cout << nc[0] << " " << nc[1] << " " << nc[2] << " " << endl; 
	//cout << "Atom number = " << na << endl; 
	fprintf(fout,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS:\n");
	fprintf(fout,"%llu\n",na);
	fprintf(fout,"ITEM: BOX BOUNDS pp pp pp\n"); 
	for (int i = 0; i < 3; i++)
		fprintf(fout,"0 %f\n",(nc[i])*a);
	fprintf(fout,"ITEM: ATOMS type x y z\n");
	
	double xx, yy, zz;
	for (int i = 0; i < nc[0]; i++)
		for (int j = 0; j < nc[1]; j++)
			for (int k = 0; k < nc[2]; k++)
			{
				for (int l = 0; l < 2; l++)
				{
					xx = (i + x[l])*a;
					yy = (j + y[l])*a;
					zz = (k + z[l])*a;
					fprintf(fout,"%i %f %f %f\n",symbol[l],xx,yy,zz);
				}
			}
}

void displayUsage(char* bin)
{
	printf("Usage: %s nx ny nz a [structure] -o output\n", bin);
	printf("Example: %s 50 50 50 3.615 -f -o singleCu.custom\n",bin);
	puts( "    -f, FCC " );
	puts( "    -b, BCC " );
	puts( "    -n, Nacl " );
	puts( "    -o, outputfile" );
	puts( "    -h, print help" );
}

struct args_t
{
	char* output; // -o output filename
	int mode;
}args={NULL,0};

static const char *optString = "fbno:h?";


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
    		case 'o': args.output = optarg; break;
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
