#include "embed.h"
#include "navier-stokes/centered.h"
#include "vof.h"
// add your parameter here
double rho1 = 1., mu1 = 0.01, rho2 = 0.85, mu2 = 0.255, radius = 0.005*0.95;

#include "tension_palas.h"
//#include "navier-stokes/perfs.h" //  to check performance in parallel -- Palas
#include "output_vtu_foreach.h" // new file which works for adaptive grid (vtk only works for uniform
#include <sys/stat.h> //  to add directories -- Palas

scalar f[], * interfaces = {f};
face vector alphav[], muv[], av[];
scalar rhov[];


int main(int argc, char * argv[]){
	// Create directories to avoid seg faults -- Palas
	struct stat st = {0};
	char name[80];
	sprintf (name, "data");
	run();
}

//INITIAL COLOR FUNCTION FIELD//
event init (t = 0){
	char filename[80];
	int pt = 0;
	
	while(1){
		//sprintf (filename, "data_%d/dump_%g", resolution,pt*0.01);
	  sprintf (filename, "snap-%d",pt);

		printf("extracting %s\n",filename);
		fflush(stdout);

		if (!restore (file = filename)){

			printf("files ended");
			exit(0);
		}
		
		scalar *list;
		bool linear;
		char name2[80];
		sprintf(name2, "vtufile-%d.vtu",pt);
		FILE *ptr2 = fopen(name2,"w");

		//~ output_vtu_bin_foreach ((scalar *) {f,p,cm}, (vector *) {u}, N, ptr2, false); // cm for embed boundary
		output_vtu_bin_foreach ((scalar *) {f,p}, (vector *) {u}, N, ptr2, false);
		fclose(ptr2);

		pt = pt + 1;
	}
}
