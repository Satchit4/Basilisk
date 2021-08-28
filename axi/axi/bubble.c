#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
double radius = 0.008; // check this value again (CGS)
#include "tension.h"
#include "output_vtu_foreach.h"

int MAXLEVEL = 10;

// Inlet
u.n[right] = dirichlet(-0.1);
u.t[right] = dirichlet(0.);
p[right]   = neumann(0.);

//outlet
p[left]  = dirichlet(0.);
u.n[left] = neumann(0.);
u.t[left] = neumann(0.);

// No flow through bottom wall, reference, http://basilisk.fr/src/test/rising.c
uf.n[bottom] = 0.;

int main(int argc, char* argv[]) {
  size (0.08); // 800um channel
  origin(-0.06, 0.0); // origin is shiefted to centre of droplet.
	
  MAXLEVEL = atoi(argv[1]);
  N = 1 << MAXLEVEL;
	
	// CGS units, put your parameters... please follow one unit system
	rho2 = 1.,  mu2 = 0.001, rho1 = 0.001, mu1 = 0.00001;
	f.sigma = 1.6;
	
  TOLERANCE = 1e-4;
  run();
}

// to give boundary condition on masked region we are defining a variable, no-slip, no-penetration
// reference: http://basilisk.fr/sandbox/Antoonvh/bathtub.c

bid tub;
u.n[tub] = dirichlet(0.);
u.t[tub] = dirichlet(0.);
p[tub]   = neumann(0.);

event init (t = 0) {		
 // we are using mask function because axi and embed files are not compatible together.	
  mask((y>0.018 || (x<-0.02 && y>0.008)) ? tub : none);
	
 //define bubble shape	
  fraction (f, sqrt(sq(x)+sq(y)) < radius);
}

/*
// Uncomment this function if you want to add gravity in simulation. (check direction of gravity)
event acceleration (i++) {
  foreach_face(y)
    av.x[] = 981.0;
  boundary ((scalar *){av});
}
*/

event logfile (i++) {
  if (i == 0)
   fprintf (ferr,"t dt \n");
   fprintf (ferr, "%g %g \n", t, dt);
}

event snapshot (t = 0; t += 0.01; t <= 200.) {	
	
  scalar *list;		
  bool linear;
  char name2[80];
  static int pt = 0;
  sprintf(name2, "bessel-%d.vtu",pt);
  FILE *ptr2 = fopen(name2,"w");	
  output_vtu_bin_foreach ((scalar *) {f, p}, (vector *) {u}, N, ptr2, false);	
  fclose(ptr2);		

  pt = pt + 1;
}

//~ event adapt (i++) {
  //~ double uemax = 0.2*normf(u.x).avg;
  //~ adapt_wavelet ({u,f}, (double[]){uemax,0.001}, minlevel = 7, maxlevel = MAXLEVEL);
//~ }

/*
#define IN_REGION (((y>0.125 && y<0.175)  && x<0.4) || (x>0.4 && y<0.3) )
event adapt (i++) {
  scalar f1[];
  scalar region[];
  foreach() 
  f1[] = f[]; 	
  foreach()
  region[] = IN_REGION*noise();
  adapt_wavelet ({region, f}, (double[]){0.01,1e-3}, minlevel = 4, maxlevel = MAXLEVEL);
  }
*/