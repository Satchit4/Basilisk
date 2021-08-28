#include "embed.h"
#include "navier-stokes/centered.h"
#include "vof.h"
double rho1 = 1.225, mu1 = 0.0000181, rho2 = 1000.0, mu2 = 0.0089, radius = 0.00008; // check these values again
#include "tension_palas.h"
#include "output_vtu_foreach.h"

scalar f[], * interfaces = {f};
face vector alphav[], muv[], av[];
scalar rhov[];

int MAXLEVEL = 10;

//~ u.n[left]  = dirichlet(y<0.005 && y>-0.005 ? 0.0423858 : 0.0);
//p[left]  = dirichlet(0.);
//u.n[left] = neumann(0.);

//u.n[right] = neumann(y<0.00038 ? 0.1:0.); // check actual value of inlet velocity
//p[right]   = dirichlet(y<0.00038 ? 300.:0.);

//u.n[right]  = dirichlet(y<0.00038 ? 0.02 : 0.0);
//p[right]    = neumann(0.);
//pf[right]   = neumann(0.);

//u.n[left] = neumann(0.);
//p[left]   = dirichlet(0.);
//pf[left]  = dirichlet(0.);

u.n[left] = neumann(0.);
u.t[left] = neumann(0.);
p[left] = dirichlet(1.);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right] = dirichlet(525);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
f[embed]   = 0.;

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
f[bottom]   = 0.;  

int main(int argc, char* argv[]) {
  size (0.0006);

  MAXLEVEL=atoi(argv[1]);
  //N = 1 << MAXLEVEL;
  origin (0.0,0.0);

  a = av;
  alpha = alphav;
  rho = rhov;
  mu = muv;
	
  f.sigma = 0.073; // check it
  TOLERANCE = 1e-4;
  run();
}

event init (t = 0) {
 /* for (scalar s in {f, u, g})
  s.prolongation = refine_injection;
  restore ("snap-141");
  for (scalar s in {f, u, g})
  s.prolongation = refine_embed_linear;	
*/		
  vertex scalar phi[];    
  foreach_vertex(){ phi[] = (((y>0.00015 && y<0.00023)  && x<0.000254) || (x>0.000254 && y<0.00038) ); }
  boundary ({phi});
  fractions (phi, cs, fs);
  refine (cs[]>0.0 && level<MAXLEVEL);
  
  fraction (f, sqrt(sq(x-0.000334)+sq(y-0.00019)) < radius);
 //boundary(all);
}

#define rho(f) (clamp(f,0,1)*(rho1 - rho2) + rho2)
#define mu(f)  (clamp(f,0,1)*(mu1 - mu2) + mu2)

event properties (i++) {
  foreach()
    rhov[] = rho(f[])*cm[];
  boundary ({rhov});
  foreach_face () {
    double ff = (f[] + f[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    muv.x[] = fm.x[]*mu(ff);
  }
  boundary ((scalar *){muv});
}
/*
event acceleration (i++) {
  foreach_face(y)
    av.y[] -= 981.0;
  boundary ((scalar *){av});
}
*/
event logfile (i++) {
  if (i == 0)
   fprintf (ferr,"t dt \n");
   fprintf (ferr, "%g %g \n", t, dt);
}

int pt = 0;
event snapshot (t = 0; t += 0.00005; t <= 0.0035) {
	// you can use these files for extracting vtu output.
  char name1[80];
  sprintf (name1, "snap-%d", pt);
  p.nodump = false;
  dump(name1);

  //~ scalar *list;		
  //~ bool linear;
  //~ char name2[80];
  //~ sprintf(name2, "bessel-%d.vtu",pt);
  //~ FILE *ptr2 = fopen(name2,"w");	
  //~ output_vtu_bin_foreach ((scalar *) {f, p, cm}, (vector *) {u}, N, ptr2, false);	
  //~ fclose(ptr2);		

  pt = pt + 1;
}



//~ event snapshot (t = 0.; t += 0.1; t <= 20.0) {

  //~ char name[80];
  //~ sprintf (name, "snapshot-%g", t);
  //~ scalar pid[];
  //~ foreach()
    //~ pid[] = fmod(pid()*(npe() + 37), npe());
  //~ boundary ({pid});
  //~ dump (name);
  //~ dump ("restart");
 
  //~ //dump("restart");
//~ }


//~ For parallel output

void backup_fields (scalar f, vector u, int nf)
{
  char name[80], subname[80];
  FILE * fp ;
	nf > 0 ? sprintf(name, "sol_%4.4d_n%3.3d.vtu", nf,pid()) : sprintf(name, "sol_n%3.3d.vtu",pid());
	fp = fopen(name, "w"); output_vtu_ascii_foreach ((scalar *) {f, p, cm}, (vector *) {u}, N, fp, false); fclose (fp);

#if _MPI
	if (pid()==0){
		nf > 0 ? sprintf(name, "sol-%d.pvtu", nf) : sprintf(name, "sol-0.pvtu");
		nf > 0 ? sprintf(subname, "sol_%4.4d", nf) : sprintf(subname, "sol");
		fp = fopen(name, "w"); output_pvtu_ascii ((scalar *) {f, p, cm}, (vector *) {u}, N, fp, subname); fclose (fp);
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	
}

event logfile2(t = 0.; t += 0.00005; t <= 0.0035) {
	
	backup_fields(f,u,t*100000);
	
}

#define IN_REGION (((y>0.00015 && y<0.00023)  && x<0.000254) || (x>0.000254 && y<0.00038))
event adapt (i++) {
  scalar f1[];
  scalar region[];
  foreach() 
  f1[] = f[]; 	
  foreach()
  region[] = IN_REGION*noise();
  //adapt_wavelet ({region,f}, (double[]){0.01,0.001}, minlevel = 3, maxlevel = MAXLEVEL);
  adapt_wavelet ({f,u}, (double[]){0.01,0.001,0.001}, minlevel = 7, maxlevel = MAXLEVEL);
}
