#include <stdio.h> //processor directive input output
/*Function,Purpose,Example
printf(),Output data (like text or variable values) to the console.,"printf(""Hello World!\n"");"
scanf(),Input data from the user and store it in a variable.,"scanf(""%d"", &my_number);"
fopen(),Open a file for reading or writing.,Used for file I/O operations.
fclose(),Close an open file.,Used for file I/O operations.*/
#include <stdlib.h>//stnumeric conv, rand,andard lib memory
/*malloc(), calloc(), realloc(), free() atoi() (ASCII to integer), atof() (ASCII to float)
rand(), srand()exit(), system() qsort(), bsearch()*/
#include <math.h>

// particle structure definition
typedef struct {
    double x, y, z;      // Position
    double vx, vy, vz;   // Velocity
    double q;            // Charge
    double m;            // Mass
} Particle;

// sim settings
const double DT = 0.1;      // Time step (delta t) 
const double T_MAX = 100.0;    // duration
const double BX = 0.0;       // Magnetic Field X component
const double BY = 0.0;       // Y component
const double BZ = 1.0;       // Z component 
/*update particles position
pointer--> particle *p
*/
/*(*p).q  is p->q access the q charge member inside the Particle structure that 
   the pointer p is pointing to */
/*void update_particle(Particle *p, double dt) {
    // EULER method
    //lorentz force (magnetic component F=q(v x B))
    double fx = p->q * (p->vy * BZ - p->vz * BY); // Fx=q(vyBz-vzBy) cross product
    double fy = p->q * (p->vy * BX - p->vx * BZ);
    double fz = p->q * (p->vx * BY - p->vy * BX);

    // acceleration (a=F/m)
    double ax = fx / p->m;
    double ay = fy / p->m;
    double az = fz / p->m;

    // update velocity (vf=vi+a*dt) and position (xf=xi+vf*dt)
    p->vx += ax * dt;
    p->vy += ay * dt;
    p->vz += az * dt;

    p->x += p->vx * dt;
    p->y += p->vy * dt;
    p->z += p->vz * dt;
}*/
/*boris algorithm*/
void update_particle_boris(Particle *p, double dt){

    double half_dt=0.5*dt;
/*T=q*B*dt/2m 
w_c=q*B/m --> cyclotron frequency
T=w_c*dt/2 */ 
    double tx=p->q*BX*half_dt/p->m;
    double ty=p->q*BY*half_dt/p->m;
    double tz=p->q*BZ*half_dt/p->m;
//v_minus initial velocity before magnetic rot applied, E=0 so v_minus is p's current velocity

    double vx_minus=p->vx;
    double vy_minus=p->vy;
    double vz_minus=p->vz;
//v_prime velocity after the first partial rotation
//v'=v_minus+(v_minus x T)
    double vx_prime=vx_minus+(vy_minus*tz-vz_minus*ty);
    double vy_prime=vy_minus+(vz_minus*tx-vx_minus*tz);
    double vz_prime=vz_minus+(vx_minus*ty-vy_minus*tx);
// S is needed to complete the rotation symetrically
//S=2T/(1+T^2)
    double t_mag_sq=tx*tx+ty*ty+tz*tz;
    double sx=2*tx/(1+t_mag_sq);
    double sy=2*ty/(1+t_mag_sq);
    double sz=2*tz/(1+t_mag_sq);
//v=v_minus+(v'x S)
    p->vx=vx_minus+(vy_prime*sz-vz_prime*sy);
    p->vy=vy_minus+(vz_prime*sx-vx_prime*sz);
    p->vz=vz_minus+(vx_prime*sy-vy_prime*sx);
// updating the position
    p->x+=p->vx*dt;
    p->y+=p->vy*dt;
    p->z+=p->vz*dt;

}

int main(){
    //output file
    FILE *output_file = fopen("trajectory_boris.dat","w");
    if (output_file== NULL ){
        printf("Couldn't open the file \n");
        return 1;
    }
    //initialize
    Particle electron ={
        .x = 0.0, .y = 0.0, .z = 0.0,
        .vx = 0.1, .vy = 1.0, .vz = 0.0, // movement in y direction
        .q = -1.0, .m = 1.0 // set to 1
    };
    // for python plot
    fprintf(output_file, "t\tx\ty\tz\n");
    double time=0.0;
    while ( time<= T_MAX){
        fprintf(output_file,"%.4f\t%.4f\t%.4f\t%.4f\n", time, electron.x, electron.y, electron.z);
        //update
        update_particle_boris(&electron,DT);
        time+=DT;
    }   
    fclose(output_file);
    printf("sim completed");
    return 0;
}

/* euler methods makes small approximation errors at every step, the velocity vector gets slightly longer
longer in each step, the particle gains artificial energy and spirals outwards instead of staying in a stable
orbit. ENERGY DRIFT*/
/*accelerate the particle using E field for half a time step, rotate the velocity vector around the B field
vector, half push again*/