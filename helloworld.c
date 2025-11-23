#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Define the Particle Structure (Struct)
typedef struct {
    double x, y, z;      // Position
    double vx, vy, vz;   // Velocity
    double q;            // Charge
    double m;            // Mass
} Particle;

// Simulation Settings
const double DT = 0.01;      // Time step (delta t)
const double T_MAX = 5.0;    // Simulation duration
const double BX = 0.0;       // Magnetic Field X component
const double BY = 0.0;       // Magnetic Field Y component
const double BZ = 1.0;       // Magnetic Field Z component (Field is currently along Z axis)

// **********************************************
// * This is where the core Physics Engine lies
// **********************************************
void update_particle(Particle *p, double dt) {
    // Basic Euler Method framework for velocity and position update
    
    // 1. Calculate the Lorentz Force (Magnetic component only: F = q(v x B))
    double fx = p->q * (p->vy * BZ - p->vz * BY); // Fx = q(Vy*Bz - Vz*By)
    double fy = p->q * (p->vz * BX - p->vx * BZ); // Fy = q(Vz*Bx - Vx*Bz)
    double fz = p->q * (p->vx * BY - p->vy * BX); // Fz = q(Vx*By - Vy*Bx)
    
    // 2. Calculate Acceleration (a = F/m)
    double ax = fx / p->m;
    double ay = fy / p->m;
    double az = fz / p->m;

    // 3. Update Velocity (v_new = v_old + a * dt)
    p->vx += ax * dt;
    p->vy += ay * dt;
    p->vz += az * dt;

    // 4. Update Position (x_new = x_old + v_new * dt)
    p->x += p->vx * dt;
    p->y += p->vy * dt;
    p->z += p->vz * dt;
}

int main() {
    // Open the Output File
    FILE *output_file = fopen("particle_trajectory.dat", "w");
    if (output_file == NULL) {
        printf("Error: Could not open the file.\n");
        return 1;
    }

    // Initialize the Particle
    Particle electron = {
        .x = 0.0, .y = 0.0, .z = 0.0, // Initial Position
        .vx = 0.1, .vy = 1.0, .vz = 0.0, // Initial Velocity (Given movement in y direction)
        .q = -1.0, .m = 1.0 // Charge and Mass (simply set to 1)
    };

    // Write the Header Row (Needed for Python plotting)
    fprintf(output_file, "t\tx\ty\tz\n");

    double time = 0.0;
    while (time <= T_MAX) {
        // Write data to the file
        fprintf(output_file, "%.4f\t%.4f\t%.4f\t%.4f\n", time, electron.x, electron.y, electron.z);

        // Update the particle's state
        update_particle(&electron, DT);

        // Advance time
        time += DT;
    }

    // Close the file and end the program
    fclose(output_file);
    printf("Simulation completed successfully. Data written to 'particle_trajectory.dat'.\n");
    return 0;
}