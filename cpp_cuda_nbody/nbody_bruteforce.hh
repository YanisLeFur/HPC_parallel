#ifndef NBODYBRUTEFORCE_H_
#define NBODYBRUTEFORCE_H_

#include <stdio.h> 
#include "parameters.hh"
#include <cuda_runtime.h>

void nbodybruteforce (particle_t * array, int nbr_particles, int nbr_iterations,const dim3 block_size) ;
void update_particles(particle_t * array, int nbr_particles, double step);


#endif /*NBODYBRUTEFORCE_H_*/
