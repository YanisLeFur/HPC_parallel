#include "nbody_bruteforce.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <cuda_runtime.h>
/*--------------------------------------------------------------------------*/
__global__ void compute_force_one_particle(particle_t * array, int nbr_particles, double step){
    
    int id = blockDim.x*blockIdx.x+threadIdx.x;
	if(id<nbr_particles){
    particle_t * p1 = &array[id];


    double x_sep, y_sep, z_sep, dist_sq, grav_base;
	double F_x=0.;
	double F_y=0.;
	double F_z=0.;
	particle_t tmp;

	for (int i = 0 ; i  < nbr_particles ; i++){
		tmp = array[i];
		if (i!=p1->id){
			x_sep = p1->x - tmp.x;
			y_sep = p1->y - tmp.y;
			z_sep = p1->z - tmp.z;
            float dist_tmp = (x_sep*x_sep + y_sep*y_sep + z_sep*z_sep);
			dist_sq = ((dist_tmp) > (0.01) ? (dist_tmp) : (0.01));

            grav_base = GRAV_CONSTANT*(p1->m)*(tmp.m)/dist_sq / sqrt(dist_sq);
			F_x += grav_base*x_sep;
			F_y += grav_base*y_sep;
			F_z += grav_base*z_sep;
		}
	}
	p1->fx = F_x;
	p1->fy = F_y;
	p1->fz = F_z;
	}
}

/*--------------------------------------------------------------------------*/
__global__ void compute_force_between_two_particles(particle_t * array, int nbr_particles, double step){    
	int idx = blockDim.x*blockIdx.x+threadIdx.x;
	int idy = blockDim.y*blockIdx.y+threadIdx.y;
   
	if((idx<nbr_particles)&&(idy<nbr_particles)){
    particle_t * p1 = &array[idx];
    double x_sep, y_sep, z_sep, dist_sq, grav_base;
	double F_x=0.;
	double F_y=0.;
	double F_z=0.;
	particle_t tmp;
	tmp = array[idy];
	if (idy!=p1->id){
		x_sep = p1->x - tmp.x;
		y_sep = p1->y - tmp.y;
		z_sep = p1->z - tmp.z;
        float dist_tmp = (x_sep*x_sep + y_sep*y_sep + z_sep*z_sep);
		dist_sq = ((dist_tmp) > (0.01) ? (dist_tmp) : (0.01));

        grav_base = GRAV_CONSTANT*(p1->m)*(tmp.m)/dist_sq / sqrt(dist_sq);
		F_x += grav_base*x_sep;
		F_y += grav_base*y_sep;
		F_z += grav_base*z_sep;
	}

	p1->fx = F_x;
	p1->fy = F_y;
	p1->fz = F_z;
	}

}
/*--------------------------------------------------------------------------*/

void nbodybruteforce (particle_t * array, int nbr_particles, int nbr_iterations,const dim3 block_size) {
	double step = 1.;
    dim3 grid_size; 
	grid_size.x = ceil(nbr_particles/float(block_size.x));
	grid_size.y = ceil(nbr_particles/float(block_size.y));

	particle_t *cuda_array{nullptr};
	cudaMallocManaged(&cuda_array, sizeof(particle_t)*nbr_particles);
	cudaMemcpy(cuda_array, array, sizeof(particle_t)*nbr_particles, cudaMemcpyHostToDevice);

    static bool first{true};
    if (first) {
		printf("Block size: %d : %d\n",block_size.x,block_size.y);
		printf("Grid size: %d : %d\n",grid_size.x,grid_size.y);
        first = false;
    }

	for (int n = 0 ; n  < nbr_iterations ; n++){
		//printf("ITERATION %d \n",n);
            #if defined(PER_ROW)
               compute_force_one_particle <<<grid_size.x, block_size.x>>> (cuda_array, nbr_particles, step);
            #else
                compute_force_between_two_particles <<<grid_size, block_size>>> (cuda_array, nbr_particles, step);
            #endif
            //TODO: did you forget to synchronize ?
            cudaDeviceSynchronize();
			update_particles(cuda_array,nbr_particles,step);
	}

cudaFree(cuda_array);

cudaError_t err = cudaGetLastError();  // add
if (err != cudaSuccess) printf("CUDA error: %s\n",cudaGetErrorString(err)); // add


}








void update_particles(particle_t * array, int nbr_particles, double step) {
	for (int i = 0 ; i  < nbr_particles ; i++){
		particle_t *p1 = &array[i];

		p1->vx += p1->fx/p1->m * step;
	  	p1->vy += p1->fy/p1->m * step;
	  	p1->vz += p1->fz/p1->m * step;

    	p1->x += p1->vx * step;
	  	p1->y += p1->vy * step;
	  	p1->z += p1->vz * step;
	}
}