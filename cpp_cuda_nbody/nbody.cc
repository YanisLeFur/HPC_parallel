#include <cuda_runtime.h>
#include <iostream>
#include <math.h>

#include <chrono>

#include "parameters.hh"
#include "nbody_bruteforce.hh"
#include "reader.hh"

void print_parameters(){
	printf("====================================================\n");
	printf("N-Body 3D simulation code for MATH-454 course EPFL  \n");
	printf("Parameters for the Barnes-Hut algorithm:\n");
	printf("\n");
	printf("Gravitational constant : %f\n",GRAV_CONSTANT);
	printf("Theta                  : %f\n",THETA);
	printf("Time step              : %f\n",TIMESTEP);
	printf("Space multiplicator    : %f\n",SIZEOFSPACE);
	printf("Number of iterations   : %d\n",NBRITERATIONS);
	printf("\n");
	printf("These parameters can be modified in \"parameters.h\"\n");
	printf("\n");
	printf("(c) 2020, Vincent Keller (Vincent.Keller@epfl.ch)\n");
	printf("====================================================\n");
	printf("\n");
}

typedef std::chrono::high_resolution_clock clk;
typedef std::chrono::duration<double> second;

static void usage(const std::string & prog_name) {
  std::cerr << prog_name << " <grid_size> <block_size [default: 32]>" << std::endl;
  exit(0);
}

int main ( int argc, char **argv ) {



	particle_t * array;
	int nbr_iterations;
	int nbr_particles;
	nbr_iterations = NBRITERATIONS;

	print_parameters();
	dim3 block_size = {32, 1, 1};

	if (argc < 3)
	{
		usage(argv[0]);
		exit(1);
	}
	if(argc>=3)   
	{ 
		std::cout<<"Read data from file"<<std::endl;
		nbr_particles = get_nbr_particles(argv[1]);
		array = read_test_case(argv[1]);
		block_size.x = 	std::stoi(argv[2]);
		std::cout<<"Number of particles : "<<nbr_particles<<std::endl;
	}
	if(argc==4){
		block_size.y = 	std::stoi(argv[2]);
	}

	  // By default, we use device 0,
  int dev_id = 0;

  cudaDeviceProp device_prop;
  cudaGetDevice(&dev_id);
  cudaGetDeviceProperties(&device_prop, dev_id);
  if (device_prop.computeMode == cudaComputeModeProhibited) {
    std::cerr << "Error: device is running in <Compute Mode Prohibited>, no "
                 "threads can use ::cudaSetDevice()"
              << std::endl;
    return -1;
  }

  auto error = cudaGetLastError();
  if (error != cudaSuccess) {
    std::cout << "cudaGetDeviceProperties returned error code " << error
              << ", line(" << __LINE__ << ")" << std::endl;
    return error;
  } else {
    std::cout << "GPU Device " << dev_id << ": \"" << device_prop.name
              << "\" with compute capability " << device_prop.major << "."
              << device_prop.minor << std::endl;
  }



	std::cout << "BRUTE FORCE simulation starting" << std::endl;
	auto t1 = clk::now();
	nbodybruteforce(array, nbr_particles, nbr_iterations,block_size);
	auto t2 = clk::now();
	second time = t2 - t1;
	std::cout << "N-Body brute force for " << nbr_particles << " particles : " << time.count() << " [s]" << std::endl;
	


	free(array);

	printf("Simulation finished \n");
	return 0;
}


