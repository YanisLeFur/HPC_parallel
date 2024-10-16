#ifndef READER_H_
#define READER_H_

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "nbody_bruteforce.hh"

#include "parameters.hh"

particle_t * read_test_case(const char *  fn);
int get_nbr_particles(const char *  fn);

#endif /*READER_H_*/
