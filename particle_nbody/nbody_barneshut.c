#include "nbody_barneshut.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stddef.h>


void distribute_array(particle_t * new_array, particle_t * to_distribute,int nbr_particles, int prank, int psize){
	int start = prank*(nbr_particles/psize);
	int end = (prank+1)*(nbr_particles/psize);
	if (prank==psize-1){end = nbr_particles;}
	int count = 0;

	for(int i = start; i<end; i++){
		new_array[count]=to_distribute[i];
		count++;
	}



}

/*
Implementation of a barnes-hut algorithm for the N-Body problem.
*/
void nbodybarneshut (particle_t * array, int nbr_particles, int nbr_iterations, int psize, int prank) 
{	
	int n;
	double step = TIMESTEP;
	node * root1;
	particle_t tmp;
	particle_t * recv_particles;
	
	//defining the MPI datatype for particle_t ------------------------------------------------------------------------------------------
	MPI_Datatype particle_type;
	int lengths[12]={1,1,1,1,1,1,1,1,1,1,1,1};
	MPI_Aint displacement[12];
	particle_t dummy_part={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0,0.};
	MPI_Aint base_address;
  	MPI_Get_address(&dummy_part,&base_address);
	MPI_Get_address(&dummy_part.x,&displacement[0]);
  	MPI_Get_address(&dummy_part.y,&displacement[1]);
	MPI_Get_address(&dummy_part.z,&displacement[2]);
  	MPI_Get_address(&dummy_part.vx,&displacement[3]);
	MPI_Get_address(&dummy_part.vy,&displacement[4]);
  	MPI_Get_address(&dummy_part.vz,&displacement[5]);
	MPI_Get_address(&dummy_part.fx,&displacement[6]);
  	MPI_Get_address(&dummy_part.fy,&displacement[7]);
	MPI_Get_address(&dummy_part.fz,&displacement[8]);
  	MPI_Get_address(&dummy_part.m,&displacement[9]);
	MPI_Get_address(&dummy_part.id,&displacement[10]);
  	MPI_Get_address(&dummy_part.V,&displacement[11]);

	displacement[0] = MPI_Aint_diff(displacement[0], base_address);
  	displacement[1] = MPI_Aint_diff(displacement[1], base_address);
	displacement[2] = MPI_Aint_diff(displacement[2], base_address);
  	displacement[3] = MPI_Aint_diff(displacement[3], base_address);
	displacement[4] = MPI_Aint_diff(displacement[4], base_address);
  	displacement[5] = MPI_Aint_diff(displacement[5], base_address);
	displacement[6] = MPI_Aint_diff(displacement[6], base_address);
  	displacement[7] = MPI_Aint_diff(displacement[7], base_address);
	displacement[8] = MPI_Aint_diff(displacement[8], base_address);
  	displacement[9] = MPI_Aint_diff(displacement[9], base_address);
	displacement[10] = MPI_Aint_diff(displacement[10], base_address);
  	displacement[11] = MPI_Aint_diff(displacement[11], base_address);

	MPI_Datatype types[12]={MPI_DOUBLE,
							MPI_DOUBLE,
							MPI_DOUBLE,
							MPI_DOUBLE,
							MPI_DOUBLE,
							MPI_DOUBLE,
							MPI_DOUBLE,
							MPI_DOUBLE,
							MPI_DOUBLE,
							MPI_DOUBLE,
							MPI_INT,
							MPI_DOUBLE};
	
	MPI_Type_create_struct(12,lengths,displacement,types,&particle_type);
  	MPI_Type_commit(&particle_type);

	//--------------------------------------------------------------------------------------------------------------------------------------------


	recv_particles = array;
	root1 = malloc(sizeof(node));	
	tmp = getMinMax(array, nbr_particles);
	init_tree(&tmp, root1);	
	construct_bh_tree(array,root1, nbr_particles);
			
	for (n = 0 ; n  < nbr_iterations ; n++){

		if(nbr_particles>1){
		nbr_particles = root1->sub_nbr_particles;
		int Nblock=nbr_particles/psize;
		if (prank==psize-1){Nblock =nbr_particles-(psize-1)*Nblock;}
		
		int displacement[psize];
		int particles_per_process[psize];

		for(int i = 0 ;i<psize;i++){
			if(i ==0){
				displacement[i]=0;
				particles_per_process[i] = nbr_particles/psize;
			}if(i==psize-1){
				displacement[i]=i*(nbr_particles/psize);
				particles_per_process[i]=nbr_particles-(i*(nbr_particles/psize));
			}
			else{
				displacement[i]=i*(nbr_particles/psize);
				particles_per_process[i]=nbr_particles/psize;
			}	
		}
		particle_t * send_particles;
		send_particles  = malloc(Nblock*sizeof(particle_t));
		distribute_array(send_particles,recv_particles,nbr_particles,prank,psize);
		
		for(int i = 0; i<Nblock; i++){
			compute_force_particle(root1,&send_particles[i]);
		}

		update_particles(send_particles,step,Nblock);

		MPI_Allgatherv(send_particles,Nblock,particle_type,recv_particles,particles_per_process,displacement,particle_type, MPI_COMM_WORLD);
			
		particle_t * tmp;
		tmp = malloc(nbr_particles*sizeof(particle_t));
		int new_nbr_particles = 0;

		for(int i = 0; i < nbr_particles;i++){
			if(!is_particle_out_of_scope(&recv_particles[i],root1)){
				tmp[new_nbr_particles]=recv_particles[i];
				new_nbr_particles++;
			}
		}

		nbr_particles = new_nbr_particles;
		recv_particles = malloc(new_nbr_particles*sizeof(particle_t));
		recv_particles = tmp;
		clean_tree(root1);
		construct_bh_tree(recv_particles,root1,nbr_particles);
		}
	}
	if(prank==0){
	printf("It remains %d particles in space \n",root1->sub_nbr_particles);	
	}
	clean_tree(root1);
	free(root1);
	MPI_Type_free(&particle_type);
}


/*

1. If the current node is an external node (and it is not body b), calculate the force exerced by the current node on b, and add this amount to b’s net force.
    
2. Otherwise, calculate the ratio s/d. If s/d < θ, treat this internal node as a single body, and calculate the force it exerts on body b, and add this amount to b’s net force.

3. Otherwise, run the procedure recursively on each of the current node’s n.

Once the computation of the force applied to the particles is complete, the new position of the particles is computed, and a new tree corresponding to the new position is created. 

*/

/*
Move all the particles from node n to new_root
*/

void update_particles(particle_t * new_particles,double step,int size_array){
	for(int i = 0; i<size_array;i++){
		double ax,ay,az;
		new_particles[i].x += new_particles[i].vx * step;	 
		new_particles[i].y += new_particles[i].vy * step;	 
		new_particles[i].z += new_particles[i].vz * step;	 
		ax = new_particles[i].fx/new_particles[i].m;
		ay = new_particles[i].fy/new_particles[i].m;
		az = new_particles[i].fz/new_particles[i].m;
		new_particles[i].vx += ax*step;
		new_particles[i].vy += ay*step;
		new_particles[i].vz += az*step;
		
	}
}



void move_all_particles(node * new_root, node * n, double step) {
	int i;
	if(n->children != NULL){
		for (i = 0; i < 8; i++){
			move_all_particles(new_root, &n->children[i], step);
		}
	}else{
		particle_t * p = n->particle;
		move_particle(new_root, n, p,step);
	}
}

/*
Compute new position/velocity of the particle
*/

void move_particle(node * root, node * n, particle_t * p, double step) {
	double ax,ay,az;

	if ((p==NULL)||(n==NULL)) return;

	p->x += p->vx * step;	 
	p->y += p->vy * step;	 
	p->z += p->vz * step;	 
	ax = p->fx/p->m;
	ay = p->fy/p->m;
	az = p->fz/p->m;
	p->vx += ax*step;
	p->vy += ay*step;
	p->vz += az*step;
		
	if (! is_particle_out_of_scope(p,root)) {
		insert_particle(p,root);
	}else{
		printf("\tparticle %d is out of scope. It will be destroyed at next iteration (%f,%f,%f)\n",p->id,p->x,p->y,p->z);
		n->particle = NULL;
	}
	
}

/*
Check if a particle is out of scope (lost body in space)
*/

bool is_particle_out_of_scope(particle_t * p, node * root){
	bool ret = false;
	if ((p->x < root->minx)||(p->y < root->miny)||(p->z < root->minz)) ret = true;
	if ((p->x > root->maxx)||(p->y > root->maxy)||(p->z > root->maxz)) ret = true;	
	return ret;
}


/*
Clean tree root
*/
void clean_tree(node * root) {
	int i;
	if (root == NULL) {return;}

	if(root->children != NULL){
		for (i = 0; i < 8; i++){
			clean_tree(&root->children[i]);
		}
		free(root->children);
		root->children = NULL;
		root->sub_nbr_particles=0;
	}
}



/*
compute the forces on the BH tree
*/

void compute_bh_force(node * n) {
	int i;
	if(n->children != NULL){
		for (i = 0; i < 8; i++){
			compute_bh_force(&n->children[i]);
		}
	}else{
		particle_t * p = n->particle;
		compute_force_particle(n,p);
	}
}

/*
Compute force of node n on particle p
*/

void compute_force_particle(node * n, particle_t * p){
	int i;
	double diffx,diffy,diffz,distance;
	double size;

	if ((n==NULL)||(n->sub_nbr_particles==0)){ return;}

	if ((n->particle != NULL)&&(n->children==NULL)) {
		compute_force(p, n->centerx, n->centery,  n->centerz, n->mass) ;
	}
	else{
		size = n->maxx - n->minx;
		diffx = n->centerx - p->x;
		diffy = n->centery - p->y;
		diffz = n->centerz - p->z;
		distance = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

//	The particle is far away. Use an approximation of the force
		if(size / distance < THETA) {
			compute_force(p, n->centerx, n->centery, n->centerz, n->mass);
		} else {

//      Otherwise, run the procedure recursively on each of the current node's children.
			for(i=0; i<8; i++) {
				compute_force_particle(&n->children[i], p);
			}
		}
	}
}



/*
Compute force 
*/

void compute_force(particle_t *p, double xpos, double ypos,  double zpos, double mass) {
	double xsep, ysep, zsep, dist_sq, gravity;

	xsep = xpos - p->x;
	ysep = ypos - p->y;
	zsep = zpos - p->z;
	dist_sq = max((xsep*xsep) + (ysep*ysep) + (zsep*zsep), 0.01);

	gravity = GRAV_CONSTANT*(p->m)*(mass)/ dist_sq / sqrt(dist_sq);

	p->fx += gravity*xsep;
	p->fy += gravity*ysep;
	p->fz += gravity*zsep;
}

/*
Compute all the forces in the particles
*/
void compute_force_in_node(node *n,node *root) {
	int i;
	if(n==NULL) return;

	if((n->particle != NULL)&&(n->children == NULL)) {
		particle_t*p = n->particle;
		p->fx = 0;
		p->fy = 0;
		p->fz = 0;
		compute_force_particle(root, p);
	}
	if(n->children != NULL) {
		for(i=0; i<8; i++) {
			compute_force_in_node(&n->children[i], root);
		}
	}
}




/*
Construction of the barnes-hut tree

Reminder:
Construction of the Barnes-Hut Tree in 2D
http://arborjs.org/docs/barnes-hut
*/

void construct_bh_tree(particle_t * array, node *  root, int nbr_particles){
	int i;
	for (i=0;i < nbr_particles; i++){
			insert_particle(&array[i],root);

	}
}



/*
Add particle p in node n or one of its children
*/

void insert_particle(particle_t * p, node * n){
	int octrant ;
	double totalmass = 0.;
	double totalx = 0.;
	double totaly = 0.;
	double totalz = 0.;
	int i;
// there is no particle
	if ((n->sub_nbr_particles == 0)&&(n->children==NULL)) {
		n->particle = p;
		n->centerx = p->x;
		n->centery = p->y;
		n->centerz = p->z;
		n->mass = p->m;
		n->sub_nbr_particles++;
		//p->node = n;
// There is already a particle
	}else{
		if (n->children==NULL){
			create_children(n);
			particle_t * particule_parent = n->particle;
// Insert the particle in the correct children
			octrant = get_octrant(particule_parent,n);
			n->particle = NULL;
			insert_particle(particule_parent,&n->children[octrant]);
		}
// insert the particle p
		octrant = get_octrant(p,n);
		insert_particle(p,&n->children[octrant]);

// Update mass and barycenter (sum of momentums / total mass)
		for(i=0; i<8; i++) {
			totalmass += n->children[i].mass;
			totalx += n->children[i].centerx*n->children[i].mass;
			totaly += n->children[i].centery*n->children[i].mass;
			totalz += n->children[i].centerz*n->children[i].mass;
		}
		n->mass = totalmass;
		n->centerx = totalx/totalmass;
		n->centery = totaly/totalmass;
		n->centerz = totalz/totalmass;
		//p->node = n;
		n->sub_nbr_particles++;
	}
}

/*
create 8 children from 1 node
*/

void create_children(node * n){
	n->children = malloc(8*sizeof(node));

	double x12 = n->minx+(n->maxx-n->minx)/2.;
	double y12 = n->miny+(n->maxy-n->miny)/2.;
	double z12 = n->minz+(n->maxz-n->minz)/2.;

	init_node(&n->children[SW_DOWN], n, n->minx, x12, n->miny, y12, n->minz, z12 );
	init_node(&n->children[NW_DOWN], n, n->minx, x12, n->miny, y12, z12, n->maxz );

	init_node(&n->children[SE_DOWN], n, n->minx, x12, y12, n->maxy, n->minz, z12 );
	init_node(&n->children[NE_DOWN], n, n->minx, x12, y12, n->maxy, z12, n->maxz );

	init_node(&n->children[SW_UP], n, x12, n->maxx, n->miny, y12, n->minz, z12 );
	init_node(&n->children[NW_UP], n, x12, n->maxx, n->miny, y12, z12, n->maxz );

	init_node(&n->children[SE_UP], n, x12, n->maxx, y12, n->maxy, n->minz, z12 );
	init_node(&n->children[NE_UP], n, x12, n->maxx, y12, n->maxy, z12, n->maxz );	
}

/*
Init a node n attached to parent parent. 
*/

void init_node(node * n, node * parent,  double minx, double maxx, double miny, double maxy, double minz, double maxz ){
	n->parent=parent;
	n->children = NULL;
	n->minx = minx;
	n->maxx = maxx;
	n->miny = miny;
	n->maxy = maxy;
	n->minz = minz;
	n->maxz = maxz;
	n->depth = parent->depth + 1;
	n->particle = NULL;
	n->sub_nbr_particles = 0.;
	n->centerx = 0.;
	n->centery = 0.;
	n->centerz = 0.;
	n->mass = 0.;
}



/*
get the "octrant" where the particle resides (octrant is a generalization in 3D of a 2D quadrant)
*/

int get_octrant(particle_t * p, node * n){
	int octrant=-1;
	double xmin = n->minx;
	double xmax = n->maxx;
	double x_center = xmin+(xmax-xmin)/2;

	double ymin = n->miny;
	double ymax = n->maxy;
	double y_center = ymin+(ymax-ymin)/2;

	double zmin = n->minz;
	double zmax = n->maxz;
	double z_center = zmin+(zmax-zmin)/2;
	if (n==NULL) printf("ERROR: node is NULL \n");
	if (p==NULL) printf("ERROR: particle is NULL \n");

	// order : x -> y -> z
	if(p->x <= x_center) {
		if(p->y <= y_center) {
			if(p->z <= z_center) {
				octrant = SW_DOWN;
			}else{
				octrant = NW_DOWN;
			}
		} else {
			if(p->z <= z_center) {
				octrant = SE_DOWN;
			}else{
				octrant = NE_DOWN;
			}
		}
	} else {
		if(p->y <= y_center) {
			if(p->z <= z_center) {
				octrant = SW_UP;
			}else{
				octrant = NW_UP;
			}
		} else {
			if(p->z <= z_center) {
				octrant = SE_UP;
			}else{
				octrant = NE_UP;
			}
		}
	}
	return octrant;
}

/*
Init the tree

Remark :We use a particle struct to transfer min and max values from main
*/

void init_tree(particle_t * particle, node * root){
	root->minx = particle->x;
	root->maxx = particle->vx;
	root->miny = particle->y;
	root->maxy = particle->vy;
	root->minz = particle->z;
	root->maxz = particle->vz;
	root->particle = NULL;
	root->sub_nbr_particles = 0;
	root->parent = NULL;
	root->children = NULL;
	root->centerx = 0.;
	root->centery = 0.;
	root->centerz = 0.;
	root->mass = 0.;
	root->depth = 0;
}

/*
============================================
Utilities for testing
============================================
*/


/* print the tree */
void print_tree(node * root){
	node * tmp;
	int i;
	if (root->children!=NULL){
		for (i =0;i<8;i++){
			tmp = &root->children[i];
			print_tree(tmp);
		}
	}
	print_node(root);

}


/* 
print a node 
*/
void print_node(node * n){
	int d = n->depth;
	int i;
	for (i=0;i<d;i++){
		printf("\t");
	}
	printf("[level %d]",d);
	/*printf(" ([%f:%f:%f])",n->centerx, n->centery,n->centerz);
	printf(" Node ");
	printf(" M = %f", n->mass);*/
	printf(" has %d particles ", n->sub_nbr_particles);
	if (n->particle!=NULL){
		particle_t * p = n->particle;
		printf(". Particle ID = %d",p->id);
	}
	printf("\n");
}


/*
print a particle 
*/
void print_particle(particle_t * p){
	printf("[Particle %d]",p->id);
	printf(" position ([%f:%f:%f])",p->x, p->y, p->z);
	printf(" M = %f", p->m);
	printf("\n");
}



