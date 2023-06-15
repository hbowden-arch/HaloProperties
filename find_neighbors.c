#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <assert.h>
#include "check_syscalls.h"
#include "make_sf_catalog.h"
#include "config_vars.h"
#include "config.h"
#include "version.h"

#include <stdbool.h>

#define FAST3TREE_TYPE struct catalog_halo
#define FAST3TREE_DIM 2
#include "fast3tree.c"

//Defining some constants
int64_t NEIGHBOR_COUNT = 50;
double MIN_STELLAR_MASS=1e9;
double MAX_NEIGHBOR_MASS_OFFSET = 0.0316227766016838; //10^-1.5
double MAX_NEIGHBOR_VELOCITY_OFFSET=1000; // Changed back to 1000 on 5/17/23
double MAX_TARGET_STELLAR_MASS=1e15;

#define CYLINDER_BINS 8
#define CYLINDER_NUM_RADII 4
double CYLINDER_BIN_SPACING = 250; //km/s
double CYLINDER_RADII[CYLINDER_NUM_RADII] = {0.5, 1.0, 2.0, 5.0}; //Make array of cylinder radii

double SCALE_FACTOR=1;
double H=100; //in km/s / Mpc/h 

//Create a structure called neighbor
struct neighbor {
  struct catalog_halo *h;
  float r; //projected distance to the neighbor
  float z; //redshift separation
};

struct fast3tree *tree=NULL;
struct fast3tree_results *results=NULL;

void find_neighbors(struct catalog_halo *h, struct fast3tree *t, struct fast3tree_results *r);
void count_cylinders(struct catalog_halo *h, struct fast3tree *t, struct fast3tree_results *r);

int main(int argc, char **argv) {
  if (argc<3) {
    fprintf(stderr, "%s", INFO_STRING);
    fprintf(stderr, "Usage: %s config_file.cfg sfr_catalog.bin [min_stellar_mass max_neighbor_mass_offset max_target_stellar_mass]\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  do_config(argv[1]);
  
  if (argc>3) {
    MIN_STELLAR_MASS = atof(argv[3]);
    if (MIN_STELLAR_MASS < 20) MIN_STELLAR_MASS = pow(10, MIN_STELLAR_MASS);
  }

  if (argc>4) {
    MAX_NEIGHBOR_MASS_OFFSET = atof(argv[4]);
    if (MAX_NEIGHBOR_MASS_OFFSET < 0)
      MAX_NEIGHBOR_MASS_OFFSET = pow(10, MAX_NEIGHBOR_MASS_OFFSET);
  }

  if (argc>5) {
    MAX_TARGET_STELLAR_MASS = atof(argv[5]);
    if (MAX_TARGET_STELLAR_MASS<20)
      MAX_TARGET_STELLAR_MASS = pow(10, MAX_TARGET_STELLAR_MASS);
  }
  
  //Set scale factor if sfr_catalog has a standard name
  char *catalog_filename = strstr(argv[2], "sfr_catalog_");
  if (catalog_filename) SCALE_FACTOR = atof(catalog_filename+strlen("sfr_catalog_"));
  H = 100.0*sqrt(Om*pow(1.0/SCALE_FACTOR, 3.0)+Ol);
  
  //Load halos from supplied file
  struct catalog_halo *halos = NULL;
  int64_t file_length = 0, num_halos = 0;
  check_slurp_file(argv[2], (void **)(&halos), &file_length);
  if (file_length % sizeof(struct catalog_halo)) {
    fprintf(stderr, "[Error] Length of file %s is not divisible by halo structure size.\n", argv[2]);
    exit(EXIT_FAILURE);
  }
  num_halos = file_length / sizeof(struct catalog_halo);


  //Filter halos to remove those that should be ignored and those below the stellar mass cut
  for (int64_t i=0; i<num_halos; i++) {

    if ((halos[i].flags & IGNORE_FLAG) || (halos[i].obs_sm < MIN_STELLAR_MASS)) {
      halos[i] = halos[num_halos-1];
      i--;
      num_halos--;
      continue;
    }

    halos[i].pos[2] += halos[i].pos[5]/H;
    if (halos[i].pos[2] < 0) halos[i].pos[2] += BOX_SIZE;
    if (halos[i].pos[2] > BOX_SIZE) halos[i].pos[2] -= BOX_SIZE;

    if (halos[i].pos[0] < 0) halos[i].pos[0] += BOX_SIZE;
    if (halos[i].pos[0] > BOX_SIZE) halos[i].pos[0] -= BOX_SIZE;

    if (halos[i].pos[1] < 0) halos[i].pos[1] += BOX_SIZE;
    if (halos[i].pos[1] > BOX_SIZE) halos[i].pos[1] -= BOX_SIZE;
  }




  //build kdtree
  struct fast3tree *tree = fast3tree_init(num_halos, halos);
  struct fast3tree_results *results = fast3tree_results_init();
  _fast3tree_set_minmax(tree,0,BOX_SIZE);

  //print header
  printf("#Neighbors and cylinder counts for %s\n", argv[2]);
  printf("#Header:");
  printf("Halo_ID(0) X Y Z Halo_Mass Stellar_Mass(0) ");
  for (int j = 1; j <= NEIGHBOR_COUNT; j++) {
	 // printf("Halo_ID(%d) UPID(%d) X(%d) Y(%d) Z(%d) Halo_Mass(%d) Stellar_Mass(%d) Distance(%d) ", j,j,j,j,j,j,j,j);
    printf("Halo_ID(%d) Stellar_Mass(%d) Distance(%d) Z_Diff(%d) ", j,j,j,j); 
  }
  for (int64_t i=0; i<CYLINDER_NUM_RADII; i++) {
	  for (int64_t j=0; j<CYLINDER_BINS; j++) {
		  printf("Cylinder_%gMpc_%g_%gkms ",CYLINDER_RADII[i],CYLINDER_BIN_SPACING*j,CYLINDER_BIN_SPACING*(j+1));
	  }
  }
  printf("\n");


  //find neighbors / counts in cylinders
  for (int64_t i=0; i<num_halos; i++) {  //iterate through all the halos
    if (halos[i].obs_sm > MAX_TARGET_STELLAR_MASS) continue; //Start at high mass end, then adjust.
    find_neighbors(halos+i, tree, results); 
    printf(" ");
    count_cylinders(halos+i, tree, results);
    printf("\n");
  }

  return 0;
}

//This is just adjusting distances between two points for periodic boundary
double periodic_distance(struct catalog_halo *h1, struct catalog_halo *h2) {  
  double s=0;
  for (int64_t i=0; i<2; i++) {
    double x = fabs(h1->pos[i]-h2->pos[i]);
    if (x>BOX_SIZE/2.0) x = BOX_SIZE-x;
    s += x*x;
  }
  return sqrt(s);
}

//Same as above but for redshift space 
double redshift_space_distance(struct catalog_halo *h1, struct catalog_halo *h2) {
  double s = fabs(h1->pos[2]-h2->pos[2]);
  if (s>BOX_SIZE/2.0) s = BOX_SIZE-s;
  return (H*s);
}


int64_t results_to_neighbor_list(struct catalog_halo *h, struct fast3tree_results *r, struct neighbor **n) {
  int64_t i, num_neighbors = 0;
  for (i=0; i<r->num_points; i++) { //Iterate through all gals returned by search
    if (redshift_space_distance(h, r->points[i])<MAX_NEIGHBOR_VELOCITY_OFFSET &&
	r->points[i]->obs_sm > h->obs_sm*MAX_NEIGHBOR_MASS_OFFSET) {
      check_realloc_every(*n, sizeof(struct neighbor), num_neighbors, 1000);  //Find address with enough space for all the neighbors
      n[0][num_neighbors].h = r->points[i]; 
      n[0][num_neighbors].r = periodic_distance(h, r->points[i]);
      n[0][num_neighbors].z = redshift_space_distance(h,r->points[i]);
      num_neighbors++;
    }
  }
  return num_neighbors;
}

//See which neighbors is closer  (*Update we need to sort them by distance)
int neighbor_distance_compare(const void *a, const void *b) {
  const struct neighbor *c = a;
  const struct neighbor *d = b;
  if (c->r < d->r) return -1;
  if (c->r > d->r) return 1;
  return 0; //If at the same distance
}

double center_mass_dist(struct neighbor *n) {
  double total_mass = 0;
  double mx = 0;
  double my = 0;
  for (int64_t i=0; i<NEIGHBOR_COUNT+1; i++) {
    total_mass += n[i].h->obs_sm;
    double dx = (n[i].h->pos[0]-n[0].h->pos[0]);
    if (dx>BOX_SIZE/2.0) dx -= BOX_SIZE;
    if (dx<-BOX_SIZE/2.0) dx += BOX_SIZE;
    mx += n[i].h->obs_sm * dx;
    double dy = (n[i].h->pos[1]-n[0].h->pos[1]);
    if (dy>BOX_SIZE/2.0) dy -= BOX_SIZE;
    if (dy<-BOX_SIZE/2.0) dy += BOX_SIZE;
    my += n[i].h->obs_sm * dy;
  }
  double Mx = mx / total_mass;
  double My = my / total_mass;
  return sqrt(Mx*Mx + My*My);
}


//Find neighbors and put them in the right order
void find_neighbors(struct catalog_halo *h, struct fast3tree *t, struct fast3tree_results *r) {
  int64_t num_neighbors = 0;
  double radius=0.125;
  struct neighbor *n = NULL;
  while (num_neighbors<NEIGHBOR_COUNT+1 && radius<BOX_SIZE/2.0) {
    fast3tree_find_sphere_periodic(t, r, h->pos, radius);
    num_neighbors = results_to_neighbor_list(h, r, &n);
    radius *= 2.0; //expand search radius till you get 50 neighbors (or you've used the whole box)
  }
  qsort(n, num_neighbors, sizeof(struct neighbor), neighbor_distance_compare); //Sort neighbors by distance
  //double r_cm  = center_mass_dist(n);
  printf("%"PRId64" %e %e %e %e %e ", n[0].h->id, n[0].h->pos[0], n[0].h->pos[1], n[0].h->pos[2], n[0].h->mp, n[0].h->obs_sm);
  for (int64_t i=1; i<NEIGHBOR_COUNT+1; i++) {
    if (i<num_neighbors)
      printf("%"PRId64" %e %e %e ", n[i].h->id, n[i].h->obs_sm, n[i].r, n[i].z);
    else
      printf("-1 -1 -1 -1");
    if (i<=NEIGHBOR_COUNT) printf(" ");
  }
}


//Put the cylinder data into redshift bins
void results_to_cylinder(struct catalog_halo *h, struct fast3tree_results *r, int64_t *counts) {
  for (int64_t i=0; i<CYLINDER_BINS; i++) counts[i] = 0;
  for (int64_t i=0; i<r->num_points; i++) {
    if (r->points[i]->obs_sm > h->obs_sm*MAX_NEIGHBOR_MASS_OFFSET) {
      double s = redshift_space_distance(h, r->points[i]);
      int64_t bin = s/CYLINDER_BIN_SPACING;
      if (bin < CYLINDER_BINS) counts[bin]++; //add count to bin if neighbor falls within that bin
    }
  }
}

void count_cylinders(struct catalog_halo *h, struct fast3tree *t, struct fast3tree_results *r) {
  int64_t counts[CYLINDER_BINS];
  for (int64_t i=0; i<CYLINDER_NUM_RADII; i++) {
    fast3tree_find_sphere_periodic(t, r, h->pos, CYLINDER_RADII[i]);
    results_to_cylinder(h, r, counts);
    //h->r_cm = center_mass_dist(n,counts);
    for (int64_t j=0; j<CYLINDER_BINS; j++) {
      printf("%" PRId64, counts[j]);
      if ((j<(CYLINDER_BINS-1)) || (i<(CYLINDER_NUM_RADII-1)))
	printf(" ");
    }
  }
}


