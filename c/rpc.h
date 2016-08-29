// rational polynomial coefficient stuff

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>


// rational polynomial coefficients (and those of the inverse model)
struct rpc {
	double numx[20];
	double denx[20];
	double numy[20];
	double deny[20];
	double scale[3], offset[3];

	double inumx[20];
	double idenx[20];
	double inumy[20];
	double ideny[20];
	double iscale[3], ioffset[3];

	double dmval[4];
	double imval[4];
};


// read an XML file specifying an RPC model
void read_rpc_file_xml(struct rpc *p, char *filename);

void print_rpc(FILE *f, struct rpc *p, char *n);

// evaluate the direct rpc model
void eval_rpc(double *result,
		struct rpc *p, double x, double y, double z);

// evaluate the inverse rpc model
void eval_rpci(double *result,
		struct rpc *p, double x, double y, double z);

// evaluate an epipolar correspondence
static void eval_rpc_pair(double xprime[2],
		struct rpc *a, struct rpc *b,
		double x, double y, double z);


// convert (long,lat,h) to ECEF coord sys. (X,Y,Z)
void geotedic_to_ECEF(double lg, double lt, double h, 
			    double *X, double *Y, double *Z);

// given (X,Y,Z) ECEF coord, computes the alt above the WGS 84 ellipsoid
double get_altitude_from_ECEF(double X, double Y, double Z);

// define the line passing by point 's' with direction vector 'v'
typedef struct 
{
	double s[3]; // a point in space
	double v[3]; // a direction vector
} SV;

// distance between a 3D point P and a line (defined by its
// normalized direction vector V and passing through a 3D point S) 
double dist_line_point3D(double *V,double *S,double *P);

// find closest 3D point from from a set of 3D lines
void find_point_opt(SV *sv_tab, int N, bool *take,
		double *point_opt,double *outerr);

// compute the height of a point given its location inside two images
// geometric solution
double rpc_height_geo(struct rpc *rpc_list, 
		double ** q_list, int *N, 
		bool findConsensus, double thres, double *outerr);
