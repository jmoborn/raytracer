#ifndef __RAYTRACER_H__
#define __RAYTRACER_H__

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <limits>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mat4.h"
#include "vec4.h"
#include "material.h"
#include "sphere.h"
#include "mesh.h"
#include "pixelmap.h"

class raytracer {
public:
	raytracer(std::string scenefile);
	~raytracer();

	void render_image(std::string& outfile);
	int intersect_scene(ray &v);
	int intersect_shadow(ray &v);
	vec4 shade(ray& v);
	vec4 trace_path(ray& v, int depth);
	void load_scene(std::string& scenefile);
	double randd();
	double randd_negative();
	int randi(int lo, int hi);
	void init_rand(int);

	static const double ray_tolerance;


private:
	vec4 read_vector(std::stringstream& ss);

	std::vector<material*> mtls;
	std::vector<object*> objs;
	std::vector<sphere*> lights;

	int image_width;
	int image_height;
	int samples1D;
	double fov;	

	int max_depth;
	double ambience;
	unsigned int *rand_seed;
};

#endif