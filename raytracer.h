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
#include <time.h>

#include "mat4.h"
#include "vec4.h"
#include "material.h"
#include "sphere.h"
#include "mesh.h"
#include "pixelmap.h"

class raytracer {
public:
	raytracer();
	~raytracer();

	vec4 shade(ray& v, int depth, int self=-1);
	vec4 trace_ray(ray& v, int depth, int self);
	double trace_shadow(ray& v, int self, int light);
	void load_scene(std::string& scenefile);
	double randd();
	double randd_negative();

	std::vector<material*> mtls;
	std::vector<object*> objs;
	std::vector<sphere*> lights;

	int max_depth;
	int shadow_samples;
	double ambience;
	double ray_tolerance;
	double samples;

	//constants
	static const double PI = 3.141592653589;

private:
	vec4 read_vector(std::stringstream& ss);
};

#endif