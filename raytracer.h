#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <limits>

#include "mat4.h"
#include "vec4.h"
#include "sphere.h"
#include "mesh.h"
#include "pixelmap.h"

class raytracer {
public:
	raytracer();
	~raytracer();

	vec4 shade(ray& v, int self=-1);
	vec4 trace_ray(ray& v, int depth, int self);
	bool trace_shadow(ray& v, int self);

	std::vector<object*> objs;
	std::vector<sphere*> lights;

	int max_depth;

	//constants
	static const double PI = 3.141592653589;
};