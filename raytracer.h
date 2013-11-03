#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>

#include "mat4.h"
#include "vec4.h"
#include "sphere.h"
#include "mesh.h"
#include "pixelmap.h"

class raytracer {
public:
	raytracer();
	~raytracer();

	std::vector<object*> objs;

	//constants
	static const double PI = 3.141592653589;
private:
	vec4 v;
	mat4 m;
	mesh g;
};