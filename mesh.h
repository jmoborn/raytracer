#ifndef __MESH_H__
#define __MESH_H__

#include <iostream>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include "vec4.h"
#include "material.h"
#include "object.h"

class face {
public:
	face(std::vector<int> pts, std::vector<int> txs, std::vector<int> nms);
	~face(){};

	std::vector<int> pnts;
	std::vector<int> nrms;
	std::vector<int> txts;

private:
	bool isQuad;
};

class mesh : public object {
public:
	mesh();
	mesh(std::string filepath, material m=material());
	~mesh(){};
	bool intersect(ray& r);
	vec4 get_normal(const vec4& p);

	std::vector<vec4> verts;
	std::vector<vec4> norms;
	std::vector<face> faces;
	std::vector<vec2> texts;

private:
	void readobj(std::string& filepath);

	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;

};

#endif