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
#include "bvh.h"

class face {
public:
	face(std::vector<int> pts, std::vector<int> txs, std::vector<int> nms);
	~face(){};

	std::vector<int> pnts;
	std::vector<int> nrms;
	std::vector<int> txts;

	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;

private:
	bool isQuad;
};

class prim_hit {
public:
	prim_hit(){};
	face* prim;
	double dist;
	double u;
	double v;
};

class mesh : public object {
public:
	mesh();
	mesh(std::string filepath, int material_index);
	~mesh(){};
	bool intersect(ray& r);
	prim_hit intersect_bvh(ray& r, bvh_node<face>* node, int level);

	std::vector<vec4> verts;
	std::vector<vec4> norms;
	std::vector<face> faces;
	std::vector<vec2> texts;

	bvh<face> hierarchy;

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