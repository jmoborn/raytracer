#include "mesh.h"

face::face(std::vector<int> pts, std::vector<int> txs, std::vector<int> nms)
{
	this->pnts = pts;
	this->txts = txs;
	this->nrms = nms;

	// xmin = ymin = zmin = std::numeric_limits<double>::infinity();
	// xmax = ymax = zmax = -std::numeric_limits<double>::infinity();
	
	// for(int i=0; i<pnts.size(); i++)
	// {
	// 	if(pnts[0]>xmax) xmax = pnts[0];
	// 	if(pnts[1]>ymax) ymax = pnts[1];
	// 	if(pnts[2]>zmax) zmax = pnts[2];
	// 	if(pnts[0]<xmin) xmin = pnts[0];
	// 	if(pnts[1]<ymin) ymin = pnts[1];
	// 	if(pnts[2]<zmin) zmin = pnts[2];
	// }

}

mesh::mesh(){}

mesh::mesh(std::string filepath, int material_index)
{
	xmin = ymin = zmin = std::numeric_limits<double>::infinity();
	xmax = ymax = zmax = -std::numeric_limits<double>::infinity();
	readobj(filepath);
	// this->shader = m;
	this->mtl_idx = material_index;
	std::vector<face*> prims;
	for(int i=0; i<faces.size(); i++)
	{
		prims.push_back(&faces[i]);
	}
	this->hierarchy.build(prims, xmin, xmax, ymin, ymax, zmin, zmax);
}

bool mesh::intersect(ray& r)
{
	// if (r.debug) std::cout << "intersecting new mesh: " << std::endl;
	prim_hit face_hit = intersect_bvh(r, hierarchy.root, 0);

	if (face_hit.prim!=NULL)
	{
		vec4 N = norms[face_hit.prim->nrms[1]]*face_hit.u + norms[face_hit.prim->nrms[2]]*face_hit.v + norms[face_hit.prim->nrms[0]]*(1 - face_hit.u - face_hit.v);
		vec2 U = texts[face_hit.prim->txts[1]]*face_hit.u + texts[face_hit.prim->txts[2]]*face_hit.v + texts[face_hit.prim->txts[0]]*(1 - face_hit.u - face_hit.v);
		// std::cout << texts[faces[i].txts[1]].u << ", " << texts[faces[i].txts[1]].v << std::endl;
		// std::cout << U.u << ", " << U.v << std::endl; 
		r.t = face_hit.dist;
		r.hit_norm = N;
		r.hit_uv = U;
		r.hit_mtl = mtl_idx;
		// r.hit_color = diffuse(U);
		// closest_face = i;
		return true;
	}
	else
	{
		return false;
	}

	// int closest_face = -1;
	// //check bounding box first
	// //for x
	// double t0 = 0;
	// double t1 = std::numeric_limits<double>::infinity();
	// double invRayDir = 1.0/r.d.x;
	// double tNear = (xmin - r.o.x)*invRayDir;
	// double tFar  = (xmax - r.o.x)*invRayDir;
	// if(tNear > tFar) std::swap(tNear, tFar);
	// t0 = tNear > t0 ? tNear : t0;
	// t1 = tFar < t1 ? tFar : t1;
	// if(t0 > t1) return false;
	// //for y
	// invRayDir = 1.0/r.d.y;
	// tNear = (ymin - r.o.y)*invRayDir;
	// tFar  = (ymax - r.o.y)*invRayDir;
	// if(tNear > tFar) std::swap(tNear, tFar);
	// t0 = tNear > t0 ? tNear : t0;
	// t1 = tFar < t1 ? tFar : t1;
	// if(t0 > t1) return false;
	// //for z
	// invRayDir = 1.0/r.d.z;
	// tNear = (zmin - r.o.z)*invRayDir;
	// tFar  = (zmax - r.o.z)*invRayDir;
	// if(tNear > tFar) std::swap(tNear, tFar);
	// t0 = tNear > t0 ? tNear : t0;
	// t1 = tFar < t1 ? tFar : t1;
	// if(t0 > t1) return false;

	// double min_dist = std::numeric_limits<double>::infinity();
	// for(int i=0; i<faces.size(); i++)
	// {
	// 	vec4 A = verts[faces[i].pnts[0]];
	// 	vec4 B = verts[faces[i].pnts[1]];
	// 	vec4 C = verts[faces[i].pnts[2]];

	// 	vec4 edge1 = B - A;
	// 	vec4 edge2 = C - A;
	// 	vec4 pvec = r.d.cross(edge2);
	// 	double det = edge1.dot(pvec);
	// 	if(det == 0) continue;
	// 	double invDet = 1/det;
	// 	vec4 tvec = r.o - A;
	// 	double u = tvec.dot(pvec) * invDet;
	// 	if(u<0 || u>1) continue;
	// 	vec4 qvec = tvec.cross(edge1);
	// 	double v = r.d.dot(qvec) * invDet;
	// 	if(v<0 || (u + v) >1) continue;
	// 	double hit = (edge2.dot(qvec) * invDet);
	// 	if(hit<0) continue;

	// 	if(hit<r.t)
	// 	{
	// 		vec4 N = norms[faces[i].nrms[1]]*u + norms[faces[i].nrms[2]]*v + norms[faces[i].nrms[0]]*(1 - u - v);
	// 		vec2 U = texts[faces[i].txts[1]]*u + texts[faces[i].txts[2]]*v + texts[faces[i].txts[0]]*(1 - u - v);
	// 		// std::cout << texts[faces[i].txts[1]].u << ", " << texts[faces[i].txts[1]].v << std::endl;
	// 		// std::cout << U.u << ", " << U.v << std::endl; 
	// 		min_dist = hit;
	// 		r.t = hit;
	// 		r.hit_norm = N;
	// 		r.hit_uv = U;
	// 		r.hit_color = diffuse(U);
	// 		closest_face = i;
	// 	}
	// }
	// if(closest_face!=-1) return true;
	// else return false;
}

prim_hit mesh::intersect_bvh(ray& r, bvh_node<face>* node, int level)
{
	// if (r.debug)
	// {
	// 	std::cout << level << std::endl;
	// 	std::cout << "  min: " << node->min[0] << " " << node->min[1] << " " << node->min[2] << std::endl;
	// 	std::cout << "  max: " << node->max[0] << " " << node->max[1] << " " << node->max[2] << std::endl;
	// 	if(level==2&&node->max[1]==0.15)
	// 	{
	// 		std::cout << "WHAT THE FUCK IS HAPPENING" << std::endl;
	// 		std::cout << node->prims.size() << std::endl;
	// 	}
	// }
	prim_hit min_hit;
	min_hit.dist = std::numeric_limits<double>::infinity();
	min_hit.prim = NULL;
	//for x
	double t0 = 0;
	double t1 = std::numeric_limits<double>::infinity();
	double invRayDir = 1.0/r.d.x;
	double tNear = (node->min[0] - r.o.x)*invRayDir;
	double tFar  = (node->max[0] - r.o.x)*invRayDir;
	if(tNear > tFar) std::swap(tNear, tFar);
	t0 = tNear > t0 ? tNear : t0;
	t1 = tFar < t1 ? tFar : t1;
	if(t0 > t1) return min_hit;
	//for y
	invRayDir = 1.0/r.d.y;
	tNear = (node->min[1] - r.o.y)*invRayDir;
	tFar  = (node->max[1] - r.o.y)*invRayDir;
	if(tNear > tFar) std::swap(tNear, tFar);
	t0 = tNear > t0 ? tNear : t0;
	t1 = tFar < t1 ? tFar : t1;
	if(t0 > t1) return min_hit;
	//for z
	invRayDir = 1.0/r.d.z;
	tNear = (node->min[2] - r.o.z)*invRayDir;
	tFar  = (node->max[2] - r.o.z)*invRayDir;
	if(tNear > tFar) std::swap(tNear, tFar);
	t0 = tNear > t0 ? tNear : t0;
	t1 = tFar < t1 ? tFar : t1;
	if(t0 > t1) return min_hit;

	// if(r.debug&&level==2&&node->max[1]==0.15) std::cout << "PASSED THE BBOX INTERSECTION BITCHES" << std::endl;

	if(node->left==NULL||node->right==NULL)
	{
		// if (r.debug) 
		// {
		// 	std::cout << "now intersecting actual primitives" << std::endl;
		// 	std::cout << "        " << node->prims.size() << std::endl;
		// 	// std::cout << "        " << 
		// }
		for(int i=0; i<node->prims.size(); i++)
		{
			vec4 A = verts[node->prims[i]->pnts[0]];
			vec4 B = verts[node->prims[i]->pnts[1]];
			vec4 C = verts[node->prims[i]->pnts[2]];

			vec4 edge1 = B - A;
			vec4 edge2 = C - A;
			vec4 pvec = r.d.cross(edge2);
			double det = edge1.dot(pvec);
			if(det == 0) continue;
			double invDet = 1/det;
			vec4 tvec = r.o - A;
			double u = tvec.dot(pvec) * invDet;
			if(u<0 || u>1) continue;
			vec4 qvec = tvec.cross(edge1);
			double v = r.d.dot(qvec) * invDet;
			if(v<0 || (u + v) >1) continue;
			double hit = (edge2.dot(qvec) * invDet);
			if(hit<0) continue;

			if(hit<r.t && hit < min_hit.dist)
			{
				min_hit.dist = hit;
				min_hit.prim = node->prims[i];
				min_hit.u = u;
				min_hit.v = v;
			}
		}
		return min_hit;
		// if(min_hit.prim != NULL)
		// {
		// 	return min_hit;
		// }
		// else
		// {
		// 	return NULL;
		//return closest hit
	}
	else
	{
		prim_hit left = intersect_bvh(r, node->left, level+1);
		prim_hit right = intersect_bvh(r, node->right, level+1);
		if (left.prim==NULL) return right;
		else if(right.prim==NULL) return left;
		else if (left.dist < right.dist) return left;
		else return right;
	}
}

vec4 mesh::get_normal(const vec4& p)
{
	return vec4(0, 1, 0);
}

void mesh::readobj(std::string& filepath)
{
	FILE *objfile = fopen(filepath.c_str(), "r");
	if(objfile==NULL)
	{
		std::cout << "Error loading file " << filepath << std::endl;
	}
	char line[128];
	while (fgets(line, sizeof(line), objfile))
	{
		std::stringstream ss;
		ss << line;
		std::string tok;
		ss >> tok;
		if (tok=="v")
		{
			std::string sx, sy, sz;
			ss >> sx >> sy >> sz;
			double px = atof(sx.c_str());
			double py = atof(sy.c_str());
			double pz = atof(sz.c_str());
			if(px>xmax) xmax = px;
			if(py>ymax) ymax = py;
			if(pz>zmax) zmax = pz;
			if(px<xmin) xmin = px;
			if(py<ymin) ymin = py;
			if(pz<zmin) zmin = pz;
			vec4 v(px, py, pz);
			verts.push_back(v);
		}
		else if (tok=="vt")
		{
			std::string sx, sy;
			ss >> sx >> sy;
			vec2 v(atof(sx.c_str()), atof(sy.c_str()));
			texts.push_back(v);
		}
		else if (tok=="vn")
		{
			std::string nx, ny, nz;
			ss >> nx >> ny >> nz;
			vec4 v(atof(nx.c_str()), atof(ny.c_str()), atof(nz.c_str()));
			norms.push_back(v);
		}
		else if (tok=="f")
		{
			std::vector<int> face_pts;
			std::vector<int> face_txs;
			std::vector<int> face_nms;
			std::string face_str;
			ss >> face_str;
			while(!ss.eof())
			{
				int pt, vt, nm;
				sscanf(face_str.c_str(), "%d/%d/%d", &pt, &vt, &nm);
				face_pts.push_back(pt-1);
				face_txs.push_back(vt-1);
				face_nms.push_back(nm-1);
				ss >> face_str;
			}
			face f(face_pts, face_txs, face_nms);
			faces.push_back(f);
		}
	}

	std::cout << "verts: " << verts.size() << std::endl;
	std::cout << "norms: " << norms.size() << std::endl;
	fclose(objfile);

	//set bbox for each polygon
	for(int i=0; i<faces.size(); i++)
	{
		// std::cout << "setting bbox for face " << i << std::endl;
		faces[i].xmin = faces[i].ymin = faces[i].zmin = std::numeric_limits<double>::infinity();
		faces[i].xmax = faces[i].ymax = faces[i].zmax = -std::numeric_limits<double>::infinity();
		double epsilon = 0.0001;
		for(int j=0; j<faces[i].pnts.size(); j++)
		{
			if(verts[faces[i].pnts[j]].x>faces[i].xmax) faces[i].xmax = verts[faces[i].pnts[j]].x;
			if(verts[faces[i].pnts[j]].y>faces[i].ymax) faces[i].ymax = verts[faces[i].pnts[j]].y;
			if(verts[faces[i].pnts[j]].z>faces[i].zmax) faces[i].zmax = verts[faces[i].pnts[j]].z;
			if(verts[faces[i].pnts[j]].x<faces[i].xmin) faces[i].xmin = verts[faces[i].pnts[j]].x;
			if(verts[faces[i].pnts[j]].y<faces[i].ymin) faces[i].ymin = verts[faces[i].pnts[j]].y;
			if(verts[faces[i].pnts[j]].z<faces[i].zmin) faces[i].zmin = verts[faces[i].pnts[j]].z;
		}
		// faces[i].xmax += epsilon;
		// faces[i].ymax += epsilon;
		// faces[i].zmax += epsilon;
		// faces[i].xmin -= epsilon;
		// faces[i].ymin -= epsilon;
		// faces[i].zmin -= epsilon;
	}
}