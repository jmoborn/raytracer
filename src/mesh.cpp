#include "mesh.h"

face::face(std::vector<int> pts, std::vector<int> txs, std::vector<int> nms)
{
	this->pnts = pts;
	this->txts = txs;
	this->nrms = nms;
}

mesh::mesh(){}

mesh::mesh(std::string filepath, int material_index)
{
	xmin = ymin = zmin = std::numeric_limits<double>::infinity();
	xmax = ymax = zmax = -std::numeric_limits<double>::infinity();
	readobj(filepath);
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
	prim_hit face_hit = intersect_bvh(r, hierarchy.root, 0);

	if (face_hit.prim!=NULL)
	{
		vec4 N = norms[face_hit.prim->nrms[1]]*face_hit.u + norms[face_hit.prim->nrms[2]]*face_hit.v + norms[face_hit.prim->nrms[0]]*(1 - face_hit.u - face_hit.v);
		vec2 U = texts[face_hit.prim->txts[1]]*face_hit.u + texts[face_hit.prim->txts[2]]*face_hit.v + texts[face_hit.prim->txts[0]]*(1 - face_hit.u - face_hit.v);
		r.t = face_hit.dist;
		r.hit_norm = N;
		r.hit_uv = U;
		r.hit_mtl = mtl_idx;
		return true;
	}
	else
	{
		return false;
	}
}

prim_hit mesh::intersect_bvh(ray& r, bvh_node<face>* node, int level)
{
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

	if(node->left==NULL||node->right==NULL)
	{

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

void mesh::readobj(std::string& filepath)
{
	FILE *objfile = fopen(filepath.c_str(), "r");
	if(objfile==NULL)
	{
		std::cout << "Error loading file " << filepath << std::endl;
		return;
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

	fclose(objfile);

	//set bbox for each polygon
	for(int i=0; i<faces.size(); i++)
	{
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
	}
}

void mesh::print_info(std::ostream &out)
{
	out << "mesh" << std::endl;
	out << "  vertices: " << verts.size() << std::endl;
	out << "  faces: " << faces.size() << std::endl;
	out << "  material index: " << mtl_idx << std::endl;
}
