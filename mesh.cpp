#include "mesh.h"

face::face(std::vector<int> pts, std::vector<int> txs, std::vector<int> nms)
{
	this->pnts = pts;
	this->txts = txs;
	this->nrms = nms;

	// for(int i=0; i<pts.size(); i++)
	// {
	// 	this->pnts[i] = pts[i];
	// }
	// for(int i=0; i<txs.size(); i++)
	// {
	// 	this->txts[i] = txs[i];
	// }
}

mesh::mesh(){}

mesh::mesh(std::string filepath, material m)
{
	xmin = ymin = zmin = std::numeric_limits<double>::infinity();
	xmax = ymax = zmax = -std::numeric_limits<double>::infinity();
	readobj(filepath);
	this->shader = m;
}

bool mesh::intersect(ray& r)
{
	
	int closest_face = -1;
	//check bounding box first
	//for x
	double t0 = 0;
	double t1 = std::numeric_limits<double>::infinity();
	double invRayDir = 1.0/r.d.x;
	double tNear = (xmin - r.o.x)*invRayDir;
	double tFar  = (xmax - r.o.x)*invRayDir;
	if(tNear > tFar) std::swap(tNear, tFar);
	t0 = tNear > t0 ? tNear : t0;
	t1 = tFar < t1 ? tFar : t1;
	if(t0 > t1) return false;
	//for y
	invRayDir = 1.0/r.d.y;
	tNear = (ymin - r.o.y)*invRayDir;
	tFar  = (ymax - r.o.y)*invRayDir;
	if(tNear > tFar) std::swap(tNear, tFar);
	t0 = tNear > t0 ? tNear : t0;
	t1 = tFar < t1 ? tFar : t1;
	if(t0 > t1) return false;
	//for z
	invRayDir = 1.0/r.d.z;
	tNear = (zmin - r.o.z)*invRayDir;
	tFar  = (zmax - r.o.z)*invRayDir;
	if(tNear > tFar) std::swap(tNear, tFar);
	t0 = tNear > t0 ? tNear : t0;
	t1 = tFar < t1 ? tFar : t1;
	if(t0 > t1) return false;

	double min_dist = std::numeric_limits<double>::infinity();
	for(int i=0; i<faces.size(); i++)
	{
		vec4 A = verts[faces[i].pnts[0]];
		vec4 B = verts[faces[i].pnts[1]];
		vec4 C = verts[faces[i].pnts[2]];

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

		if(hit<r.t)
		{
			// vec4 N = edge1.cross(edge2);
			// vec4 N = (norms[faces[i].nrms[0]] + norms[faces[i].nrms[1]] + norms[faces[i].nrms[2]]) * (1.0/3.0);
			// std::cout << "origin " << A.x << " " << A.y << " " << A.z << std::endl;
			// std::cout << "u axis " << B.x << " " << B.y << " " << B.z << std::endl;
			// std::cout << "v axis " << C.x << " " << C.y << " " << C.z << std::endl;
			// std::cout << "barycentric u: " << u << std::endl;
			// std::cout << "barycentric v: " << v << std::endl;
			// std::cout << "normal x" << norms[faces[i].nrms[0]].x << " " << norms[faces[i].nrms[0]].y << " " << norms[faces[i].nrms[0]].z << std::endl;
			// std::cout << "normal y" << norms[faces[i].nrms[1]].x << " " << norms[faces[i].nrms[1]].y << " " << norms[faces[i].nrms[1]].z << std::endl;
			// std::cout << "normal x" << norms[faces[i].nrms[2]].x << " " << norms[faces[i].nrms[2]].y << " " << norms[faces[i].nrms[2]].z << std::endl;
			// vec4 N = norms[faces[i].nrms[0]]*u + norms[faces[i].nrms[1]]*v + norms[faces[i].nrms[2]]*hit;
			vec4 N = norms[faces[i].nrms[1]]*u + norms[faces[i].nrms[2]]*v + norms[faces[i].nrms[0]]*(1 - u - v);
			min_dist = hit;
			r.t = hit;
			r.hit_norm = N;
			r.hit_color = diffuse();
			closest_face = i;
		}
	}
	if(closest_face!=-1) return true;
	else return false;
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
		printf("Error loading file\n");
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
}