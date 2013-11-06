#include "mesh.h"

face::face(std::vector<int> pts, std::vector<int> txs)
{
	this->pnts = pts;
	this->txts = txs;

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

mesh::mesh(std::string filepath, vec4 c)
{
	readobj(filepath);
	color = c;
}

bool mesh::intersect(ray& r)
{
	double min_dist = std::numeric_limits<double>::infinity();
	int closest_face = -1;
	for(int i=0; i<faces.size(); i++)
	{
		vec4 A = verts[faces[i].pnts[0]];
		vec4 B = verts[faces[i].pnts[1]];
		vec4 C = verts[faces[i].pnts[2]];
		// vec4 A(-.5, 0, -1.9);
		// vec4 B(.25, -.4, -2);
		// vec4 C(-.5, 0, -2);
		// vec4 B(.25, -.4, -2);
		// vec4 C(.5, .25, -2);
		// vec4 AB = B - A;
		// vec4 AC = C - A;
		// vec4 N = AB.cross(AC);
		// double n_dot_r = N.dot(r.d);
		// if(n_dot_r==0) continue;
		// double D = N.dot(A);
		// double hit = -(N.dot(r.o) + D)/n_dot_r;
		// if(hit<0) continue;
		// double tmp = r.t;
		// r.t = hit;
		// vec4 P = r.end();
		// r.t = tmp;
		// // N.normalize();

		// vec4 AP = P - A;
		// vec4 ABxAP = AB.cross(AP);
		// double v_num = N.dot(ABxAP);
		// if(v_num<0) continue;

		// vec4 BP = P - B;
		// vec4 BC = C - B;
		// vec4 BCxBP = BC.cross(BP);
		// if(N.dot(BCxBP)<0) continue;

		// vec4 CP = P - C;
		// vec4 CA = A - C;
		// vec4 CAxCP = CA.cross(CP);
		// float u_num = N.dot(CAxCP);
		// if(u_num < 0) continue;

		// double nom = N.dot(N);
		// double u = u_num/nom;
		// double v = v_num/nom;

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
		// vec4 N = edge1.cross(edge2);
		// N.normalize();
		// double D = N.dot(A);
		// if(D>=0) continue;
		// double hit = -(N.dot(r.o) + D)/(N.dot(r.d));

		if(hit<r.t)
		{
			vec4 N = edge1.cross(edge2);
			N.normalize();
			// std::cout << "we have a hit!" << std::endl;
			min_dist = hit;
			r.t = hit;
			N.normalize();
			r.hit_norm = N;
			r.hit_color = diffuse();
			closest_face = i;
		}
	}
	// std::cout << min_dist << std::endl;
	if(closest_face!=-1) return true;
	else return false;
}

vec4 mesh::get_normal(const vec4& p)
{
	return vec4(0, 1, 0);
}

vec4 mesh::diffuse()
{
	return color;
}

vec4 mesh::reflect()
{
	return vec4(1, 1, 1);
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
			vec4 v(vec4(atof(sx.c_str()), atof(sy.c_str()), atof(sz.c_str())));
			verts.push_back(v);
		}
		else if (tok=="vt")
		{
			std::string sx, sy;
			ss >> sx >> sy;
			vec2 v(vec2(atof(sx.c_str()), atof(sy.c_str())));
			texts.push_back(v);
		}
		else if (tok=="n")
		{
			// std::string num;
			// ss >> num;
			// double x = atof(num.c_str());
			// ss >> num;
			// double y = atof(num.c_str());
			// ss >> num;
			// double z = atof(num.c_str());
		}
		else if (tok=="f")
		{
			std::vector<int> face_pts;
			std::vector<int> face_txs;
			//add normals
			std::string face_str;
			ss >> face_str;
			while(!ss.eof())
			{
				int pt, vt, nm;
				sscanf(face_str.c_str(), "%d/%d/%d", &pt, &vt, &nm);
				face_pts.push_back(pt-1);
				face_txs.push_back(vt-1);
				//add normal
				ss >> face_str;
			}
			face f(face_pts, face_txs);
			faces.push_back(f);
		}
	}
	fclose(objfile);
}

// int main()
// {
// 	mesh m("crayon_box.obj");
// 	for(int i=0; i<m.faces.size(); i++)
// 	{
// 		for(int j=0; j<m.faces[i].pnts.size(); j++)
// 		{
// 			std::cout << m.faces[i].pnts[j] << " ";
// 		}
// 		std::cout << std::endl;
// 	}
// }