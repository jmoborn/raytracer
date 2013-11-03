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

mesh::mesh(std::string filepath)
{
	readobj(filepath);
}

bool mesh::intersect(ray& r)
{
	return false;
}

vec4 mesh::get_normal(const vec4& p)
{
	return vec4(0, 1, 0);
}

vec4 mesh::get_color()
{
	return vec4(1, 0, 0);
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