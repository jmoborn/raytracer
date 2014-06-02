#include "mat4.h"

mat4::mat4()
{
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{
			if(i==j)
			{
				data[i][j] = 1.0;
			}
			else
			{
				data[i][j] = 0.0;
			}
		}
	}
}

void mat4::translate(double tx, double ty, double tz)
{
	mat4 tr;
	tr.data[0][3] = tx;
	tr.data[1][3] = ty;
	tr.data[2][3] = tz;
	(*this) *= tr;
}
void mat4::scale(double sx, double sy, double sz)
{
	mat4 ts;
	ts.data[0][0] = sx;
	ts.data[1][1] = sy;
	ts.data[2][2] = sz;
	(*this) *= ts;
}
void mat4::rotateX(double theta)
{
	double rad = theta * M_PI / 180.0;
	mat4 rx;
	rx.data[1][1] = cos(rad);
	rx.data[1][2] = -sin(rad);
	rx.data[2][1] = sin(rad);
	rx.data[2][2] = cos(rad);
	(*this) *= rx;
}
void mat4::rotateY(double theta)
{
	double rad = theta * M_PI / 180.0;
	mat4 ry;
	ry.data[0][0] = cos(rad);
	ry.data[0][2] = sin(rad);
	ry.data[2][0] = -sin(rad);
	ry.data[2][2] = cos(rad);
	(*this) *= ry;
}
void mat4::rotateZ(double theta)
{
	double rad = theta * M_PI / 180.0;
	mat4 rz;
	rz.data[0][0] = cos(rad);
	rz.data[0][1] = -sin(rad);
	rz.data[1][0] = sin(rad);
	rz.data[1][1] = cos(rad);
	(*this) *= rz;
}

double & mat4::operator()(const int i, const int j)
{
	return data[i][j];
}

double mat4::operator()(const int i, const int j) const
{
	return data[i][j];
}

mat4 & mat4::operator*=(const mat4& m)
{
	double newdata[4][4];
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{
			double sum = 0;
			for(int k=0; k<4; k++)
			{
				sum += data[i][k]*m.data[k][j];
			}
			newdata[i][j] = sum;
		}
	}
	// std::copy(newdata, newdata + 4, data);
	std::copy(&newdata[0][0], &newdata[0][0] + 4*4, &data[0][0]);
	return *this;
}

mat4 & mat4::operator*=(const double c)
{
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{
			data[i][j] *= c;
		}
	}
	return *this;
}

mat4 operator*(const mat4& m1, mat4& m2)
{
	mat4 a;
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{
			double sum = 0;
			for(int k =0; k<4; k++)
			{
				sum += m1.data[i][k]*m2.data[k][j];
			}
			a.data[i][j] = sum;
		}
	}
	return a;
}
mat4 operator*(const mat4&m, const double c)
{
	mat4 a;
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{
			a.data[i][j] = m.data[i][j] * c;
		}
	}
	return a;
}

// int main()
// {
// 	mat4 m;
// 	m(2,0) = 1;
// 	m(3,1) = 1;
// 	m(0,1) = 1;
// 	mat4 w;
// 	w(0,0) = 3;
// 	w(1,1) = 2;
// 	w(3,3) = 0;
// 	mat4 a;
// 	a = m * w;
// 	//m *= 2;
// 	for(int i=0; i<4; i++)
// 	{
// 		for(int j=0; j<4; j++)
// 		{
// 			std::cout << a(i,j) << " ";
// 		}
// 		std::cout << std::endl;
// 	}
// }
