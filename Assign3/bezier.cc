#include "bezier.h"
#include <vector>
#include <iostream>
#include <string>
using namespace std;

void coord::init()
{
	//set everything to zero
	x = 0;
	y = 0;
	z = 0;
	
	r = 0;
	g = 0;
	b = 0;
	
	nx = 0;
	ny = 0;
	nz = 0;
	
	t1 = 0;
	t2 = 0;
}

//set the boolean flags
void setAttributeIndices(string vType, bool &g, bool &c, bool &n, bool &t)
{
	g = false;
	c = g;
	n = g;
	t = g;
	
	for(int i = 0; i < vType.size(); i++)
	{
		if(vType[i] == 'P')
			g = true;
		if(vType[i] == 'N' || vType[i] == 'D')
			n = true;
		if(vType[i] == 'C')
			c = true;
		if(vType[i] == 'T')
			t = true;
	}
}

//set the indices
void setAttributeIndices(string vType,
						 int &vSize,
						 int &g,
						 int &c,
						 int &n,
						 int &t)
{
	int j = 0;
	
	for(int i = 0; i < vType.size(); i++)
	{
		if(vType[i] == 'P')
		{
			g = j;
			vSize += 3;
			j += 3;
		}
		if(vType[i] == 'N' || vType[i] == 'D')
		{
			n = j;
			vSize += 3;
			j += 3;
		}
		if(vType[i] == 'C')
		{
			c = j;
			vSize += 3;
			j += 3;
		}
		if(vType[i] == 'T')
		{
			t = j;
			vSize += 2;
			j += 2;
		}
	}
}

vector<coord> arr2vect(const float arr[], int degree, string vType)
{
	//coordinates to return
	vector<coord> coords;
	//indices for the geometry, color, normal, texture
	//vSize is the vertex size
	int g, c, n, t, vSize;
	vSize = 0;
	//booleans to see if geometric, color, normals, and texture is used
	bool geo = false; bool col = false; bool nor = false; bool tex = false;
	
	//call set_attribute _indices to fill in above variables
	setAttributeIndices(vType, vSize, g, c, n, t);
	setAttributeIndices(vType, geo, col, nor, tex);
	
	//go in and fill information in
	for(int i = 0; i < degree + 1; i++)
	{
		//get to the correct coordinate in array
		int j = vSize*i;
		//create new coordinate object and initialize it
		coord v;
		v.init();
		
		//if geometry is set, add it in
		//do the same for color, normals, and texture
		if(geo)
		{
			v.x = arr[j+g];
			v.y = arr[j+g+1];
			v.z = arr[j+g+2];
		}
		if(col)
		{
			v.r = arr[j+c];
			v.g = arr[j+c+1];
			v.b = arr[j+c+2];
		}
		if(nor)
		{
			v.nx = arr[j+n];
			v.ny = arr[j+n+1];
			v.nz = arr[j+n+2];
		}
		if(tex)
		{
			v.t1 = arr[j+t];
			v.t2 = arr[j+t+1];
		}
		
		//push back c
		coords.push_back(v);
	}
	
	//return coords
	return coords;
}

vector<coord> arr2vect(const float arr[], 
					   int u_degree, 
					   int v_degree, 
					   string vType)
{
	//vector to be returned
	vector<coord> coords;
	
	//insert using above function.  Our limits have changed though
	coords = arr2vect(arr, (u_degree+1)*(v_degree+1), vType);
	
	return coords;
}

coord interpolate(coord a, coord b, float t)
{
	//coordinate to be returned
	coord c;
	c.init();
	
	//calculate the interpolation for each value
	c.x = a.x + t*(b.x - a.x);
	c.y = a.y + t*(b.y - a.y);
	c.z = a.z + t*(b.z - a.z);
	
	//color gets interpolated too
	c.r = a.r + t*(b.r - a.r);
	c.g = a.g + t*(b.g - a.g);
	c.b = a.b + t*(b.b - a.b);
	
	//same with normals
	c.nx = a.nx + t*(b.nx - a.nx);
	c.ny = a.ny + t*(b.ny - a.ny);
	c.nz = a.nz + t*(b.nz - a.nz);
	
	//same with texture
	c.t1 = a.t1 + t*(b.t1 - a.t1);
	c.t2 = a.t2 + t*(b.t2 - a.t2);
	
	//return the coordinate
	return c;
}

coord mult(coord c, float m)
{
	//multiply all values by m
	c.x *= m;
	c.y *= m;
	c.z *= m;
	
	c.r *= m;
	c.g *= m;
	c.b *= m;
	
	c.nx *= m;
	c.ny *= m;
	c.nz *= m;
	
	c.t1 *= m;
	c.t2 *= m;
	
	//return c
	return c;
}

coord subtract(coord b, coord a)
{
	//subtract values in a from b
	b.x -= a.x;
	b.y -= a.y;
	b.z -= a.z;
	
	b.r -= a.r;
	b.g -= a.g;
	b.b -= a.b;
	
	b.nx -= a.nx;
	b.ny -= a.ny;
	b.nz -= a.nz;
	
	b.t1 -= a.t1;
	b.t2 -= a.t2;
	
	//return b
	return b;
}

coord calcBezierCurve(vector<coord> cPts, int degree, float t)
{
	//DeCastelJau's algorithm at time t
	for(int i = 0; i < degree; i++)
	{
		for(int j = 0; j < degree - i; j++)
		{
			cPts[j] = interpolate(cPts[j], cPts[j+1], t);
		}
	}
	
	return cPts[0];
}

vector<coord> calcBezierCurve(vector<coord> cPts, int degree, int n)
{
	//pts to be returned
	vector<coord> pts;
	
	//for each t time, calculate point on curve
	for(float t = 0; t <= 1; t += 1.0 / (n+1))
	{
		coord c = calcBezierCurve(cPts, degree, t);
		
		pts.push_back(c);
	}
	
	return pts;
}

vector<coord> calcBezierPatch(vector<coord> cPts,
							  int u_degree,
							  int v_degree,
							  int n)
{
	//points to be returned
	vector<coord> pts;
	//step for time based on n
	float tStep = 1.0 / n;
	//time variables to step with
	float t1, t2;
	
	//go through each time step and calculate the point on curve
	for(int i = 0; i <= n; i++)
	{
		t2 = i * tStep;
		for(int j = 0; j <= n; j++)
		{
			t1 = j*tStep;
			
			//calc the point on curve
			coord c = calcBezierPatch(cPts, u_degree, v_degree, t1, t2);
			//push it back
			pts.push_back(c);
		}
	}
	
	//return the vector of points
	return pts;
}

coord calcBezierPatch(vector<coord> cPts,
					  int u_degree,
					  int v_degree,
					  float t1,
					  float t2)
{
	//coordinate to be returned
	coord c;
	
	//vectors of u coordinates (rows) to calculate with
	vector<coord> u;
	//vectors of v coords (columns) to calculate with
	vector<coord> v; 
	
	for(int i = 0; i < v_degree + 1; i++)
	{
		//get the u coordinates to calculate with
		for(int j = 0; j < u_degree + 1; j++)
		{
			int ij = i * (u_degree + 1) + j;
			u.push_back(cPts[ij]);
		}
		
		//calculate at time t1 and push back to v
		c = calcBezierCurve(u, u_degree, t1);
		v.push_back(c);
		//clear u vector as we'll have multiple
		u.clear();
	}
	
	//calculate the u points at time t2
	c = calcBezierCurve(v, v_degree, t2);
	
	//return that point
	return c;
}

vector<coord> getDerBezierCurve(vector<coord> cPts, int degree)
{
	//control points for derivative
	vector<coord> derPts;
	
	//for each of the points, get new derivative point
	for(int i = 0; i < degree - 1; i++)
	{
		//new coordinate to be pusehd back
		coord c;
		c.init();
		
		//subtract adjacent points from one another
		c = subtract(cPts[i+1], cPts[i]);
		//multiply by n
		c = mult(c, degree);
		
		//push back c
		derPts.push_back(c);
	}
	
	//return control points for derivative
	return derPts;
}

void getDerBezierPatch(vector<coord> cPts,
					   vector<coord> &du,
					   vector<coord> &dv, 
					   int u_degree,
					   int v_degree,
					   int n)
{
	//get d/du
	for(int i = 0; i < (u_degree+1)*(v_degree+1);)
	{
		coord c;
		c.init();
		c = subtract(cPts[i], cPts[i+1]);
		//c = mult(c, u_degree);
		
		du.push_back(c);
		
		i++;
		if((i+1) % (u_degree+1) == 0)
			i++;
	}
	
	//get d/dv
	for(int i = 0; i < (u_degree+1)*(v_degree); i++)
	{
		coord c;
		c.init();
		c = subtract(cPts[i], cPts[i+u_degree+1]);
		//c = mult(c, v_degree);
		
		dv.push_back(c);
	}
	
	return;
}
