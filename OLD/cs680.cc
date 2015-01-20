#include "cs680.h"
#include "polyread.h"
#include "bezier.h"
#include <cmath>
#include <vector>
#include <iostream>
using namespace std;
using namespace render_direct;

/*
	Notes to self
	
	n_divisions => number of divisions to use
	
	Before each superquadric use:
	obj_normal_flag = true;
	render_m_attr.clear();
    index = render_m_attr.add_geometry();	//index = geometry info start
    index = render_m_attr.add_shading_offset();	//index = norm vect info 
    
    attr_point p;
    ...
    p.coord[index+0] = ...
    p.coord[index+1] = ...
    
*/

//constants for conversions to radians and pi
#define PI 3.141592654
const float TORAD = PI/180;

float sgn(float a)
{
	if(a > 0)
	{
		return 1;
	}
	else if(a < 0)
	{
		return -1;
	}
	else if(a == 0)
	{
		return 0;
	}
	
	return 0;
}

float cos(float a,float b)
{
	float sol = 0;
	
	sol = sgn(cos(a))*abs(pow(abs(cos(a)),b));
	
	return sol;
}

float sin(float a, float b)
{
	float sol = 0;
	
	sol = sgn(sin(a))*abs(pow(abs(sin(a)),b));
	
	return sol;
}

void outputPolySet(vector<coordinate> p, 
					vector<coordinate> n, 
					int vSize, 
					int nSize)
{
	int u0 = 0;
	int u1 = 1;
	int v0 = nSize;
	int v1 = v0 + 1;
	cout << "PolySet \"PN\"\n";
	cout << p.size() << " " << (vSize - 1)*(nSize - 1) << endl;
	for(int i = 0; i < vSize * nSize; i++)
	{
		cout << p[i].x << " " << p[i].y << " " << p[i].z << " \n";
		cout << n[i].x << " " << n[i].y << " " << n[i].z << endl;
	}
	for(int i = 0; i < (vSize - 1)*(nSize - 1); i++)
	{
		u1 = u0 + 1;
		v1 = v0 + 1;
		cout << u0 << " " << u1 << " " << v1 << " " << v0 << " -1\n";
		
		u0++;
		v0++;
		
		if((u0 + 1) % nSize == 0)
		{
			u0++;
			v0++;
		}
	}
}

namespace render_direct{

vector<coordinate> getNormVect(vector<coordinate> points, 
							   int vMult)
{
	/*
		vector structure (assuming 6 n_divisions)
				     u
		    | 0  1  2  3  4  5 
		  --|------------------
		  0 | 0  1  2  3  4  5
		  1 | 6  7  8  9  10 11
		v 2 | 12 13 14 15 16 17
		  3 | 18 19 20 21 22 23
		  4 | 24 25 26 27 28 29
		  5 | 30 31 32 33 34 35
		  
		(u,v) = v*n_divisions + u;
	*/
	
	vector<coordinate> normals;
	
	/*
	cout << "v  : " << vMult << endl;
	cout << "n  : " << n_divisions << endl;
	cout << "v*n: " << vMult*n_divisions << endl;
	*/
	
	int vSize = n_divisions;
	int nSize = vMult*n_divisions;
	int u0 = 0;
	int v0 = nSize+1;
	
	for(int i = 0; i < (vSize+1)*(nSize+1); i++)
	{
		int u1 = u0 + 1;
		int v1 = v0 + 1;
		
			//get current vertex and adjacent ones to cross with
			coordinate c0 = points[u0];
			coordinate c1 = points[u1];
			coordinate c2 = points[v1];
			coordinate c3 = points[v0];
			
			//compute the normal vectors
			coordinate n0;
			coordinate n1;
			coordinate n2;
			
			//crossing diagonals
			n0.x = c2.x - c0.x;
			n0.y = c2.y - c0.y;
			n0.z = c2.z - c0.z;
			
			n1.x = c3.x - c1.x;
			n1.y = c3.y - c1.y;
			n1.z = c3.z - c1.z;
			
			//compute the cross product
			n2.x = n0.y*n1.z - n0.z*n1.y;
			n2.y = n0.z*n1.x - n0.x*n1.z;
			n2.z = n0.x*n1.y - n0.y*n1.x;
			
			//n2.x += c0.x;
			//n2.y += c0.y;
			//n2.z += c0.z;
			
			//cout << "Pushing back vertex " << count << " / \n";
			//	cout << points.max_size() << endl;
			
			//add to normals vector
			normals.push_back(n2);
		
		u0++;
		v0++;
		/*
		if((u0) % nSize == 0)
		{
			u0++;
			v0++;
		}
		*/
	}
	
	return normals;
}

int render_bezier_curve(const string & vertex_type, 
						int degree, 
						const float * vertex)
{
	vector<vert> vertList;
	int vertSize = get_vertex_size(vertex_type);
	int geometry, color, normal, texture;
	
	set_attribute_indices(vertex_type, vertSize, geometry, color, normal
						  texture);
	
	//fill the vector vertList with vertices 
	for(int i = 0; i <= degree; i++)
	{
		//set j to the appropriate location
		int j = vertSize*i;
		
		//vert to be pushed to vertList
		vert v;
		
		//if geometry is included, we will insert geometry
		if(geometry != -1)
		{
			int g = j+geometry;
			v.g.x = vertex[g+0];
			v.g.y = vertex[g+1];
			v.g.z = vertex[g+2];
		}
		
		//same for color, normal, etc.
		if(color != -1)
		{
			int c = j+color;
			v.c.x = vertex[c+0];
			v.c.y = vertex[c+1];
			v.c.z = vertex[c+2];
		}
		
		if(normal != -1)
		{
			int n = j+normal;
			v.n.x = vertex[n+0];
			v.n.y = vertex[n+1];
			v.n.z = vertex[n+2];
		}
		
		if(texture != -1)
		{
			int t = j+texture;
			v.t.x = vertex[t+0];
			v.t.y = vertex[t+1];
			v.y.z = 0;
		}
		
		//pushback to vector
		vertList.push_back(v);
	}
	
	int n = n_divisions;
	float step = n / 1.0;
	vector<coordinate> outList;
	outList.push_back(vertList[0].g);
	
	for(float t = step; t < 1; t+= step)
	{
		vector<coordinate> vList = copyVertList(vertList);
		
		outList.push_back(calcBezier(vList, degree, t));
	}
	
	outList.push_back(vertList[vertList.size() - 1].g);
	
	/*
		*****************************************************************
		*****************************************************************
		OUTLIST NOW HAS THE VERTICES THAT WILL DRAW THE CURVE IN ORDER!!!
		*****************************************************************
		*****************************************************************
	*/
	
	for(int i = 0; i < vertList.size() - 1; i++)
	{
		coordinate g[2];
		g[0] = outList[i];
		g[1] = outList[i+1];
		
		//draw line outList[i] -> outList[i+1]
		for(int j = 0; j < 2; j++)
		{
			attr_point p;
				
			p.coord[0] = g[j].x;
			p.coord[1] = g[j].y;
			p.coord[2] = g[j].z;
			p.coord[3] = 1;
			
			p.coord[gDex+0] = g[j].x;
			p.coord[gDex+1] = g[j].y;
			p.coord[gDex+2] = g[j].z;
			
			int mode = (j+1) / 2;
			
			int mode = (j + 1) / 4;
			
			if(mode == 1)
			{
				coordinate n;
				n.x = g[1].x - g[0].x;
				n.y = g[1].y - g[0].y;
				n.z = g[1].z - g[0].z;
				poly_normal[0] = 1;
				poly_normal[1] = 1;
				poly_normal[2] = (n.x + n.y) / n.z;
			}
			
			poly_pipeline(p, mode);
		}
	}
	
	return 0;
}

int render_bezier_patch(const string &vertex_type,
						int u_degree, 
						int v_degree,
						const float * vertex)
{
	vector<vert> vertList;
	int vertSize = get_vertex_size(vertex_type);
	int geometry, color, normal, texture;
	
	set_attribute_indices(vertex_type, vertSize, geometry, color, normal
						  texture);
	
	vector<coordinate> outList;
	
	//fill the vector vertList with vertices 
	for(int i = 0; i < (u_degree+1) * (v_degree+1); i++)
	{
		//set j to the appropriate location
		int j = vertSize*i;
		
		//vert to be pushed to vertList
		vert v;
		
		//if geometry is included, we will insert geometry
		if(geometry != -1)
		{
			int g = j+geometry;
			v.g.x = vertex[g+0];
			v.g.y = vertex[g+1];
			v.g.z = vertex[g+2];
		}
		
		//same for color, normal, etc.
		if(color != -1)
		{
			int c = j+color;
			v.c.x = vertex[c+0];
			v.c.y = vertex[c+1];
			v.c.z = vertex[c+2];
		}
		
		if(normal != -1)
		{
			int n = j+normal;
			v.n.x = vertex[n+0];
			v.n.y = vertex[n+1];
			v.n.z = vertex[n+2];
		}
		
		if(texture != -1)
		{
			int t = j+texture;
			v.t.x = vertex[t+0];
			v.t.y = vertex[t+1];
			v.y.z = 0;
		}
		
		//pushback to vector
		vertList.push_back(v);
	}
	
	int n = n_divisions;
	float step = n / 1.0;
	vector<vert> vList = copyVertList(vertList);
	//int u = 0;
	//int v = 0;
	
	for(float t = 0; t <= 1; t += step)
	{
		for(int v = 0; v <= v_degree; v++)
		{
			for(int u = 0; u <= u_degree; u++)
			{
				outList.push_back(calcuv(u,v,u_degree,v_degree, vList,t));
			}
		}
		
		outList.push_back(o);
	}
	
	//geometry index and normal index
	render_m_attr.clear();
	render_m_attr.add_geometry();
	//render_m_attr.add_normal();
	render_m_attr.add_shading_offset();
	int gDex = render_m_attr.geometry; 
	//int nDex = render_m_attr.normal; 
	
	//cerr << "gDex: " << gDex << endl;
	//cerr << "nDex: " << nDex << endl;
	
	int vSize = n_divisions;
	int nSize = n_divisions;
	int u0 = 0;
	int v0 = nSize+1;
	//int faceNum = 0;
	
	for(int i = 0; i < (vSize)*(nSize); i++)
	{
		int u1 = u0 + 1;
		int v1 = v0 + 1;
			
			//cerr << "Starting face: " << faceNum << endl;
			
			vector<coordinate> g;
			vector<coordinate> n;
			g.clear();
			n.clear();
			
			//cerr << u0 << endl;
			
			g.push_back(outList[u0]);
			g.push_back(outList[u1]);
			g.push_back(outList[v1]);
			g.push_back(outList[v0]);
			
			for(int j = 0; j < 4; j++)
			{
					
				//cerr << "Setting up attr_point geo\n";
				
				attr_point p;
				
				p.coord[0] = g[j].x;
				p.coord[1] = g[j].y;
				p.coord[2] = g[j].z;
				p.coord[3] = 1;
				
				p.coord[gDex+0] = g[j].x;
				p.coord[gDex+1] = g[j].y;
				p.coord[gDex+2] = g[j].z;
				
				//cerr << "Setting up attr_point nor\n";
				
				//p.coord[nDex+0] = n[j].x;
				//p.coord[nDex+1] = n[j].y;
				//p.coord[nDex+2] = n[j].z;
				
				int mode = (j + 1) / 4;
				if(mode == 1)
				{
					coordinate d1 = g[2] - g[0];
					coordinate d2 = g[3] - g[1];
					
					coordinate n = calcNormal(d1, d2);
					
					poly_normal[0] = n.x;
					poly_normal[1] = n.y;
					poly_normal[2] = n.z;
				}
				
				//cerr << "Set mode: " << mode << endl;
				//cerr << "Starting polypipeline " << j+1 << endl;
				
				poly_pipeline(p, mode);
				
				//cerr << "Finished polypipeline " << j+1 << endl;
				
			}
			
			//cerr << "Finished polypipeline " << faceNum++ << "\n\n";
			
		u0++;
		v0++;
		
		if((u0+1) % (nSize+1) == 0)
		{
			u0++;
			v0++;
		}
	}
	
	return 0;
}

} //close namespace renderdirect

/*
int REDirect::rd_bezier_curve(const string & vertex_type, 
						  	  int degree, 
						  	  const float * vertex)
{
	return 0;
}

int REDirect::rd_bezier_patch(const string &vertex_type,
						  	  int u_degree, 
						  	  int v_degree,
						  	  const float * vertex)
{
	return 0;
}
*/

int REDirect::rd_sqsphere(float radius,
						  float north,
						  float east,
						  float zmin,
						  float zmax,
						  float theatamax)
{
	int vMult = 2;
	float step_phi = PI / n_divisions;
	float step_theta = 2*PI / (vMult*n_divisions);
	float phi = -PI / 2;
	float theta = 0;
	
	//vector of points
	vector<coordinate> points;
	//vector of normals
	vector<coordinate> normals;
	
	/*
		vector structure (assuming 6 n_divisions)
				     u
		    | 0  1  2  3  4  5 
		  --|------------------
		  0 | 0  1  2  3  4  5
		  1 | 6  7  8  9  10 11
		v 2 | 12 13 14 15 16 17
		  3 | 18 19 20 21 22 23
		  4 | 24 25 26 27 28 29
		  5 | 30 31 32 33 34 35
		  
		(u,v) = v*n_divisions + u;
	*/
	
	//fill in the geometry points
	
	int count = 0;
	
	//v controls the phi values
	for(int v = 0; v < n_divisions + 1; v++)
	{
		//u controls the theta values
		for(int u = 0; u < vMult*n_divisions + 1; u++)
		{
			//angle to compute
			float aPhi = v*step_phi+phi;
			float aTheta = u*step_theta+theta;
			
			//cout << "aPhi " << aPhi << endl;
			//cout << "aTheta " << aTheta << endl;
			
			//get geometry and add it to the vector
			coordinate c;
			c.x = radius*cos(aPhi,north)*cos(aTheta,east);
			c.y = radius*cos(aPhi,north)*sin(aTheta,east);
			c.z = radius*sin(aPhi,north);
			
			/*
			cout << "Pushing back vertex " << count++ << " / ";
				cout << points.max_size() << endl;
			*/
			
			points.push_back(c);
			//cout << "point " << count << " added: ";
			//cout << points[count].x << " " << points[count].y << " ";
			//cout << points[count].z << endl;
			count++;
		}
	}
	
	//get the normal vectors
	normals = getNormVect(points, vMult);
	
	/*
	 * FOR USE LATER!
	 *
	 
	//geometry index and normal index
	render_m_attr.clear();
	g = render_m_attr.add_geometry();
	n = render_m_attr.add_shading_offset();
	p.coord[g+0] = x;
	p.coord[g+1] = y;
	p.coord[g+2] = z;
	*
	*
	*/
	
	/*
		vector structure (assuming 6 n_divisions)
				     u
		    | 0  1  2  3  4  5 
		  --|------------------
		  0 | 0  1  2  3  4  5
		  1 | 6  7  8  9  10 11
		v 2 | 12 13 14 15 16 17
		  3 | 18 19 20 21 22 23
		  4 | 24 25 26 27 28 29
		  5 | 30 31 32 33 34 35
		  
		(u,v) = v*n_divisions + u;
	*/
	
	//outputPolySet(points, normals, vMult*n_divisions+1, n_divisions+1);
	
	//geometry index and normal index
	render_m_attr.clear();
	render_m_attr.add_geometry();
	//render_m_attr.add_normal();
	render_m_attr.add_shading_offset();
	int gDex = render_m_attr.geometry; 
	//int nDex = render_m_attr.normal; 
	
	//cerr << "gDex: " << gDex << endl;
	//cerr << "nDex: " << nDex << endl;
	
	int vSize = n_divisions;
	int nSize = vMult*n_divisions;
	int u0 = 0;
	int v0 = nSize+1;
	//int faceNum = 0;
	
	for(int i = 0; i < (vSize)*(nSize); i++)
	{
		int u1 = u0 + 1;
		int v1 = v0 + 1;
			
			//cerr << "Starting face: " << faceNum << endl;
			
			vector<coordinate> g;
			vector<coordinate> n;
			g.clear();
			n.clear();
			
			//cerr << u0 << endl;
			
			g.push_back(points[u0]);
			g.push_back(points[u1]);
			g.push_back(points[v1]);
			g.push_back(points[v0]);
			
			n.push_back(normals[u0]);
			n.push_back(normals[u1]);
			n.push_back(normals[v1]);
			n.push_back(normals[v0]);
			
			//nDex = 0;
			//gDex = 3;
			
			for(int j = 0; j < 4; j++)
			{
					
				//cerr << "Setting up attr_point geo\n";
				
				attr_point p;
				
				p.coord[0] = g[j].x;
				p.coord[1] = g[j].y;
				p.coord[2] = g[j].z;
				p.coord[3] = 1;
				
				p.coord[gDex+0] = g[j].x;
				p.coord[gDex+1] = g[j].y;
				p.coord[gDex+2] = g[j].z;
				
				//cerr << "Setting up attr_point nor\n";
				
				//p.coord[nDex+0] = n[j].x;
				//p.coord[nDex+1] = n[j].y;
				//p.coord[nDex+2] = n[j].z;
				
				int mode = (j + 1) / 4;
				if(mode == 1)
				{
					poly_normal[0] = n[0].x;
					poly_normal[1] = n[0].y;
					poly_normal[2] = n[0].z;
				}
				
				//cerr << "Set mode: " << mode << endl;
				//cerr << "Starting polypipeline " << j+1 << endl;
				
				poly_pipeline(p, mode);
				
				//cerr << "Finished polypipeline " << j+1 << endl;
				
			}
			
			//cerr << "Finished polypipeline " << faceNum++ << "\n\n";
			
		u0++;
		v0++;
		
		if((u0+1) % (nSize+1) == 0)
		{
			u0++;
			v0++;
		}
	}
	
	return 0;
}

int REDirect::rd_sqtorus(float radius1,
						 float radius2,
						 float north,
						 float east,
						 float phimin,
						 float phimax,
						 float thetamax)
{
	float step = 2*PI / n_divisions;
	float phi = -PI;
	float theta = 0;
	int vMult = 1;
	
	//vector of points
	vector<coordinate> points;
	//vector of normals
	vector<coordinate> normals;
	
	for(int v = 0; v < vMult*n_divisions+1; v++)
	{
		for(int u = 0; u < n_divisions+1; u++)
		{
			//angles to compute with
			float aPhi = v*step+phi;
			float aTheta = u*step+theta;
			
			//set coordinate values
			coordinate c;
			//cos(u)(R + r cos(v))
			c.x = cos(aTheta,east)*(radius1 + radius2*cos(aPhi,north));
			//sin(u)(R + r cos(v))
			c.y = sin(aTheta,east)*(radius1 + radius2*cos(aPhi,north));
			//r sin(v)
			c.z = radius2*sin(aPhi,north);
			
			//add to vector
			points.push_back(c);
		}
	}
	
	//get the normal vectors
	normals = getNormVect(points,vMult);
	
	//outputPolySet(points, normals, vMult*n_divisions+1, n_divisions+1);
	
	//geometry index and normal index
	render_m_attr.clear();
	render_m_attr.add_geometry();
	//render_m_attr.add_normal();
	render_m_attr.add_shading_offset();
	int gDex = render_m_attr.geometry; 
	//int nDex = render_m_attr.normal; 
	
	//cerr << "gDex: " << gDex << endl;
	//cerr << "nDex: " << nDex << endl;
	
	int vSize = n_divisions;
	int nSize = vMult*n_divisions;
	int u0 = 0;
	int v0 = nSize+1;
	//int faceNum = 0;
	
	for(int i = 0; i < (vSize)*(nSize); i++)
	{
		int u1 = u0 + 1;
		int v1 = v0 + 1;
			
			//cerr << "Starting face: " << faceNum << endl;
			
			vector<coordinate> g;
			vector<coordinate> n;
			g.clear();
			n.clear();
			
			//cerr << u0 << endl;
			
			g.push_back(points[u0]);
			g.push_back(points[u1]);
			g.push_back(points[v1]);
			g.push_back(points[v0]);
			
			n.push_back(normals[u0]);
			n.push_back(normals[u1]);
			n.push_back(normals[v1]);
			n.push_back(normals[v0]);
			
			//nDex = 0;
			//gDex = 3;
			
			for(int j = 0; j < 4; j++)
			{
					
				//cerr << "Setting up attr_point geo\n";
				
				attr_point p;
				
				p.coord[0] = g[j].x;
				p.coord[1] = g[j].y;
				p.coord[2] = g[j].z;
				p.coord[3] = 1;
				
				p.coord[gDex+0] = g[j].x;
				p.coord[gDex+1] = g[j].y;
				p.coord[gDex+2] = g[j].z;
				
				//cerr << "Setting up attr_point nor\n";
				
				//p.coord[nDex+0] = n[j].x;
				//p.coord[nDex+1] = n[j].y;
				//p.coord[nDex+2] = n[j].z;
				
				int mode = (j + 1) / 4;
				if(mode == 1)
				{
					poly_normal[0] = n[0].x;
					poly_normal[1] = n[0].y;
					poly_normal[2] = n[0].z;
				}
				
				//cerr << "Set mode: " << mode << endl;
				//cerr << "Starting polypipeline " << j+1 << endl;
				
				poly_pipeline(p, mode);
				
				//cerr << "Finished polypipeline " << j+1 << endl;
				
			}
			
			//cerr << "Finished polypipeline " << faceNum++ << "\n\n";
			
		u0++;
		v0++;
		
		if((u0+1) % (nSize+1) == 0)
		{
			u0++;
			v0++;
		}
	}
	
	return 0;
}