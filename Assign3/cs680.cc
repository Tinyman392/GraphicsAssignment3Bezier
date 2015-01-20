#include "cs680.h"
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
//const float TORAD = PI/180;

//coordinate structure for vectors
struct coordinate
{
	float x;
	float y;
	float z;
};

//calculate the cross product of vectors n0 and n1
coord cross(coord n0, coord n1)
{
	coord n2;
	
	n2.x = n0.y*n1.z - n0.z*n1.y;
	n2.y = n0.z*n1.x - n0.x*n1.z;
	n2.z = n0.x*n1.y - n0.y*n1.x;
	
	return n2;
}

//get booleans based on the vType passed in
void getBools(string vType, bool &g, bool &c, bool &n, bool &t)
{
	//calls setAttributeIndices
	setAttributeIndices(vType, g, c, n, t);
}

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

//calculates normals for superquadrics
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

//render the bezier curve
int render_bezier_curve(const string & vertex_type, 
						int degree, 
						const float * vertex)
{
	//points to be drawn
	vector<coord> pts;
	//control points
	vector<coord> cPts;
	//booleans for vertex info
	bool geo, col, nor, tex;
	//indices for vertex info
	int gDex, cDex, nDex, tDex;
	
	getBools(vertex_type, geo, col, nor, tex);
	cPts = arr2vect(vertex, degree, vertex_type);
	pts = calcBezierCurve(cPts, degree, n_divisions);
	
	//needs to be called to run
	render_m_attr.clear();
	render_m_attr.add_shading_offset();
	
	//check if booleans are set
	//of they are, get the indices for it and let the 
		//renderer know we'll be using them
	if(geo)
	{
		render_m_attr.add_geometry();
		gDex = render_m_attr.geometry;
	}
	if(col)
	{
		render_m_attr.add_color();
		cDex = render_m_attr.color;
	}
	if(nor)
	{
		render_m_attr.add_normal();
		nDex = render_m_attr.normal;
	}
	if(tex)
	{
		render_m_attr.add_texture();
		tDex = render_m_attr.texture;
	}
	
	//mode for drawing
	//set to 0 to initially move
	int mode = 0;
	//number of divisions
	int n = n_divisions;
	for(int i = 0; i < n+1; i++)
	{
		//attribute point to draw
		attr_point p;
		
		//set first 4 coordinates
		//needed to run properly
		p.coord[0] = pts[i].x;
		p.coord[1] = pts[i].y;
		p.coord[2] = pts[i].z;
		p.coord[3] = 1;
		
		//if boolean is set, we will set attribute point appropriately
		if(geo)
		{
			p.coord[gDex + 0] = pts[i].x;
			p.coord[gDex + 1] = pts[i].y;
			p.coord[gDex + 2] = pts[i].z;
		}
		if(col)
		{
			p.coord[cDex + 0] = pts[i].r;
			p.coord[cDex + 1] = pts[i].g;
			p.coord[cDex + 2] = pts[i].b;
		}
		if(nor)
		{
			p.coord[nDex + 0] = pts[i].nx;
			p.coord[nDex + 1] = pts[i].ny;
			p.coord[nDex + 2] = pts[i].nz;
		}
		if(tex)
		{
			p.coord[tDex + 0] = pts[i].t1;
			p.coord[tDex + 1] = pts[i].t2;
		}
		
		//send to line_pipeline to draw
		line_pipeline(p, mode);
		//set mode to 1
		mode = 1;
	}
	
	return 0;
}

int render_bezier_patch(const string &vertex_type,
						int u_degree, 
						int v_degree,
						const float * vertex)
{
	//vector of
		//points to draw
		//control points
		//derivative control points
		//derivative points on curve
	vector<coord> pts;
	vector<coord> cPts;
	vector<coord> du;
	vector<coord> dv;
	vector<coord> duPts;
	vector<coord> dvPts;
	//booleans for vertex information
	bool geo, col, nor, tex;
	//indices for vertex information
	int gDex, cDex, nDex, tDex;
	
	//get the booleans
	getBools(vertex_type, geo, col, nor, tex);
	//convert the vertex array to points
	cPts = arr2vect(vertex, u_degree, v_degree, vertex_type);
	//get the points to plot
	pts = calcBezierPatch(cPts, u_degree, v_degree, n_divisions);
	//get the derivative control points
	getDerBezierPatch(cPts, du, dv, u_degree, v_degree, n_divisions);
	
	//get the derivative points
	duPts = calcBezierPatch(du, u_degree-1, v_degree, n_divisions);
	dvPts = calcBezierPatch(dv, u_degree, v_degree-1, n_divisions);
	
	//needed to run properly
	render_m_attr.clear();
	render_m_attr.add_shading_offset();
	
	//set normals since we are calculating them
	render_m_attr.add_normal();
	nDex = render_m_attr.normal;
	
	//check booleans to see what's needed
	//if something is set, we will get index and let renderer know
		//we'll be using it
	if(geo)
	{
		render_m_attr.add_geometry();
		gDex = render_m_attr.geometry;
	}
	if(col)
	{
		render_m_attr.add_color();
		cDex = render_m_attr.color;
	}
	if(tex)
	{
		render_m_attr.add_texture();
		tDex = render_m_attr.texture;
	}
	
	//indices for faces
	int u1, u2, v1, v2;
	//number of divisions
	int n = n_divisions;
	for(int i = 0; i < (n)*(n+1);)
	{
		//vectors for the face
			//and derivative points of the face
		vector<coord> face;
		vector<coord> derv;
		vector<coord> deru;
		//make sure they are clean
		face.clear();
		derv.clear();
		deru.clear();
		
		//set face corners appropriately
		u1 = i;
		u2 = u1 + 1;
		v1 = u1 + n + 1;
		v2 = v1 + 1;
		
		//push back points to each vector
		face.push_back(pts[u1]);
		face.push_back(pts[u2]);
		face.push_back(pts[v2]);
		face.push_back(pts[v1]);
		
		derv.push_back(dvPts[u1]);
		derv.push_back(dvPts[u2]);
		derv.push_back(dvPts[v2]);
		derv.push_back(dvPts[v1]);
		
		deru.push_back(duPts[u1]);
		deru.push_back(duPts[u2]);
		deru.push_back(duPts[v2]);
		deru.push_back(duPts[v1]);
		
		//draw the face
		for(int j = 0; j < 4; j++)
		{
			//attribute point to be passed
			attr_point p;
			
			//set first four coords appropriately
				//needed to avoid seg faults
			p.coord[0] = face[j].x;
			p.coord[1] = face[j].y;
			p.coord[2] = face[j].z;
			p.coord[3] = 1;
			
			//if flags are set, set those points too
			if(geo)
			{
				p.coord[gDex + 0] = face[j].x;
				p.coord[gDex + 1] = face[j].y;
				p.coord[gDex + 2] = face[j].z;
			}
			if(col)
			{
				p.coord[cDex + 0] = face[j].r;
				p.coord[cDex + 1] = face[j].g;
				p.coord[cDex + 2] = face[j].b;
			}
			if(nor)
			{
				p.coord[nDex + 0] = face[j].nx;
				p.coord[nDex + 1] = face[j].ny;
				p.coord[nDex + 2] = face[j].nz;
			}
			else
			{
				coord n = cross(deru[j],derv[j]);
				p.coord[nDex + 0] = n.x;
				p.coord[nDex + 1] = n.y;
				p.coord[nDex + 2] = n.z;
			}
			if(tex)
			{
				p.coord[tDex + 0] = face[j].t1;
				p.coord[tDex + 1] = face[j].t2;
			}
			
			//set mode to draw on last vertex
			int mode = (j + 1) / 4;
			
			//if mode is 1, we need a normal as well
			if(mode == 1)
			{
				coord n1, n2;
				n1 = subtract(face[3], face[1]);
				n2 = subtract(face[2], face[0]);
				
				coord n = cross(n2, n1);
				//n = cross(dvPts[j],duPts[i]);
				
				poly_normal[0] = n.x;
				poly_normal[1] = n.y;
				poly_normal[2] = n.z;
			}
			
			//send to pipeline
			poly_pipeline(p, mode);
		}
		
		//increment i
		i++;
		//check to see if at edge
		if((i+1) % (n+1) == 0)
			i++;
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
	
	int vSize = n_divisions;
	int nSize = vMult*n_divisions;
	int u0 = 0;
	int v0 = nSize+1;
	//int faceNum = 0;
	
	for(int i = 0; i < (vSize)*(nSize); i++)
	{
		int u1 = u0 + 1;
		int v1 = v0 + 1;
			
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