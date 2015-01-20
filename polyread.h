#ifndef polyread_h
#define polyread_h

#include <iostream>
#include <string>
using std::istream;
using std::string;


/****************************************************************************/
// Function declarations

static void EatSpace(istream& );
// Removes extraneous spaces and comments from input

int  get_vertex_size(const string & vertex_type);
  // Given a vertex type string returns the number of array elements
  // taken up by a single vertex.


int parse_polyset(istream& istr, string & vtype, 
		  int & nvertex, int & nface, 
		  float coord [], int vertex_list[]);

  // Given an input stream returns 
  //     the vertex type string (vtype), 
  //     the number of vertices (nvertex),
  //     the number of faces (nface),
  //     and the values of the vertex coordinate (coord)
  //     and the values of the polygon vertex indices (vertex_list)
  // The arrays must already have memory allocated to them.

int  set_attribute_indices(const string & vertex_type, 
			   int & vertex_size,
			   int & geometry,
			   int & color,
			   int & normal,
			   int & texture);

  // Given a vertex type string returns
  //    the number of elements taken by a vertex (vertex_size),
  //    the index positions of the start of various attributes within a vertex.
  //    Current attributes include geometry, color, normal, and texture.
  //    The sizes of each attribute are as follows:
  //    geometry - 3 values
  //    color - 3 values
  //    normal - 3 values
  //    texture - 2 values

/****************************************************************************/

#endif