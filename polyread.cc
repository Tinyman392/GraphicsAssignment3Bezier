#include "polyread.h"

const int ILLEGAL_TOKEN = -1;
const int UNEXPECTED_EOF = -2;
const int ILLEGAL_VERTEX_INDEX = -3;
const int ILLEGAL_VERTEX_LIST_MAX_SIZE = -4;
const int ILLEGAL_VERTEX_LIST_OVERFLOW = -5;
const int UNKNOWN_VERTEX_TYPE = -6;
const int OK = 1;
const int VERTEX_LIST_MAX_SIZE = 1000000;
const int VERTEX_LIST_OVERFLOW = -7;


// Global variables

int rd_input_line = 0;
// Useful for error messages

/****************************************************************************/

static void	EatSpace(istream& file)
{
  // Eat white space and comments
  
  int dummychar;
  
  int	valid = 0;
  
  // eat whitespace and comments
  
  do{
    dummychar = file.get();
    if(file.eof()) return;               // Read past end of file
    switch(dummychar){
    case ' ':                            // Eat white space
    case '\t':
    case '\r':   // Needed for DOS files
      break;
    case '\n':
      rd_input_line++;
      break;
    case '#':                            // Eat comments
      while('\n' != file.get() && !file.eof());
      if(file.eof())  return;
      rd_input_line++;
      break;
    default:
      file.putback(dummychar);
      valid = 1;
      break;
    }
  }while(!valid);
}

int parse_polyset(istream& istr, string & vtype, 
			 int & nvertex, int & nface, 
			 float coord [], int vertex_list[])
{
  string token;

  int vsize;
  int size;

  int err = 0;

  EatSpace(istr);

  istr >> token;

  if(token != "PolySet")
    return ILLEGAL_TOKEN;

  EatSpace(istr);

  if(istr.eof())
    return UNEXPECTED_EOF;

  istr >> vtype;
  vsize = get_vertex_size(vtype);
  if(vsize < 0) // error code
    return vsize;

  EatSpace(istr);
  istr >> nvertex;
  EatSpace(istr);
  istr >> nface;

  size = vsize * nvertex;

  for(int i = 0; i < size; i++)
    {
      EatSpace(istr);
      if(istr.eof())
	{
	  return UNEXPECTED_EOF;
	}
      double dvalue;

      EatSpace(istr);
      istr >> dvalue;

      coord[i] = dvalue;
    }

  int index_cnt = 0;
  int index;

  for(int i = 0; i < nface;)
    {
      EatSpace(istr);
      if(istr.eof())
	{
	  return UNEXPECTED_EOF;
	}

      EatSpace(istr);
      istr >> index;

      if(index >= nvertex)
	{
	  return ILLEGAL_VERTEX_INDEX;
	}
      vertex_list[index_cnt++] = index;
      if(index_cnt > VERTEX_LIST_MAX_SIZE)
	{
	  return VERTEX_LIST_OVERFLOW;
	}
      if(index == -1) i++;  // New face
    }

  
  return err;
}


int  get_vertex_size(const string & vertex_type)
{
  const char * attribute = &vertex_type[0];
  int i = vertex_type.length();
  int vertex_size = 0;
  
  for(; i; attribute++, i--)
    {
      switch(*attribute)
        {
        case '\"': // Ignore \" 
          break;
        case 'P':
          vertex_size += 3;
          break;
        case 'N': case 'D':
          vertex_size += 3;
          break;
        case 'C':
          vertex_size += 3;
          break;
        case 'T':
          vertex_size += 2;
          break;
        default:
          break;
        }
    }

  return vertex_size;
}


int  set_attribute_indices(const string & vertex_type, 
                                  int & vertex_size,
                                  int & geometry,
                                  int & color,
                                  int & normal,
                                  int & texture)
{
  int err = OK;
  const char * attribute = &vertex_type[0];
  int i = vertex_type.length();
  vertex_size = 0;
  
  geometry = color = normal = texture = -1;  // Not used

  for(; i; attribute++, i--)
    {
      switch(*attribute)
        {
        case '\"': // Ignore \" 
          break;
        case 'P':
          geometry = vertex_size;
          vertex_size += 3;
          break;
        case 'N': case 'D':
          normal = vertex_size;
          vertex_size += 3;
          break;
        case 'C':
          color = vertex_size;
          vertex_size += 3;
          break;
        case 'T':
          texture = vertex_size;
          vertex_size += 2;
          break;
        default:
          err = UNKNOWN_VERTEX_TYPE;
          break;
        }
    }

  return err;
}