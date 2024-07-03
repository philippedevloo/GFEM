// ---- buildcrack_cgal.h ---- Sat Jun 18 12:26:21 CST 2004
// (c) Copyright: C. Armando Duarte, 2002-2004. All rights reserved
// See README.CopyRight file for further details.
//

#ifndef BUILDCRACK_CGALH
#define BUILDCRACK_CGALH

#include <iosfwd>
#include <typeinfo>
#include <vector>
#include <algorithm> 

#include "functors.h"
#include "crackel.h"
#include "cracknod.h"
#include "config.h"

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/tags.h>

/*
 * CGAL configuration, like kernel selection, customization of Polyhedron, and
 * typedefs for Polyhedron surface
 *
 * NOTE: All other CGAL related typedefs should be added to config_cgal.h
 */

// ***************************************************
// ***************************************************

// CGAL Kernel selection
// =====================

// typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt  Kernel;
//typedef CGAL::Exact_predicates_exact_constructions_kernel<realType> Kernel;

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel       CGAL_EP_EC_Kernel;


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel     CGAL_EP_IC_Kernel;

#if USE_CGAL_EXACT_ARITHMETIC

// NOTE: To pass all QA tests, we may need to use this predicate:
typedef CGAL_EP_EC_Kernel     CGAL_Kernel;

#else

// Recommended for Delaunay_triangulation_3
// The performance is close to CGAL::Simple_cartesian<realType> and
// order of magnitude better than 
// CGAL::Exact_predicates_exact_constructions_kernel CGAL_Kernel; 

typedef CGAL_EP_IC_Kernel     CGAL_Kernel;

//#include <CGAL/Simple_cartesian.h>
//typedef CGAL::Simple_cartesian<realType>        CGAL_Kernel;
#endif

// ***************************************************
// ***************************************************

// Customized facets and vertices definition
// ( see http://www.cgal.org/Tutorials/Polyhedron/index.html pdf tutorials  
//  and source code files therein - HINT: In the source files, look for 
// enriched_polyhedron.h   )
//
// 

// CAD see:
// http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Polyhedron/Chapter_main.html#Section_25.5
//
// http://www.cgal.org/Manual/latest/doc_html/cgal_manual/HalfedgeDS/Chapter_main.html#sectionHdsExamples

// This is used on EdgeIntersectData to set what is the status of the source or target node
enum ESourceTargetStatus {ENotOnSurf, EOnTrueSurfFacet, EOnTempSurfFacet};

// Use this enumerator to set the crack vertex status
enum CrackVertexType { InsideCrSurface,  OnCrFrontInterior,
                       OnCrFrontBegin, OnCrFrontEnd, OnCrFrontStart, 
                       OnCrSurfaceContour /* at bnd of domain */ };

// Use this enumerator to set the crack facet status
enum CrackFacetType { TrueCrSurfaceFacet, TempCrSurfaceFacet };

std::ostream& 
operator<<(std::ostream& out, const CrackVertexType& neType);

std::ostream& 
operator<<(std::ostream& out, const CrackFacetType& feType);

class PolySurfFacetUserData;
class PolySurfVertexUserData;

// customized facet used for crack surface representation
// See documentation for CGAL::HalfedgeDS_face_base< Refs > Class
template <class Refs, class T = CGAL::Tag_true /*face with a reference to an incident halfedge by default*/>
class Customized_face : 
  public CGAL::HalfedgeDS_face_base<Refs, T> {

public:
  // default constructor is mandatory
  Customized_face() : facetStatus(TrueCrSurfaceFacet), 
                      facetID(-1), is_flipped(false), p_user_data(0) {}

  // Custom copy constructor. Implicit one will not work since we want
  // a deep copy of PolySurfFacetUserData object, not just a copy of
  // pointers.
  Customized_face(const Customized_face& cust_face_in);

  Customized_face<Refs, T>& operator=(const Customized_face& cust_face_in);
    
  ~Customized_face();

  // access/set data members

  void setFacetStatus( const CrackFacetType fStatus ) 
  { facetStatus = fStatus; }

  CrackFacetType getFacetStatus() const { return facetStatus; }

  void setFacetID(const int ID) { facetID = ID; }

  int getFacetID() const { return facetID;}

  // Get/set PolySurfFacetUserData object stored in this facet
  const PolySurfFacetUserData* userData() const { return p_user_data;}
  PolySurfFacetUserData*& userData() { return p_user_data;}
    
  // Get/set if a facet is flipped
  const bool& isFlipped() const { return is_flipped;}
  bool& isFlipped() { return is_flipped;}

protected:
  // TrueCrSurfaceFacet -> facet represents the actual crack surface (default).
  //                            
  // TempCrSurfaceFacet -> facet is used for virtual advance only (.i.e just 
  //                       for cutting the elements that have BranchFn 
  //                       enrichment in at least one onde)
  //
  CrackFacetType facetStatus;

  int facetID; // facet ID (e.g., set by user in .crf file)
    
  // tag if facet is flipped. Can only happen in temp facets
  bool is_flipped;

  // Pointer to object containing user data associated with this facet.
  // Actual object is application dependent.
  PolySurfFacetUserData* p_user_data;

};

// customized vertex. 
// If T = CGAL::Tag_true, vertex with a reference to an incident halfedge: Use this ALWAYS
// P is the type of Point stored
// See documentation for CGAL::HalfedgeDS_vertex_base< Refs >
template < class Refs, class T, class P >
class Customized_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P> {

public:
  // by default, all vertices are set as InsideCrSurface
  Customized_vertex() : vertexStatus(InsideCrSurface),
    p_user_data(0) { }

  // repeat mandatory constructors
  explicit Customized_vertex(const P& pt)
    : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt) ,
      vertexStatus(InsideCrSurface),
      p_user_data(0) {}

  // Added by Piyush Gupta
  // Custom copy constructor. Implicit one will not work since we want
  // a deep copy of PolySurfVertexUserData object, not just a copy of
  // pointers.
  Customized_vertex(const Customized_vertex& cust_vertex_in); 

  Customized_vertex<Refs, T, P>& operator=(const Customized_vertex& cust_vertex_in);

  ~Customized_vertex();

  // access/set data members

  // set vertex status
  void setVertexStatus( const CrackVertexType vStatus ) 
  { vertexStatus = vStatus; }

  // get vertex status
  CrackVertexType getVertexStatus( void ) const { return vertexStatus; }
  
  // check if crack vertex is on Front
  // To replace use of, e.g., getVertexStatus() == OnCrFront
  bool isVertexOnCrFront() const { return ( vertexStatus == OnCrFrontInterior || 
                                            vertexStatus == OnCrFrontBegin || 
                                            vertexStatus == OnCrFrontEnd ||
                                            vertexStatus == OnCrFrontStart );  }
 
  // Get/set PolySurfVertexUserData object stored in this vertex
  PolySurfVertexUserData* userData() const { return p_user_data;}
  PolySurfVertexUserData*& userData() { return p_user_data;}
 
protected:
  // InsideCrSurface -> vertex is inside of the crack surface. 
  //
  // If not, vertex is on the boundary of the crack surface
  //
  // OnCrFrontInterior          -> vertex is on the crackfront interior
  // OnCrFrontBegin             -> vertex is at the Begining of the crackfront
  // OnCrFrontEnd               -> vertex is at the End of the crackfront
  // OnCrFrontStart             -> vertex used as starting point for continuos crack front
  // OnCrSurfaceContour         -> vertex is on the boundary but it is not a 
  //                               crack front vertex. In this case, vertex is 
  //                               on the boundary of the domain of analysis
  CrackVertexType vertexStatus;

  // Added by Piyush Gupta
  // Pointer to the object containing user data associated with this vertex
  // Object can be application dependent
  PolySurfVertexUserData* p_user_data;

};

//
// defintion of the items that were set for the crgeoeng3d
// See documentation for PolyhedronItems_3 Concept R and CGAL::Polyhedron_items_3 Class
struct Customized_items : public CGAL::Polyhedron_items_3  {
  
  // face wrapper ( definition of the face using our customized face )
  template <class Refs, class Traits>
  struct Face_wrapper {
    typedef Customized_face< Refs, CGAL::Tag_true> Face;
  };
    
  // vertex wrapper ( definition of the vertex using our customized vertex )
  template < class Refs, class Traits>  
  struct Vertex_wrapper {    
    typedef typename Traits::Point_3 Point;
    typedef Customized_vertex< Refs, CGAL::Tag_true, Point > Vertex;
  };
  
};

// ****************************************************************
// ****************************************************************

// Helper class to associate user info with a vertex of a triagulation
// representing a crack surface. This data is passed at surface
// remeshing function.
// 
// NOTE that a CGAL triangulation is NOT the same thing as a
// Polyhedron surface! 
//
// Assignment of data to vertex or facets of a Polygedron is done by
// customizing the vertex and facet used to build the polyhedron (see
// Customized_vertex and Customized_face above).
//

class PolySurfVertexInfo {
public:
  PolySurfVertexInfo() : vertex_status(InsideCrSurface) {};

  PolySurfVertexInfo(const CrackVertexType vertexStatus_in) : 
    vertex_status(vertexStatus_in) {};

  // Added by Piyush Gupta
  // Need it for printing the crack status vertex in remesh function
  CrackVertexType getVertexInfo() {return vertex_status;}
  
  // Default copy constructor and assigment operator are fine

  // Type of surface vertex
  CrackVertexType vertex_status;
};

// ***************************************************
// ***************************************************

// typedefs for Polyhedron surface
// ===============================

// NOTE: All other CGAL related typedefs should be added to config_cgal.h

// TBD: Change this to format: CGAL_TypeName

typedef CGAL::Polyhedron_3<CGAL_Kernel, Customized_items>  Polyhedron;
typedef Polyhedron::HalfedgeDS                        HalfedgeDS;

//typedef CGAL::Polyhedron_3<CGAL_Kernel, Customized_items>  Polyhedron<CGAL_Kernel>;
// type alias used to hide a template parameter
template <class Kernel_Type>
using CGAL_Polyhedron_3 = CGAL::Polyhedron_3<Kernel_Type, Customized_items>;

typedef Polyhedron::Halfedge               Halfedge;
typedef Polyhedron::Halfedge_handle        Halfedge_handle;
typedef Polyhedron::Halfedge_const_handle  Halfedge_const_handle;
typedef Polyhedron::Halfedge::Facet_handle Facet_handle;
typedef Polyhedron::Facet_iterator         Facet_iterator;
typedef Polyhedron::Facet_const_iterator   Facet_const_iterator;
typedef Polyhedron::Facet_const_iterator   facetConstIter;
typedef Polyhedron::Edge_iterator          Edge_iterator;
typedef Polyhedron::Edge_const_iterator    Edge_const_iterator;

typedef Polyhedron::Vertex                                 Vertex;
typedef Polyhedron::Vertex_handle                          Vertex_handle;
typedef Polyhedron::Vertex_const_handle                    VertexConstHandle;
typedef Polyhedron::Vertex_const_iterator                  verConstIter;
typedef Polyhedron::Vertex_iterator                        vertexIter;
typedef Polyhedron::Halfedge_iterator                      halfIter;
typedef Polyhedron::Halfedge_const_iterator                halfConstIter;
typedef Polyhedron::Halfedge_around_facet_const_circulator halfFacetConstCirc;
typedef Polyhedron::Halfedge_around_facet_circulator       halfFacetCirc;
typedef Polyhedron::Halfedge_around_vertex_const_circulator halfVertConstCirc; 
typedef Polyhedron::Halfedge_around_vertex_circulator       halfVertCirc; 

typedef CGAL::Inverse_index< verConstIter >        Index;
typedef CGAL::Inverse_index< facetConstIter >      FacetIndex;
typedef CGAL::Inverse_index< halfIter >            Index_h;
typedef CGAL::Inverse_index< halfFacetConstCirc >  Index_circ;

// ***************************************************
// ***************************************************
// ***************************************************

// Class used to create a Polyhedron surface representing a crack surface.
//
// Vertices and facets of the Polyhedron are given to this class using methods
// addNode and addElem. They are stored as CrackNode and CrackEl objects.
//
// Member CrackBuilder::operator() creates a Polyhedron surface using them.

template <class HalfedgeDS>
class CrackBuilder : public CGAL::Modifier_base<HalfedgeDS> {

public:
  // constructor
  CrackBuilder() {}

  // method that is going to be used by the polyhedral surface object
  void operator()( HalfedgeDS& hds );

  void addElem( const int elemID, const int nNodes, 
                const ConcreteRigidArray1d< int, MAXNNODESCRACKELM > &uConnect,
                const bool isCrackFront);
  // we don't need to set cgalIndex when we call this function,
  // just when we create a CGAL vertex ( see operator () ) 
  void addNode( int nodeID,  
                const ConcreteRigidArray1d < realType, 3 >& coordinates );

  // Moved by Piyush Gupta
  // This function maps from userID to CGAL index
  // See Josuttis' Book pg. 413
  std::vector< CrackNode* >::const_iterator getCgalIndexWithID( const int ID );
  
  //  void addCrackFront( ConcreteRigidArray1D<int, MAXNNODESCRACKELM> &uConnect );

  void freeAll();

  // print
  void print(std::ostream & out = std::cout) const;

  void setCrackFront( Polyhedron &crSurface );

  void setCrackFrontEnds( Polyhedron &crSurface );

  //
  // function to set up crack front vertices when reading info from FOM Library
  void setCrackFront( Polyhedron &crSurface,
                      const Subscript* p_bordID2VertID,
                      const Subscript* p_bordVertLabel,
                      const Subscript& nBordVerts );

  // destructor
  virtual ~CrackBuilder();

private:
  // These containers are used to map from userID to cgalIdex.
  // To improve searching performance we must sort the containers first using
  // Binary Predicates (See Josuttis' Book pg. 123 and 
  // also the file SetSolver/mesh/functors.h ) 
  //
  // See Josuttis' Book pg. 397
  void  sortNodes();

  typedef std::vector <CrackNode *> crackNodeContainerType;
  typedef std::vector <CrackEl *>   crackElContainerType; 
  
  crackNodeContainerType  crackVertices;
  crackElContainerType  crackFacets;
  crackElContainerType  crackFrontEdges;

};
// 
// the only way to compile this was using inline on the definition.
// CAD: That's bcz you are defining a function in a header file which
//      is included multiple times --> will lead to multiple definitions
//      if inline is not used.
//
inline std::ostream& operator<<(std::ostream& os, 
                         const CrackBuilder<HalfedgeDS>& crackBuilder) {
  crackBuilder.print(os);
  return os;
}

inline std::ostream& 
operator<<(std::ostream& out, const CrackVertexType& neType) {

  switch(neType) {
  case InsideCrSurface : {
    out << "InsideCrSurface" ;
    break;
  }
  case OnCrFrontInterior : {
    out << "OnCrFrontInterior" ;
    break;
  }
  case OnCrFrontBegin : {
    out << "OnCrFrontBegin" ;
    break;
  }
  case OnCrFrontEnd : {
    out << "OnCrFrontEnd" ;
    break;
  }  
  case OnCrFrontStart : {
    out << "OnCrFrontStart" ;
    break;
  }    
  case OnCrSurfaceContour : {
    out << "OnCrSurfaceContour" ;
    break;
  }
  default: {
    std::cerr << "WARNING: invalid CrackVertexType to be printed!!"
         << " Vertex Type = " << static_cast<int>(neType) << std::endl;
    break; }
  }
  return out;
}

// enum CrackFacetType { TrueCrSurfaceFacet, TempCrSurfaceFacet };
inline std::ostream& 
operator<<(std::ostream& out, const CrackFacetType& feType) {

  switch(feType) {
  case TrueCrSurfaceFacet : {
    out << "TrueCrSurfaceFacet";
    break;
  }
  case TempCrSurfaceFacet : {
    out << "TempCrSurfaceFacet";
    break;
  }
  default: {
    std::cerr << "WARNING: invalid CrackFacetType to be printed!!"
         << " Facet Type = " << static_cast<int>(feType) << std::endl;
    break; }
  }
  return out;
}


// **************************************************************
// **************************************************************
// **************************************************************

// Creation of a Polyhedron from another one: Copy-constructor like function but
// the kernels used by the two Polyhedrons can be different

template <class Polyhedron_Type_in, class Polyhedron_Type_out>
class CopyPolyhedron : public CGAL::Modifier_base<typename Polyhedron_Type_out::HalfedgeDS> {

private:
  typedef CGAL::Modifier_base<typename Polyhedron_Type_out::HalfedgeDS> Base;

public:
  typedef Polyhedron_Type_in   Poly_Type_in;
  typedef Polyhedron_Type_out  Poly_Type_out;

  // constructor
  CopyPolyhedron(const Poly_Type_in& poly_surf_in) 
    : Base(), poly_surf_in(poly_surf_in) { }

  // Will copy to hds_out (of type Polyhedron_Type_out) the polyhedron
  // surface of type Poly_Type_in stored in poly_surf_in (which
  // was passed to the constructor)
  void operator()(typename Poly_Type_out::HalfedgeDS& hds_out);

private:
  const Poly_Type_in& poly_surf_in;
}; 
  
#ifdef XLC_QNOTEMPINC
#include "buildcrack_cgal.C"
#endif


#endif
