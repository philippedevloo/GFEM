// ---- config_cgal.h ---- Fri Oct 28 16:40:07 CDT 2011
// (c) Copyright: C. Armando Duarte, 2002-2011. All rights reserved
// See README.CopyRight file for further details.
//

#ifndef CONFIG_CGALH
#define CONFIG_CGALH

/*
 * typedefs for CGAL classes
 *
 * Some of the typedefs were originally in buildcrack_cgal.h and crgeoeng3d_cgal.h
 */

#include <CGAL/Inverse_index.h>

//
// includes for tetrahedralization
#include <CGAL/Triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>
//#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_hierarchy_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
// CGAL Delaunay triangulation with vertex and cell info provided by user
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

// includes for 2D triangulations
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

// includes for typedefs for CrackGeoEng3D
#include <CGAL/enum.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>

// includes for AABB tree: typedefs used by class CrackManager3D_2
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>

// includes for crack remesh
#include <CGAL/AABB_polyhedral_oracle.h>

#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Timer.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>

//#include "Complex_2_in_triangulation_3_wInfo_polyhedron_builder.h"

#include <CGAL/Cartesian_converter.h>

// Kernel selection and Polyhedron customization
#include "buildcrack_cgal.h"

// ***************************************************
// ***************************************************

namespace iset {
  // This function allows us to switch kernels without modifying the code
  // Just overload them with the data type used by the kernel

#if USE_CGAL_EXACT_ARITHMETIC

  // NOTE: May be we only need this one!!
  inline realType to_realType(const CGAL_Kernel::FT& value) {
    return CGAL::to_double( value ); }


  //#ifdef USE_CGAL_WITH_GMP

  // NOTE: This is used with CGAL 3.3.1
  //       Requires libgmp and that CGAL was compiled with it!
  //  inline realType 
  //  to_realType(const CGAL::Lazy_exact_nt<CGAL::Gmpq>& value) {
  //  return CGAL::to_double( value ); }

  // NOTE: This is used with CGAL 3.2.1
  // inline realType 
  // to_realType(const CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> >& value) {
  //  return CGAL::to_double( value ); }

#else

  inline realType to_realType(const realType& value) {
    return value; }

#endif

}/* end of namespace iset */

// ***************************************************
// ***************************************************

// typedefs for Polyhedron surface: SEE buildcrack_cgal.h

// ***************************************************
// ***************************************************

// Miscellaneous typedefs

// CAD: Use of typedefs like these is NO longer acceptable.
//
// (*) All CGAL types used in ISET must have a template argument for
// the CGAL kernel, so we know which kernel is used with each type.

//typedef CGAL::Point_3<CGAL_Kernel>                 CGALPoint3;
//typedef CGAL::Point_2<CGAL_Kernel>                 CGALPoint2;
//typedef CGAL::Segment_3<CGAL_Kernel>               CGALSegment3;
//typedef CGAL::Segment_2<CGAL_Kernel>               CGALSegment2;
//typedef CGAL::Line_3<CGAL_Kernel>                  CGALLine3;
//typedef CGAL::Vector_3<CGAL_Kernel>                CGALVector3;
//typedef CGAL::Direction_3<CGAL_Kernel>             CGALDirection3;
//typedef CGAL::Plane_3<CGAL_Kernel>                 CGALPlane3;
//typedef CGAL::Sphere_3<CGAL_Kernel>                CGALSphere3;

//typedef CGAL::Object                          CGALObject;

//typedef CGAL::Aff_transformation_3<CGAL_Kernel>    CGALTransform;

// ***************************************************
// ***************************************************

// Type definitions used by class CrackGeoEng3D and class CrackManager

// No longer used typedefs

//typedef CGAL::Tetrahedron_3<CGAL_Kernel>       CGALTet;
//typedef CGAL::Triangle_2<CGAL_Kernel>          CGALTria2;
//typedef CGAL::Triangle_3<CGAL_Kernel>          CGALTria3;
//typedef CGAL::Triangle_3<CGAL_Kernel>          CGALTriangle3;

//typedef std::list<CGALPoint2>             Container;
//typedef CGAL::Polygon_2<CGAL_Kernel,Container> CGALPolygon2;

// NOTE: CGAL::Orientation is a typedef to CGAL::Sign
typedef CGAL::Sign                        Orientation;


// ***************************************************
// ***************************************************
// Type definitions for CGAL tetrahedralization

// No longer used typedefs

// Used by CrackManager3D

//typedef CGAL::Triangulation_cell_base_3<CGAL_Kernel>        Cb;
//typedef CGAL::Triangulation_vertex_base_3<CGAL_Kernel>      Vb;
//typedef CGAL::Triangulation_data_structure_3<Vb,Cb>    TDS;
//typedef CGAL::Triangulation_data_structure_3<Vb,Cb>    CGAL_Triangulation_data_structure_3;


//typedef CGAL::Delaunay_triangulation_3<CGAL_Kernel>     Triangulation;
//typedef CGAL::Delaunay_triangulation_3<CGAL_Kernel>     CGAL_Triangulation;

//typedef Triangulation::Point                       PointDel3D; //JG: maybe erase not necessary we already have CGALPoint3
//typedef CGAL::Delaunay_triangulation_3<CGAL_Kernel>::Point    CGAL_PointDel3D;

// ***************************************************
// ***************************************************
// Delaunay tetrahedralization with cell and vertex info: Used by CrackManager3D2
//
// NOTE: Use CGAL_Kernel_for_Del instead of CGAL_Kernel in the declarations below.
//
// ALWAYS use CGAL_EP_IC_Kernel for Delaunay. Elsewhere in ISET, CGAL_EP_EC_Kernel can be used.
// CGAL_EP_EC_Kernel significatly slows down Delaunay tetratization and may lead to exceptions
// being thrown by some CGAL functions used by Delaunay tetratization.
//
// NOTE: CGAL_Kernel (the CGAL kernel seleted through Makefile variable
// USE_CGAL_EXACT_ARITHMETIC) can be CGAL_EP_EC_Kernel even if
// CGAL_EP_IC_Kernel is used here

typedef CGAL_EP_IC_Kernel   CGAL_Kernel_for_Del;

class DelVertexInfo;
class DelCellInfo;

// DONE: Update these typedefs following guideline (*) above.

//typedef CGAL::Point_3<CGAL_Kernel_for_Del>                 CGAL_Del_Point3;
//typedef CGAL::Triangle_3<CGAL_Kernel_for_Del>              CGAL_Del_Tria3;

// Typedef only used here
//typedef CGAL::Triangulation_cell_base_3<CGAL_Kernel_for_Del>        CGAL_Del_Cb;
template <class Kernel_Type>
using CGAL_Del_Cb = CGAL::Triangulation_cell_base_3<Kernel_Type>;

//typedef CGAL::Triangulation_vertex_base_3<CGAL_Kernel_for_Del>      CGAL_Del_Vb;
template <class Kernel_Type>
using CGAL_Del_Vb = CGAL::Triangulation_vertex_base_3<Kernel_Type>;

//typedef CGAL::Triangulation_cell_base_with_info_3<DelCellInfo, CGAL_Kernel_for_Del, CGAL_Del_Cb>       CGAL_CellBaseWInfo;
template <class Kernel_Type>
using CGAL_CellBaseWInfo = CGAL::Triangulation_cell_base_with_info_3<DelCellInfo, Kernel_Type, CGAL_Del_Cb<Kernel_Type> >; 

//typedef CGAL::Triangulation_vertex_base_with_info_3<DelVertexInfo, CGAL_Kernel_for_Del, CGAL_Del_Vb>   CGAL_VertexBaseWInfo;
template <class Kernel_Type>
using CGAL_VertexBaseWInfo = CGAL::Triangulation_vertex_base_with_info_3<DelVertexInfo, Kernel_Type, CGAL_Del_Vb<Kernel_Type> >;

//typedef CGAL::Triangulation_data_structure_3<CGAL_VertexBaseWInfo, CGAL_CellBaseWInfo>  CGAL_Triangulation_data_structure_3_WInfo;
//typedef CGAL::Triangulation_data_structure_3<CGAL_VertexBaseWInfo, Cb>  CGAL_Triangulation_data_structure_3_WInfo;
template <class Kernel_Type>
using CGAL_Triangulation_data_structure_3_WInfo = CGAL::Triangulation_data_structure_3<CGAL_VertexBaseWInfo<Kernel_Type>, CGAL_CellBaseWInfo<Kernel_Type> >;

//typedef CGAL::Delaunay_triangulation_3<CGAL_Kernel_for_Del, CGAL_Triangulation_data_structure_3_WInfo > CGAL_Delaunay_Triangulation_WInfo;
template <class Kernel_Type>
using CGAL_Delaunay_Triangulation_WInfo = CGAL::Delaunay_triangulation_3<Kernel_Type, CGAL_Triangulation_data_structure_3_WInfo<Kernel_Type> >;

// Used to convert objects like Point_3<CGAL_Kernel> to Point_3<CGAL_Kernel_for_Del> and vice-versa.
// If both kernels are the same, these converters do nothing.
//typedef CGAL::Cartesian_converter<CGAL_Kernel, CGAL_Kernel_for_Del> CGAL_converter_Kernel_to_Del_Kernel;
//typedef CGAL::Cartesian_converter<CGAL_Kernel_for_Del, CGAL_Kernel> CGAL_converter_Del_Kernel_to_Kernel;
//
// C++11 style for consistency with previous template declarations
using CGAL_converter_Kernel_to_Del_Kernel = CGAL::Cartesian_converter<CGAL_Kernel, CGAL_Kernel_for_Del>;
using CGAL_converter_Del_Kernel_to_Kernel = CGAL::Cartesian_converter<CGAL_Kernel_for_Del, CGAL_Kernel>;

// ***************************************************
// ***************************************************
// Delaunay triangulation with cell and vertex info: Used by CrackManager2D
//
// NOTE: Use CGAL_Kernel_for_Del instead of CGAL_Kernel in the declarations below.

// Typedef only used here
template <class Kernel_Type>
using CGAL_Del_Fb_2D = CGAL:: Triangulation_face_base_2<Kernel_Type>;

template <class Kernel_Type>
using CGAL_Del_Vb_2D = CGAL::Triangulation_vertex_base_2<Kernel_Type>;

template <class Kernel_Type>
using CGAL_FaceBaseWInfo_2D = CGAL::Triangulation_face_base_with_info_2<DelCellInfo, Kernel_Type, CGAL_Del_Fb_2D<Kernel_Type> >;

template <class Kernel_Type>
using CGAL_VertexBaseWInfo_2D = CGAL::Triangulation_vertex_base_with_info_2<DelVertexInfo, Kernel_Type, CGAL_Del_Vb_2D<Kernel_Type> >;

template <class Kernel_Type>
using CGAL_Triangulation_data_structure_2D_WInfo = CGAL::Triangulation_data_structure_2<CGAL_VertexBaseWInfo_2D<Kernel_Type>, CGAL_FaceBaseWInfo_2D<Kernel_Type> >;

template <class Kernel_Type>
using CGAL_Delaunay_Triangulation_2D_WInfo = CGAL::Delaunay_triangulation_2<Kernel_Type, CGAL_Triangulation_data_structure_2D_WInfo<Kernel_Type> >;


// ***************************************************
// ***************************************************

// Type definitions used by class CrackManager3D2

// NOTE: Use style: CGAL_TypeName or CGAL_Type_Name

// NOTE: Must use exact constructions for AABB intersection computation.
// Otherwise round-off may lead to _large_ errors!
//
// However, it is NOT straightforward to use CGAL_EP_EC_Kernel for AABB tree
// and CGAL_EP_IC_Kernel for CGAL_Kernel. AABB tree uses a Polyhedron surface
// (which uses CGAL_Kernel (the CGAL kernel seleted through Makefile variable
// USE_CGAL_EXACT_ARITHMETIC)).
//
// Thus, we currently MUST use CGAL_Kernel_for_AABB == CGAL_Kernel
//
// Therefore, in order to use CGAL_EP_EC_Kernel for AABB tree, USE_CGAL_EXACT_ARITHMETIC = YES 
// MUST be set at file Solver.mkdefs

typedef CGAL_Kernel   CGAL_Kernel_for_AABB;

// Replace in code: 
// CGAL_AABB_Polyhedron  with  CGAL_Polyhedron_3<CGAL_Kernel_for_AABB>
//typedef CGAL::Polyhedron_3<CGAL_Kernel_for_AABB, Customized_items>  CGAL_AABB_Polyhedron;

//typedef CGAL::AABB_face_graph_triangle_primitive<CGAL_Polyhedron_3<CGAL_Kernel_for_AABB>>  CGAL_AABB_Primitive;
template <class Kernel_Type>
using CGAL_AABB_Primitive = CGAL::AABB_face_graph_triangle_primitive< CGAL_Polyhedron_3<Kernel_Type> >;

//typedef CGAL::AABB_traits<CGAL_Kernel_for_AABB, CGAL_AABB_Primitive>                       CGAL_AABB_Traits;
template <class Kernel_Type>
using CGAL_AABB_Traits = CGAL::AABB_traits<Kernel_Type, CGAL_AABB_Primitive<Kernel_Type> >;

//typedef CGAL::AABB_tree<CGAL_AABB_Traits>                            CGAL_AABB_Tree;
template <class Kernel_Type>
using CGAL_AABB_Tree = CGAL::AABB_tree< CGAL_AABB_Traits<Kernel_Type> >; 

//typedef CGAL_AABB_Tree::Object_and_primitive_id                      CGAL_AABB_Tree_Object_and_primitive_id;
template <class Kernel_Type>
using CGAL_AABB_Tree_Object_and_primitive_id = typename CGAL_AABB_Tree<Kernel_Type>::Object_and_primitive_id;

//typedef CGAL_AABB_Tree::Point_and_primitive_id                       CGAL_AABB_Tree_Point_and_primitive_id;
template <class Kernel_Type>
using CGAL_AABB_Tree_Point_and_primitive_id = typename CGAL_AABB_Tree<Kernel_Type>::Point_and_primitive_id;

//typedef CGAL_AABB_Tree::Primitive_id                                 CGAL_AABB_Tree_Primitive_id;
template <class Kernel_Type>
using CGAL_AABB_Tree_Primitive_id = typename CGAL_AABB_Tree<Kernel_Type>::Primitive_id;

//typedef CGAL::Point_3<CGAL_Kernel_for_AABB>                 CGAL_AABB_Point3;
// Replace in code: CGAL_AABB_Point3 with CGAL::Point_3<CGAL_Kernel_for_AABB>

//typedef CGAL::Segment_3<CGAL_Kernel_for_AABB>               CGAL_AABB_Segment3;
// Replace in code: CGAL_AABB_Segment3 with CGAL::Segment_3<CGAL_Kernel_for_AABB>

// Used to convert objects like Point_3<CGAL_Kernel> to Point_3<CGAL_Kernel_for_AABB> and vice-versa.
// If both kernels are the same, these converters do nothing.
typedef CGAL::Cartesian_converter<CGAL_Kernel, CGAL_Kernel_for_AABB> CGAL_converter_Kernel_to_AABB_Kernel;
typedef CGAL::Cartesian_converter<CGAL_Kernel_for_AABB, CGAL_Kernel> CGAL_converter_AABB_Kernel_to_Kernel;


// ***************************************************
// ***************************************************
// Added by Piyush Gupta
//typedef std::list< CGAL::Triangle_3<CGAL_Kernel_for_AABB> >::iterator                      CGAL_AABB_Tria_Facet_Iterator;
template <class Kernel_Type>
using CGAL_AABB_Tria_Facet_Iterator = typename std::list< CGAL::Triangle_3<Kernel_Type> >::iterator;

//typedef CGAL::AABB_triangle_primitive<CGAL_Kernel_for_AABB, CGAL_AABB_Tria_Facet_Iterator> CGAL_AABB_Tria_Facets_Primitive;
template <class Kernel_Type>
using CGAL_AABB_Tria_Facets_Primitive = CGAL::AABB_triangle_primitive<Kernel_Type, CGAL_AABB_Tria_Facet_Iterator<Kernel_Type> >;

//typedef CGAL::AABB_traits<CGAL_Kernel_for_AABB, CGAL_AABB_Tria_Facets_Primitive>           CGAL_AABB_Tria_Facets_Traits;
template <class Kernel_Type>
using CGAL_AABB_Tria_Facets_Traits = CGAL::AABB_traits<Kernel_Type, CGAL_AABB_Tria_Facets_Primitive<Kernel_Type> >;

//typedef CGAL::AABB_tree<CGAL_AABB_Tria_Facets_Traits>                                      CGAL_AABB_Tria_Facets_Tree;
template <class Kernel_Type>
using CGAL_AABB_Tria_Facets_Tree = CGAL::AABB_tree<CGAL_AABB_Tria_Facets_Traits<Kernel_Type> >;


// ***************************************************
// Used by class CrackFrontGeoEng3D
//
//typedef CGAL::AABB_halfedge_graph_segment_primitive< CGAL_Polyhedron_3<CGAL_Kernel_for_AABB> >  CGAL_AABB_Poly_Segment_Primitive;
template <class Kernel_Type>
using CGAL_AABB_Poly_Segment_Primitive = CGAL::AABB_halfedge_graph_segment_primitive< CGAL_Polyhedron_3<Kernel_Type> >;

//typedef CGAL::AABB_traits<CGAL_Kernel_for_AABB, CGAL_AABB_Poly_Segment_Primitive>           CGAL_AABB_Poly_Segment_Traits;
template <class Kernel_Type>
using CGAL_AABB_Poly_Segment_Traits = CGAL::AABB_traits< Kernel_Type, CGAL_AABB_Poly_Segment_Primitive<Kernel_Type> >; 

//typedef CGAL::AABB_tree<CGAL_AABB_Poly_Segment_Traits>                                      CGAL_AABB_Poly_Segment_Tree;
template <class Kernel_Type>
using CGAL_AABB_Poly_Segment_Tree = CGAL::AABB_tree< CGAL_AABB_Poly_Segment_Traits<Kernel_Type> >;


// No longer used by class CrackFrontGeoEng3D
//typedef std::vector<CGAL_AABB_Segment3>::iterator                                      CGAL_AABB_Segment_Iterator;
//typedef CGAL::AABB_segment_primitive<CGAL_Kernel_for_AABB, CGAL_AABB_Segment_Iterator> CGAL_AABB_Segment_Primitive;
//typedef CGAL::AABB_traits<CGAL_Kernel_for_AABB, CGAL_AABB_Segment_Primitive>           CGAL_AABB_Segment_Traits;
//typedef CGAL::AABB_tree<CGAL_AABB_Segment_Traits>                                      CGAL_AABB_Segment_Tree;

// ***************************************************
// ***************************************************

// NOTE: Algorithm used in mesh simplication does not work with
// CGAL::Exact_predicates_exact_constructions_kernel 
//
// CGAL functions used for crack surface remesh may not work with
// CGAL::Exact_predicates_exact_constructions_kernel

//#if USE_CGAL_EXACT_ARITHMETIC
//typedef CGAL_EP_EC_Kernel    CGAL_Kernel_for_remesh;
//#else
//typedef CGAL_EP_IC_Kernel    CGAL_Kernel_for_remesh;
//#endif

typedef CGAL_EP_IC_Kernel   CGAL_Kernel_for_remesh;

// Typedefs used for Crack Surface remesh
//
// NOTE: Should use style: CGAL_TypeName or CGAL_Type_Name

// TBD: Update these typedefs following guideline (*) above.

//typedef CGAL::Point_3<CGAL_Kernel_for_remesh>                          CGAL_remesh_Point3;

//typedef CGAL::Polyhedron_3<CGAL_Kernel_for_remesh, Customized_items>   CGAL_remesh_Polyhedron;
// Use instead: CGAL_Polyhedron_3<CGAL_Kernel_for_remesh>

// traits class. Typedef only used here
typedef CGAL::Robust_circumcenter_traits_3<CGAL_Kernel_for_remesh>     CGAL_remesh_RC_Kernel;

// vertex and cell types. Typedef only used here
typedef CGAL::Triangulation_vertex_base_3<CGAL_Kernel_for_remesh>      CGAL_remesh_Vb;

// Typedefs only used here
typedef CGAL::Triangulation_vertex_base_with_info_3<PolySurfVertexInfo, CGAL_Kernel_for_remesh, CGAL_remesh_Vb>  CGAL_VertexBaseWInfo2;
typedef CGAL::Surface_mesh_vertex_base_3<CGAL_remesh_RC_Kernel, CGAL_VertexBaseWInfo2>              CGAL_remesh_Vb_WInfo;
typedef CGAL::Surface_mesh_cell_base_3<CGAL_remesh_RC_Kernel>                                       CGAL_remesh_SCb;
//
#ifdef CGAL_VERSION_4_4
typedef CGAL::Triangulation_cell_base_with_circumcenter_3<CGAL_remesh_RC_Kernel, CGAL_remesh_SCb>   CGAL_remesh_Cb_with_circumcenter;
#else
typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<CGAL_remesh_RC_Kernel, CGAL_remesh_SCb>   CGAL_remesh_Cb_with_circumcenter;
#endif

// data structure for triangulation with info for surf mesh. Typedef only used here
typedef CGAL::Triangulation_data_structure_3<CGAL_remesh_Vb_WInfo, CGAL_remesh_Cb_with_circumcenter>  CGAL_remesh_Tds_WInfo;

// Delaunay triangulation with vertex info for surf mesh
typedef CGAL::Delaunay_triangulation_3<CGAL_remesh_RC_Kernel, CGAL_remesh_Tds_WInfo>      CGAL_remesh_Tr;
//
// Can use instead:  CGAL_remesh_triangulation_3<CGAL_remesh_RC_Kernel>
template <class Kernel_Type>
using CGAL_remesh_triangulation_3 = CGAL::Delaunay_triangulation_3<Kernel_Type, CGAL_remesh_Tds_WInfo>;

// container complex 2 in triangulation 3
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<CGAL_remesh_Tr>                   CGAL_remesh_C2t3;

//typedef CGAL::Complex_2_in_triangulation_3<CGAL_remesh_Tr> CGAL_remesh_C2t3;

// criteria for surf mesh
typedef CGAL::Surface_mesh_default_criteria_3<CGAL_remesh_Tr>                             CGAL_remesh_Criteria;

// Typedef only used here
typedef CGAL::Mesh_3::Robust_intersection_traits_3<CGAL_Kernel_for_remesh>                CGAL_remesh_IGT;

// DONE: Does not work: Try using CGAL_Kernel equal to kernel of Polyhedral surface
//typedef CGAL::AABB_polyhedral_oracle<Polyhedron, CGAL_Kernel, CGAL_remesh_IGT>   CGAL_remesh_Input_RSurface;
typedef CGAL::AABB_polyhedral_oracle<CGAL_Polyhedron_3<CGAL_Kernel_for_remesh>, CGAL_Kernel_for_remesh, CGAL_remesh_IGT>   CGAL_remesh_Input_RSurface;

// Used to convert objects like Point_3<CGAL_Kernel> to Point_3<CGAL_Kernel_for_remesh> and vice-versa.
// If both kernels are the same, these converters do nothing.
typedef CGAL::Cartesian_converter<CGAL_Kernel, CGAL_Kernel_for_remesh>           CGAL_converter_Kernel_to_remesh_Kernel;
typedef CGAL::Cartesian_converter<CGAL_Kernel_for_remesh, CGAL_Kernel>           CGAL_converter_remesh_Kernel_to_Kernel;
typedef CGAL::Cartesian_converter<CGAL_Kernel_for_remesh, CGAL_remesh_IGT>       CGAL_remesh_Converter;

// ***************************************************
// ***************************************************

// CGAL Kernel used when computing intersections: Must be CGAL_EP_EC_Kernel even if user sets USE_CGAL_EXACT_ARITHMETIC == false
typedef CGAL_EP_EC_Kernel   CGAL_Kernel_for_intersect;

typedef CGAL::Cartesian_converter<CGAL_Kernel, CGAL_Kernel_for_intersect>        CGAL_converter_Kernel_to_intersect_Kernel;
typedef CGAL::Cartesian_converter<CGAL_Kernel_for_intersect, CGAL_Kernel>        CGAL_converter_intersect_Kernel_to_Kernel;

// ***************************************************
// ***************************************************

#endif
