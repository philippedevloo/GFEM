// ---- testCrackMg2.C ---- Sat Oct 29 14:52:32 CDT 2011
// (c) Copyright: C. Armando Duarte, 2002-2004. All rights reserved
// See README.CopyRight file for further details.
//


#define DBUG_ON 2

#include <cstdlib>
#include <string>
#include <cstdio>
#include <typeinfo>
#include <algorithm>
#include <set>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "ConcreteFormedArray1d.h"

#include "config.h"
#include "config_cgal.h"
#include "error.h"
#include "geomesh.h"
#include "geonod.h"
#include "compmesh.h"
#include "linanalysis.h"
#include "graphmesh.h"
#include "phdatafilesol.h"
#include "matelas3d.h"
#include "compel.h"
#include "compnod.h"
#include "crackmg3d2_cgal.h"
#include "crgeoeng3d2_cgal.h"
#include "compgeom.h"
#include "geoeltet10.h"

#include "edgeintersectdata_cgal.h"
#include "edge.h"
#include "crcutelinfo2_cgal.h"
#include "gnodcrackeleminfo_cgal.h"
#include "crfrontgeoeng3d_cgal.h"

using std::cerr;
using std::cout;
using std::endl;

/*
Triangulation_data_structure_3< Triangulation_vertex_base_3<DelaunayTriangulationTraits_3>, 
                                Triangulation_cell_base_3<DelaunayTriangulationTraits_3> >

CGAL::Delaunay_triangulation_3<Kernel, TriangulationDataStructure_3, LocationPolicy>
*/

// *************************************************************
// *************************************************************

typedef CGAL::Point_3<CGAL_Kernel_for_AABB>                 CGAL_AABB_Point3;

namespace{

  // Define (scoped) pointer to objects.
  // The object pointed by a smart pointer is destroyed when 
  // the smart pointer is, i.e., at the end of execution in this case.

  boost::scoped_ptr<Analysis> scp_analysis(0);

#if 0 // not being used
  std::string createNewName( char* datafile, 
                             const std::string suffix = std::string("dx"), 
                             const std::string altsuffix = std::string("xx") );

  void test_orientation_calc(CGAL_AABB_Tree<CGAL_Kernel_for_AABB>& aabb_tree_cr_surf);

  void test_intersection_w_segm_calc(CGAL_AABB_Tree<CGAL_Kernel_for_AABB>& aabb_tree_cr_surf);

  void test_interserction_w_edge_segments(const CGAL_AABB_Tree<CGAL_Kernel_for_AABB>& aabb_tree_cr_surf, 
                                          const boost::shared_ptr<GeoMesh>& shp_gMesh,
                                          std::set<EdgeIntersectData*, 
                                          EdgeUnsortedID<EdgeIntersectData*> >&
                                          tet_edges_w_intersections);

  int create_CrackCutElInfo( const boost::shared_ptr<GeoMesh>& shp_gMesh, 
                             const std::set<EdgeIntersectData*, EdgeUnsortedID<EdgeIntersectData*> >&
                             tet_edges_w_intersections,
                             std::vector<CrackCutElInfo2 *>& crack_elems_info,
                             std::vector<GNodCrackedElemsInfo *>& gnod_cracked_elems);

  void create_GNodCrackedElemsInfo(const EdgeIntersectData& edge,
                                   const GeoNod *p_gnode_0, const GeoNod *p_gnode_1,
                                   std::vector<GNodCrackedElemsInfo *>& gnod_cracked_elems);
  
  void drawToDX_CrackCutElInfo2(const CrackCutElInfo2* p_cut_el_info2,
                                const boost::shared_ptr<CrackGeoEng3D2>& shp_crack_geoeng,
                                const boost::shared_ptr<GeoMesh>& shp_gMesh, 
                                const std::vector<GNodCrackedElemsInfo *>& gnod_cracked_elems,
                                const std::string& file_name,
                                const bool draw_elem_vertex_intersections = true);

  
  void test_intersection_w_crack_edges(const boost::shared_ptr<GeoMesh>& shp_gMesh,
                                       const boost::shared_ptr<CrackGeoEng3D2>& shp_crack_geoeng,
                                       std::vector<CrackCutElInfo2 *>& crack_elems_info);

  void create_Polyhedron_from_GeoEl(const GeoEl* p_geoel,
                                    Polyhedron& geo_polyhedron);

  void test_intersection_w_crack_vertices(const boost::shared_ptr<GeoMesh>& shp_gMesh,
                                          const boost::shared_ptr<CrackGeoEng3D2>& shp_crack_geoeng,
                                          std::vector<CrackCutElInfo2 *>& crack_elems_info);


  void test_Delaunay_Triangulation(const std::vector<GNodCrackedElemsInfo *>& gnod_cracked_elems,
                                   const std::vector<CrackCutElInfo2 *>& crack_elems_info);
#endif
} // end of namespace

  // ****************************************************
  // ****************************************************

int main(int argc, char** argv) {

  if (argc != 3) {
    cerr << "\n" << argv[0] 
         << " requires a .grf and a .crf given on the command line" << endl;
    return EXIT_FAILURE;
  }

  // GeoMesh
  boost::shared_ptr<GeoMesh> shp_gMesh(new GeoMesh);
  shp_gMesh->setName("Geomesh");

  // Computational stuff
  boost::shared_ptr<CompMesh> shp_cMesh( new CompMesh(*shp_gMesh) );

  // Read mesh:
  PhDataFileSol dataFile( argv[1] );

  if ( !dataFile.read( *shp_gMesh, *shp_cMesh) ) {
    cerr << "\n*** invalid data file " << argv[1] << " ***" << endl;
    return EXIT_FAILURE;
  } 

  std::string mesh_name = std::string( argv[1] ) + "_" + std::string( argv[2] );
  shp_gMesh->setName( mesh_name ); 

  cout << "\nNumber of GeoNods = " << shp_gMesh->numNodes() << endl;
  cout << "Number of GeoEls  = " << shp_gMesh->numElems() << endl;
  cout << "Number of CompEls = " << shp_cMesh->numElems() << endl;
  cout << "Number of Materials = " << shp_cMesh->numMats() << endl;

  //  cout << "gMesh after reading data file: " << endl;
  //  cout << *shp_gMesh << endl;

  //  cout << "cMesh after reading data file: " << endl;
  //  cout << *shp_cMesh << endl;

  // ------------------------------------

  // CGAL stuff
  // **********
 
  std::cout << "\n\n***Tests with AABB Tree and CrackSurface **\n\n";

  // CrackManager3D2 creates a CrackGeoEng3D2 for crack surface provided at
  // command line

  // CAD NOTE: This surface has NO temporary advancement. It may be a good
  // idea to create another surface with the temporary advancement and use it
  // to compute intersections with elements. Thus, we would have two surfaces:
  // One with temporary advance and another one without it!

  const int cm_ID = 1;
  const std::string geo_eng_crf_file( argv[2] );
  CrackManager3D2  crackmg_v2(cm_ID, shp_gMesh /* problem geomesh */, geo_eng_crf_file);


#if 0
  CrackManager3D  crackmg_v1(cm_ID, shp_gMesh /* problem geomesh */, geo_eng_crf_file,
                             false /* isLocal */, true /* sliceTetBranchFnc */);

  boost::shared_ptr<CrackGeoEng3D2> shp_crack_geoeng_v1( crackmg_v1.crackGeoEng() );

  // Temporary advancement of crack front to fully cut elements partially cut
  // by actual crack surface

  ConcreteRigidArray1d<realType,3> displ_v1;
  displ_v1(0) = 0.6; displ_v1(1) = 0.0; displ_v1(2) = 0.0;

  shp_crack_geoeng_v1->advanceCrackFront(displ_v1, 
                                         true /*tempAdvance*/);

  // openDX output with crack surface
  std::string outFileName_v1 = argv[2];
  outFileName_v1 += "_wTempAdvance_v1.dx";
  shp_crack_geoeng_v1->openOutFile( outFileName_v1  );
  shp_crack_geoeng_v1->drawCrack();

#endif

  
  boost::shared_ptr<CrackGeoEng3D2> shp_crack_geoeng( crackmg_v2.crackGeoEng() );

  // openDX output with crack surface
  std::string outFileName = argv[2];
  outFileName += ".vtk";
  shp_crack_geoeng->drawCrackSurface( outFileName, EVTKStyle );

#if 0
  //
  // constructs CGAL_AABB_Segment_Tree AABB tree for crack front and
  // the internal search tree for efficient distance computations.
  //
  std::vector<CGAL::Segment_3<CGAL_Kernel>> frontSegments;
  shp_crack_geoeng->reportCrackFront(frontSegments);

  PRINT_ELEMENTS(frontSegments,"\nALL crack front segments:\n");

  //  CGAL_AABB_Segment_Tree cr_front_tree( frontSegments.begin(), frontSegments.end() );
  CGAL_AABB_Segment_Tree cr_front_tree;
  cr_front_tree.rebuild( frontSegments.begin(), frontSegments.end() );

  cr_front_tree.accelerate_distance_queries();

  // computes the closest point from a point query 
  //
  // Using one_tet.grf planar_crack.crf for these tests
  //
  CGAL::Point_3<CGAL_Kernel_for_AABB> point_query(2.5, 2.0, 0.25);
  CGAL::Point_3<CGAL_Kernel_for_AABB> closest_pt = cr_front_tree.closest_point(point_query);

  std::cout << "\nNumber of crack fronts in tree: " << cr_front_tree.size();

  std::cout << "\nQuery point is: " << point_query
            << "\nClosest point is: " << closest_pt;

  CGAL_Kernel_for_AABB::FT squared_distance =  cr_front_tree.squared_distance(point_query);
  std::cout << "\nDistance from query point to front: " 
            << std::sqrt( iset::to_realType(squared_distance) ) << std::endl;

  CGAL_AABB_Segment_Tree::Point_and_primitive_id  pt_and_primit = 
    cr_front_tree.closest_point_and_primitive (point_query);
  std::cout << "\nClosest point is: " << pt_and_primit.first;

  // closest primitive id (== a segment handle)
  CGAL_AABB_Segment_Primitive::Id  hd_segment = pt_and_primit.second;

  CGAL::Segment_3<CGAL_Kernel>* pt_segment = &(*pt_and_primit.second);
  
  std::cout << "\nClosest front edge: " << *hd_segment << endl;
  std::cout << "\nClosest front edge: " << *pt_segment << endl;

  // Checking rebuild:
  cr_front_tree.rebuild( frontSegments.begin(), frontSegments.end() );
  cr_front_tree.accelerate_distance_queries();
#endif

#if 1
  //
  // constructs CGAL_AABB_Poly_Segment_Tree<CGAL_Kernel_for_AABB> AABB tree for crack front
  // and the internal search tree for efficient distance computations.
  //
  std::vector<Polyhedron::Halfedge_handle> frontHalfEdges;
  shp_crack_geoeng->reportCrackFront(frontHalfEdges);

  CGAL_AABB_Poly_Segment_Tree<CGAL_Kernel_for_AABB> cr_front_poly_tree;
  cr_front_poly_tree.rebuild( frontHalfEdges.begin(), frontHalfEdges.end(),
                              *(shp_crack_geoeng->crackSurface()) );

  cr_front_poly_tree.accelerate_distance_queries();

  // computes the closest point from a point query 
  //
  // Using one_tet.grf planar_crack.crf for these tests
  //
  CGAL_AABB_Point3 point_query(2.5, 2.0, 0.25);
  CGAL_AABB_Point3 closest_pt = cr_front_poly_tree.closest_point(point_query);

  std::cout << "\n\ntestCrackMg2: Number of crack fronts in tree: " << cr_front_poly_tree.size();

  std::cout << "\nQuery point is: " << point_query
            << "\nClosest point is: " << closest_pt;

  CGAL_Kernel_for_AABB::FT squared_distance =  cr_front_poly_tree.squared_distance(point_query);
  std::cout << "\nDistance from query point to front: " 
            << std::sqrt( iset::to_realType(squared_distance) ) << std::endl;

  CGAL_AABB_Poly_Segment_Tree<CGAL_Kernel_for_AABB>::Point_and_primitive_id  pt_and_primit = 
    cr_front_poly_tree.closest_point_and_primitive (point_query);
  std::cout << "\nClosest point is: " << pt_and_primit.first;

  // closest primitive id (== a halhedge handle)
  Polyhedron::Halfedge_handle hd_halfedge = pt_and_primit.second.halfedge();
  CGAL_AABB_Poly_Segment_Primitive<CGAL_Kernel_for_AABB>::Id  hd_halfedge_a = pt_and_primit.second;

  CGAL::Segment_3<CGAL_Kernel> segment( hd_halfedge->vertex()->point(), 
                        hd_halfedge->opposite()->vertex()->point() );
  std::cout << "\nClosest front edge: " << segment << endl;

  // Checking rebuild:
  cr_front_poly_tree.rebuild( frontHalfEdges.begin(), frontHalfEdges.end(),
                              *(shp_crack_geoeng->crackSurface()) );
  cr_front_poly_tree.accelerate_distance_queries();
  

#endif

  //
  // Same as above but using CrackFrontGeoEng3D object in GeoEng3D
  //
  boost::shared_ptr<CrackFrontGeoEng3D> shp_crack_front_geoeng = 
    shp_crack_geoeng->crackFrontGeoEng();

  std::cout << "\n\ntestCrackMg2: Using CrackFrontGeoEng3D:"
            << "\nNumber of crack fronts in tree: "
            << shp_crack_front_geoeng->numFrontEdges();

  // computes the closest point from a point query
  Polyhedron::Halfedge_handle hd_hedge_w_close_pt = 0;

  // Compute closest crack front point to query_pt.
  shp_crack_front_geoeng->getClosestPtCrackFront(point_query,
                                                 closest_pt,
                                                 hd_hedge_w_close_pt);

  std::cout << "\nQuery point is: " << point_query
            << "\nClosest point is: " << closest_pt;

  CGAL::Segment_3<CGAL_Kernel> hedge( hd_hedge_w_close_pt->vertex()->point(), 
                      hd_hedge_w_close_pt->opposite()->vertex()->point() );
  std::cout << "\nClosest front edge: " << hedge << endl;

  // Compute squared distance of query_pt w.r.t. crack front(s).
  realType squared_dist = 0.;
  shp_crack_front_geoeng->getSquaredDistCrackFront(point_query,
                                                   squared_dist);

  std::cout << "\nDistance from query point to front: " 
            << std::sqrt( iset::to_realType(squared_distance) ) << std::endl;

#if 1
  // Temporary advancement of crack front to fully cut elements partially cut
  // by actual crack surface

  ConcreteRigidArray1d<realType,3> displ;
  displ(0) = 0.5; displ(1) = 0.0; displ(2) = 0.0;

  shp_crack_geoeng->advanceCrackFront(displ, 
                                      true /*tempAdvance*/);
  outFileName = argv[2];
  outFileName += "_wTempAdvance.vtk";
  shp_crack_geoeng->drawCrackSurface( outFileName, EVTKStyle );

#endif


#define TEST_BASIC_FNS 0
#if TEST_BASIC_FNS


  // Create AABB tree for crack surface

  // DONE: Add this functionality to crack Geo Eng

  // See:
  // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/AABB_tree/Chapter_main.html#Subsection_63.3.2

  Polyhedron* p_crack_surface = shp_crack_geoeng->crackSurface();

  Polyhedron::Facet_iterator 
    cr_facet_pos = p_crack_surface->facets_begin(),
    cr_facet_end = p_crack_surface->facets_end();

  CGAL_AABB_Tree<CGAL_Kernel_for_AABB> aabb_tree_cr_surf( cr_facet_pos, cr_facet_end );

  // Tests with hardwired query points and segments:
  // -----------------------------------------------

  //  test_intersection_w_segm_calc( aabb_tree_cr_surf );

  //  test_orientation_calc( aabb_tree_cr_surf );


  // Tests for Tets read from .grf
  // ----------------------------


  // DONE: Add this functionality to crack mg


  // Container with pointers to all EdgeIntersectData objects for this mesh
  std::set<EdgeIntersectData*, EdgeUnsortedID<EdgeIntersectData*> > tet_edges_w_intersections;

  test_interserction_w_edge_segments(aabb_tree_cr_surf, shp_gMesh,
                                     tet_edges_w_intersections);

  // - DONE: Update scaling factor used later on to check equality of intersections
  // points on facet/volume with existing intersection points.  We use this
  // representative value to set the scaling factor for intersection point
  // objects associated with intersections at facets/volume.  This allow us to
  // have a single class CompareIntersectPoints to be used with Edges, facets
  // and volume!
  //
  // - 1) DONE: Return number of new intersection stored  <-- EdgeIntersectData::storeIntersection
  //
  // - 2) VOID (should insert all): Insert in tet_edges_w_intersections only edges that have new intersections

  // DONE: Add this functionality to crack mg

  // DONE: Create CrackCutElInfo2 for each Tet cut by crack surface
  //
  // Collect info from EdgeIntersectData into CrackCutElInfo objects. Check if
  // edges of element have intersections with crack surface

  // CrackCutElInfo2 * are stored by index of geoEl from global problem.
  // This is needed to provide direct access (no searches) to a
  // CrackCutElInfo for a given global problem or local integration element.
  std::vector<CrackCutElInfo2 *> crack_elems_info;

  // DONE: Create GNodCrackedElemsInfo for each node belonging to elements that
  // are cut by the crack surface
  // Collect info from EdgeIntersectData into GNodCrackedElemsInfo objects.
  //
  // GNodCrackedElemsInfo* are stored by the index of GeoNods from global problem.
  std::vector<GNodCrackedElemsInfo *> gnod_cracked_elems;

  int num_elems_edge_intersect = create_CrackCutElInfo( shp_gMesh, tet_edges_w_intersections,
                                                        crack_elems_info, gnod_cracked_elems );
  std::cout << "\nnum_elems_edge_intersect = " << num_elems_edge_intersect << std::endl;


  // DONE for CM2 (below): Compute distance and orientation of GeoNods
  // belonging to tets that are cut by crack surface
  //
  // NOTE: May need also this info for GeoEls that are neighbor to elements cut by crack surface(?)


  // DONE: Check if this is done at CM2 

  iset::cleanUpContainer(tet_edges_w_intersections, "Cleaning up Container with EdgeIntersectData");


  // DONE: Add this functionality to CM2

  // Intersection of crack surface edges with _faces_ of tets: 
  // see: FiberGeoEng3d::storeIntersectionWithTet, Aditya's compute_facet_intersections
  //

  // DONE: - Consistency check: There can be no new intersection with tet vertices when processing 
  // intersections with facets/volume ! 

  // DONE: Create tetrahedron from each GeoElTets, use their triangular faces to define AABB tree.
  // Are triangles closed or open? Are facets in a polyhedron closed or open?
  //
  // see http://www.cgal.org/Manual/latest/doc_html/cgal_manual/AABB_tree/Chapter_main.html#Subsection_63.3.1


  // VOID (done together with other edges) : Intersection of crack front edges with Tets

  test_intersection_w_crack_edges(shp_gMesh, shp_crack_geoeng, crack_elems_info);

  // DONE: Compute intersections of interior of tets with vertex of crack surface
  // see  Aditya's compute_tet_space_intersections, Aditya's check_internal_point
  //
  // - Consistency check: There can be no new intersection with tet vertices when processing 
  // intersections with facets/volume ! 

  test_intersection_w_crack_vertices(shp_gMesh, shp_crack_geoeng,
                                     crack_elems_info);

  PRINT_ELEMENTS_VALUE( crack_elems_info, 
                        "\nALL CrackCutElInfo2 Objects after test_intersection_w_crack_vertices:" );


  // Create a Delaunay triangulation for each CrackCutElInfo2 object
  test_Delaunay_Triangulation(gnod_cracked_elems,
                              crack_elems_info);



  // TBD: Setup enrichment type flags for all tets with CrackCutElInfo2


#if 1
  // DONE: Add this functionality to CM2
  // Dump CrackCutElInfo2 to dx files for visualization
  
  std::vector<CrackCutElInfo2 *>::const_iterator
    cut_info_pos = crack_elems_info.begin(),
    cut_info_end = crack_elems_info.end();

  char elementIndex_chr [10];
  for( ; cut_info_pos != cut_info_end; ++cut_info_pos ) {

    CrackCutElInfo2 *p_crack_cut_info2 = *cut_info_pos;

    if ( !p_crack_cut_info2 )
      continue; // no CrackCutElInfo2 for this element
  
    // transform number to string
    sprintf(elementIndex_chr, "%d", p_crack_cut_info2->geoElProb()->ID() );
    // file name suffix
    std::string file_name = "CrackCutElInfo2_elem_ID_" + std::string(elementIndex_chr) + ".dx"; 
    file_name = createNewName( argv[2], file_name );
    file_name = std::string( argv[1] ) + "_" + file_name;

    drawToDX_CrackCutElInfo2(p_crack_cut_info2,shp_crack_geoeng, shp_gMesh, gnod_cracked_elems,
                             file_name);
  }
#endif


  iset::cleanUpContainer(crack_elems_info, "Cleaning up container with CrackCutElInfo2");

  iset::cleanUpContainer(gnod_cracked_elems, "Cleaning up container with GNodCrackedElemsInfo");

  // -----------------------------------
#endif // TEST_BASIC_FNS

#define REFINE_MESH 0
#if REFINE_MESH

  // Crack front for problem: edge_cracked_crack_propag_3D_hpGFEM_mesh_11x11x3.tcl
  std::vector<Segment> segmentList;
  Segment A_B;

  A_B.A(0) = 2.1;
  A_B.A(1) = 2.0;
  A_B.A(2) = 0.0;

  A_B.B(0) = 2.1;
  A_B.B(1) = 2.0;
  A_B.B(2) = 1.0;

  segmentList.push_back( A_B );

  int n_ref = 16;
  for(int iref=0; iref<n_ref; ++iref) {
    shp_gMesh->refineIfIntersect(segmentList);
  }

  realType maxEdgeLen, minEdgeLen, maxEdgeLen_ratio, maxRadiusLen;
  shp_gMesh->elemSizes( maxEdgeLen, minEdgeLen, 
                        maxEdgeLen_ratio, maxRadiusLen );

  cout << "\nRange of element sizes for the entire mesh:\n" 
       << "Maximum Edge Size = " << maxEdgeLen << endl
       << "Minimum Edge Size = " << minEdgeLen << endl
       << "Ratio = " << maxEdgeLen/minEdgeLen << endl
       << "Maximum Edge lenght ratio in an element = " << maxEdgeLen_ratio
       << endl
       << "Max. radius-length (Q) = "<< maxRadiusLen << endl 
       << "this value should be bounded as small as possible. \n" 
       << "Q = 0.612 is optimal, if we do not consider sliver elems.\n" 
       << endl;
  
#endif


  // Repeat above tests now using CrackManager3D2 object
  // ===================================================

  // This was used for early tests: Call functions one-by-one
  // SEE below another test where we call processCrack instead.

#define TEST_CM2_INDIVIDUAL_FNS 0
#if TEST_CM2_INDIVIDUAL_FNS

  bool calc_crack_edge_intersec_all_elems  = true;
  bool calc_crack_vertices_intersect_elems = true;

  // Compute intersection between all element edges in GeoMeh shp_gMesh_prob
  // and the crack surface.
  crackmg_v2.calcElemEdgeIntersecSurface();

  // Create CrackCutElInfo2 for each GeoEl cut by crack surface.
  // Create GNodCrackedElemsInfo for nodes that are on the crack surface (only).
  // Return value = number of elements with edges cut by or on crack surface
  int num_elems_edge_inters = crackmg_v2.createCrackCutElInfo();

  std::cout << "\nNum of elems with edge intersections = " << num_elems_edge_inters << std::endl;

  // TBD: Instead of not calling this method at all, pass a flag for it to use
  // crack front edges only!

  if (calc_crack_edge_intersec_all_elems) {
    // For all elements in GeoMesh, compute intersections between edges of crack
    // surface and element faces
    crackmg_v2.calcCrackEdgeIntersecAllElems();
  }

  // TBD: Instead of not calling this method at all, pass a flag for it to use
  // crack front vertices only!

  if (calc_crack_vertices_intersect_elems) {
    // For all elements in GeoMesh, compute intersections between element
    // interior and crack surface vertices
    crackmg_v2.calcCrackVerticesIntersectElems();
  }

  // Loop over GeoEls of problem GeoMesh that HAVE CrackCutElInfo2 and set,
  // for their GeoNods, the orientation and unsigned distance w.r.t. crack
  // surface.

  // NOTE: May need also this info for GeoEls that are neighbor to elements
  // cut by crack surface(?)
  crackmg_v2.setNodalOrientAndDistToCrackSurf();


  // Create a Delaunay triangulation for each primitive element with a
  // CrackCutElInfo2 object
  crackmg_v2.cutElements();


  // TBD: Setup enrichment type flags for all tets with CrackCutElInfo2
  

  // Draws all CompNod enrichments marked (but not necessarily used) in vector
  // p_crack_mg->compNod_prob_enrichment[.] and all intersections (between
  // crack surface and GeoMesh) stored in all CrackCutElInfo2 objetcs of
  // p_crack_mg. If draw_elem_vertex_intersections = false, intersections at
  // GeoEl vertices are NOT drawn
  // If draw_to_single_file = true, data from all CrackCutElInfo2 are combined
  // in a single file. Otherwise, one file per CrackCutElInfo2 will be created.

  bool draw_elem_vertex_intersections = true;
  bool draw_to_single_file = false;
  std::string file_name_prefix = std::string( argv[1] ) + "_" + std::string( argv[2] );
  crackmg_v2.drawAllCrackCutElInfo(file_name_prefix, 
                                   draw_elem_vertex_intersections, draw_to_single_file);


  draw_to_single_file = true;
  std::string full_file_name = file_name_prefix + "_CrackCutElInfo2_ALL_elems.dx";  
  crackmg_v2.drawAllCrackCutElInfo(full_file_name, 
                                   draw_elem_vertex_intersections, draw_to_single_file);

  // Clean up containers with intersection data.
  // Use this to same memory after these data are no longer needed.

  crackmg_v2.cleanUpIntersectData();

#endif

  // ------------------------------------

#define USE_CM2_a  1
#if USE_CM2_a

  // Another CM2 object to perform test of higher level methods (like processCrack)
  // The temp advancement of crack surface and most everything is done at processCrack

  CrackManager3D2*  p_crackmg_v2_a = 
    new CrackManager3D2( 2 /*cm_ID*/, shp_gMesh /* problem geomesh */, geo_eng_crf_file);

  // Apply temporary advance at crack front
  // Find elements cut by crack surface in shp_gMesh_prob, cut elements,
  // create integration elements, set-up CompNod enrichments, etc.
  p_crackmg_v2_a->processCrack();

#if 0
  // Create list of points representing the crack surface. The points are
  // intersections stored in CrackCutElInfo2 and GNodCrackedElemsInfo objects.
  std::vector< std::pair<CGAL::Point_3<CGAL_Kernel>, DelSurfVertexInfo> > surf_vertexes;
  p_crackmg_v2_a->vertexesForCrackSurface(surf_vertexes);
#endif

  // Draw all CompNod enrichments marked (but not necessarily used) in vector
  // p_crack_mg->compNod_prob_enrichment[.] and all intersections (between
  // crack surface and GeoMesh) stored in all CrackCutElInfo2 objetcs of
  // p_crack_mg. If draw_elem_vertex_intersections = false, intersections at
  // GeoEl vertices are NOT drawn
  // If draw_to_single_file = true, data from all CrackCutElInfo2 are combined
  // in a single file. Otherwise, one file per CrackCutElInfo2 will be created.

  bool draw_elem_vertex_intersections = true;
  bool draw_to_single_file = false;
  std::string file_name_prefix = std::string( argv[1] ) + "_" + std::string( argv[2] );
  p_crackmg_v2_a->drawAllCrackCutElInfo(file_name_prefix, 
                                   draw_elem_vertex_intersections, draw_to_single_file);


  draw_to_single_file = true;
  std::string full_file_name = file_name_prefix + "_CrackCutElInfo2_ALL_elems.dx";  
  p_crackmg_v2_a->drawAllCrackCutElInfo(full_file_name, 
                                     draw_elem_vertex_intersections, draw_to_single_file);

  // Draw Computational Crack Surface

  std::string comp_surface_file = file_name_prefix + "_CompSurface_above.dx";
  CrackMn3DGraphDX2::EsurfaceToDrawType surf_to_draw = CrackMn3DGraphDX2::EAboveSurface;
  p_crackmg_v2_a->drawComputationalCrackSurface( comp_surface_file, surf_to_draw);

  comp_surface_file = file_name_prefix + "_CompSurface_below.dx";
  surf_to_draw = CrackMn3DGraphDX2::EBelowSurface;
  p_crackmg_v2_a->drawComputationalCrackSurface( comp_surface_file, surf_to_draw);


  // Clean up containers with intersection data.
  // Use this to same memory after these data are no longer needed.

  // p_crackmg_v2_a->cleanUpIntersectData();

#endif
 // ------------------------------------

#define USE_CM1_a  0
#if USE_CM1_a

  // Another CM1 object to perform test of higher level methods (like processCrack)
  // The temp advancement of crack surface and most everything is done at processCrack

  CrackManager3D*  p_crackmg_v1_a = 
    new CrackManager3D(cm_ID, shp_gMesh /* problem geomesh */, geo_eng_crf_file,
                       false /* isLocal */, true /* sliceTetBranchFnc */);

  p_crackmg_v1_a->setOptStraightFront(true);

  // Apply temporary advance at crack front
  // Find elements cut by crack surface in shp_gMesh_prob, cut elements,
  // create integration elements, set-up CompNod enrichments, etc.
  p_crackmg_v1_a->processCrack( 1 /*max_num_terms_stepfn*/,
                                true /*useBranchFunctions*/,
                                0 /*p_localProblem*/,
                                false /*locProbCrackMgWithAllElems*/);
#endif


#if 0
  // linanalysis.
  scp_analysis.reset( new LinAnalysis( shp_cMesh, shp_gMesh ) );
  
#if USE_CM1_a
  scp_analysis->store( p_crackmg_v1_a );

#elif USE_CM2_a 
  scp_analysis->store( p_crackmg_v2_a );

#endif

  ConcreteFormedArray1d<std::string> scalNames(1), vecNames(1);
  scalNames(0) = "Solution";
  vecNames(0) = "Derivative";
  const std::string dx_filename = createNewName( argv[1] );
  short dimension = shp_gMesh->dimension();

  // Graph mesh
  scp_analysis->createGraphMesh(scalNames, vecNames, dx_filename, 
                                dimension);
  //  cout << *(scp_analysis->graphMesh()) << endl;


  scp_analysis->assemble();
  scp_analysis->factAndSolve();

  scp_analysis->postProcess( 0 /*resolution*/ );
#endif
  // ------------------------------------


  return EXIT_SUCCESS;
}
// **************************************************************
// **************************************************************

namespace{
#if 0 // not being used
std::string createNewName( char* datafile, const std::string suffix, 
                           const std::string altsuffix ) {

  // creates a new name for a file based on the name of datafile

  std::string filename, basename, extname, dxname;

  filename = datafile;

  // search period in file name
  std::string::size_type idx = filename.find('.');
  if( idx == std::string::npos ) {
    // file name does not contain any period
    dxname = filename + '.' + suffix;
  }
  else {
    /* split file name into base name and extension
     * base name contains all characters before the period
     * extension contains all characters after the period
     */
    
    basename = filename.substr(0, idx);
    extname = filename.substr(idx+1);
    if ( extname.empty() ) {
      // contains period but no extension: append suffix
      dxname = filename;
      dxname += suffix;
    }
    else if ( extname == suffix ) {
      // replace extension suffix with altsuffix
      dxname = filename;
      dxname.replace( idx+1, extname.size(), altsuffix );
    }
    else {
      // replace any extension with suffix
      dxname = filename;
      dxname.replace( idx+1, std::string::npos, suffix);
    }
  }

  DBGPRNT1(cout << "\nfilename = " << filename << " new name = " << dxname 
           << endl);

  return dxname;
}
#endif
// **************************************************************
// **************************************************************
#if 0 // not used
 void test_orientation_calc(CGAL_AABB_Tree<CGAL_Kernel_for_AABB>& aabb_tree_cr_surf) {

   // Tests of algorithm to find:
   //   - minimum distance of a point w.r.t. crack surface
   //   - orientation of a point w.r.t. crack surface

   std::cout << "\n\n***Tests with Orientation of Points w.r.t. CrackSurface **" 
             << std::endl;

   // computes internal KD-tree data structure to accelerate distance queries
   aabb_tree_cr_surf.accelerate_distance_queries();

   std::vector<CGAL::Point_3<CGAL_Kernel>> all_query_points;

   // query points
   all_query_points.push_back( CGAL::Point_3<CGAL_Kernel>(0.5, 3.0, 0.3) );

   all_query_points.push_back( CGAL::Point_3<CGAL_Kernel>(2.0, 1.9, 0.3) );
   all_query_points.push_back( CGAL::Point_3<CGAL_Kernel>(2.0, 1.9, 0.0) );

   all_query_points.push_back( CGAL::Point_3<CGAL_Kernel>(2.0, 2.0, 0.3) );
   all_query_points.push_back( CGAL::Point_3<CGAL_Kernel>(2.0, 2.0, 0.0) );

   all_query_points.push_back( CGAL::Point_3<CGAL_Kernel>(2.0, 2.0, -1.0e-14) );

   all_query_points.push_back( CGAL::Point_3<CGAL_Kernel>(2.0, 2.0, -1.0) );

   std::vector<CGAL::Point_3<CGAL_Kernel>>::iterator 
     query_pt_pos = all_query_points.begin(),
     query_pt_end = all_query_points.end();
   
   //   CGAL_converter_Kernel_to_AABB_Kernel converter_K_to_AABB_K;
   //   CGAL_converter_AABB_Kernel_to_Kernel converter_AABB_K_to_K;

   for( ; query_pt_pos != query_pt_end; ++query_pt_pos ) {

     CGAL::Point_3<CGAL_Kernel> query_pt = *query_pt_pos;
     
     std::cout << "\nQuery Point: " << query_pt << std::endl;

     // computes squared distance from query_pt
     Kernel::FT sqd = aabb_tree_cr_surf.squared_distance(query_pt);
     std::cout << "squared distance from query point: " << sqd << std::endl;

     // computes closest point
     CGAL::Point_3<CGAL_Kernel> closest = aabb_tree_cr_surf.closest_point(query_pt);
     std::cout << "closest point on surface: " << closest << std::endl;

     // computes closest point and primitive id (== a facet handle)
     CGAL_AABB_Tree_Point_and_primitive_id<CGAL_Kernel_for_AABB>  pt_and_primit = 
       aabb_tree_cr_surf.closest_point_and_primitive(query_pt);
     std::cout << "closest point on surface: " << pt_and_primit.first << std::endl;

     Polyhedron::Facet_handle hd_facet = pt_and_primit.second; // closest primitive id (== a facet handle)

     std::cout << "\nID of facet with closest point: " 
               << hd_facet->getFacetID() << std::endl;
     
     // Find orientation of query point w.r.t. crack surface

     // get orientation based on closest facet 

     // WARNING: This test may give false CGAL::COPLANAR on surfaces with
     // angles between facets smaller or equal to 90 degrees. This is not
     // physically consistent though. Furthermore, the code can detect this
     // inconsistency when computing orientation of GeoEl nodes.  Those nodes
     // that are COPLANAR should have been flag as such by the cutting
     // algorithm.  Thus, the orientation computed here should never return
     // COPLANAR for nodes not found to be on the crack surface by the cutting
     // algorithm.

     // How to FIX: If an invalid COPLANAR result is found, build aa_bb tree
     // for element using only the facets that the element intersect. It is
     // less likey we will get COPLANAR again (but possible).
     //
     // Better: Apply small random perturbation on the position of the query
     // point and search again for facet with closest point. This seems quite
     // a robust solution!
     //
     // We can write a function to check the orientation of nodes after
     // creating descendants: Using a face of a descendant that has only one
     // face on creack surface, use one of the facets that contain this face
     // to compute the orientation of the other node of the descendant (the
     // one that is not on the crack surface. Can check against all facets
     // that contain the descendant face. This is similar to one of the
     // functions JP wrote.

     // After computing intersections; 
     //  - loop over CrackCutElInfo and mark nodes that are COPLANAR
     //  - Loop over GeoEls with CrackCutElInfo and set orientation of their nodes 
     //    that are not COPLANAR
     //  - When computing orientation of descend, use orientation of nodes that 
     //    have it; use neighboring info to transmit this info to other descend; 
     //    Use list of faces potentially on crack to set orientation of neighbors 
     //    across crack surface

     // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Kernel_23_ref/Function_orientation.html#Index_anchor_464
     CGAL::Orientation orient_of_point = 
       CGAL::orientation( hd_facet->halfedge()->vertex()->point(),
                          hd_facet->halfedge()->next()->vertex()->point(),
                          hd_facet->halfedge()->next()->next()->vertex()->point(),
                          query_pt );

     std::cout << "\nOrientation of query point w.r.t. facet with closet point: " 
               << orient_of_point << std::endl;
   }
 }
 // **************************************************************
 // **************************************************************

 void test_intersection_w_segm_calc(CGAL_AABB_Tree<CGAL_Kernel_for_AABB>& aabb_tree_cr_surf) {

   std::cout << "\n\n***Intersection of Segments and CrackSurface **\n\n";

#if 0
   // test creation of AABB tree

   CGAL::Point_3<CGAL_Kernel> p(1.0, 0.0, 0.0);
   CGAL::Point_3<CGAL_Kernel> q(0.0, 1.0, 0.0);
   CGAL::Point_3<CGAL_Kernel> r(0.0, 0.0, 1.0);
   CGAL::Point_3<CGAL_Kernel> s(0.0, 0.0, 0.0);
   
   Polyhedron polyhedron;
   polyhedron.make_tetrahedron(p, q, r, s);
   
   // constructs AABB tree
   CGAL_AABB_Tree<CGAL_Kernel_for_AABB> aabb_tree(polyhedron.facets_begin(),polyhedron.facets_end());
#endif

   // Define some query segments and test intersection
   // Crack surface from: planar_crack.crf

   // Construct several segments to test special cases. Put them in a
   // vectors and loop over this vector below

   // NOTE: CGAL::Segment_3<CGAL_Kernel> is topologically closed
   std::vector<CGAL::Segment_3<CGAL_Kernel>> all_query_segments;

   // constructs segment query perpendicular to crack surface, intersection at
   // single facet
   all_query_segments.push_back( CGAL::Segment_3<CGAL_Kernel>(CGAL::Point_3<CGAL_Kernel>(0.5, 3.0, 0.3), 
                                              CGAL::Point_3<CGAL_Kernel>(0.5, 0.0, 0.3)) );

   // possibly intersect two facets
   all_query_segments.push_back( CGAL::Segment_3<CGAL_Kernel>(CGAL::Point_3<CGAL_Kernel>(0.5, 3.0, 0.5), 
                                              CGAL::Point_3<CGAL_Kernel>(0.5, 0.0, 0.5)) );

   // possibly intersect at a vertex with four facets
   all_query_segments.push_back( CGAL::Segment_3<CGAL_Kernel>(CGAL::Point_3<CGAL_Kernel>(1.0, 3.0, 0.0), 
                                              CGAL::Point_3<CGAL_Kernel>(1.0, 0.0, 0.0)) );

   // possibly segment intersection with facet 2, source is intersection
   all_query_segments.push_back( CGAL::Segment_3<CGAL_Kernel>(CGAL::Point_3<CGAL_Kernel>(0.9, 2.0, 0.3), 
                                              CGAL::Point_3<CGAL_Kernel>(0.9, 2.0, 2.0)) );

   // possibly segment intersections with facets 2 and 3
   all_query_segments.push_back( CGAL::Segment_3<CGAL_Kernel>(CGAL::Point_3<CGAL_Kernel>(0.5, 2.0, -0.3), 
                                              CGAL::Point_3<CGAL_Kernel>(0.5, 2.0, 2.0)) );

   // possibly segment and point intersections with facets 
   all_query_segments.push_back( CGAL::Segment_3<CGAL_Kernel>(CGAL::Point_3<CGAL_Kernel>(1.0, 2.0, -0.3), 
                                              CGAL::Point_3<CGAL_Kernel>(1.0, 2.0, 2.0)) );

   // possibly segment and point intersections with facets. source and target
   // are intersections
   all_query_segments.push_back( CGAL::Segment_3<CGAL_Kernel>(CGAL::Point_3<CGAL_Kernel>(1.0, 2.0, 0.3), 
                                              CGAL::Point_3<CGAL_Kernel>(1.0, 2.0, 0.9)) );

   std::vector<CGAL::Segment_3<CGAL_Kernel>>::iterator 
     query_seg_pos = all_query_segments.begin(),
     query_seg_end = all_query_segments.end();
   
   for( ; query_seg_pos != query_seg_end; ++query_seg_pos ) {
     
     CGAL::Segment_3<CGAL_Kernel> segment_query = *query_seg_pos;

     std::cout << "\nSegment query: " << segment_query << std::endl;

     // tests intersections with segment query
     if( aabb_tree_cr_surf.do_intersect(segment_query) )
       std::cout << "intersection(s) found" << std::endl;
     else
       std::cout << "NO intersection found" << std::endl;

     // computes #intersections with segment query
     int num_intersect = aabb_tree_cr_surf.number_of_intersected_primitives(segment_query);
     std::cout << "number of intersection(s) found = " << num_intersect << std::endl;

#if 0
     // computes first encountered intersection with segment query
     // (generally a point)
     boost::optional< CGAL_AABB_Tree_Object_and_primitive_id<CGAL_Kernel_for_AABB> > intersection =
       aabb_tree_cr_surf.any_intersection(segment_query);

     if( intersection ) {
       // gets intersection object
       CGAL_AABB_Tree_Object_and_primitive_id<CGAL_Kernel_for_AABB> op = *intersection;
       CGAL::Object object = op.first;

       // check if it is a point
       const CGAL::Point_3<CGAL_Kernel>* p_point;
       if(  p_point = CGAL::object_cast<CGAL::Point_3<CGAL_Kernel>>(&object) ) {
         std::cout << "first encountered intersection object is a point" << std::endl;
       std::cout << (*p_point) << std::endl;
       }
       // TBD: Check if it is a segment otherwise ERROR
     }

     // computes all intersected primitives with segment query as primitive ids
     std::list<CGAL_AABB_Tree_Primitive_id> all_primitives;
     aabb_tree_cr_surf.all_intersected_primitives(segment_query, std::back_inserter(all_primitives));

#endif

     // computes all intersections with segment query as pairs <object, primitive_id>
     // primitive id == a facet handle

     // Object used to process intersections and discard repeated ones
     EdgeIntersectData edge( segment_query );

     std::list< CGAL_AABB_Tree_Object_and_primitive_id<CGAL_Kernel_for_AABB> > all_intersections;
     aabb_tree_cr_surf.all_intersections( edge.segment(), std::back_inserter(all_intersections) );

     std::list< CGAL_AABB_Tree_Object_and_primitive_id<CGAL_Kernel_for_AABB> >::iterator
       inters_pos = all_intersections.begin(),
       inters_end = all_intersections.end();
     
     for( ; inters_pos != inters_end; ++inters_pos) {
       
       // Object, which can be a point or a segment
       CGAL::Object object = inters_pos->first;
       // check if it is a point
       const CGAL::Point_3<CGAL_Kernel>* p_point=0;
       const CGAL::Segment_3<CGAL_Kernel>* p_segm=0;
       if( (p_point = CGAL::object_cast<CGAL::Point_3<CGAL_Kernel>>(&object)) ) {
         std::cout << "\nEncountered intersection object is a point" << std::endl;
         std::cout << (*p_point) << std::endl;
       }
       else if ( (p_segm = CGAL::object_cast<CGAL::Segment_3<CGAL_Kernel>>(&object)) ) {
         // Check if it is a segment otherwise error out
         
         std::cout << "\nEncountered intersection object is a segment" << std::endl;
         std::cout << (*p_segm) << std::endl;
       }
       else {
         std::cerr << "\nInvalid object returned by aaa_bb tree" << std::endl;
       }

       // primitive id (== a facet handle) of intersected facet
       Polyhedron::Facet_handle hd_facet = inters_pos->second;
       std::cout << "ID of facet with intersection object: " 
                 << hd_facet->getFacetID() << std::endl;

       // Process intersections: 
       //  - Discard repeated ones; 
       //  - Discard those that are at end points of segm (store this info);
       //  - set flags if end-points are on the surface
       //  - process end points of segments as two points
       
       edge.storeIntersection( *inters_pos );

     } // next intersecion


     // Print Edge: Number of intersections etc.
     cout << "\nEdge after storing all intersections: " << edge << endl;


   } // next query segment
 }
 // **************************************************************
 // **************************************************************

 void test_interserction_w_edge_segments(const CGAL_AABB_Tree<CGAL_Kernel_for_AABB>& aabb_tree_cr_surf, 
                                         const boost::shared_ptr<GeoMesh>& shp_gMesh,
                                         std::set<EdgeIntersectData*, EdgeUnsortedID<EdgeIntersectData*> >&
                                         tet_edges_w_intersections){

   // Compute intersection between all edges in the GeoMeh and the crack surface
   // store output in tet_edges_w_intersections

   cout << "\n***test_interserction_w_edge_segments***" << endl;

   ERRCHK( tet_edges_w_intersections.size(), 
           "Container for EdgeIntersectData is not empty on entry!");

   const GeoMesh::elemContainerType* p_geoElContainer = shp_gMesh->elemContn();
   GeoMesh::elemContainerType::const_iterator 
     epos = p_geoElContainer->begin(),
     epos_end = p_geoElContainer->end();
   GeoEl* p_gel_prob;
   GeoNod *p_gnode_0, *p_gnode_1;

   std::list< CGAL_AABB_Tree_Object_and_primitive_id<CGAL_Kernel_for_AABB> > all_intersections;
   EdgeIntersectData *p_edgeData = new EdgeIntersectData(0, 0);

   for(; epos != epos_end; ++epos) {

     p_gel_prob = (*epos);
     if ( !p_gel_prob->active() )
       continue;
     
     // Loop over edges of this GelEl
     int num_edge = p_gel_prob->numEdges();
     for(Subscript iedg=0; iedg<num_edge; ++iedg) {

       p_gnode_0 = p_gel_prob->nodePtr( p_gel_prob->edgeNode(iedg,0) );
       p_gnode_1 = p_gel_prob->nodePtr( p_gel_prob->edgeNode(iedg,1) );
       
       DBGPRNT(cout << "\nElem ID = " << p_gel_prob->ID() << "\nEdge Index: " << iedg);
       DBGPRNT(cout << "\nGeoNods on edge: [" << p_gnode_0->ID() << ", " 
               << p_gnode_1->ID() << "]" << endl);

       if ( !p_edgeData )
         p_edgeData = new EdgeIntersectData(p_gnode_0->ID(), p_gnode_1->ID());
       else {
         // Re-cycle available Edge that was not inserted in container bcz it
         // did not intersect crack surface
         p_edgeData->setNodeID(0, p_gnode_0->ID());
         p_edgeData->setNodeID(1, p_gnode_1->ID());
         // Edge to be re-cycled should be empty otherwise it should be in container!
         ERRCHK( p_edgeData->numIntersections(), "\nRe-cycling Edge with intersections!");
       }

       // check if edge [ p_gnode_0->ID() -- p_gnode_1->ID() ] has been processed:
       // check if it is in tet_edges_w_intersections container
       //
       // NOTE: edge [ p_gnode_1->ID() -- p_gnode_0->ID() ] IS equal to 
       //       edge [ p_gnode_0->ID() -- p_gnode_1->ID() ]
       // according to sorting criteria EdgeUnsortedID !!
       if ( tet_edges_w_intersections.find( p_edgeData ) != 
            tet_edges_w_intersections.end() ) {
         DBGPRNT(cout << "\nEdge already in container" << endl);
         continue;
       }

       // Process Edge:
       DBGPRNT(cout << "\nNew Edge: Compute Intersections" << endl);
       
       // Set end-points of EdgeIntersectData (only IDs were set above)
       CGAL::Segment_3<CGAL_Kernel> edge( iset::ConcRigArr1dToCGAL::Point_3<CGAL_Kernel>( p_gnode_0->coor() ),
                          iset::ConcRigArr1dToCGAL::Point_3<CGAL_Kernel>( p_gnode_1->coor() ) );
       p_edgeData->setSegment( edge );

       // computes all intersections of edge with crack surface as pairs <object, primitive_id>
       // primitive id == a facet handle
       
       all_intersections.clear();
       aabb_tree_cr_surf.all_intersections( edge, std::back_inserter(all_intersections) );

       std::list< CGAL_AABB_Tree_Object_and_primitive_id<CGAL_Kernel_for_AABB> >::iterator
         inters_pos = all_intersections.begin(),
         inters_end = all_intersections.end();
     
       int num_new_intersect = 0;
       for( ; inters_pos != inters_end; ++inters_pos) {

         // Process intersections: 
         //  - Discard repeated ones; 
         //  - Discard those that are at end points of segm (store this info);
         //  - set flags if end-points are on the surface
         //  - process end points of segments as two points
       
         num_new_intersect += p_edgeData->storeIntersection( *inters_pos );

       } // next intersection
       DBGPRNT(cout << "\nnum_new_intersect = " << num_new_intersect << endl);

       // sanity check
       ERRCHK( num_new_intersect != p_edgeData->numIntersections(),
               "num_new_intersect != p_edgeData->numIntersections() ??");

       // Insert this new Edge in container even if no intersections found
       // (otherwise, another element using this edge will try to find
       // intersections again).
       // 
       if ( (tet_edges_w_intersections.insert( p_edgeData ) ).second ) {
         // a new edge was inserted in the container.
         p_edgeData = NULL;
       } else {
         // There is a matching edge in the container. It should have been
         // found above!
         ERREXIT("\nThere is a matching edge in the container!:");
       }

     } // next edge of this element
   } // next GeoEl

   // Last created Edge may have not been inserted in Container
   delete p_edgeData;

   // Print all objects in container
   PRINT_ELEMENTS_VALUE(tet_edges_w_intersections, "\nEdgeIntersectData for all GelEls:");

 }
 // **************************************************************
 // **************************************************************

  int create_CrackCutElInfo( const boost::shared_ptr<GeoMesh>& shp_gMesh, 
                             const std::set<EdgeIntersectData*, EdgeUnsortedID<EdgeIntersectData*> >&
                             tet_edges_w_intersections,
                             std::vector<CrackCutElInfo2 *>& crack_elems_info,
                             std::vector<GNodCrackedElemsInfo *>& gnod_cracked_elems) {

    // Create CrackCutElInfo2 for each GeoEl cut by crack surface
    //
    // Collect info from EdgeIntersectData into CrackCutElInfo2 objects.
    //
    // CrackCutElInfo2 pointers are stored by index of GeoEls from global
    // problem.  This is needed to provide direct access (no searches) to a
    // CrackCutElInfo for a given global problem or local integration element.

    // Create GNodCrackedElemsInfo for each node belonging to elements that
    // are cut by the crack surface
    //
    // Collect info from EdgeIntersectData into GNodCrackedElemsInfo objects.
    //
    // GNodCrackedElemsInfo pointers are stored by the index of GeoNods from
    // global problem.

    // Return value = number of elements with edges cut by crack surface

    cout << "\n***create_CrackCutElInfo**" << endl;

    // GeoMesh must be consolidated since we use element index below
    ERRCHK( !shp_gMesh->consolidated(),
            "create_CrackCutElInfo: GeoMesh must be consolidated");

    int num_elem_intersect_crack = 0;

    ERRCHK( crack_elems_info.size(), 
          "create_CrackCutElInfo: crack_elems_info not empty");
    crack_elems_info.resize(shp_gMesh->numElems(), NULL);

    ERRCHK( gnod_cracked_elems.size(),
            "create_CrackCutElInfo: gnod_cracked_elems not empty");
    gnod_cracked_elems.resize(shp_gMesh->numNodes(), NULL);


    const GeoMesh::elemContainerType* p_geoElContainer = shp_gMesh->elemContn();
    GeoMesh::elemContainerType::const_iterator 
      epos = p_geoElContainer->begin(),
      epos_end = p_geoElContainer->end();
    GeoEl* p_gel_prob =  NULL;
    GeoNod *p_gnode_0, *p_gnode_1;

    // Object used for searches
    EdgeIntersectData edge_for_searches(0, 0);

    CrackCutElInfo2* p_crack_cut_info = NULL;

    for(; epos != epos_end; ++epos) {

      p_gel_prob = (*epos);
      if ( !p_gel_prob->active() )
        continue;
     
      // Loop over edges of this GelEl
      p_crack_cut_info = NULL;
      int num_edge = p_gel_prob->numEdges();
      int num_intersect = 0;

      for(Subscript iedg=0; iedg<num_edge; ++iedg) {

        p_gnode_0 = p_gel_prob->nodePtr( p_gel_prob->edgeNode(iedg,0) );
        p_gnode_1 = p_gel_prob->nodePtr( p_gel_prob->edgeNode(iedg,1) );
       
        DBGPRNT(cout << "\nElem ID = " << p_gel_prob->ID() << "\nEdge Index: " << iedg);
        DBGPRNT(cout << "\nGeoNods on edge: [" << p_gnode_0->ID() << ", " 
                << p_gnode_1->ID() << "]" << endl);

        // check if an edge with end nodes p_gnode_0->ID() -- p_gnode_1->ID()
        // is in container. This is the case if it intersects the crack surface
       
        edge_for_searches.setNodeID(0, p_gnode_0->ID());
        edge_for_searches.setNodeID(1, p_gnode_1->ID());

        std::set<EdgeIntersectData*, EdgeUnsortedID<EdgeIntersectData*> >::const_iterator edge_pos;
        edge_pos = tet_edges_w_intersections.find( &edge_for_searches );
        
        if ( edge_pos == tet_edges_w_intersections.end() ) {
          // Even edges that do NOT intersect crack surface are in the container
          ERREXIT("create_CrackCutElInfo: Could not find this edge in container!");
        }

        // Number of intersections with this edge (including with end points!)
        int num_intersect = (*edge_pos)->numIntersections();

        DBGPRNT(cout << "\nNum times this edge intersects crack surface: " 
                << num_intersect << "\n");
        if ( !num_intersect )
          continue; // goto next edge of this element

        // Collect intersections from this edge into CrackCutElInfo2 object
        // for p_gel_prob.
        // Mark GeoNods that are on crack surface as well
        if ( !p_crack_cut_info )
          // element may already have CrackCutElInfo2
          p_crack_cut_info = crack_elems_info[ p_gel_prob->index() ];
        if ( !p_crack_cut_info ) {
          // first edge with intersections for this element: Need CrackCutElInfo2
          CrackManager3D2* p_crack_mg = NULL;
          p_crack_cut_info = new CrackCutElInfo2( p_gel_prob, p_crack_mg );
          crack_elems_info[ p_gel_prob->index() ] = p_crack_cut_info;
          ++num_elem_intersect_crack;
        }
        // sanity check
        ERRCHK( p_gel_prob != p_crack_cut_info->geoElProb(),
                "Corrupted crack_elems_info vector!");
          
        // Store data (intersections, crack surface facets) from this edge into
        // CrackCutElInfo2 object for this GeoEl
        int num_new_intersect = p_crack_cut_info->addEdgeIntersectionData( iedg, **edge_pos );
        num_intersect += num_new_intersect;
        DBGPRNT(cout <<"\nNum intersects on edge " << iedg << " : " 
                << num_new_intersect << endl);

        // Mark GeoNods that are on crack surface as well, using info from *edge_pos

        // NOTE: The returned edge may have IDs = [ p_gnode_1->ID(), p_gnode_0->ID() ]
        // since this edge IS equal to [ p_gnode_0->ID(), p_gnode_1->ID() ] !!!

        create_GNodCrackedElemsInfo( **edge_pos, p_gnode_0, p_gnode_1,
                                     gnod_cracked_elems);

      }// next edge of GeoEl
      DBGPRNT(cout <<"\nNum intersects for this element (not counting vertex nodes intersect): "
              << num_intersect << endl);

    } // next GeoEl

    // Print containers  crack_elems_info and gnod_cracked_elems

    //    PRINT_ELEMENTS_VALUE( crack_elems_info, "\nCrackCutElInfo2 Objects:" );

    //    PRINT_ELEMENTS_VALUE( gnod_cracked_elems, "\nGNodCrackedElemsInfo Objects:" );

    return num_elem_intersect_crack;

  }
  // **************************************************************
  // **************************************************************

  void create_GNodCrackedElemsInfo(const EdgeIntersectData& edge,
                                   const GeoNod *p_gnode_0, const GeoNod *p_gnode_1,
                                   std::vector<GNodCrackedElemsInfo *>& gnod_cracked_elems) {

    // Update GNodCrackedElemsInfo objects for GeoNods p_gnode_0 and p_gnode_1
    // using data in edge:
    //  - Set orientation of GeoNods if they are CGAL::COPLANAR
    //  - Store crack facet that has node
    //
    // GNodCrackedElemsInfo* are stored by the index of GeoNods from global problem.
    // 
    // NOTE: p_gnode_0 may correspond to edge.source or edge.target since
    //      edge [ID_0, ID_1] IS equal to edge [ID_1, ID_0]
    // same thing applies to p_gnode_1

    
    // Get GeoNods corresponding to edge.source and edge.target using
    // available IDs in edge object. This avoids searching for geoNods based
    // on their IDs
    const GeoNod 
      *p_gnode_source = p_gnode_0,
      *p_gnode_target = p_gnode_1;

    if ( p_gnode_0->ID() != edge.nodeID(0) ) {
      // (0,1) and (1,0) 
      p_gnode_source = p_gnode_1;
      p_gnode_target = p_gnode_0;

      ERRCHK( p_gnode_0->ID() != edge.nodeID(1), "create_GNodCrackedElemsInfo: Invalid arguments");
      ERRCHK( p_gnode_1->ID() != edge.nodeID(0), "create_GNodCrackedElemsInfo: Invalid arguments");
    }
    // (0,1) and (0,1)
    else if ( p_gnode_1->ID() != edge.nodeID(1) ) {
      ERREXIT("create_GNodCrackedElemsInfo: Invalid arguments");
    }

    // Check if source/target are on crack surface and which facets they belong to

    GNodCrackedElemsInfo *p_gnode_info = NULL;

    if ( edge.isSourceOnSurface() ) {
      p_gnode_info = gnod_cracked_elems[ p_gnode_source->index() ];
      if (!p_gnode_info) {
        p_gnode_info =  new GNodCrackedElemsInfo( p_gnode_source, 0 /*dim_side_w_node*/ );
        gnod_cracked_elems[ p_gnode_source->index() ] = p_gnode_info;
      }

      p_gnode_info->orientation() = CGAL::COPLANAR;
      // Object GNodCrackedElemsInfo will have facet containing node of last
      // edge passed to this method
      p_gnode_info->facetWNode() = edge.facetWSource();

      DBGPRNT(cout << "GeoNod w/ ID : " << p_gnode_source->ID() << " is on crack surf" << endl);
      DBGPRNT(cout << "Facet w/ Node : " << p_gnode_info->facetWNode()->getFacetID() << endl);
    }

    if ( edge.isTargetOnSurface() ) {
      p_gnode_info = gnod_cracked_elems[ p_gnode_target->index() ];
      if (!p_gnode_info) {
        p_gnode_info =  new GNodCrackedElemsInfo( p_gnode_target,  0 /*dim_side_w_node*/ );
        gnod_cracked_elems[ p_gnode_target->index() ] = p_gnode_info;
      }

      p_gnode_info->orientation() = CGAL::COPLANAR;
      p_gnode_info->facetWNode() = edge.facetWTarget();

      DBGPRNT(cout << "GeoNod w/ ID : " << p_gnode_target->ID() << " is on crack surf" << endl);
      DBGPRNT(cout << "Facet w/ Node : " << p_gnode_info->facetWNode()->getFacetID() << endl);
    }     
  }

  // **************************************************************
  // **************************************************************


  void drawToDX_CrackCutElInfo2(const CrackCutElInfo2* p_cut_el_info2,
                                const boost::shared_ptr<CrackGeoEng3D2>& shp_crack_geoeng,
                                const boost::shared_ptr<GeoMesh>& shp_gMesh, 
                                const std::vector<GNodCrackedElemsInfo *>& gnod_cracked_elems,
                                const std::string& file_name,
                                const bool draw_elem_vertex_intersections) {

    // DONE: CAD: TBD: Move this method to CrackMn3DGraphDX2

    // The purpose of this method is to generate a DX output for debugging
    // the element cutting process.
    //
    // This implementation is based on CrackMn3DGraphDX::drawAllCrackCutElInfo
    // Here we draw a single CrackCutElInfo2. The same dx network used for
    // output from drawAllCrackCutElInfo can be used for the output from this method.
    // Try crack_analysis.net

    // TBD: It draws CompNod enrichments marked for the nodes of
    // GeoEl->reference().  NOTE: These enrichments may NOT necessarily be used.
    // 
    // Draw intersections between crack surface and GeoEl. If
    // draw_elem_vertex_intersections = false, intersections at GeoEl vertices
    // are NOT drawn
    //
    // TBD: It also draws the base vectors used by Branch Functions.

    // GNodCrackedElemsInfo* are stored in gnod_cracked_elems by the index of
    // GeoNods from global problem.
    //


    DBGPRNT(cout << "\ndrawToDX_CrackCutElInfo2:" << endl);

    DataFileO odataFile(file_name);
    std::ofstream& outf = odataFile.outFile;

    CompMesh* p_cMesh = shp_gMesh->reference();
    const CompMesh::nodContainerType*  p_compNods  = p_cMesh->nodContn(); 

    CompMesh::nodContainerType::const_iterator  
      now_compNod = p_compNods->begin();
    CompMesh::nodContainerType::const_iterator  
      compNod_end = p_compNods->end();

    //  Subscript nowIndex = 0;

    // Mark CompNods with Branch Fns, Step Fns and Local solutions assigned in
    // CrackCutElInfo but may actually NOT be used!
    std::vector<Subscript> compNodWithStepFnct(0), compNodWithBranchFnct(0), 
      compNodWithLocalSol(0);

    /*
      
      TBD
      if (! shp_crack_geoeng->has_tempCrackSurfaceFacet() ) {
      
      for (; now_compNod != compNod_end ; ++now_compNod ){
      nowIndex = (*now_compNod)->index();
      if( p_crack_mg->compNod_prob_enrichment[nowIndex] == CrackCutElInfo::StepFunction )
      compNodWithStepFnct.push_back( nowIndex );
      else if( p_crack_mg->compNod_prob_enrichment[nowIndex] == 
      CrackCutElInfo::BranchFunction )
      compNodWithBranchFnct.push_back( nowIndex );
      else if( p_crack_mg->compNod_prob_enrichment[nowIndex] == 
      CrackCutElInfo::LocalSolution )
      compNodWithLocalSol.push_back( nowIndex ); // added by DJK
      }

      }
    */

  
    std::vector<Subscript>::const_iterator  compNodIndex, compNodIndex_end;
    CompNod* p_compNod;
    Subscript now_index;

    // print coordinates for those nodes with step function

    Subscript nCompNodsWithStepFnct = compNodWithStepFnct.size();
    
    // Work around to avoid DX Field with empty object: added by DJK
    if ( nCompNodsWithStepFnct == 0 ) {
      nCompNodsWithStepFnct = 1;
      compNodWithStepFnct.push_back((*(p_compNods->begin()))->index());
    }
    
    outf << "# The irregular positions, which are "
         << nCompNodsWithStepFnct << " three-dimensional points \n";
    outf << "object \"compNodsWithStepFnct\" class array type float rank 1 shape 3 items " 
         << nCompNodsWithStepFnct << " data follows \n";
    
    compNodIndex = compNodWithStepFnct.begin();
    compNodIndex_end = compNodWithStepFnct.end();
    
    for( ; compNodIndex != compNodIndex_end; ++compNodIndex ){
      
      now_index = (*compNodIndex);
      p_compNod = p_cMesh->compNodFromArray( now_index );
      
      for(Subscript icoor=0; icoor<3; ++icoor)
        outf << "\t" << p_compNod->geoNod()->coor( icoor );
      outf << endl;
    }
    outf << "attribute \"dep\" string \"positions\" \n\n";
    
    
    // print data to plot glyphs: DX field needs a data object
    outf << "object \"dataForStepFnct\" class array type float rank 1 shape 1 items " 
         << nCompNodsWithStepFnct << " data follows\n";
    for( Subscript i = 0; i < nCompNodsWithStepFnct; i++ ){
      outf << "\t" << 1.0 << "\n";
    }
    outf << "attribute \"dep\" string \"positions\"\n\n";


    // print positions for branchFunctions

    Subscript nCompNodsWithBranchFnct = compNodWithBranchFnct.size();

    // Work around to avoid DX Field with empty object: added by DJK
    if ( nCompNodsWithBranchFnct == 0 ) {
      nCompNodsWithBranchFnct = 1;
      compNodWithBranchFnct.push_back((*(p_compNods->begin()))->index());
    }
    
    outf << "# The irregular positions, which are "
         << nCompNodsWithBranchFnct << " three-dimensional points \n";
    outf << "object \"compNodsWithBranchFnct\" class array type float rank 1 shape 3 items " 
         << nCompNodsWithBranchFnct << " data follows \n";
    
    compNodIndex = compNodWithBranchFnct.begin();
    compNodIndex_end = compNodWithBranchFnct.end();

    for( ; compNodIndex != compNodIndex_end; ++compNodIndex ){
      
      now_index = (*compNodIndex);
      p_compNod = p_cMesh->compNodFromArray( now_index ) ;
      
      for(Subscript icoor=0; icoor<3; ++icoor)
        outf << "\t" << p_compNod->geoNod()->coor( icoor );
      outf << endl;
    }
    outf << "attribute \"dep\" string \"positions\" \n\n";
    
    
    // print data to plot glyphs
    outf << "object \"dataForBranchFnct\" class array type float rank 1 shape 1 items " 
         << nCompNodsWithBranchFnct << " data follows\n";
    for( Subscript i = 0; i < nCompNodsWithBranchFnct; i++ ){
      outf << "\t" << 1.0 << "\n";
    }
    outf << "attribute \"dep\" string \"positions\"\n\n";
    

    // print positions for LocalSolution

    Subscript nCompNodsWithLocalSol = compNodWithLocalSol.size();
    
    // Work around to avoid DX Field with empty object: added by DJK
    if ( nCompNodsWithLocalSol == 0 ) {
      nCompNodsWithLocalSol = 1;
      compNodWithLocalSol.push_back((*(p_compNods->begin()))->index());
    }
    
    outf << "# The irregular positions, which are "
         << nCompNodsWithLocalSol << " three-dimensional points \n";
    outf << "object \"compNodsWithLocalSol\" class array type float rank 1 shape 3 items " 
         << nCompNodsWithLocalSol << " data follows \n";

    compNodIndex = compNodWithLocalSol.begin();
    compNodIndex_end = compNodWithLocalSol.end();
    
    for( ; compNodIndex != compNodIndex_end; ++compNodIndex ){
      
      now_index = (*compNodIndex);
      p_compNod = p_cMesh->compNodFromArray( now_index ) ;
      
      for(Subscript icoor=0; icoor<3; ++icoor)
        outf << "\t" << p_compNod->geoNod()->coor( icoor );
      outf << endl;
    }
    outf << "attribute \"dep\" string \"positions\" \n\n";
    
    
    // print data to plot glyphs
    outf << "object \"dataForLocalSol\" class array type float rank 1 shape 1 items " 
         << nCompNodsWithLocalSol << " data follows\n";
    for( Subscript i = 0; i < nCompNodsWithLocalSol; i++ ){
      outf << "\t" << 1.0 << "\n";
    }
    outf << "attribute \"dep\" string \"positions\"\n\n";
    
    
    // print intersections between GeoEl and crack surface
    
    const std::set<IntersectPointData, CompareIntersectPoints >& 
      intersect_points = *(p_cut_el_info2->intersectPoints());
    
    // whether or not dump intersections between crack surface and element
    // vertices
    
    Subscript n_intersect = intersect_points.size();

    if (draw_elem_vertex_intersections) {
      GNodCrackedElemsInfo *p_gnode_info = NULL;
      const GeoEl* p_geoel = p_cut_el_info2->geoElProb();
      int num_nodes = p_geoel->numNodes();
      
      for(int inod=0; inod<num_nodes; ++inod) {
        p_gnode_info = gnod_cracked_elems[ p_geoel->nodePtr(inod)->index() ];
        if ( p_gnode_info && p_gnode_info->orientation() == CGAL::COPLANAR ) {
          // This GeoNod is on the crack surface
          ++n_intersect;
        }
      }
    }
    
    outf << "# The irregular positions, which are "
         << n_intersect  << " three-dimensional points \n";
    outf << "object \"intersectionPoints\" class array type float rank 1 shape 3 items " 
         << n_intersect << " data follows \n";
    
    std::set<IntersectPointData, CompareIntersectPoints >::const_iterator
      intersect_pos = intersect_points.begin(),
      intersect_end = intersect_points.end();
    
    for( ; intersect_pos != intersect_end; ++intersect_pos ){
      
      const CGAL::Point_3<CGAL_Kernel>& intersectPoint = intersect_pos->intersectPoint();
      for(Subscript icoor=0; icoor<3; ++icoor) {
        outf << "\t" << intersectPoint[icoor];
      }
      outf << endl;
    }
    
    if (draw_elem_vertex_intersections) {
      GNodCrackedElemsInfo *p_gnode_info = NULL;
      const GeoEl* p_geoel = p_cut_el_info2->geoElProb();
      int num_nodes = p_geoel->numNodes();
      
      for(int inod=0; inod<num_nodes; ++inod) {
        p_gnode_info = gnod_cracked_elems[ p_geoel->nodePtr(inod)->index() ];
        if ( p_gnode_info && p_gnode_info->orientation() == CGAL::COPLANAR ) {
          // This GeoNod is on the crack surface
          for(Subscript icoor=0; icoor<3; ++icoor) {
            outf << "\t" << p_geoel->nodePtr(inod)->coor(icoor);
          }
          outf << endl;
        }
      }
    }
    
    outf << "attribute \"dep\" string \"positions\" \n\n";
    
    outf << "object \"dataForIntersections\" class array type float rank 1 shape 1 items " 
         << n_intersect << " data follows\n";
    
    for( Subscript i = 0; i < n_intersect; i++ ){
      outf << "\t" << 1.0 << "\n";
    }
    outf << "attribute \"dep\" string \"positions\"\n\n";

    // plot coordinate axes (base vectors) for the branch function enrichment
    // the CoorSys is plotted at the node instead of the average point position
    // to avoid repeated CoorSys at the same location

    ConcreteRigidArray1d<realType,3> point, closePt;
    ConcreteRigidArray1d<realType,3> closePt_tangent, closePt_normal, closePt_binormal;

    //Added by Piyush Gupta
    // NOTE: Can use either true or false since origin of coord. sys. returned 
    // by crack manager is not dumped to file
    bool averageCrackFrontCoordSystem = false;  

    if ( shp_crack_geoeng->has_tempCrackSurfaceFacet() ){
      // Can not compute base vectors if surface has temporary advance since
      // front is not a border edge in this case (methods used by
      // shp_crack_geoeng->closestPointOnCrackFront will fail)
      nCompNodsWithBranchFnct = 1;
      int first_comp_nod_index = compNodWithBranchFnct[0];
      compNodWithBranchFnct.resize(1);
      compNodWithBranchFnct[0] = first_comp_nod_index;

      closePt = 0.;
      closePt_tangent = 0.;
      closePt_normal = 0.;
      closePt_binormal = 0.;
    }
    
    outf << "object \"dataForBranchFncTan\" class array type double rank 1 shape 3 items " 
         << nCompNodsWithBranchFnct << " data follows\n";
    
    compNodIndex = compNodWithBranchFnct.begin();
    compNodIndex_end = compNodWithBranchFnct.end();
    for( ; compNodIndex != compNodIndex_end; ++compNodIndex ){
      
      now_index = (*compNodIndex);
      p_compNod = p_cMesh->compNodFromArray( now_index ) ;
      for(Subscript icoor=0; icoor<3; ++icoor )
        point(icoor) =  p_compNod->geoNod()->coor(icoor);
      
      if ( !shp_crack_geoeng->has_tempCrackSurfaceFacet() ) {
        shp_crack_geoeng->closestPointOnCrackFront(point, 
                                                   closePt, closePt_tangent, 
                                                   closePt_normal, closePt_binormal,
                                                   averageCrackFrontCoordSystem);
      }
      
      for(Subscript icoor=0; icoor<3; ++icoor)
        outf << "\t" << closePt_tangent( icoor );
      outf << endl;
    }
    outf << "attribute \"dep\" string \"positions\"\n\n";
    
    outf << "object \"dataForBranchFncBiN\" class array type double rank 1 shape 3 items " 
         << nCompNodsWithBranchFnct << " data follows\n";

    compNodIndex = compNodWithBranchFnct.begin();
    compNodIndex_end = compNodWithBranchFnct.end();
    for( ; compNodIndex != compNodIndex_end; ++compNodIndex ){
      
      now_index = (*compNodIndex);
      p_compNod = p_cMesh->compNodFromArray( now_index ) ;
      for(Subscript icoor=0; icoor<3; ++icoor )
        point(icoor) =  p_compNod->geoNod()->coor(icoor);
      
      if ( !shp_crack_geoeng->has_tempCrackSurfaceFacet() ) {
        shp_crack_geoeng->closestPointOnCrackFront(point, closePt, closePt_tangent, 
                                                   closePt_normal, closePt_binormal,
                                                   averageCrackFrontCoordSystem);
      }
      
      for(Subscript icoor=0; icoor<3; ++icoor)
        outf << "\t" << closePt_binormal( icoor );
      outf << endl;
    }
    outf << "attribute \"dep\" string \"positions\"\n\n";
    
    outf << "object \"dataForBranchFncNorm\" class array type double rank 1 shape 3 items " 
         << nCompNodsWithBranchFnct << " data follows\n";
    
    compNodIndex = compNodWithBranchFnct.begin();
    compNodIndex_end = compNodWithBranchFnct.end();
    for( ; compNodIndex != compNodIndex_end; ++compNodIndex ){
      
      now_index = (*compNodIndex);
      p_compNod = p_cMesh->compNodFromArray( now_index ) ;
      for(Subscript icoor=0; icoor<3; ++icoor )
        point(icoor) =  p_compNod->geoNod()->coor(icoor);
      
      if ( !shp_crack_geoeng->has_tempCrackSurfaceFacet() ) {
        shp_crack_geoeng->closestPointOnCrackFront(point, closePt, closePt_tangent, 
                                                   closePt_normal, closePt_binormal,
                                                   averageCrackFrontCoordSystem);
      }
      
      for(Subscript icoor=0; icoor<3; ++icoor)
        outf << "\t" << closePt_normal( icoor );
      outf << endl;
    }
    outf << "attribute \"dep\" string \"positions\"\n\n";

    
    // Create DX fields with the objects defined above
    
    outf << "object \"stepFunction\" class field\n"
         << "component \"positions\" value \"compNodsWithStepFnct\" \n"
         << "component \"data\" value \"dataForStepFnct\"\n\n";
    
    outf << "object \"branchFunction\" class field\n"
         << "component \"positions\" value \"compNodsWithBranchFnct\" \n"
         << "component \"data\" value \"dataForBranchFnct\"\n\n";
    
    // added by DJK
    outf << "object \"LocalSolution\" class field\n"
         << "component \"positions\" value \"compNodsWithLocalSol\" \n"
         << "component \"data\" value \"dataForLocalSol\"\n\n";
    
    outf << "object \"interPoints\" class field\n"
         << "component \"positions\" value  \"intersectionPoints\" \n"
         << "component \"data\" value \"dataForIntersections\" \n\n";
    
    outf << "object \"BranchFncTan\" class field\n"
         << "component \"positions\" value \"compNodsWithBranchFnct\" \n"
         << "component \"data\" value \"dataForBranchFncTan\" \n\n";
    
    outf << "object \"BranchFncBiN\" class field\n"
         << "component \"positions\" value \"compNodsWithBranchFnct\" \n"
         << "component \"data\" value \"dataForBranchFncBiN\" \n\n";
    
    outf << "object \"BranchFncNorm\" class field\n"
         << "component \"positions\" value \"compNodsWithBranchFnct\" \n"
         << "component \"data\" value \"dataForBranchFncNorm\" \n\n";
    
    // collect all fields in a group
    
    outf << "object \"all\" class group\n"
         << "member 0 value \"stepFunction\"\n"
         << "member 1 value \"branchFunction\"\n"
         << "member 2 value \"interPoints\"\n"
         << "member 3 value \"BranchFncTan\"\n"
         << "member 4 value \"BranchFncBiN\"\n"
         << "member 5 value \"BranchFncNorm\"\n"
         << "member 6 value \"LocalSolution\"\n" // DJK
         << "end\n" << endl;
    
    // close file
    odataFile.close();
    
  }
  // **************************************************************
  // **************************************************************


  void test_intersection_w_crack_edges(const boost::shared_ptr<GeoMesh>& shp_gMesh,
                                       const boost::shared_ptr<CrackGeoEng3D2>& shp_crack_geoeng,
                                       std::vector<CrackCutElInfo2 *>& crack_elems_info) {


    std::cout << "\n**test_intersection_w_crack_edges**" << std::endl;

    Polyhedron* p_crack_surface = shp_crack_geoeng->crackSurface();

    CGAL_AABB_Tree<CGAL_Kernel_for_AABB> aabb_tree_for_geoel;

    std::list< CGAL_AABB_Tree_Object_and_primitive_id<CGAL_Kernel_for_AABB> > all_intersections;

    // This is our customized polyhedron
    Polyhedron tet_polyhedron;

    // Loop over ALL GeoEls of problem GeoMesh (NOT just those that have
    // CrackCutElInfo2 since there can be elements that are cut by the crack
    // surface along one of its faces but has no intersections along element
    // edges. CrackCutElInfo2 has only intersection with element edges on
    // entry.

    const GeoMesh::elemContainerType* p_geoElContainer = shp_gMesh->elemContn();
    GeoMesh::elemContainerType::const_iterator 
      epos = p_geoElContainer->begin(),
      epos_end = p_geoElContainer->end();

    CrackCutElInfo2* p_crack_cut_info = NULL;
    for(; epos != epos_end; ++epos) {
      GeoEl* p_geoEl_prob = *epos;

      // skip inactive elements
      if ( !p_geoEl_prob->active() )
        continue;

      DBGPRNT0(std::cout << "\nGeoEl ID = " <<  p_geoEl_prob->ID() << std::endl);
      DBGPRNT0(std::cout << "GeoEl: " << (*p_geoEl_prob) << std::endl);


      // Create a Polyhedron out of this GeoEl

      create_Polyhedron_from_GeoEl(p_geoEl_prob,
                                   tet_polyhedron);

 
      // Check intersections between edges of crack surface and this tet
      
      // AABB tree for GeoEl
      aabb_tree_for_geoel.clear();
      aabb_tree_for_geoel.rebuild(tet_polyhedron.facets_begin(), tet_polyhedron.facets_end(),
                                  tet_polyhedron );

      //
      // Loop over edges of crack surface 

      // CAD to JG: Is this the right way to loop over edges ?
      // It seems that this just loop over every other halfedge
      // iterator over all edges (every other halfedge). 
      Polyhedron::Edge_iterator 
        edge_pos = p_crack_surface->edges_begin(),
        edge_end = p_crack_surface->edges_end();

      p_crack_cut_info = NULL;
      int iedge = 0;
      int num_new_intersect_w_elem = 0;
      for( ; edge_pos != edge_end; ++edge_pos, ++iedge) {
        
        std::cout << "\ncrack surface edge # " << iedge << std::endl;

        // Create a query segment for this Edge of crack surface
        // CAD to JG: Is this correct?? 
        CGAL::Segment_3<CGAL_Kernel> segment_query( edge_pos->opposite()->vertex()->point(), 
                                    edge_pos->vertex()->point() );

        // source vertex
        DBGPRNT0(cout << "Source: " << segment_query.source() << endl);
        // target vertex 
        DBGPRNT0(cout << "Target: " << segment_query.target() << endl);
        
        // Get intersections

        // computes all intersections with segment query as pairs <object, primitive_id>
        // primitive id == a facet handle
        
        // Object used to process intersections and discard repeated ones
        //      EdgeIntersectData edge( segment_query );

        all_intersections.clear();

        aabb_tree_for_geoel.all_intersections( segment_query, 
                                               std::back_inserter(all_intersections) );

        DBGPRNT0(cout << "Number of intersection between this edge and GeoEl : " 
                 << all_intersections.size() << endl);

        // If intersections found, need CrackCutElInfo2 associated with
        // p_geoEl_prob to store them

        // Get CrackCutElInfo2 object for this element

        if ( all_intersections.size() == 0 )
          continue; // try next crack surface edge

        // NOTE: May need to create CrackCutElInfo2 for this element if it does
        // not have one Need to set scalling factor for it.
        //
        // Give warning: Crack surface intersect facet but does not intersect
        // ANY edge of this element
        
        // CrackCutElInfo2 pointers are stored by index of GeoEls from global
        // problem.
        p_crack_cut_info = crack_elems_info[ p_geoEl_prob->index() ];
        
        if ( !p_crack_cut_info ) {
          // first intersections for this element with crack surface: Need
          // CrackCutElInfo2
          
          // FIXME when porting to CM2!!
          CrackManager3D2* p_crack_mg = NULL;
          p_crack_cut_info = new CrackCutElInfo2( p_geoEl_prob, p_crack_mg );
          crack_elems_info[ p_geoEl_prob->index() ] = p_crack_cut_info;
          
          cout << "\nNOTE: GeoEl with ID " << p_geoEl_prob->ID() 
               << "\nHas face intersections but NO edge intersections " << endl;
        }
        // sanity check
        ERRCHK( p_geoEl_prob != p_crack_cut_info->geoElProb(),
                "Corrupted crack_elems_info vector!");
        

        std::list< CGAL_AABB_Tree_Object_and_primitive_id<CGAL_Kernel_for_AABB> >::iterator
          inters_pos = all_intersections.begin(),
          inters_end = all_intersections.end();
        
        int itersect = 0;
        for( ; inters_pos != inters_end; ++inters_pos, itersect++) {
          DBGPRNT0(cout << "\nProcessing intersection #: " << itersect << endl);

#if 0       
          // Object, which can be a point or a segment
          CGAL::Object object = inters_pos->first;
          // check if it is a point
          const CGAL::Point_3<CGAL_Kernel>* p_point=0;
          const CGAL::Segment_3<CGAL_Kernel>* p_segm=0;
          if( (p_point = CGAL::object_cast<CGAL::Point_3<CGAL_Kernel>>(&object)) ) {
            std::cout << "\nEncountered intersection object is a point" << std::endl;
            std::cout << (*p_point) << std::endl;
          }
          else if ( (p_segm = CGAL::object_cast<CGAL::Segment_3<CGAL_Kernel>>(&object)) ) {
            // Check if it is a segment otherwise error out
         
            std::cout << "\nEncountered intersection object is a segment" << std::endl;
            std::cout << (*p_segm) << std::endl;
          }
          else {
            std::cerr << "\nInvalid object returned by aaa_bb tree" << std::endl;
          }
#endif

          // primitive id (== a facet handle) of intersected facet
          Polyhedron::Facet_handle hd_facet = inters_pos->second;
          std::cout << "Index of tet face with intersection object: " 
                    << hd_facet->getFacetID() << std::endl;

          // Process intersections: 
          //  - Discard repeated ones; 
          //  - process end points of segments as two points
       
          num_new_intersect_w_elem += p_crack_cut_info->
            processIntersectWElemFacet( *inters_pos, edge_pos);

        } // next intersecion with this element


        // Print Edge: Number of intersections etc.
        //cout << "\nEdge after storing all intersections: " << iedge << endl;

      } // next crack surface edge

      cout << "\nGeoEl with ID " << p_geoEl_prob->ID() <<
        " num_new_intersect_w_elem = " << num_new_intersect_w_elem << endl;


      // DONE TBD: Check if intersections with facets are close to element
      // nodes. Discard if this is the case and mark nodes as at crack surface.
      // - Consistency check: There can be no NEW intersection with tet vertices when processing 
      // intersections with facets/volume ! 
      // Update status of GeoNods if this happens (at crack front, on temp surface, etc)


    } // next element in GeoMesh

  }
  // **************************************************************
  // **************************************************************

  void create_Polyhedron_from_GeoEl(const GeoEl* p_geoel,
                                    Polyhedron& geo_polyhedron) {

    // Created a Polyhedron out of the faces of a GeoEl
    // This is our customized Polyhedron which can have ids assigned to facets
    //
    // Implemented for tet4 and tet10 elements only

    
    if ( (typeid( *p_geoel ) != typeid( GeoElTet )) &&
         (typeid( *p_geoel ) != typeid( GeoElTet10 )) ) {
      ERREXIT("create_Polyhedron_from_GeoEl implemented for tet4 and tet10 only");
    }

    // cgal_to_iset_tet_facet_index[ifacet_cgal] = 
    // iset tet facet index corresponding to cgal tet facet index ifacet_cgal
    int cgal_to_iset_tet_facet_index[] = {3,1,2,0};

    // Create a Polyhedron out of this tet
    geo_polyhedron.erase_all();

    //   Halfedge_const_handle handle_hedege =
    geo_polyhedron.make_tetrahedron(iset::ConcRigArr1dToCGALPoint3( p_geoel->nodePtr(0)->coor() ),
                                    iset::ConcRigArr1dToCGALPoint3( p_geoel->nodePtr(1)->coor() ),
                                    iset::ConcRigArr1dToCGALPoint3( p_geoel->nodePtr(2)->coor() ),
                                    iset::ConcRigArr1dToCGALPoint3( p_geoel->nodePtr(3)->coor() ) );

    DBGPRNT0(std::cout << "geo_polyhedron: " << geo_polyhedron << std::endl);

    // Set IDs for facets so we can find out which facets are intersected by
    // crack surface edges

    Facet_iterator 
      facet_pos = geo_polyhedron.facets_begin(),
      facet_end = geo_polyhedron.facets_end();
    int ifacet_cgal = 0;
    for( ; facet_pos != facet_end; ++facet_pos, ++ifacet_cgal) {
      // Set index for facets according to iset convention for face index of a tet
      facet_pos->setFacetID( cgal_to_iset_tet_facet_index[ifacet_cgal] );

      DBGPRNT2(std::cout << "\nFacet : " << ifacet_cgal << endl;
               Polyhedron::Halfedge_around_facet_circulator circ_hedge = facet_pos->facet_begin();
               do {
                 std::cout << "\nFacet Vertex : " << circ_hedge->vertex()->point() << std::endl;
                 ++circ_hedge; // next hedge of this facet
               } while( circ_hedge != facet_pos->facet_begin() ) );
    }

    // Another option to build tet_polyhedron using our customized polyhedron
    // CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> polyhedronBuilder( tet_polyhedron, true /*verbose*/);

  }
  // **************************************************************
  // **************************************************************

  void test_intersection_w_crack_vertices(const boost::shared_ptr<GeoMesh>& shp_gMesh,
                                          const boost::shared_ptr<CrackGeoEng3D2>& shp_crack_geoeng,
                                          std::vector<CrackCutElInfo2 *>& crack_elems_info) {


    std::cout << "\n**test_intersection_w_crack_vertices**" << std::endl;

    // Find crack surface vertices that are in the interior of GeoEls

    // Loop over GeoEls of problem GeoMesh that HAVE CrackCutElInfo2 (there
    // can be NO element that has a crack vertex but no intersections with
    // crack edges, thus, a CrackCutElInfo2 must have already been created for
    // elements that have crack vertices).

    // ASSUME intersections between element edges and crack surface and
    // element faces and crack edges have already been computed and stored in
    // CrackCutElInfo2.

    Polyhedron* p_crack_surface = shp_crack_geoeng->crackSurface();

    const GeoMesh::elemContainerType* p_geoElContainer = shp_gMesh->elemContn();
    GeoMesh::elemContainerType::const_iterator 
      epos = p_geoElContainer->begin(),
      epos_end = p_geoElContainer->end();

    CrackCutElInfo2* p_crack_cut_info = NULL;
    for(; epos != epos_end; ++epos) {
      GeoEl* p_geoEl_prob = *epos;

      // skip inactive elements
      if ( !p_geoEl_prob->active() )
        continue;

      // Check if this element has a CrackCutElInfo2 object
      //
      // CrackCutElInfo2 pointers are stored by index of GeoEls from global
      // problem.
      p_crack_cut_info = crack_elems_info[ p_geoEl_prob->index() ];
      if ( !p_crack_cut_info )
        continue; // see why above


      DBGPRNT0(std::cout << "\nGeoEl ID = " <<  p_geoEl_prob->ID() << std::endl);
      DBGPRNT0(std::cout << "GeoEl: " << (*p_geoEl_prob) << std::endl);

      // This element has intersections with crack surface: May have vertices of
      // surface inside it
      //
      // Create a CGAL::Tetrahedron_3 out of this GeoEl

      CGAL::Tetrahedron_3<Kernel> 
        cgal_tet( iset::ConcRigArr1dToCGALPoint3( p_geoEl_prob->nodePtr(0)->coor() ),
                  iset::ConcRigArr1dToCGALPoint3( p_geoEl_prob->nodePtr(1)->coor() ),
                  iset::ConcRigArr1dToCGALPoint3( p_geoEl_prob->nodePtr(2)->coor() ),
                  iset::ConcRigArr1dToCGALPoint3( p_geoEl_prob->nodePtr(3)->coor() ) );

      // Check which crack vertices are inside cgal_tet

      // Loop over vertices of crack surface 

      Polyhedron::Vertex_iterator 
        vertex_pos = p_crack_surface->vertices_begin(),
        vertex_end = p_crack_surface->vertices_end();

      int ivertex = 0;
      int num_new_vertex_inside_elem = 0;
      for( ; vertex_pos != vertex_end; ++vertex_pos, ++ivertex) {
        
        // NOTE: Could use flag to skip vertices that are not at the crack front.

        CGAL::Point_3<CGAL_Kernel> vertex_point = vertex_pos->point();

        std::cout << "\ncrack surface vertex # " << ivertex << std::endl;
        std::cout << "\n Point : " << vertex_point << std::endl;
        
        if ( !cgal_tet.has_on_bounded_side( vertex_point ) )
          continue;

        std::cout << "\nThis point IS on the bounded side of element" << std::endl;

        // If intersection is found, need CrackCutElInfo2 associated with
        // p_geoEl_prob to store it
        ERRCHK(!p_crack_cut_info, 
               "Element has a crack vertex but NO CrackCutElInfo2");
        
        // sanity check
        ERRCHK( p_geoEl_prob != p_crack_cut_info->geoElProb(),
                "Corrupted crack_elems_info vector!");
        
        
        // Process intersection: 
        //  - Discard if repeated, etc;
#if 1   
        if ( p_crack_cut_info->
             processIntersectWCrackVertex( vertex_pos, vertex_point) )
          ++num_new_vertex_inside_elem;
#endif
        
      } // next crack surface vertex

      cout << "\nGeoEl with ID " << p_geoEl_prob->ID() <<
        " num_new_vertex_inside_elem = " << num_new_vertex_inside_elem << endl;
      
      
      // DONE: Check if intersections with vertex are close to element
      // nodes. Discard if this is the case and mark nodes as at crack surface.
      // - Consistency check: There can be no NEW intersection with tet vertices when processing 
      // intersections with facets/volume ! 
      // Update status of GeoNods if this happens (at crack front, on temp surface, etc)
      
      
    } // next element in GeoMesh
    
  }
#endif
  // **************************************************************
  // **************************************************************

#if 0
  void create_Tetrahedron_from_GeoEl(const GeoEl* p_geoel,
                                     Polyhedron& geo_polyhedron) {

    CGAL::Tetrahedron_3<Kernel> 
      t( Point_3<Kernel> p0, Point_3<Kernel> p1, Point_3<Kernel> p2, Point_3<Kernel> p3); 

  }

#endif
  // **************************************************************
  // **************************************************************

#if 0
  struct Plane_equation {
    template <class Facet>
    typename Facet::Plane_3 operator()( Facet& f) {
      typename Facet::Halfedge_handle h = f.halfedge();
      typedef typename Facet::Plane_3  Plane;
      return Plane( h->vertex()->point(),
                    h->next()->vertex()->point(),
                    h->next()->next()->vertex()->point());
    }
  };
  
  // **************************************************************
  // **************************************************************


  bool hasPointOnBoundedSide(Polyhedron& geo_polyhedron) {

    /*

    // check our Polydedron supports plane equations for facets

    //    Polyhedron::Facet::plane() == CGAL::Tag_true)


    cout << "\n**Polyhedron::Supports_facet_plane :" 
         << Polyhedron::Supports_facet_plane.value << endl

         << (Polyhedron::Supports_facet_plane == CGAL::Tag_true) << endl;
    */

    // Sample point
    CGAL::Point_3<CGAL_Kernel> p(1.0, 0.0, 0.0);
    CGAL::Point_3<CGAL_Kernel> q(0.0, 1.0, 0.0);
    CGAL::Point_3<CGAL_Kernel> r(0.0, 0.0, 1.0);

    Polyhedron::Facet::Plane_3 plane_0(p, q, r);

    CGAL::Polyhedron_3<Kernel>::Facet::Plane_3 plane_1(p, q, r);


    /*
    // Compute a plane for each facet of geo_polyhedron and assign to facet.plane

    std::transform( geo_polyhedron.facets_begin(), geo_polyhedron.facets_end(), 
                    geo_polyhedron.planes_begin(),
                    Plane_equation() );

    // Sample point
    CGAL::Point_3<CGAL_Kernel> p(1.0, 0.0, 0.0);

    */

    /*    
    // loop over crack facets
    Polyhedron::Facet_const_iterator
      facet_pos = geo_polyhedron.facets_begin(),
      facet_end = geo_polyhedron.facets_end();

    bool point_is_on_bounded_side = true;
    for( ; facet_pos != facet_end; ++facet_pos) {

      if( facet_pos->plane().has_on_negative_side( p ) ) {
        point_is_on_bounded_side = false;
        break;
      }
    }
    */

    /*
    //  iterator over all plane equations.
    Polyhedron::Plane_iterator  
      plane_pos = geo_polyhedron.planes_begin(),
      plane_end = geo_polyhedron.planes_end();

    */

    /*
    Polyhedron::Facet::Plane_3 plane_0 = *( geo_polyhedron.planes_begin() );
      
    cout << "\n** Plane.a :" << plane_0.a() << endl;
    */

    return false;
  }
#endif

  // **************************************************************
  // **************************************************************

#if 0 // not being used
  void test_Delaunay_Triangulation(const std::vector<GNodCrackedElemsInfo *>& gnod_cracked_elems,
                                   const std::vector<CrackCutElInfo2 *>& crack_elems_info) {

    // Create a Delaunay Triangulation for each CrackCutElInfo2 in crack_elems_info

    DBGPRNT0( cout << "\n****** computeDelaunay3D ******\n"; );

    CGAL_Delaunay_Triangulation_WInfo delaunayTetElemMesh;

    std::list< std::pair<CGAL::Point_3<CGAL_Kernel_for_Del>, DelVertexInfo> > inputVertexWInfoList;

    std::vector<CrackCutElInfo2 *>::const_iterator 
      cut_info2_pos = crack_elems_info.begin(),
      cut_info2_end = crack_elems_info.end();

    CrackCutElInfo2* p_cut_el_info2 = NULL;
    for( ; cut_info2_pos != cut_info2_end; ++cut_info2_pos ) {

      p_cut_el_info2 = *cut_info2_pos;

      if ( !p_cut_el_info2 )
        continue; // no CrackCutElInfo2 for this element

      inputVertexWInfoList.clear();

      // Insert a std::pair<CGAL::Point_3<CGAL_Kernel>, DelVertexInfo> for each GeoNod and for
      // each intersection between GeoEl and crack surface

      // Insert vertices with GNodCrackedElemsInfo, one for each GeoNod of GeoEl
      // These are always inserted, even if intersections are on a plane (see below)
      
      const GeoEl* p_geoel = p_cut_el_info2->geoElProb();

      DBGPRNT0( cout << "\n***** Insert tetrahedron nodes *****\n" << endl
                << *p_geoel << endl; );
      DBGPRNT0( cout << "\nList of primitive nodes:" << endl; );
      for( int inod=0; inod<p_geoel->numNodes(); ++inod ) {

        const GeoNod* p_geoNod = p_geoel->nodePtr( inod );
        const GNodCrackedElemsInfo* p_gnode_info = gnod_cracked_elems[ p_geoNod->index() ];

        DBGPRNT0( cout << "\n node " << inod << " : " << std::setprecision( 15 )
              << " coord. = " <<  p_geoNod->coor(); );

        // store nodal coordinates for generating Delaunay3D tetratization
        if (p_gnode_info) {
          inputVertexWInfoList.push_back( std::make_pair( iset::ConcRigArr1dToCGALPoint3<CGAL_Kernel_for_Del>( p_geoNod->coor() ),  
                                                          DelVertexInfo(p_gnode_info) ) );
        }
        else { // FIXME: This should be an ERROR!
          inputVertexWInfoList.push_back( std::make_pair( iset::ConcRigArr1dToCGALPoint3<CGAL_Kernel_for_Del>( p_geoNod->coor() ),  
                                                          DelVertexInfo() ) );
        }
      }
      DBGPRNT0( cout << "\n\n***** Insert tetrahedron nodes done *****\n"<< endl;);

      // loop over intersection points to store them in the list

      // TBD (ONLY if element has no intersection at the crack front !): 
      //      Check if intersection points (including those at GeoNods) are on
      //      a plane. If that is the case, use below only (at most 4) intersections at element
      //      edges.
      //      
      //      NOTE: We use intersection at geonods to check planarity of
      //            intersections (ONLY). 
      //      Write a function that returns a list of edge intersections (only) if
      //      above planarity condition is satisfied.
      //
      // Check if (selected) edge + node intersection define a
      // plane. "Selected" edge means that we pick at most one intersection from
      // each edge. 
      // If we find only 1 or 2 points --> crack is at bnd or single face/interior of element.
      // Handle these cases appropriately.  

      // If yes, check if remaining edge intersections, pure face and pure
      // interior intersections belong to this plane.  "Pure" means that it is
      // not an edge intersection too.  For simplicity, we can just check all
      // intersections.
      //
      // DONE: Look at JP's calculation of plane for possible useful CGAL function.


      const std::set<IntersectPointData, CompareIntersectPoints >& 
        intersect_points = *(p_cut_el_info2->intersectPoints());

      std::set<IntersectPointData, CompareIntersectPoints >::const_iterator
        intersect_pos = intersect_points.begin(),
        intersect_end = intersect_points.end();
    
      DBGPRNT0( cout << "\n***** Insert intersection points *****\n" << endl
            << "number of intersection points = "  
                << intersect_points.size() << endl
                << "\n IntersectPointData coordinates: " << endl; );

      CGAL_converter_Kernel_to_Del_Kernel converter_K_to_Del_K;
      int i_intersect_point = 0;
      for( ; intersect_pos != intersect_end; ++intersect_pos, ++i_intersect_point ){
        
        inputVertexWInfoList.push_back( std::make_pair(  converter_K_to_Del_K( intersect_pos->intersectPoint() ),  
                                                         DelVertexInfo( &(*intersect_pos) ) ) );

        DBGPRNT0( cout << "\n intersect point index " << i_intersect_point << endl 
                  << std::setprecision( 15 )
                  << " coord. = " << intersect_pos->intersectPoint() );
      }
      DBGPRNT0(cout << "\n\n***** Insert intersection points done *****\n" << endl;);



      //
      // Generate CGAL Delaunay tetrahedralization
      //
      DBGPRNT0( cout << "\n**** Generating tetrahedralization: **** \n" << endl; );

      // CAD WARNING: The insert function is not guaranteed to insert the points
      // following the order in the list, as spatial_sort is used to improve
      // efficiency. Thus, the index of primitive nodes is NOT [1:4], etc.

      delaunayTetElemMesh.clear();
      size_t num_insert_points = 
        delaunayTetElemMesh.insert( inputVertexWInfoList.begin(), inputVertexWInfoList.end() );

      DBGPRNT0( cout<<"\n Number of inserted points in Delaunay = " 
                << num_insert_points << endl);

      ERRCHK( num_insert_points != (intersect_points.size() + p_geoel->numNodes()),
              " Delaunay rejected some vertices: Are there repeated vertices?");

      DBGPRNT0( cout<<"\n **** Tetrahedralization DONE **** \n" << endl;);


      // Print tetratization info

#define CHECK_DELAUNAY_TRIANGULATION 1
#if CHECK_DELAUNAY_TRIANGULATION
      bool val_del_triang = delaunayTetElemMesh.is_valid( true /*verbose */);
      DBGPRNT0(cout <<"\n Is Delaunay Triangulation valid (0/1)? :  " << std::flush << 
               val_del_triang << endl);
      ERRCHK( !val_del_triang, "computeDelaunay3D: Invalid Del. Triang.!!");
#endif

      DBGPRNT0( cout <<"\n *** CGAL 3D Delaunay tetrahedralization INFO *** \n" 
                << "\n number of finite vertices = "<< delaunayTetElemMesh.number_of_vertices()
                << "\n number of vertices        = "<< delaunayTetElemMesh.tds().number_of_vertices()
                << "\n number of finite edges    = "<< delaunayTetElemMesh.number_of_finite_edges()
                << "\n number of finite facets   = "<< delaunayTetElemMesh.number_of_finite_facets()
                << "\n number of finite cells    = "<< delaunayTetElemMesh.number_of_finite_cells()
                << "\n number of cells           = "<< delaunayTetElemMesh.number_of_cells()
                << endl;);

      // print Delaunay vertex index and coordinates

#define PRINT_DELAUNAY_VERTEX_AND_CELLS 1
#if PRINT_DELAUNAY_VERTEX_AND_CELLS

      cout << "\n\n *** Delaunay finite vertex coordinates after insertion*** \n";
  
      CGAL::Inverse_index< CGAL_Delaunay_Triangulation_WInfo::Finite_vertices_iterator >  
        del_fin_vert_index( delaunayTetElemMesh.finite_vertices_begin(), 
                            delaunayTetElemMesh.finite_vertices_end() );

      CGAL_Delaunay_Triangulation_WInfo::Finite_vertices_iterator 
        del_vert_pos = delaunayTetElemMesh.finite_vertices_begin(),
        del_vert_end = delaunayTetElemMesh.finite_vertices_end();

      for ( ; del_vert_pos != del_vert_end; ++del_vert_pos) {

        cout << "\n index: " << del_fin_vert_index[ del_vert_pos ] << "\n\tcoord:";
        for( Subscript iCoord = 0;  iCoord < 3;  ++iCoord ) {
          cout <<" \t" << std::setw(12) << std::setprecision( 15 )
               << del_vert_pos->point().cartesian( iCoord );
        }

        const DelVertexInfo& vertex_info = del_vert_pos->info();
        if ( vertex_info.p_gnodCrackedElm ) {
          cout << "\n\tVertex is at GeoNod with ID = " 
               << vertex_info.p_gnodCrackedElm->geoNod()->ID() << endl;
        }
        else if ( vertex_info.p_intersectData ) {
          cout << "\n\tVertex is at an intersection point " << endl;
        }
        // FIXME: one of obove pointers must be not NULL
        
      }
#endif

      // DONE: Put this code in another function that:
      //   compute volume of cells, check if theur sum equal volume primitive element,
      //   create SubElemInfo for each cell

      // compute volume of the primitive element
      //





    } // next CrackCutElInfo2



  }
#endif
  // **************************************************************
  // **************************************************************

} /* end of unamed namespace */
