/*!
*      \file omt_detri2.h
*      \brief Interface between OMT adn Detri2
*
*		Exchange (input and output) data between OMTMesh and Detri2
*
*	   \author Hang Si and David Gu
*      \date 05/07/2020
*
*/

#ifndef _OMT_DERTRI2_H_
#define _OMT_DERTRI2_H_

#include <vector>

#include "Geometry/Point.h"
#include "OTMesh.h"

#include "detri2.h"

#ifdef _WIN32
#define EXPORTIT __declspec(dllexport)
#else
#define EXPORTIT
#endif

// Generate initial weighted DT.
EXPORTIT
bool generate_wdt(/* input data */
                  std::vector<MeshLib::CPoint> &ptlist,
                  std::vector<double> &input_weights,
                  std::vector<std::pair<int, int>> &boundary_edge_list,
                  bool refine_mesh_flag,
                  /* output data */
                  detri2::Triangulation **outputTr,
                  std::vector<int> &missing_point_list);

// Update a given weighted DT whose vertex weights are changed.
EXPORTIT
bool remesh_wdt(/* input data */
                detri2::Triangulation *Tr,  /* Input/Output */
                std::vector<double> *new_vertex_weights,
                bool insert_missing_point,  /* option: true or false */
                /* output data */
                std::vector<int> *missing_point_list,
                std::vector<double> *updated_vertex_weights);

EXPORTIT
bool get_voronoi_vertices(detri2::Triangulation *Tr,
                          detri2::Triangulation *OMT_domain);

EXPORTIT
bool get_voronoi_cell(detri2::Triangulation *Tr,
                      int vertex_index,
                      std::vector<MeshLib::CPoint> &ptlist,
                      double &area,
                      MeshLib::CPoint &masspt);

EXPORTIT
bool get_voronoi_edge(detri2::Triangulation *Tr,
                      int vertex_index1,
                      int vertex_index2,
                      double &length);

// For visuliazation
EXPORTIT
bool export_Detri2_to_OMTmesh(detri2::Triangulation *Tr,   /* Input */
                              MeshLib::COMTMesh *OMTmesh); /* Output */

#endif //_OMT_DERTRI2_H_ defined
