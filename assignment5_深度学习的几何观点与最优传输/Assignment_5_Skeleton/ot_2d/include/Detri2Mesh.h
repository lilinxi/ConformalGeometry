
/*!  \file   Detri2OMTMesh
 *   \brief  converting Detri2::triangulation to OMTMesh
 *   \author Hang Si and David Gu
 *   \date   documented on 07/15/2020
 *
 *   Base Class for Detri2 and COMTMesh conversion
 */

#ifndef _DETRI2_OMTMESH_
#define _DETRI2_OMTMESH_

#include <vector>

#include "omt_detri2.h"

namespace MeshLib
{

/*-------------------------------------------------------------------------------------------------------------------------------------

        Detri2 Mesh Interface

--------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CDetri2Mesh Class
 *
 *  class for Detri2 Mesh Conversion
 *
 */
class CDetri2Mesh
{
  public:
    CDetri2Mesh(){};
    ~CDetri2Mesh(){};

  protected:
    /*! compute Weighted Delaunay and Power Voronoi
     */
    bool __detri2_WDT(COMTMesh* mesh, detri2::Triangulation** outputTr);
    
    /*! generate background triangulation
     */
    void __detri2_generate_disk(detri2::Triangulation*& domainTr, double& total_target_area);
    
    /*! convert Detri2 triangulation to Mesh
     */
    void __detri2_to_mesh(detri2::Triangulation* outputTr, detri2::Triangulation* domainTr, COMTMesh*& pMesh);
};

/*
        generate power delaunay triangulation, if there are missing points return false; if there is no missing point, return true
*/
inline bool CDetri2Mesh::__detri2_WDT(COMTMesh* pMesh,                 // input mesh, the vertex weight is set
                                      detri2::Triangulation** outputTr // output triangulation,
)
{
    std::vector<MeshLib::CPoint> ptlist;
    std::vector<double> input_weights;
    std::vector<int> missing_point_list;
    std::vector<std::pair<int, int>> boundary_edge_list;

    for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); viter++)
    {
        COMTMesh::CVertex* pv = *viter;
        // CPoint p = pv->point();
        CPoint p = CPoint(pv->uv()[0], pv->uv()[1], 0);
        ptlist.push_back(p);
        double weight = pv->weight();
        input_weights.push_back(weight);
    }

    generate_wdt(ptlist, input_weights, boundary_edge_list, false, outputTr, missing_point_list);

    if (!missing_point_list.empty())
    {
        std::cout << "Missing Point" << std::endl;
        return false;
    }

    return true;
};

/*
        generate the background weighted Delaunay triangulation,
        generate a disk with raidus one
*/
inline void CDetri2Mesh::__detri2_generate_disk(detri2::Triangulation*& domainTr, double& total_target_area)
{
    std::vector<MeshLib::CPoint> ptlist;
    std::vector<double> input_weights;
    std::vector<int> missing_point_list;
    std::vector<std::pair<int, int>> boundary_edge_list;
    int n = 32;
    double R = 1.0;

    for (int i = 0; i < n; i++)
    {
        double ang = i * 2.0 * PI / n;
        CPoint p(R * cos(ang), R * sin(ang), 0);
        ptlist.push_back(p);
    }

    CPolygon3D plygon;
    for (size_t i = 0; i < ptlist.size(); i++)
    {
        CSegment3D s(ptlist[i], ptlist[(i + 1) % ptlist.size()]);
        plygon.add(s);
    }
    total_target_area = plygon.area();

    for (int i = 0; i < n - 1; i++)
    {
        boundary_edge_list.push_back(std::pair<int, int>(i + 1, i + 2));
    }
    boundary_edge_list.push_back(std::pair<int, int>(n, 1));

    for (int i = 0; i < n; i++)
    {
        input_weights.push_back(0);
    }

    generate_wdt(ptlist, input_weights, boundary_edge_list, true, &domainTr, missing_point_list);

    // for debugging
    /*
    COMTMesh* pM = new COMTMesh;
    export_Detri2_to_OMTmesh(domainTr, pM);
    for (COMTMesh::MeshVertexIterator viter(pM); !viter.end(); viter++)
    {
        COMTMesh::CVertex* pv = *viter;
        pv->point() = CPoint(pv->uv()[0], pv->uv()[1], 0);
    }
    pM->labelBoundary();
    pM->write_m("Circle.m");
    */
};

/*
        convert weighted Delaunay triangulation, power diagram to a mesh
        with edge length, dual length
        vertex dual cell, dual center, and dual celll area
*/
inline void CDetri2Mesh::__detri2_to_mesh(detri2::Triangulation* outputTr, detri2::Triangulation* domainTr, COMTMesh*& pWDT)
{
    // if the output mesh pointer is non-empty, delete the old mesh
    if (pWDT != NULL)
        delete pWDT;

    // generate a new mesh
    pWDT = new COMTMesh;
    // convert WDT to the mesh
    export_Detri2_to_OMTmesh(outputTr, pWDT);
    pWDT->labelBoundary();
    // pWDT->write_m("Debug.m");

    get_voronoi_vertices(outputTr, domainTr);
    // get vertex dual cell, dual centerand dual area
    for (COMTMesh::MeshVertexIterator viter(pWDT); !viter.end(); viter++)
    {
        COMTMesh::CVertex* pv = *viter;
        double area;
        CPoint center;
        std::vector<MeshLib::CPoint> ptlist;
        get_voronoi_cell(outputTr, pv->id(), ptlist, area, center);
        pv->dual_area() = area;
        pv->dual_center() = center;

        for (size_t i = 0; i < ptlist.size(); i++)
        {
            CPoint start = ptlist[i];
            CPoint end = ptlist[(i + 1) % ptlist.size()];
            CSegment3D s(start, end);
            pv->dual_cell().add(s);
        }
        // break;
        // std::cout << "Voronoi cell " << pv->index() << " " << ptlist.size() << " " << area << std::endl;
    }
    // get edge length
    for (COMTMesh::MeshEdgeIterator eiter(pWDT); !eiter.end(); eiter++)
    {
        COMTMesh::CEdge* pe = *eiter;
        COMTMesh::CVertex* pv1 = pWDT->edgeVertex1(pe);
        COMTMesh::CVertex* pv2 = pWDT->edgeVertex2(pe);
        // pe->length() = pWDT->edgeLength(pe);
        pe->length() = (pv1->uv() - pv2->uv()).norm();
    }
    // get dual edge length
    for (COMTMesh::MeshEdgeIterator eiter(pWDT); !eiter.end(); eiter++)
    {
        COMTMesh::CEdge* pe = *eiter;
        COMTMesh::CVertex* pv1 = pWDT->edgeVertex1(pe);
        COMTMesh::CVertex* pv2 = pWDT->edgeVertex2(pe);
        double length;
        get_voronoi_edge(outputTr, pv1->id(), pv2->id(), length);
        pe->dual_length() = length;
    }
};
} // namespace MeshLib
#endif //! DETRI2_MESH_H_
