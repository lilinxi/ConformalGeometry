/*!  \file OMTDynamicMesh.h
 *   \brief Optimal Transport DynamicMesh
 *   \author David Gu
 *   \date   documented on 05/01/2020
 *
 *   Dynamic Mesh for viewer
 */
#ifndef _OMT_DYNAMIC_MESH_H_
#define _OMT_DYNAMIC_MESH_H_

#include "Mesh/Boundary.h"
#include "Mesh/DynamicMesh.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Iterators.h"
#include "Mesh/Vertex.h"
#include "Parser/parser.h"
#include "Parser/traits_io.h"
#include "Geometry/Polygon3D.h"

namespace MeshLib
{

/*! \brief CMOTVertex class
 *
 *   Vertex class for optimal transport
 *   Trait: vertex rgb color, target area, dual area, index,
 *          dual cell, dual center, update direction, weight
 */
class COMTVertex : public CVertex
{
  protected:
    /*! vertex rgb color */
    CPoint m_rgb;

    /*! vertex target area */
    double m_target_area;

    /*! verte dual area */
    double m_dual_area;

    /*! vertex index */
    int m_index;

    /*! dual cell in 3D*/
    CPolygon3D m_dual_cell;

    /*! duall cell center */
    CPoint m_dual_center;

    /*! update direction */
    double m_update_direction;

    /*! weight */
    double m_weight;

  public:
    /*! CViewerVertex Constructor */
    COMTVertex()
    {
        m_rgb = CPoint(1, 1, 1); // default color is white
        m_index = 0;
        m_weight = 0;
        m_normal = CPoint(0, 0, 1);
        m_target_area = 0;
        m_dual_area = 0;
        m_update_direction = 0;
    }

    /*! vertex rgb color */
    CPoint& rgb() { return m_rgb; };

    /*! vertex target area */
    double& area() { return m_target_area; };

    /*! vertex target area */
    double& target_area() { return m_target_area; };

    /*! vertex dual area */
    double& dual_area() { return m_dual_area; };

    /*! read vertex rgb, uv from vertex string */
    void _from_string();

    /*! vertex index*/
    int& index() { return m_index; };

    /*! vertex conjudate cell */
    CPolygon3D& dual_cell() { return m_dual_cell; };

    /*! dual cell center */
    CPoint& dual_center() { return m_dual_center; };

    /*! update direction */
    double& update_direction() { return m_update_direction; };

    /*! weight */
    double& weight() { return m_weight; };
};

// read vertex rgb, uv from vertex string
inline void COMTVertex::_from_string()
{
    CParser parser(m_string);

    for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
    {
        CToken* token = *iter;
        if (token->m_key == "uv")
        {
            token->m_value >> m_uv;
        }
        else if (token->m_key == "rgb")
        {
            token->m_value >> m_rgb;
        }
        else if (token->m_key == "normal")
        {
            token->m_value >> m_normal;
        }
    }
};

/*! \brief COMTEdge class
 *
 *   Edge class for optimal transport
 *   Trait : Edge length, dual length
 */
class COMTEdge : public CEdge
{
  protected:
    /*! edge length */
    double m_length;
    
    /*! dual length*/
    double m_dual_length;

  public:
    /*! COTVertex Constructor */
    COMTEdge() { }

    /*! edge length */
    double& length() { return m_length; };

    /*! dual length */
    double& dual_length() { return m_dual_length; };
};

/*! \brief CMOTFace class
 *
 *  Face class for optimal transport
 *  Trait: face normal, area, dual point
 */
class COMTFace : public CFace
{
  protected:
    /*! face normal */
    CPoint m_normal;

    /*! face area  */
    double m_area;

    /*! face dual point */
    CPoint m_dual_point;

  public:
    /*! face normal */
    CPoint& normal() { return m_normal; };

    /*! face area */
    double& area() { return m_area; };

    /*! face dual point*/
    CPoint& dual_point() { return m_dual_point; };
};

/*! \brief COMTHalfEdge class
 *
 *   HalfEdge class for optimal transport
 *   Trait : corner uv, corner id, corner angle
 */
class COMTHalfEdge : public CHalfEdge
{
  public:
    /*! corner uv */
    CPoint2& uv() { return m_uv; };

    /*! corner id*/
    int& cid() { return m_cid; };

    /*! corner angle */
    double& angle() { return m_angle; };

  protected:
    /*! edge sharp */
    CPoint2 m_uv;

    /*! halfedge corner id*/
    int m_cid;

    /*! corner angle */
    double m_angle;
};

/*-------------------------------------------------------------------------------------------------------------------------------------

        Optimal Transport Mesh

--------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief COTMesh class
 *
 *	mesh class for optimal transport
 *
 */
template <typename V, typename E, typename F, typename H>
class COMTDynamicMesh : public CDynamicMesh<V, E, F, H>
{
  public:
    typedef V CVertex;
    typedef E CEdge;
    typedef F CFace;
    typedef H CHalfEdge;

    typedef CBoundary<V, E, F, H>                   CBoundary;
    typedef CLoop<V, E, F, H>                       CLoop;

    typedef MeshVertexIterator<V, E, F, H>          MeshVertexIterator;
    typedef MeshFaceIterator<V, E, F, H>            MeshFaceIterator;
    typedef MeshEdgeIterator<V, E, F, H>            MeshEdgeIterator;
    typedef VertexVertexIterator<V, E, F, H>        VertexVertexIterator;
    typedef FaceVertexIterator<V, E, F, H>          FaceVertexIterator;
    typedef VertexEdgeIterator<V, E, F, H>          VertexEdgeIterator;
    typedef VertexFaceIterator<V, E, F, H>          VertexFaceIterator;
    typedef FaceHalfedgeIterator<V, E, F, H>        FaceHalfedgeIterator;
    typedef VertexInHalfedgeIterator<V, E, F, H>    VertexInHalfedgeIterator;
    typedef VertexOutHalfedgeIterator<V, E, F, H>   VertexOutHalfedgeIterator;
};

typedef COMTDynamicMesh<COMTVertex, COMTEdge, COMTFace, COMTHalfEdge> COMTMesh;

} // namespace MeshLib
#endif //! OMT_DYNAMIC_MESH_H_
