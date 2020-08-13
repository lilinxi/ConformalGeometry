#ifndef _SPHERICAL_HARMONIC_MAP_MESH_H_
#define _SPHERICAL_HARMONIC_MAP_MESH_H_

#include "Mesh/BaseMesh.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Vertex.h"

#include "Mesh/Boundary.h"
#include "Mesh/Iterators.h"
#include "Parser/parser.h"

namespace MeshLib {
    class CSHMVertex;

    class CSHMEdge;

    class CSHMFace;

    class CSHMEdge;

/*! \brief CSHMVertex class
 *
 *   Vertex class for spherical harmonic map
 */
    class CSHMVertex : public CVertex {
    public:
        /*! Constructor */
        CSHMVertex() : m_u(0, 0, 0) {};

        /*!	Vertex spherical harmonic map image coordinates
         */
        CPoint &u() { return m_u; };

        /*! vertex area
        */
        double &area() { return m_area; };

    protected:

        /*! Vertex spherical harmonic map image coordinates */
        CPoint m_u;
        /*! vertex area */
        double m_area;
    };

/*! \brief CSHMEdge class
 *
 *   Edge class for spherical harmonic map
 */
    class CSHMEdge : public CEdge {
    public:
        /*! Constructor */
        CSHMEdge() : m_length(0), m_weight(0) {};

        /*!	Edge weight */
        double &weight() { return m_weight; };

        /*! Edge length */
        double &length() { return m_length; };

    protected:
        /*!	Edge weight */
        double m_weight;

        /*! Edge length */
        double m_length;
    };

/*! \brief CSHMFace class
 *
 *   Face class for spherical harmonic map
 */
    class CSHMFace : public CFace {
    public:
        /*! Constructor */
        CSHMFace() {};

        /*! face normal */
        CPoint &normal() { return m_normal; };

        /*! face area */
        double &area() { return m_area; };

    protected:
        /*! face normal */
        CPoint m_normal;
        /*! face area */
        double m_area;
    };

/*! \brief CSHMHalfEdge class
 *
 *   HalfEdge class for spherical harmonic map
 */
    class CSHMHalfEdge : public CHalfEdge {
    public:
        /*!	CSHMHalfEdge constructor */
        CSHMHalfEdge() : m_angle(0) {};

        /*!	Corner angle */
        double &angle() { return m_angle; };

    protected:
        /*! Corner angle */
        double m_angle;
    };

/*! \brief CSHMMesh class
 *
 *	Mesh class for spherical harmonic map
 *
 */
    template<typename V, typename E, typename F, typename H>
    class TSHMMesh : public CBaseMesh<V, E, F, H> {
    public:
        typedef V CVertex;
        typedef E CEdge;
        typedef F CFace;
        typedef H CHalfEdge;

        typedef CBoundary<V, E, F, H> CBoundary;
        typedef CLoop<V, E, F, H> CLoop;

        typedef MeshVertexIterator<V, E, F, H> MeshVertexIterator;
        typedef MeshEdgeIterator<V, E, F, H> MeshEdgeIterator;
        typedef MeshFaceIterator<V, E, F, H> MeshFaceIterator;
        typedef MeshHalfEdgeIterator<V, E, F, H> MeshHalfEdgeIterator;

        typedef VertexVertexIterator<V, E, F, H> VertexVertexIterator;
        typedef VertexEdgeIterator<V, E, F, H> VertexEdgeIterator;
        typedef VertexFaceIterator<V, E, F, H> VertexFaceIterator;
        typedef VertexInHalfedgeIterator<V, E, F, H> VertexInHalfedgeIterator;
        typedef VertexOutHalfedgeIterator<V, E, F, H> VertexOutHalfedgeIterator;

        typedef FaceVertexIterator<V, E, F, H> FaceVertexIterator;
        typedef FaceEdgeIterator<V, E, F, H> FaceEdgeIterator;
        typedef FaceHalfedgeIterator<V, E, F, H> FaceHalfedgeIterator;
    };

    typedef TSHMMesh<CSHMVertex, CSHMEdge, CSHMFace, CSHMHalfEdge> CSHMMesh;
}

#endif // !_SPHERICAL_HARMONIC_MAP_MESH_H_
