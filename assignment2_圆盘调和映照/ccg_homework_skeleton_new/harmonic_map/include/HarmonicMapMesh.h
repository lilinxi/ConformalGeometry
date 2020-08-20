#ifndef _HARMONIC_MAP_MESH_H_
#define _HARMONIC_MAP_MESH_H_

#include "Mesh/BaseMesh.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"
#include "Mesh/Vertex.h"

#include "Mesh/Boundary.h"
#include "Mesh/Iterators.h"
#include "Parser/parser.h"

namespace MeshLib
{
class CHarmonicMapVertex;
class CHarmonicMapEdge;
class CHarmonicMapFace;
class CHarmonicMapHalfEdge;

/*! \brief CHarmonicMapVertex class
 *
 *   Vertex class for cut graph algoritm
 *   Trait : Vertex valence
 */
class CHarmonicMapVertex : public CVertex
{
  public:
    /*! Constructor */
    CHarmonicMapVertex() : 
        m_index(0), 
        m_rgb(229.0 / 255.0, 162.0 / 255.0, 141.0 / 255.0){};
    
    /*! Vertex index */
    int& idx() { return m_index; };

    /*! Vertex color */
    CPoint& rgb() { return m_rgb; };

    /*!
     *	Read vertex traits to vertex string
     */
    void _from_string();

  protected:
    /*! Vertex index */
    int m_index;

    /*! Vertex color */
    CPoint m_rgb;
};

inline void CHarmonicMapVertex::_from_string()
{
    CParser parser(m_string);
    for (std::list<CToken*>::iterator iter = parser.tokens().begin(); 
         iter != parser.tokens().end(); ++iter)
    {
        CToken* token = *iter;

        if (token->m_key == "rgb")
        {
            token->m_value >> m_rgb;
        }
    }
}


/*! \brief CHarmonicMapEdge class
 *
 *   Edge class for cut graph algorithm
 *   Trait : Edge sharp
 */
class CHarmonicMapEdge : public CEdge
{
  public:
    /*! Constructor */
    CHarmonicMapEdge() : m_length(0), m_weight(0) {};

    /*!	Edge weight */
    double& weight() { return m_weight; };
    
    /*! Edge length */
    double& length() { return m_length; };
    
  protected:
    /*!	Edge weight */
    double m_weight;

    /*! Edge length */
    double m_length;
};

/*! \brief CHarmonicMapFace class
 *
 *   Face class for cut graph algorithm
 *   Trait : Face touched flag
 */
class CHarmonicMapFace : public CFace
{
  public:
    /*! Constructor */
    CHarmonicMapFace() {};

    /*! face normal */
    CPoint& normal() { return m_normal; };
  protected:
    /*! face normal */
    CPoint m_normal;
};

/*! \brief CHarmonicMapHalfEdge class
 *
 *   HalfEdge class for cut graph algorithm
 */
class CHarmonicMapHalfEdge : public CHalfEdge
{
  public:
    /*!	CHarmonicHalfEdge constructor */
    CHarmonicMapHalfEdge() : m_angle(0) {};
    
    /*!	Corner angle */
    double& angle() { return m_angle; };

  protected:
    /*! Corner angle */
    double m_angle;
};

/*! \brief CHarmonicMapMesh class
 *
 *	Mesh class for cut graph algorithm
 *
 */
template <typename V, typename E, typename F, typename H>
class THarmonicMapMesh : public CBaseMesh<V, E, F, H>
{
  public:
    typedef V CVertex;
    typedef E CEdge;
    typedef F CFace;
    typedef H CHalfEdge;

    typedef CBoundary<V, E, F, H>                   CBoundary;
    typedef CLoop<V, E, F, H>                       CLoop;

    typedef MeshVertexIterator<V, E, F, H>          MeshVertexIterator;
    typedef MeshEdgeIterator<V, E, F, H>            MeshEdgeIterator;
    typedef MeshFaceIterator<V, E, F, H>            MeshFaceIterator;
    typedef MeshHalfEdgeIterator<V, E, F, H>        MeshHalfEdgeIterator;

    typedef VertexVertexIterator<V, E, F, H>        VertexVertexIterator;
    typedef VertexEdgeIterator<V, E, F, H>          VertexEdgeIterator;
    typedef VertexFaceIterator<V, E, F, H>          VertexFaceIterator;
    typedef VertexInHalfedgeIterator<V, E, F, H>    VertexInHalfedgeIterator;
    typedef VertexOutHalfedgeIterator<V, E, F, H>   VertexOutHalfedgeIterator;

    typedef FaceVertexIterator<V, E, F, H>          FaceVertexIterator;
    typedef FaceEdgeIterator<V, E, F, H>            FaceEdgeIterator;
    typedef FaceHalfedgeIterator<V, E, F, H>        FaceHalfedgeIterator;
};

typedef THarmonicMapMesh<CHarmonicMapVertex, CHarmonicMapEdge, CHarmonicMapFace, CHarmonicMapHalfEdge> CHarmonicMapMesh;
}

#endif // !_HARMONIC_MAP_MESH_H_
