#ifndef _HODGE_DECOMPOSITION_MESH_H_
#define _HODGE_DECOMPOSITION_MESH_H_

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
class CHodgeDecompositionVertex;
class CHodgeDecompositionEdge;
class CHodgeDecompositionFace;
class CHodgeDecompositionHalfEdge;

/*! \brief CHodgeDecompositionVertex class
 */
class CHodgeDecompositionVertex : public CVertex
{
  public:
    /*! Constructor */
    CHodgeDecompositionVertex() : 
        m_index(0), 
        m_form(0),
        m_rgb(229.0 / 255.0, 162.0 / 255.0, 141.0 / 255.0),
        m_touched(false),
        m_father(0) {};
    
    /*! Vertex index */
    int& idx() { return m_index; };

    /*! Vertex color */
    CPoint& rgb() { return m_rgb; };

    /*!
     *	Read vertex traits to vertex string
     */
    void _from_string();

    /*! zero form */
    double& form() { return m_form; };

    /*! Wheter the vertex has been accessed */
    bool& touched() { return m_touched; };

    /*! Vertex father id */
    int& father() { return m_father; };

  protected:
    /*! Vertex index */
    int m_index;

    /*! Vertex color */
    CPoint m_rgb;

    /*! 0-form */
    double m_form;

    /*! Vertex touched trait */
    bool m_touched;

    /*! Vertex father id */
    int m_father;
};

inline void CHodgeDecompositionVertex::_from_string()
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

        if (token->m_key == "father")
        {
            std::string line = strutil::trim(token->m_value, "()");
            m_father = strutil::parseString<int>(line);
        }
    }
}


/*! \brief CHodgeDecompositionEdge class
 */
class CHodgeDecompositionEdge : public CEdge
{
  public:
    /*! Constructor */
    CHodgeDecompositionEdge() : m_length(0), m_weight(0), m_form(0) {};

    /*!	Edge weight */
    double& weight() { return m_weight; };
    
    /*! Edge length */
    double& length() { return m_length; };
    
    /*! form */
    double& du() { return m_form; };

    /*! read one-form from file */
    void _from_string();

    /*! duv */
    CPoint2& duv() { return m_duv; };

protected:
    /*!	Edge weight */
    double m_weight;

    /*! Edge length */
    double m_length;

    /*! form */
    double m_form;

    /*! duv */
    CPoint2 m_duv;
};

//read harmonic 1-form trait "du" to the trait m_du
inline void CHodgeDecompositionEdge::_from_string()
{
    CParser parser(m_string);
    for (std::list<CToken*>::iterator titer = parser.tokens().begin(); titer != parser.tokens().end(); titer++)
    {
        CToken* pT = *titer;
        if (pT->m_key == "du")
        {
            std::string line = strutil::trim(pT->m_value, "()");
            m_form = strutil::parseString<double>(line);
            break;
        }
    }
};

/*! \brief CHodgeDecompositionFace class
 */
class CHodgeDecompositionFace : public CFace
{
  public:
    /*! Constructor */
      CHodgeDecompositionFace() { m_form = 0; m_index = 0; };

    /*! face normal */
    CPoint& normal() { return m_normal; };

    /*! face 2-form */
    double& form() { return m_form; };

    /*! index */
    int& idx() { return m_index; };

  protected:
    /*! face normal */
    CPoint m_normal;

    /*! 2-form */
    double m_form;

    /*! index */
    int m_index;
};

/*! \brief CHodgeDecompositionHalfEdge class
 */
class CHodgeDecompositionHalfEdge : public CHalfEdge
{
  public:
    /*!	CHarmonicHalfEdge constructor */
      CHodgeDecompositionHalfEdge() : m_angle(0) { m_form = 0; };
    
    /*!	Corner angle */
    double& angle() { return m_angle; };

    /*! 1-form */
    double& form() { return m_form; };

  protected:
    /*! Corner angle */
    double m_angle;

    /* 1-form */
    double m_form;
};

/*! \brief CHodgeDecompositionMesh class
 *
 *	Mesh class for cut graph algorithm
 *
 */
template <typename V, typename E, typename F, typename H>
class THodgeDecompositionMesh : public CBaseMesh<V, E, F, H>
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

typedef THodgeDecompositionMesh<CHodgeDecompositionVertex, CHodgeDecompositionEdge, CHodgeDecompositionFace, CHodgeDecompositionHalfEdge> CHodgeDecompositionMesh;
}

#endif // !_HARMONIC_MAP_MESH_H_
