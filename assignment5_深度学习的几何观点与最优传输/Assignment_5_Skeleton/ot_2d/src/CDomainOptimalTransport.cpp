/*!  \file   DomainOptimalTransport.cpp
 *   \brief  OptimalTransport
 *   \author David Gu
 *   \date   documented on 05/18/2020
 *
 *   Class for Planar Optimal Transportation, Semi-Discrete Algorithm
 */

#include <Eigen/Eigen>

#include "CDomainOptimalTransport.h"

namespace MeshLib
{
CDomainOptimalTransport::CDomainOptimalTransport(COMTMesh* pMesh) : CBaseOT(pMesh)
{
    m_pWDT = NULL;
    m_domainTr = NULL;
    m_outputTr = NULL;
}

CDomainOptimalTransport ::~CDomainOptimalTransport()
{
    delete m_domainTr;
    delete m_outputTr;
}

/* initialize the potential function as quadratic
 */
void CDomainOptimalTransport::_initialize()
{
    /* Total target area */
    double total_target_area;

    /*! compute the domain triangulation, a convex polygon */
    __detri2_generate_disk(m_domainTr, total_target_area);

    /*! set the target measure */
    _set_target_measure(m_pMesh, total_target_area);

    /*! initialize the vertex index, starting from zero */
    index(m_pMesh);

    /* initialize the vertex weight */
    for (COMTMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
    {
        COMTMesh::CVertex* pv = *viter;
        pv->weight() = 0;
    }

    /*! normalize the vertex uv coordinates */
    _normalize_uv(m_pMesh);

    /*! Compute Weighted Delaunay of Base Mesh Vertices */
    __detri2_WDT(m_pMesh, &m_outputTr);
    
    /*! convert the Weighted Delaunay Triangulation to a mesh */
    __detri2_to_mesh(m_outputTr, m_domainTr, m_pWDT);

    /*! copy vertex_weight, vertex_target_area, vertex_index from base mesh to the Weighted Delaunay Triangulation Mesh */
    _copy_mesh(m_pMesh, m_pWDT);
};

/*! Gradient Descende method to compute the OT Map */
void CDomainOptimalTransport::__gradient_descend(
    COMTMesh* pInput,  // input  mesh, the vertex->target_area, dual_area need to be set
    COMTMesh*& pOutput // output mesh, the vertex->target_area, dual_area, dual_cell, uv (dual cell center) have been set, the connectivity is updated
)
{
    /*! compute the gradient, which equals to (target_area - dual_area) */
    for (COMTMesh::MeshVertexIterator viter(pInput); !viter.end(); viter++)
    {
        COMTMesh::CVertex* pv = *viter;
        //insert your code here
        //to compute the gradient
        double grad = 0.0;
        pv->update_direction() = grad;
    }

    /*! set initial step length */
    double step_length = 0.1;
    /*! pointer to the output WDT */
    detri2::Triangulation* pTr = NULL;
    /*! temparary mesh */
    COMTMesh* pMesh = NULL;

    /*! damping iteration */
    while (true)
    {
        // update the vertex weight
        for (COMTMesh::MeshVertexIterator viter(pInput); !viter.end(); viter++)
        {
            COMTMesh::CVertex* pv = *viter;
            pv->weight() -= step_length * pv->update_direction();
        }

        // compute Weighted Delaunay Triangulation
        if (!__detri2_WDT(pInput, &pTr))
        {
            // if there are missing points,
            // roll back the vertex weight
            for (COMTMesh::MeshVertexIterator viter(pInput); !viter.end(); viter++)
            {
                COMTMesh::CVertex* pv = *viter;
                pv->weight() += step_length * pv->update_direction();
            }
            delete pTr;
            pTr = NULL;
            // reduce the step length by half
            step_length /= 2.0;
            // try again
            continue;
        }
        // no missing point

        // convert the detri2::Triangulation to OMTMesh
        // this will set the vertex->dual_area, vertex->dual_cell, vertex->dual_center;
        // edge->length, edge->dual_length;
        __detri2_to_mesh(pTr, m_domainTr, pMesh);
        // copy vertex->target_carea, vertex->weight, vertex->index from the input mesh
        // to the newly generated WDT mesh
        _copy_mesh(pInput, pMesh);

        // if the old WDT pointer is nonempty, delete it
        if (m_outputTr != NULL)
        {
            delete m_outputTr;
        }
        // update the WDT pointer by the current one
        m_outputTr = pTr;
        pOutput = pMesh;

        break;
    }
 };

/*! Newton's method to compute the OT Map */
void CDomainOptimalTransport::__newton(COMTMesh* pInput,  // input mesh, with vertex->target_area, vertex->dual_area, vertex->weight, vertex->index
                                       COMTMesh*& pOutput // output mesh, with copied vertex->target_area, vertex->index,
                                                          // and updated vertex->dual_area, vertex->weight, vertex->dual_cell, vertex->dual_center, connectivity
)
{
    /*! Use Newton's method to compute the update_direction */
    __update_direction(pInput);

    /*! set initial step length */
    double step_length = 1.0;
    /*! pointer to the output WDT */
    detri2::Triangulation* pTr = NULL;
    /*! temparary mesh */
    COMTMesh* pMesh = NULL;

    /*! damping iteration */
    while (true)
    {
        // update the vertex weight
        for (COMTMesh::MeshVertexIterator viter(pInput); !viter.end(); viter++)
        {
            COMTMesh::CVertex* pv = *viter;
            pv->weight() -= step_length * pv->update_direction();
        }

        // compute Weighted Delaunay Triangulation
        bool success = __detri2_WDT(pInput, &pTr);

        if (!success )
        {
            //insert your code here
            // if there are missing points,
            // roll back the vertex weight
            // reduce the step length by half

            delete pTr;
            pTr = NULL;
            continue;
        }

        // no missing point
        
        // convert the detri2::Triangulation to OMTMesh
        // this will set the vertex->dual_area, vertex->dual_cell, vertex->dual_center;
        // edge->length, edge->dual_length;
        __detri2_to_mesh(pTr, m_domainTr, pMesh);
        // copy vertex->target_carea, vertex->weight, vertex->index from the input mesh
        // to the newly generated WDT mesh
        _copy_mesh(pInput, pMesh);

        // if the old WDT pointer is nonempty, delete it
        if (m_outputTr != NULL)
        {
            delete m_outputTr;
        }
        // update the WDT pointer by the current one
        m_outputTr = pTr;
        pOutput = pMesh;

        break;
    }
}

} // namespace MeshLib