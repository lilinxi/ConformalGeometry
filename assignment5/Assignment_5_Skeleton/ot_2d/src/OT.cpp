/*!  \file   OptimalTransport.h
 *   \brief  OptimalTransport
 *   \author David Gu
 *   \date   documented on 05/18/2020
 *
 *   Base Class for Planar Optimal Transportation, Semi-Discrete Algorithm
 */

#include "OT.h"

namespace MeshLib
{
/*! copy vertex weight, target area and index from the input mesh to the output mesh */
void CBaseOT::_copy_mesh(COMTMesh* pInput, COMTMesh* pOutput)
{
    for (COMTMesh::MeshVertexIterator viter(pInput), witer(pOutput); !viter.end(); viter++, witer++)
    {
        COMTMesh::CVertex* pv = *viter;
        COMTMesh::CVertex* pw = *witer;
        pw->weight() = pv->weight();
        pw->target_area() = pv->target_area();
        pw->index() = pv->index();
        pw->rgb() = pv->rgb();
        pw->normal() = pv->normal();
        pw->string() = pv->string();
    }
}

/*! normalize the uv coordinates to be within the unit disk */
void CBaseOT::_normalize_uv(COMTMesh* pMesh)
{
    double total_area = 0;
    // calculate the total area of the mesh
    for (COMTMesh::MeshFaceIterator fiter(pMesh); !fiter.end(); ++fiter)
    {
        COMTMesh::CFace* pf = *fiter;
        std::vector<CPoint2> uvs;
        for (COMTMesh::FaceVertexIterator fviter(pf); !fviter.end(); fviter++)
        {
            COMTMesh::CVertex* pv = *fviter;
            uvs.push_back(pv->uv());
        }
        total_area += (uvs[1] - uvs[0]) ^ (uvs[2] - uvs[0]);
    }

    std::cout << "Total Area " << total_area << std::endl;

    CPoint2 s(0, 0);
    double total_length = 0;
    for (COMTMesh::MeshEdgeIterator eiter(pMesh); !eiter.end(); ++eiter)
    {
        COMTMesh::CEdge* pe = *eiter;
        if (!pe->boundary())
            continue;

        COMTMesh::CVertex* pv1 = pMesh->edgeVertex1(pe);
        COMTMesh::CVertex* pv2 = pMesh->edgeVertex2(pe);

        (pv1->uv() + pv2->uv()) * pMesh->edgeLength(pe) / 2.0;
        total_length += pMesh->edgeLength(pe);
    }

    s = s / total_length;

    for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        COMTMesh::CVertex* v = *viter;
        CPoint2 p = v->uv();
        p = p - s;
        v->uv() = p;
    }
    double d = 0;

    for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        COMTMesh::CVertex* v = *viter;
        CPoint2 p = v->uv();
        d = (d > p.norm()) ? d : p.norm();
    }

    for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        COMTMesh::CVertex* v = *viter;
        CPoint2 p = v->uv();
        p = p / d;
        v->uv() = p;
        // shrink into the small inner disk
        v->uv() *= 0.98;
    }
};

/*! compute the maximal, relative error of the input mesh */
void CBaseOT::_compute_error(COMTMesh* pMesh)
{
    double max_error = -1e+10;
    double total_error = 0;

    for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); viter++)
    {
        COMTMesh::CVertex* pv = *viter;
        double da = fabs(pv->target_area() - pv->dual_area());
        double error = da / pv->target_area();
        if (error > max_error)
        {
            max_error = error;
        }
        total_error += da * da;
    }
    std::cout << "Max relative error is " << max_error << " Total L2 error is " << total_error << std::endl;
};

/*! solve the linear system */

bool CBaseOT::__solve(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, Eigen::VectorXd& result)
{
    // solve the linear system
    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double> >	solver;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    std::cout << "Eigen decomposition:" << std::endl;
    solver.compute(A);
    std::cout << "Eigen decomposition finished" << std::endl;

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Waring: eigen decomposition failed!!!!" << std::endl;
        return false;
    }

    result = solver.solve(b);

    if (solver.info() == Eigen::Success)
    {
        return true;
    }
    std::cerr << "Waring: eigen solver failed!!!!" << std::endl;
    return false;
}

/*! Set the update direction for every height */

void CBaseOT::__update_direction(COMTMesh* m_pMesh)
{
    /* set vertex ID, index Lookup Table */

    std::vector<int> ids;
    ids.resize(m_pMesh->numVertices());
    for (COMTMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
    {
        COMTMesh::CVertex* pv = *viter;
        ids[pv->index()] = pv->id();
    }

    /* set gradient vector */

    Eigen::VectorXd m_gradient;
    m_gradient.resize(m_pMesh->numVertices());

    for (COMTMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
    {
        COMTMesh::CVertex* pv = *viter;
        //insert your code here, 
        //compute the gradient for each vertex
        double gradient = 0;
        m_gradient[pv->index()] = gradient;
    }

    /* compute the Hessian matrix */
    Eigen::SparseMatrix<double> hessian;
    __compute_hessian_matrix(*m_pMesh, hessian);

    /* solve hessian equation */
    Eigen::VectorXd m_direction;
    if (!__solve(hessian, m_gradient, m_direction))
    {
        std::cout << "Numerical Error" << std::endl;
    }
    else
    {
        for (int i = 0; i < m_pMesh->numVertices(); i++)
        {
            int id = ids[i];
            COMTMesh::CVertex* pv = m_pMesh->idVertex(id);
            pv->update_direction() = m_direction[i];
        }
    }
};

/*! set target area, assume the vertex area has been set already */

void CBaseOT::_set_target_measure(COMTMesh*& pMesh, double total_target_area)
{
    double total_area = 0;
    /* compute the vertex area */
    for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); viter++)
    {
        COMTMesh::CVertex* pV = *viter;
        double s = 0;
        for (COMTMesh::VertexFaceIterator vfiter(pV); !vfiter.end(); vfiter++)
        {
            COMTMesh::CFace* pF = *vfiter;
            s += pF->area();
        }
        pV->target_area() = s / 3.0;
        total_area += pV->target_area();
    }
    /*! set the target area proportional to the vertex area*/
    for (COMTMesh::MeshVertexIterator viter(pMesh); !viter.end(); viter++)
    {
        COMTMesh::CVertex* pV = *viter;
        pV->target_area() *= (total_target_area / total_area);
    }
};

/* compute the hessian matrix */

void CBaseOT::__compute_hessian_matrix(COMTMesh& mesh, Eigen::SparseMatrix<double>& hessian)
{
    std::vector<Eigen::Triplet<double>> hessian_coefficients;
    /* all the off diagonal elements */
    for (COMTMesh::MeshEdgeIterator eiter(&mesh); !eiter.end(); eiter++)
    {
        COMTMesh::CEdge* pe = *eiter;
        //insert your coder here
        //compute the off diagonal element in Hessian matrix
        double weight = 1.0;
        COMTMesh::CVertex* pv1 = mesh.edgeVertex1(pe);
        COMTMesh::CVertex* pv2 = mesh.edgeVertex2(pe);
        int ids = pv1->index();
        int idt = pv2->index();
        hessian_coefficients.push_back(Eigen::Triplet<double>(ids, idt, weight));
        hessian_coefficients.push_back(Eigen::Triplet<double>(idt, ids, weight));
    }
    /* all the diagonal elements */
    for (COMTMesh::MeshVertexIterator viter(&mesh); !viter.end(); viter++)
    {
        COMTMesh::CVertex* pv = *viter;
        double sum = 0;
        for (COMTMesh::VertexEdgeIterator veiter(pv); !veiter.end(); veiter++)
        {
            COMTMesh::CEdge* pe = *veiter;
            //insert your code here
            //compute the diagonal element in Hessian matrix
            double weight = 1.0;
            sum += weight;
        }
        int id = pv->index();
        hessian_coefficients.push_back(Eigen::Triplet<double>(id, id, -sum));
    }

    hessian.resize(mesh.numVertices(), mesh.numVertices());
    hessian.setZero();
    hessian.setFromTriplets(hessian_coefficients.begin(), hessian_coefficients.end());
};

/*! reindex all the vertices, starting from zero */
void CBaseOT::index(COMTMesh* m_pMesh)
{
    int idx = 0;
    for (COMTMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
    {
        COMTMesh::CVertex* pv = *viter;
        pv->index() = idx++;
    }
};

} // namespace MeshLib
