/*!
 *  \file BaseHolomorphicForm.h
 *  \brief Algorithm for computing holomorphic 1-forms, Hodge star operator
 *  \author David Gu
 *  \date Document 10/12/2010
 *
 *  Algorithm for computing holomorphic 1-forms, Hodge Star operator
 *
 */
/********************************************************************************************************************************
 *
 *      BaseHolomorphic 1-form Class
 *
 *       Copyright (c) Stony Brook University
 *
 *    Purpose:
 *
 *       Compute the holomorphic 1-forms
 *
 *       David Gu June 27, 2008,  gu@cs.stonybrook.edu
 *
 *
 *    Input:
 *
 *           Original mesh, the mesh cut through a shortest path connecting one inner boundary and the exterior boundary
 *
 *    Output:
 *
 *           Closed non-exact Harmonic 1-form. and the mesh with UV coordinates.
 *
 *********************************************************************************************************************************/

#ifndef _BASE_HOLOMORPHIC_FORM_H_
#define _BASE_HOLOMORPHIC_FORM_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <list>
#include <math.h>
#include <queue>
#include <vector>

namespace MeshLib
{

/*! \brief CWedgeOperator class
 *
 *	Wedge Star Operator
 */
template <typename M>
class CWedgeOperator
{
  public:
    /*! CWedgeOperator constructor
     * \param solid0 harmonic 1-form \f$\omega_0\f$
     * \param solid1 harmonic 1-form \f$\omega_1\f$
     */
    CWedgeOperator(M* mesh0, M* mesh1);
    /*! CWedgeOperator destructor
     */
    ~CWedgeOperator();
    /*! wedge product
     * \return
     *	\f[
                            \int \omega_0 \wedge \omega_1
    \f]
     */
    double wedge_product();
    /*! wedge product
     * \return
     *	\f[
                            \int \omega_0 \wedge {}^*\omega_1
    \f]
     */
    double wedge_star_product();

  private:
    /*! Two input harmonic 1-forms $\f\omega_0,\omega_1\f$. */
    M* m_pMesh[2];
};

template <typename M>
double CWedgeOperator<M>::wedge_product()
{

    double p = 0;
    for (typename M::MeshFaceIterator fiter(m_pMesh[0]); !fiter.end(); ++fiter)
    {
        typename M::CFace* f0 = *fiter;
        typename M::CFace* f1 = m_pMesh[1]->idFace(f0->id());

        std::vector<typename M::CHalfEdge*> h0, h1;
        for (typename M::FaceHalfedgeIterator fhiter(f0); !fhiter.end(); ++fhiter)
        {
            typename M::CHalfEdge* h = *fhiter;
            h0.push_back(h);
        }

        for (typename M::FaceHalfedgeIterator fhiter(f1); !fhiter.end(); ++fhiter)
        {
            typename M::CHalfEdge* h = *fhiter;
            h1.push_back(h);
        }

        typename M::CVertex* s0 = m_pMesh[0]->halfedgeSource(h0.front());
        for (size_t i = 0; i < 3; i++)
        {
            typename M::CHalfEdge* h = h1.front();
            typename M::CVertex* s = m_pMesh[1]->halfedgeSource(h);
            if (s->id() != s0->id())
            {
                h1.erase(h1.begin());
                h1.push_back(h);
            }
            else
            {
                break;
            }
        }

        std::vector<double> du0, du1;
        for (size_t i = 0; i < 3; i++)
        {
            typename M::CHalfEdge* h = h0[i];
            typename M::CEdge* e = m_pMesh[0]->halfedgeEdge(h);
            double du = (h == m_pMesh[0]->edgeHalfedge(e, 0)) ? e->du() : -e->du();
            du0.push_back(du);
        }

        for (size_t i = 0; i < 3; i++)
        {
            typename M::CHalfEdge* h = h1[i];
            typename M::CEdge* e = m_pMesh[1]->halfedgeEdge(h);
            double du = (h == m_pMesh[1]->edgeHalfedge(e, 0)) ? e->du() : -e->du();
            du1.push_back(du);
        }

        //TODO insert your code here
        //to compute p from du0[3], du1[3];
        //TODO insert your code here,

    }
    return p;
};

// Hodge Star operator
template <typename M>
double CWedgeOperator<M>::wedge_star_product()
{

    double p = 0;
    for (typename M::MeshFaceIterator fiter(m_pMesh[0]); !fiter.end(); ++fiter)
    {
        typename M::CFace* f0 = *fiter;
        typename M::CFace* f1 = m_pMesh[1]->idFace(f0->id());

        std::vector<typename M::CHalfEdge*> h0, h1;
        for (typename M::FaceHalfedgeIterator fhiter(f0); !fhiter.end(); ++fhiter)
        {
            typename M::CHalfEdge* h = *fhiter;
            h0.push_back(h);
        }

        for (typename M::FaceHalfedgeIterator fhiter(f1); !fhiter.end(); ++fhiter)
        {
            typename M::CHalfEdge* h = *fhiter;
            h1.push_back(h);
        }

        typename M::CVertex* s0 = m_pMesh[0]->halfedgeSource(h0.front());
        for (size_t i = 0; i < 3; i++)
        {
            typename M::CHalfEdge* h = h1.front();
            typename M::CVertex* s = m_pMesh[1]->halfedgeSource(h);
            if (s->id() != s0->id())
            {
                h1.erase(h1.begin());
                h1.push_back(h);
            }
            else
            {
                break;
            }
        }

        std::vector<double> du0, du1, theta;
        for (size_t i = 0; i < 3; i++)
        {
            typename M::CHalfEdge* h = h0[i];
            typename M::CEdge* e = m_pMesh[0]->halfedgeEdge(h);
            double du = (h == m_pMesh[0]->edgeHalfedge(e, 0)) ? e->du() : -e->du();
            du0.push_back(du);
            theta.push_back(h0[(i + 1) % 3]->angle());
        }

        for (size_t i = 0; i < 3; i++)
        {
            typename M::CHalfEdge* h = h1[i];
            typename M::CEdge* e = m_pMesh[1]->halfedgeEdge(h);
            double du = (h == m_pMesh[1]->edgeHalfedge(e, 0)) ? e->du() : -e->du();
            du1.push_back(du);
        }

        //TODO insert your code here
        //to compute p from du0[3], du1[3];
        //TODO insert your code here,

    }
    return p;
};

// CWedgeOperator constructor
//\param mesh0, mesh1 two input harmonic 1-forms
template <typename M>
CWedgeOperator<M>::CWedgeOperator(M* mesh0, M* mesh1)
{
    m_pMesh[0] = mesh0;
    m_pMesh[1] = mesh1;
};

// CWedgeOperator destructor
template <typename M>
CWedgeOperator<M>::~CWedgeOperator(){};

/*! \brief CHolomorphicForm class
 *
 *	Compute holomorphic forms on a mesh
 */
template <typename M>
class CBaseHolomorphicForm
{
  public:
    /*! CHolomorphicForm constructor
     *	\param meshes the list of meshes with stores the harmonic 1-forms, which form the basis
     *   of the first cohomology group of the mesh
     */
    CBaseHolomorphicForm(std::vector<M*>& meshes);
    /*! CHolomorphicForm destructor
     */
    ~CBaseHolomorphicForm();
    /*!	The list of meshes storing harmonic form bases
     */
    std::vector<M*>& meshes() { return m_meshes; };

    /*! Compute the conjugate harmonic 1-forms for the base harmonic 1-forms by solving the equation,
            Assume \f$ {}^*\omega_i = \sum_j \lambda_{ij} \omega_j \f$, then
            \f[
                    \int_M \omega_k \wedge {}^*\omega_i = \sum_j \lambda_{ij} \int_M \omega_k \wedge \omega_j
            \f]
            the left hand side can be computed using face-vector represenation.
    */
    void conjugate();

  protected:
    /*!	The list of meshes storing harmonic form bases
     */
    std::vector<M*> m_meshes;
};

// CBaseHolomorphicForm constructor
//\param meshes are the basis of harmonic 1-form group
template <typename M>
CBaseHolomorphicForm<M>::CBaseHolomorphicForm(std::vector<M*>& meshes)
{
    for (typename std::vector<M*>::iterator miter = meshes.begin(); miter != meshes.end(); miter++)
    {
        M* pM = *miter;
        m_meshes.push_back(pM);
    }
};

// CBaseHolomorphicForm destructor
template <typename M>
CBaseHolomorphicForm<M>::~CBaseHolomorphicForm(){};

// Compute the conjugate harmonic 1-forms for the base harmonic 1-forms
template <typename M>
void CBaseHolomorphicForm<M>::conjugate()
{

    //_angle_structure();

    int n = m_meshes.size();

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            CWedgeOperator<M> wo(m_meshes[i], m_meshes[j]);
            A(i, j) = wo.wedge_product();
            // std::cout << "Wedge Product " << i << "," << j << " " << A(i, j) << " ";
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < n; i++)
    {
        Eigen::VectorXd b(n);

        for (int j = 0; j < n; j++)
        {
            CWedgeOperator<M> wo(m_meshes[i], m_meshes[j]);
            b(j) = wo.wedge_star_product();
        }

        // Eigen::VectorXd x = solver.solve(b);
        Eigen::VectorXd x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

        std::cout << "Hodge Star Coefficients: ";
        for (int j = 0; j < n; j++)
        {
            std::cout << x(j) << " ";
        }
        std::cout << std::endl;

        for (typename M::MeshEdgeIterator eiter(m_meshes[i]); !eiter.end(); ++eiter)
        {
            typename M::CEdge* e = *eiter;
            e->duv()[0] = e->du();
            e->duv()[1] = 0;

            int id1 = m_meshes[i]->edgeVertex1(e)->id();
            int id2 = m_meshes[i]->edgeVertex2(e)->id();

            for (int k = 0; k < n; k++)
            {
                typename M::CVertex* w1 = m_meshes[k]->idVertex(id1);
                typename M::CVertex* w2 = m_meshes[k]->idVertex(id2);

                typename M::CEdge* edge = m_meshes[k]->vertexEdge(w1, w2);
                e->duv()[1] += edge->du() * x(k);
            }
        }
    }
};

} // namespace MeshLib
#endif