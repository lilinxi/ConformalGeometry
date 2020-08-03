#include <Eigen/Sparse>
#include <float.h>
#include <math.h>
#include <time.h>

#include "HodgeDecomposition.h"
#include "WedgeProduct.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279
#endif

MeshLib::CHodgeDecomposition::CHodgeDecomposition()
{
    m_pMesh = NULL;
}

void MeshLib::CHodgeDecomposition::set_mesh(CHodgeDecompositionMesh* pMesh)
{
    using M = CHodgeDecompositionMesh;

    m_pMesh = pMesh;

    // 1. index all the vertices
    int vid = 0;
    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
    {
        M::CVertex* pV = *viter;
        pV->idx() = vid++;
    }

    // 2. index all the faces
    int fid = 0;
    for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter)
    {
        M::CFace* pf = *fiter;
        pf->idx() = fid++;
    }
}

void MeshLib::CHodgeDecomposition::_calculate_edge_weight(bool using_geometry)
{
    using M = CHodgeDecompositionMesh;

    if (!using_geometry)
    {
        for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter)
        {
            M::CEdge* pE = *eiter;
            pE->weight() = (!pE->boundary()) ? 1.0 : 0.5;
        }
        return;
    }
    // using geometry to compute the cotangent edge weight

    // 1. compute edge length
    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter)
    {
        M::CEdge* pE = *eiter;
        M::CVertex* v1 = m_pMesh->edgeVertex1(pE);
        M::CVertex* v2 = m_pMesh->edgeVertex2(pE);
        pE->length() = (v1->point() - v2->point()).norm();
    }

    // 2. compute corner angle
    for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter)
    {
        M::CFace* pF = *fiter;
        M::CHalfEdge* pH[3];
        pH[0] = m_pMesh->faceHalfedge(pF);
        for (int i = 0; i < 3; i++)
        {
            pH[(i + 1) % 3] = m_pMesh->halfedgeNext(pH[i]);
        }

        double len[3];
        for (int i = 0; i < 3; i++)
        {
            len[i] = m_pMesh->halfedgeEdge(pH[i])->length();
        }

        for (int i = 0; i < 3; i++)
        {
            double a = len[(i + 1) % 3], b = len[(i + 2) % 3], c = len[i];
            pH[(i + 1) % 3]->angle() = _inverse_cosine_law(a, b, c);
        }
    }

    // 3. compute edge weight
    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter)
    {
        M::CEdge* pE = *eiter;

        if (!pE->boundary())
        {
            double theta[2];
            theta[0] = m_pMesh->halfedgeNext(m_pMesh->edgeHalfedge(pE, 0))->angle();
            theta[1] = m_pMesh->halfedgeNext(m_pMesh->edgeHalfedge(pE, 1))->angle();
            pE->weight() = std::cos(theta[0]) / std::sin(theta[0]) + std::cos(theta[1]) / std::sin(theta[1]);
        }
        else
        {
            double theta = m_pMesh->halfedgeNext(m_pMesh->edgeHalfedge(pE, 0))->angle();
            pE->weight() = std::cos(theta) / std::sin(theta);
        }
    }
}

double MeshLib::CHodgeDecomposition::_inverse_cosine_law(double a, double b, double c)
{
    double cs = (a * a + b * b - c * c) / (2.0 * a * b);
    assert(cs <= 1.0 && cs >= -1.0);
    return std::acos(cs);
}

void MeshLib::CHodgeDecomposition::_d(int dimension)
{
    if (dimension >= 2)
        return;

    using M = CHodgeDecompositionMesh;

    if (dimension == 0)
    {
        for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
        {
            M::CFace* pf = *fiter;
            for (M::FaceHalfedgeIterator fhiter(pf); !fhiter.end(); fhiter++)
            {
                M::CHalfEdge* ph = *fhiter;
                //TODO insert your code here, 
                //convert vertex->form() to halfedge->form()
                //TODO insert your code here,
            }
        }
        return;
    }

    if (dimension == 1)
    {
        for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
        {
            M::CFace* pf = *fiter;
            pf->form() = 0;
            for (M::FaceHalfedgeIterator fhiter(pf); !fhiter.end(); fhiter++)
            {
                M::CHalfEdge* ph = *fhiter;
                //TODO insert your code here, 
                //convert halfedge->form() to face->form()
                //TODO insert your code here,


            }
        }
        return;
    }
}

void MeshLib::CHodgeDecomposition::_delta(int dimension)
{
    using M = CHodgeDecompositionMesh;

    if (dimension == 2)
    {
        for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
        {
            M::CFace* pf = *fiter;
            for (M::FaceHalfedgeIterator fhiter(pf); !fhiter.end(); fhiter++)
            {
                M::CHalfEdge* ph = *fhiter;
                //TODO insert your code here, 
                //convert face->form() to halfedge->form()
                //TODO insert your code here,


            }
        }
        return;
    }

    if (dimension == 1)
    {
        for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
        {
            M::CVertex* pv = *viter;
            //TODO insert your code here, 
            //convert halfedge->form() to vertex->form()
            //TODO insert your code here,

        }
        return;
    }
}

void MeshLib::CHodgeDecomposition::_remove_exact_form()
{
    using M = CHodgeDecompositionMesh;
    M* pM = m_pMesh;
    // remove exact-form
    for (M::MeshFaceIterator fiter(pM); !fiter.end(); fiter++)
    {
        M::CFace* pf = *fiter;
        for (M::FaceHalfedgeIterator fhiter(pf); !fhiter.end(); fhiter++)
        {
            M::CHalfEdge* ph = *fhiter;
            //TODO insert your code here,
            //remove d vertex->form() from halfedge->form()

            //TODO insert your code here,

        }
    }
}

void MeshLib::CHodgeDecomposition::_normalize()
{
    using M = CHodgeDecompositionMesh;

    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); eiter++)
    {
        M::CEdge* pe = *eiter;
        M::CHalfEdge* ph = m_pMesh->edgeHalfedge(pe, 0);
        pe->du() = ph->form();
    }
    CWedgeOperator<M> W(m_pMesh, m_pMesh);
    double p = W.wedge_star_product();
    p = sqrt(p);

    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); eiter++)
    {
        M::CEdge* pe = *eiter;
        M::CHalfEdge* ph = m_pMesh->edgeHalfedge(pe, 0);
        pe->du() /= p;
    }

    for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
    {
        M::CFace* pf = *fiter;
        double w = 0;
        for (M::FaceHalfedgeIterator fhiter(pf); !fhiter.end(); fhiter++)
        {
            M::CHalfEdge* ph = *fhiter;
            ph->form() /= p;
        }
    }
}

void MeshLib::CHodgeDecomposition::_compute_exact_form()
{

    using M = CHodgeDecompositionMesh;

    // 1. Initialize

    // 2. Set the matrix A
    std::vector<Eigen::Triplet<double>> A_coefficients;

    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
    {
        M::CVertex* pV = *viter;

        int vid = pV->idx();

        double sw = 0;
        for (M::VertexVertexIterator witer(pV); !witer.end(); ++witer)
        {
            M::CVertex* pW = *witer;
            int wid = pW->idx();

            M::CEdge* e = m_pMesh->vertexEdge(pV, pW);
            double w = e->weight();
            A_coefficients.push_back(Eigen::Triplet<double>(vid, wid, w));

            sw += w;
        }
        A_coefficients.push_back(Eigen::Triplet<double>(vid, vid, -sw));
    }

    Eigen::SparseMatrix<double> A(m_pMesh->numVertices(), m_pMesh->numVertices());
    A.setZero();
    A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());

    // 3. compute delta tau
    _delta(1);

    // 4. Solve the equations
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    std::cerr << "Eigen Decomposition" << std::endl;
    solver.compute(A);
    std::cerr << "Eigen Decomposition Finished" << std::endl;

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Waring: Eigen decomposition failed" << std::endl;
    }

    Eigen::VectorXd b(m_pMesh->numVertices());
    // set boundary constraints vector b
    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
    {
        M::CVertex* pV = *viter;
        int id = pV->idx();
        b(id) = pV->form();
    }

    Eigen::VectorXd x = solver.solve(b); // Ax=b
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Waring: Eigen decomposition failed" << std::endl;
    }

    // set the images of the harmonic map to interior vertices
    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
    {
        M::CVertex* pV = *viter;
        int id = pV->idx();
        pV->form() = x(id);
    }
}

void MeshLib::CHodgeDecomposition::_random_form()
{
    using M = CHodgeDecompositionMesh;

    // std::srand((unsigned)time(NULL));
    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); eiter++)
    {
        double r = 2.0 * ((double) std::rand() / RAND_MAX - 0.5);
        M::CEdge* pe = *eiter;
        M::CHalfEdge* ph = m_pMesh->edgeHalfedge(pe, 0);
        ph->form() = r;
        M::CHalfEdge* ps = m_pMesh->edgeHalfedge(pe, 1);
        if (ps != NULL)
        {
            ps->form() = -r;
        }
    }
}

void MeshLib::CHodgeDecomposition::_compute_coexact_form()
{

    using M = CHodgeDecompositionMesh;

    // 1. Initialize

    // 2. Set the matrix A
    std::vector<Eigen::Triplet<double>> A_coefficients;

    for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter)
    {
        M::CFace* pf = *fiter;

        int fid = pf->idx();

        double sw = 0;
        for (M::FaceHalfedgeIterator fhiter(pf); !fhiter.end(); ++fhiter)
        {
            M::CHalfEdge* ph = *fhiter;
            M::CHalfEdge* ps = m_pMesh->halfedgeSym(ph);
            if (ps != NULL)
            {
                //TODO insert your code here
                //TODO insert your code here,
            }
            else // boundary face
            {
                //TODO insert your code here
                //TODO insert your code here,
            }
        }
        A_coefficients.push_back(Eigen::Triplet<double>(fid, fid, -sw));
    }

    Eigen::SparseMatrix<double> A(m_pMesh->numFaces(), m_pMesh->numFaces());
    A.setZero();
    A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());

    // 3. compute delta tau
    _d(1);

    // 4. Solve the equations
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    std::cerr << "Eigen Decomposition" << std::endl;
    solver.compute(A);
    std::cerr << "Eigen Decomposition Finished" << std::endl;

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Waring: Eigen decomposition failed" << std::endl;
    }

    Eigen::VectorXd b(m_pMesh->numFaces());
    // set boundary constraints vector b
    for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter)
    {
        M::CFace* pF = *fiter;
        int id = pF->idx();
        b(id) = pF->form();
    }

    Eigen::VectorXd x = solver.solve(b); // Ax=b
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Waring: Eigen decomposition failed" << std::endl;
    }

    // set the images of the harmonic map to interior vertices
    for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter)
    {
        M::CFace* pF = *fiter;
        int id = pF->idx();
        pF->form() = x(id);
    }
}

void MeshLib::CHodgeDecomposition::_test_closedness()
{
    double max_error = -1e+10;
    double total_squared_error = 0;

    using M = CHodgeDecompositionMesh;
    for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); fiter++)
    {
        M::CFace* pf = *fiter;
        double w = 0;
        for (M::FaceHalfedgeIterator fhiter(pf); !fhiter.end(); fhiter++)
        {
            M::CHalfEdge* ph = *fhiter;
            w += ph->form();
        }
        max_error = (max_error > fabs(w)) ? max_error : fabs(w);
        total_squared_error += w * w;
    }
    total_squared_error /= m_pMesh->numFaces();
    total_squared_error = sqrt(total_squared_error); // root mean error

    std::cout << "Closedness: Max Error " << max_error << "  Root Mean Squared Error " << total_squared_error << std::endl;
}

void MeshLib::CHodgeDecomposition::_remove_coexact_form()
{
    using M = CHodgeDecompositionMesh;
    M* pM = m_pMesh;
    // remove coexact-form
    for (M::MeshFaceIterator fiter(pM); !fiter.end(); fiter++)
    {
        M::CFace* pf = *fiter;
        for (M::FaceHalfedgeIterator fhiter(pf); !fhiter.end(); fhiter++)
        {
            //TODO insert your code here
            //remove delta face->from() from halfedge->form()
            //TODO insert your code here,
        }
    }
}

void MeshLib::CHodgeDecomposition::random_harmonic_form()
{
    _random_form();

    _calculate_edge_weight(false);
    _compute_coexact_form();
    _remove_coexact_form();
    _test_closedness();

    _calculate_edge_weight(true);
    _compute_exact_form();
    _remove_exact_form();
    _test_coclosedness();
    _normalize();
}

void MeshLib::CHodgeDecomposition::_test_coclosedness()
{
    using M = CHodgeDecompositionMesh;
    double max_error = -1e+10;
    double total_squared_error = 0;
    int interior_vertices = 0;

    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); viter++)
    {
        M::CVertex* pv = *viter;
        double w = 0;
        for (M::VertexOutHalfedgeIterator vhiter(m_pMesh, pv); !vhiter.end(); vhiter++)
        {
            M::CHalfEdge* ph = *vhiter;
            M::CEdge* pe = m_pMesh->halfedgeEdge(ph);
            w += pe->weight() * ph->form();
        }
        if (pv->boundary())
        {
            M::CHalfEdge* ph = m_pMesh->vertexMostCcwInHalfEdge(pv);
            M::CEdge* pe = m_pMesh->halfedgeEdge(ph);
            w -= pe->weight() * ph->form();
        }
        if (pv->boundary())
            continue;
        interior_vertices++;
        max_error = (max_error > fabs(w)) ? max_error : fabs(w);
        total_squared_error += w * w;
    }
    
    total_squared_error /= interior_vertices;
    total_squared_error = sqrt(total_squared_error); // root mean error

    std::cout << "CoClosedness: Max Error " << max_error << "  Root Mean Squared Error " << total_squared_error << std::endl;
}

void MeshLib::CHodgeDecomposition::_exact_harmonic_form()
{
    
    using M = CHodgeDecompositionMesh;

    // 1. Initialize
    int vid = 0; // interior vertex id
    int bid = 0; // boundary vertex id
    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
    {
        M::CVertex* pV = *viter;

        if (pV->boundary())
            pV->idx() = bid++;
        else
            pV->idx() = vid++;
    }

    int interior_vertices = vid;
    int boundary_vertices = bid;

    // 2. Set the matrix A and B
    std::vector<Eigen::Triplet<double>> A_coefficients;
    std::vector<Eigen::Triplet<double>> B_coefficients;

    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
    {
        M::CVertex* pV = *viter;
        if (pV->boundary())
            continue;
        int vid = pV->idx();

        double sw = 0;
        for (M::VertexVertexIterator witer(pV); !witer.end(); ++witer)
        {
            M::CVertex* pW = *witer;
            int wid = pW->idx();

            M::CEdge* e = m_pMesh->vertexEdge(pV, pW);
            double w = e->weight();

            if (pW->boundary())
            {
                B_coefficients.push_back(Eigen::Triplet<double>(vid, wid, w));
            }
            else
            {
                A_coefficients.push_back(Eigen::Triplet<double>(vid, wid, -w));
            }
            sw += w;
        }
        A_coefficients.push_back(Eigen::Triplet<double>(vid, vid, sw));
    }

    Eigen::SparseMatrix<double> A(interior_vertices, interior_vertices);
    A.setZero();
    Eigen::SparseMatrix<double> B(interior_vertices, boundary_vertices);
    B.setZero();
    A.setFromTriplets(A_coefficients.begin(), A_coefficients.end());
    B.setFromTriplets(B_coefficients.begin(), B_coefficients.end());

    // 3. Solve the equations
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    std::cerr << "Eigen Decomposition" << std::endl;
    solver.compute(A);
    std::cerr << "Eigen Decomposition Finished" << std::endl;

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "Waring: Eigen decomposition failed" << std::endl;
    }

    {
        Eigen::VectorXd b(boundary_vertices);
        // set boundary constraints vector b
        for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
        {
            M::CVertex* pV = *viter;
            if (!pV->boundary())
                continue;
            int id = pV->idx();
            b(id) = pV->form();
        }

        Eigen::VectorXd c(interior_vertices);
        c = B * b;

        Eigen::VectorXd x = solver.solve(c); // Ax=c
        if (solver.info() != Eigen::Success)
        {
            std::cerr << "Waring: Eigen decomposition failed" << std::endl;
        }

        // set the images of the harmonic map to interior vertices
        for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
        {
            M::CVertex* pV = *viter;
            if (pV->boundary())
                continue;
            int id = pV->idx();
            pV->form() = x(id);
        }
    }
}

void MeshLib::CHodgeDecomposition::_set_boundary_condition(int bid)
{
    using M = CHodgeDecompositionMesh;
    M::CBoundary bnd(m_pMesh);
    std::vector<M::CLoop*>& Ls = bnd.loops();

    for (size_t i = 0; i < Ls.size(); i++)
    {
        std::list<M::CHalfEdge*>& pHs = Ls[i]->halfedges();
        for (std::list<M::CHalfEdge*>::iterator hiter = pHs.begin(); hiter != pHs.end(); hiter++)
        {
            M::CHalfEdge* ph = *hiter;
            M::CVertex* pv = m_pMesh->halfedgeTarget(ph);
            pv->form() = (i == bid) ? -1 : 0;
        }
    }
}

void MeshLib::CHodgeDecomposition::exact_harmonic_form(int bnd_id)
{
    using M = CHodgeDecompositionMesh;

    _set_boundary_condition(bnd_id);
    _calculate_edge_weight(true);
    _exact_harmonic_form();
    _d(0);
    _test_closedness();
    _test_coclosedness();
    _normalize();
}

void MeshLib::CHodgeDecomposition::integration(
    CHodgeDecompositionMesh* pForm, 
    CHodgeDecompositionMesh* pDomain)
{
    using M = CHodgeDecompositionMesh;

    M::CVertex* head = NULL;

    for (M::MeshVertexIterator viter(pDomain); !viter.end(); viter++)
    {
        M::CVertex* v = *viter;
        v->touched() = false;
        head = v;
    }

    std::queue<M::CVertex*> vqueue;
    head->uv() = CPoint2(0, 0);
    head->touched() = true;

    vqueue.push(head);

    while (!vqueue.empty())
    {
        head = vqueue.front();
        vqueue.pop();

        for (M::VertexEdgeIterator veiter(head); !veiter.end(); veiter++)
        {
            M::CEdge* e = *veiter;
            M::CVertex* v1 = pDomain->edgeVertex1(e);
            M::CVertex* v2 = pDomain->edgeVertex2(e);

            M::CVertex* tail = (head != v1) ? v1 : v2;
            if (tail->touched())
                continue;

            tail->touched() = true;
            vqueue.push(tail);

            int id1 = v1->father();
            // if there is no "father" field for the vertex, then directly use the vertex id
            if (id1 == 0)
                id1 = v1->id();

            int id2 = v2->father();
            // if there is no "father" field for the vertex, then directly use the vertex id
            if (id2 == 0)
                id2 = v2->id();

            M::CVertex* w1 = pForm->idVertex(id1);
            M::CVertex* w2 = pForm->idVertex(id2);

            M::CEdge* we = pForm->vertexEdge(w1, w2);

            if (pForm->edgeVertex1(we) == w1)
            {
                e->duv() = we->duv();
            }
            else
            {
                e->duv() = CPoint2(0, 0) - we->duv();
            }

            if (tail == v2)
            {
                tail->uv() = head->uv() + e->duv();
            }
            else
            {
                tail->uv() = head->uv() - e->duv();
            }
        }
    }
}
