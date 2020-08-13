#include <float.h>
#include <cmath>

#include "SphericalHarmonicMap.h"

void MeshLib::CSphericalHarmonicMap::set_mesh(CSHMMesh *pMesh) {
    m_pMesh = pMesh;

    // 1. compute vertex normal
    _calculate_normal();

    // 2. compute the weights of edges
    _calculate_edge_weight();

    // 3. initialize the map
    using M = CSHMMesh;
    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
        M::CVertex *pV = *viter;
        pV->u() = pV->normal(); // 高斯 Map
    }
}

void MeshLib::CSphericalHarmonicMap::_calculate_normal() {
    using M = CSHMMesh;
    for (CSHMMesh::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter) {
        CSHMMesh::CFace *pF = *fiter;
        CPoint p[3];
        CSHMMesh::CHalfEdge *he = m_pMesh->faceHalfedge(pF);

        for (int k = 0; k < 3; k++) {
            CSHMMesh::CVertex *pv = m_pMesh->halfedgeTarget(he);
            p[k] = pv->point();
            he = m_pMesh->halfedgeNext(he);
        }

        CPoint fn = (p[1] - p[0]) ^(p[2] - p[0]);
        pF->area() = fn.norm() / 2.0;
        pF->normal() = fn / fn.norm();
    }

    for (CSHMMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
        CSHMMesh::CVertex *v = *viter;
        CPoint n(0, 0, 0);
        double area = 0;
        for (CSHMMesh::VertexFaceIterator vfiter(v); !vfiter.end(); ++vfiter) {
            CSHMMesh::CFace *pF = *vfiter;
            n += pF->normal() * pF->area();
            area += pF->area();
        }
        v->area() = area / 3.0;
        n = n / n.norm();
        v->normal() = n;
    }
}

double MeshLib::CSphericalHarmonicMap::step_one(int steps, double step_length) {
    if (!m_pMesh) {
        std::cerr << "Should set mesh first!" << std::endl;
        return DBL_MAX;
    }

    using M = CSHMMesh;

    for (int i = 0; i < steps; ++i) {
        for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
            M::CVertex *pV = *viter;
            // TODO insert your code here
            // 1. compute vertex laplacian
            CPoint laplacian(0, 0, 0);
            double totalWeight = 0.0;
            for (M::VertexVertexIterator vviter(pV); !vviter.end(); ++vviter) {
                M::CVertex *pW = *vviter;
                M::CEdge *e = m_pMesh->vertexEdge(pV, pW);
                double w = e->weight();
                laplacian = laplacian + (pW->u() - pV->u()) * w;
                totalWeight += w;
            }
            laplacian = laplacian / totalWeight;
            // 2. get the normal component
            CPoint normal = pV->u() * (laplacian * pV->u());
            // 3. get the tangent_component
            CPoint tangent = laplacian - normal;

            pV->deltaU = tangent;// 先保存 deltaU，事后计算步长

            // 4. update u
//            CPoint u = pV->u() + tangent * step_length;
            // 5. normalize the vertex->u() to the unit sphere
//            pV->u() = u / u.norm();
            // TODO insert your code here
        }
        // TODO insert your code here
        // 使用线性搜索方法重新计算步长，能量导数为 0

        // 计算步长
        double upper = 0.0; // 分子
        double lower = 0.0; // 分母
        for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter) {
            M::CEdge *pE = *eiter;
            M::CVertex *pV = m_pMesh->edgeVertex1(pE);
            M::CVertex *pW = m_pMesh->edgeVertex2(pE);
            CPoint d = pV->u() - pW->u();
            CPoint dd = pV->deltaU - pW->deltaU;
            upper -= pE->weight() * (d * dd);
            lower += pE->weight() * (dd * dd);
        }
        double lambda = upper / lower;
        printf("lambda：%f\n", lambda);

        for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
            M::CVertex *pV = *viter;
            // 4. update u
            pV->u() += pV->deltaU * lambda;
            // 5. normalize the vertex->u() to the unit sphere
            pV->u() /= pV->u().norm();
        }
        // TODO insert your code here
    }
    // TODO insert your code here
    // 6. normalize the mapping, such that mass center is at the origin
    this->_normalize();
    // TODO insert your code here
    // 7. compute the harmonic energy
    double E = _calculate_harmonic_energy();
    std::cout << "After " << steps << " steps, harmonic energy is " << E << std::endl;
    return E;
}

void MeshLib::CSphericalHarmonicMap::map(double step_length, double epsilon) {
    if (!m_pMesh) {
        std::cerr << "Should set mesh first!" << std::endl;
        return;
    }

    double E_prev = _calculate_harmonic_energy();
    double E = 0;
    while (true) {
        E = step_one(10, step_length);

        if (std::fabs(E - E_prev) < epsilon) break; // 能量不在减少时终止迭代
        E_prev = E;
    };
}

void MeshLib::CSphericalHarmonicMap::_calculate_edge_weight() {
    using M = CSHMMesh;

    // 1. compute edge length
    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter) {
        M::CEdge *pE = *eiter;
        M::CVertex *v1 = m_pMesh->edgeVertex1(pE);
        M::CVertex *v2 = m_pMesh->edgeVertex2(pE);
        pE->length() = (v1->point() - v2->point()).norm();
    }

    // 2. compute corner angle
    for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter) {
        M::CFace *pF = *fiter;
        M::CHalfEdge *pH[3];
        pH[0] = m_pMesh->faceHalfedge(pF);
        for (int i = 0; i < 3; i++) {
            pH[(i + 1) % 3] = m_pMesh->halfedgeNext(pH[i]);
        }

        double len[3];
        for (int i = 0; i < 3; i++) {
            len[i] = m_pMesh->halfedgeEdge(pH[i])->length();
        }

        for (int i = 0; i < 3; i++) {
            double a = len[(i + 1) % 3], b = len[(i + 2) % 3], c = len[i];
            pH[(i + 1) % 3]->angle() = _inverse_cosine_law(a, b, c);
        }
    }

    // 3. compute edge weight
    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter) {
        M::CEdge *pE = *eiter;

        if (!pE->boundary()) {
            double theta[2];
            theta[0] = m_pMesh->halfedgeNext(m_pMesh->edgeHalfedge(pE, 0))->angle();
            theta[1] = m_pMesh->halfedgeNext(m_pMesh->edgeHalfedge(pE, 1))->angle();
            pE->weight() = 0.5 * (std::cos(theta[0]) / std::sin(theta[0]) + std::cos(theta[1]) / std::sin(theta[1]));
        } else {
            double theta = m_pMesh->halfedgeNext(m_pMesh->edgeHalfedge(pE, 0))->angle();
            pE->weight() = 0.5 * std::cos(theta) / std::sin(theta);
        }
    }
}

double MeshLib::CSphericalHarmonicMap::_calculate_harmonic_energy() {
    using M = CSHMMesh;

    double energy = 0;
    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter) {
        M::CEdge *pE = *eiter;
        M::CVertex *pV = m_pMesh->edgeVertex1(pE);
        M::CVertex *pW = m_pMesh->edgeVertex2(pE);
        CPoint d = pV->u() - pW->u();
        energy += pE->weight() * (d * d);// weight of edge * d * d
    }
    return energy;
}

double MeshLib::CSphericalHarmonicMap::_inverse_cosine_law(double a, double b, double c) {
    double cs = (a * a + b * b - c * c) / (2.0 * a * b);
    assert(cs <= 1.0 && cs >= -1.0);
    return std::acos(cs);
}

void MeshLib::CSphericalHarmonicMap::_normalize() {
    using M = CSHMMesh;

    CPoint center(0, 0, 0);
    double area = 0;

    // TODO insert your code here
    //move the mass center of vertex->u() to the origin
    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
        M::CVertex *pV = *viter;
        center = center + pV->u() * pV->area();
        area += pV->area();
    }
    center = center / area; // 计算重心

    for (CSHMMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
        CSHMMesh::CVertex *pV = *viter;
        CPoint u = pV->u();
        u = u - center;
        u = u / u.norm();
        pV->u() = u; // 减去重心，使重心为（0，0，0）
    }
    // TODO insert your code here
}