//#include <float.h>
//#include <math.h>
//#include <cmath>
//#include "SphericalHarmonicMap.h"
//
//void MeshLib::CSphericalHarmonicMap::set_mesh(CSHMMesh *pMesh) {
//    m_pMesh = pMesh;
//
//    // 1. compute vertex normal
//    _calculate_normal();
//
//    // 2. compute the weights of edges
//    _calculate_edge_weight();
//
//    // 3. initialize the map
//    using M = CSHMMesh;
//    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
//        M::CVertex *pV = *viter;
//        pV->u() = pV->normal();
//    }
//}
//
//void MeshLib::CSphericalHarmonicMap::_calculate_normal() {
//    using M = CSHMMesh;
//    for (CSHMMesh::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter) {
//        CSHMMesh::CFace *pF = *fiter;
//        CPoint p[3];
//        CSHMMesh::CHalfEdge *he = m_pMesh->faceHalfedge(pF);
//
//        for (int k = 0; k < 3; k++) {
//            CSHMMesh::CVertex *pv = m_pMesh->halfedgeTarget(he);
//            p[k] = pv->point();
//            he = m_pMesh->halfedgeNext(he);
//        }
//
//        CPoint fn = (p[1] - p[0]) ^(p[2] - p[0]);
//        pF->area() = fn.norm() / 2.0;
//        pF->normal() = fn / fn.norm();
//    }
//
//    for (CSHMMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
//        CSHMMesh::CVertex *v = *viter;
//        CPoint n(0, 0, 0);
//        double area = 0;
//        for (CSHMMesh::VertexFaceIterator vfiter(v); !vfiter.end(); ++vfiter) {
//            CSHMMesh::CFace *pF = *vfiter;
//            n += pF->normal() * pF->area();
//            area += pF->area();
//        }
//        v->area() = area / 3.0;
//        n = n / n.norm();
//        v->normal() = n;
//    }
//}
//
//double MeshLib::CSphericalHarmonicMap::step_one(int steps, double step_length) {
//    if (!m_pMesh) {
//        std::cerr << "Should set mesh first!" << std::endl;
//        return DBL_MAX;
//    }
//
//    // static int time = 0;
//    // static int div_time = 0;
//    // time++;
//    // if(time > 10) {
//    //     div_time++;
//    //     for(int i = 0; i < time; i++)
//    //         step_length /= 10.0;
//    //     time = 0;
//    // }
//
//    using M = CSHMMesh;
//    for (int i = 0; i < steps; ++i) {
//        int cnt = 0;
//        for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
//            M::CVertex *pV = *viter;
//            //insert your code here
//
//            // 1. compute vertex laplacian
//            CPoint delta_phi(0, 0, 0);
//            double sw = 0.0;
//            for (M::VertexVertexIterator vviter(pV); !vviter.end(); vviter++) {
//                M::CVertex *pW = *vviter;
//                M::CEdge *pE = m_pMesh->vertexEdge(pV, pW);
//                delta_phi += (pW->u() - pV->u()) * pE->weight();
//                sw += pE->weight();
//            }
//            delta_phi /= sw;
//
//            // 2. get the noraml component
//            CPoint n_component = pV->u() * (pV->u() * delta_phi);
//            cnt++;
//            // 3. get the tangent_component
//            delta_phi -= n_component;
//            pV->delta_u = delta_phi;
//        }
//
//        // Eigen::VectorXd direction(3 * cnt);
//        // int icnt = 0;
//        // for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
//        // {
//        //     M::CVertex* pV = *viter;
//        //     direction[icnt + 0] = pV->delta_u[0];
//        //     direction[icnt + 1] = pV->delta_u[1];
//        //     direction[icnt + 2] = pV->delta_u[2];
//        //     icnt += 3;
//        // }
//
//        // for(size_t k = 0; k < searched_direction.size(); ++k) {
//        //     Eigen::VectorXd &dir = searched_direction[k];
//        //     assert(dir.size() == direction.size());
//        //     double proj = direction.dot(dir);
//        //     direction = direction - proj * dir;
//        // }
//
//        // icnt = 0;
//        // for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
//        // {
//        //     M::CVertex* pV = *viter;
//        //     pV->delta_u[0] = direction[icnt + 0];
//        //     pV->delta_u[1] = direction[icnt + 1];
//        //     pV->delta_u[2] = direction[icnt + 2];
//        //     icnt += 3;
//        // }
//
//        // direction.normalize();
//        // searched_direction.push_back(direction);
//
//        double den = 0.0;
//        double num = 0.0;
//
//        for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter) {
//            M::CEdge *pE = *eiter;
//            M::CVertex *pV = m_pMesh->edgeVertex1(pE);
//            M::CVertex *pW = m_pMesh->edgeVertex2(pE);
//            CPoint d = pV->u() - pW->u();
//            CPoint dd = pV->delta_u - pW->delta_u;
//            //energy += pE->weight() * (d * d);
//            den += pE->weight() * (dd * dd);
//            num -= pE->weight() * (d * dd);
//        }
//
//        double lambda = num / den;
//        for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
//            M::CVertex *pV = *viter;
//            pV->u() += pV->delta_u * lambda;
//            // 5. normalize the vertex->u() to the unit sphere
//            pV->u() /= pV->u().norm();
//        }
//
//    }
//
//    // 6. normalize the mapping, such that mass center is at the origin
//    _normalize();
//    // 7. compute the harmonic energy
//    double E = _calculate_harmonic_energy();
//    std::cout << "After " << steps << " steps, harmonic energy is " << E << std::endl;
//    return E;
//}
//
//void MeshLib::CSphericalHarmonicMap::map(double step_length, double epsilon) {
//    using M = CSHMMesh;
//    if (!m_pMesh) {
//        std::cerr << "Should set mesh first!" << std::endl;
//        return;
//    }
//
//    double E_prev = _calculate_harmonic_energy();
//    double E = 0;
//    while (true) {
//        E = step_one(10, step_length);
//
//        if (std::fabs(E - E_prev) < epsilon) break;
//        E_prev = E;
//    };
//
////    CPoint north(0, 0, 1);
////    CPoint south(0, 0, -1);
////    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
////        M::CVertex *pV = *viter;
//    //insert your code here
//
////        CPoint pt = pV->u();
////        if (pt[2] < 1.0) {
////            pV->uv[0] = pt[0] / (1 - pt[2]);
////            pV->uv[1] = pt[1] / (1 - pt[2]);
////        }
////    }
//}
//
//void MeshLib::CSphericalHarmonicMap::_calculate_edge_weight() {
//    using M = CSHMMesh;
//
//    // 1. compute edge length
//    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter) {
//        M::CEdge *pE = *eiter;
//        M::CVertex *v1 = m_pMesh->edgeVertex1(pE);
//        M::CVertex *v2 = m_pMesh->edgeVertex2(pE);
//        pE->length() = (v1->point() - v2->point()).norm();
//    }
//
//    // 2. compute corner angle
//    for (M::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter) {
//        M::CFace *pF = *fiter;
//        M::CHalfEdge *pH[3];
//        pH[0] = m_pMesh->faceHalfedge(pF);
//        for (int i = 0; i < 3; i++) {
//            pH[(i + 1) % 3] = m_pMesh->halfedgeNext(pH[i]);
//        }
//
//        double len[3];
//        for (int i = 0; i < 3; i++) {
//            len[i] = m_pMesh->halfedgeEdge(pH[i])->length();
//        }
//
//        for (int i = 0; i < 3; i++) {
//            double a = len[(i + 1) % 3], b = len[(i + 2) % 3], c = len[i];
//            pH[(i + 1) % 3]->angle() = _inverse_cosine_law(a, b, c);
//        }
//    }
//
//    // 3. compute edge weight
//    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter) {
//        M::CEdge *pE = *eiter;
//
//        if (!pE->boundary()) {
//            double theta[2];
//            theta[0] = m_pMesh->halfedgeNext(m_pMesh->edgeHalfedge(pE, 0))->angle();
//            theta[1] = m_pMesh->halfedgeNext(m_pMesh->edgeHalfedge(pE, 1))->angle();
//            pE->weight() = 0.5 * (std::cos(theta[0]) / std::sin(theta[0]) + std::cos(theta[1]) / std::sin(theta[1]));
//        } else {
//            double theta = m_pMesh->halfedgeNext(m_pMesh->edgeHalfedge(pE, 0))->angle();
//            pE->weight() = 0.5 * std::cos(theta) / std::sin(theta);
//        }
//    }
//}
//
//double MeshLib::CSphericalHarmonicMap::_calculate_harmonic_energy() {
//    using M = CSHMMesh;
//
//    double energy = 0;
//    for (M::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter) {
//        M::CEdge *pE = *eiter;
//        M::CVertex *pV = m_pMesh->edgeVertex1(pE);
//        M::CVertex *pW = m_pMesh->edgeVertex2(pE);
//        CPoint d = pV->u() - pW->u();
//        energy += pE->weight() * (d * d);
//    }
//    return energy;
//}
//
//double MeshLib::CSphericalHarmonicMap::_inverse_cosine_law(double a, double b, double c) {
//    double cs = (a * a + b * b - c * c) / (2.0 * a * b);
//    assert(cs <= 1.0 && cs >= -1.0);
//    return std::acos(cs);
//}
//
//void MeshLib::CSphericalHarmonicMap::_normalize() {
//    using M = CSHMMesh;
//
//    CPoint center(0, 0, 0);
//    double area = 0;
//    //insert your code here
//    //move the mass center of vertex->u() to the origin
//
//    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
//        M::CVertex *pV = *viter;
//        center += pV->u() * pV->area();
//        area += pV->area();
//    }
//
//    center /= area;
//    for (M::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter) {
//        M::CVertex *pV = *viter;
//        pV->u() -= center;
//        pV->u() /= pV->u().norm();
//    }
//}