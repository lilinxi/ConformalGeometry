/*!  \file OptimalTransport.h
 *   \brief OptimalTransport
 *   \author David Gu
 *   \date   documented on 05/01/2020
 *
 *   Base Class for Optimal Transportation, Semi-Discrete Algorithm
 */

#ifndef _BASE_OT_H_
#define _BASE_OT_H_

#include <map>
#include <vector>
#include <Eigen/Eigen>

#include "OTMesh.h"

namespace MeshLib
{

/*-------------------------------------------------------------------------------------------------------------------------------------

        Optimal Transport

--------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CBaseOT class
 *
 *  base class for optimal transport
 *
 */
class CBaseOT
{
  public:
    CBaseOT(COMTMesh* pMesh) { m_pMesh = pMesh; };
    ~CBaseOT(){};

  public:
    /*! normalize the uv coordinates to be the unit disk */
    virtual void _normalize_uv(COMTMesh* pMesh);
    /*! compute the maximal relative error */
    virtual void _compute_error(COMTMesh* pMesh);
    /*! set target are */
    virtual void _set_target_measure(COMTMesh*& pMesh, double total_target_area = PI);
    /*! initialize the mapping, idendity*/
    virtual void _initialize(COMTMesh* pChull){};
    /*! copy different attributes from input mesh to the output mesh */
    void _copy_mesh(COMTMesh* pInput, COMTMesh* pOutput);
    
  protected:
    COMTMesh* m_pMesh;
    /*! compute the hessian matrix, from edge length and dual edge length, and vertex index */
    void __compute_hessian_matrix(COMTMesh& mesh, Eigen::SparseMatrix<double>& hessian);
    /*! compute the vertex update direction using Newton's method */
    void __update_direction(COMTMesh* m_pMesh);
    /*! use Eigen to solve a linear system */
    bool __solve(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, Eigen::VectorXd& result);
    /*! reindex all the vertices, starting from zero */
    void index(COMTMesh* pMesh);

};

}; // namespace MeshLib

#endif //! BASE_OPTIMAL_TRANSPORT_H_
