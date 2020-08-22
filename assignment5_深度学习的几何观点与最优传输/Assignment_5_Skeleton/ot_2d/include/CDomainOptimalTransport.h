/*!  \file   DomainOptimalTransport.h
 *   \brief  OptimalTransport
 *   \author David Gu
 *   \date   documented on 05/18/2020
 *
 *   Base Class for Planar Optimal Transportation, Semi-Discrete Algorithm
 */

#ifndef _DOMAIN_OPTIMAL_TRANSPORT_H_
#define _DOMAIN_OPTIMAL_TRANSPORT_H_

#include "Detri2Mesh.h"
#include "OT.h"

namespace MeshLib
{

/*-------------------------------------------------------------------------------------------------------------------------------------

         Optimal Transport Map to a Convex Planar Domain

--------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CDomainOptimalTransport class
 *
 *  class for planar polygonal domain optimal transport
 *
 *  WDT  - Detri2::Weighted Delaunay Triangulation
 *  mesh - MeshLib::COMTMesh
 *
 */
class CDomainOptimalTransport : public CBaseOT, public CDetri2Mesh
{
  public:
    CDomainOptimalTransport(COMTMesh* pMesh);

    ~CDomainOptimalTransport();

    /*! initialize the potential function as quadratic
     */
    void _initialize();
    /*! gradient descend
     */
    void __gradient_descend(COMTMesh* pInput, COMTMesh*& pOutput);

    /*! Newton's method
     */
    void __newton(COMTMesh* pInput, COMTMesh*& pOutput);

    /* data members */
  public:
    // Weighted Delaunay Mesh
    COMTMesh*& pWeightedDT() { return m_pWDT; };

    // base mesh
    COMTMesh*& PMesh() { return m_pMesh; };

  protected:
    // triangle mesh for the Weighted Delaunay Triangulation
    COMTMesh* m_pWDT = NULL;

    // triangulation for the background domain
    detri2::Triangulation* m_domainTr = NULL;

    // triangulation for the Weighted Delaunay
    detri2::Triangulation* m_outputTr = NULL;
};

} // namespace MeshLib
#endif //! DOMAIN_OPTIMAL_TRANSPORT_H_
