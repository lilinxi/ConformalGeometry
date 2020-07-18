#ifndef _HARMONIC_MAP_H_
#define _HARMONIC_MAP_H_

#include "HarmonicMapMesh.h"

namespace MeshLib
{

/*! \brief CHarmonicMap class
 *
 *   Harmonic map algorithm that maps a topological disk to 
 *   the unit disk with minimal harmonic energy.
 */
class CHarmonicMap
{
  public:
    /*!
     *  CHarmonicMap constructor
     */
    CHarmonicMap() : m_pMesh(NULL) {};

    /*!
     *  Set mesh and initialization 
     *  \param pMesh input mesh
     */
    void set_mesh(CHarmonicMapMesh* pMesh);

    /*! 
     *  Take next one step
     *  \return maximal error of current step 
     */
    double step_one();

    /*!
     *  Iterative meshod
     *  \param epsilon error threshold
     */
    void iterative_map(double epsilon = 1e-5);
    
    /*!
     *  Directly solving the harmonic map
     */
    void map();

  protected:
    /*!
     *  Compute edge weight
     */
    void _calculate_edge_weight();

    /*!	
     *  Fix the boundary vertices to the unit circle
     *  using arc length parameter
     */
    void _set_boundary();

    /*! 
     * Compute angle using cosine law
     * \param a first edge length
     * \param b second edge length 
     * \param c third edge length
     * \return angle between 1st and 2nd edges
     */
    double _inverse_cosine_law(double a, double b, double c);
  protected:
    /*!
     * The input surface mesh
     */
    CHarmonicMapMesh* m_pMesh;
};
}
#endif // !_HARMONIC_MAP_H_
