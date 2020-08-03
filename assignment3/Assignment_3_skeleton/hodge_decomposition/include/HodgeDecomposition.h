#ifndef _HODGE_DECOMPOSITION_H_
#define _HODGE_DECOMPOSITION_H_

#include "HodgeDecompositionMesh.h"

namespace MeshLib
{

/*! \brief CHodgeDecomposition class
 *
 *   Hodge Decomposition
 */
class CHodgeDecomposition
{
  public:
    /*!
     *  CHodgeDecomposition constructor
     */
    CHodgeDecomposition();

    void set_mesh(CHodgeDecompositionMesh* pMesh);

    /*!
    *   compute a random holomorphic form
    */
    void random_harmonic_form();
    /*!
     * compute exact harmonic form
     */
    void exact_harmonic_form( int bnd_id );
    
    void integration(CHodgeDecompositionMesh* pForm, CHodgeDecompositionMesh* pDomain);
  protected:

    /*! 
     * Compute angle using cosine law
     * \param a first edge length
     * \param b second edge length 
     * \param c third edge length
     * \return angle between 1st and 2nd edges
     */
    double _inverse_cosine_law(double a, double b, double c);

    /*! exterior differentiation operator */
    void _d(int dimension);

    /*! delta operator */
    void _delta(int dimension);

    /*!
     *  Compute edge weight
     */
    void _calculate_edge_weight(bool using_geometry);
    /*! find the exact form */
    void _compute_exact_form();

    /*! find the coexact form */
    void _compute_coexact_form();

    /*! compute the harmonic form */
    void _normalize();

    /*! random 1-form */
    void _random_form();

    /*! verify if the 1-form is closed */
    void _test_closedness();

    /*! verify if the 1-form is closed */
    void _test_coclosedness();

    /*! remove df from \omega */
    void _remove_exact_form();
    /*! remove \delta \eta from \omega */
    void _remove_coexact_form();

    /*! exact harmonic form */
    void _exact_harmonic_form();

    /*! set boundary conditions for exact harmonic forms */
    void _set_boundary_condition(int boundary_id);

  protected:
    /*!
     * The input surface mesh
     */
    CHodgeDecompositionMesh* m_pMesh;
};
}
#endif // !_HODGE_DECOMPOSITION_H_
