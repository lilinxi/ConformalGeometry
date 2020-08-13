#ifndef _SPHERICAL_HARMONIC_MAP_H_
#define _SPHERICAL_HARMONIC_MAP_H_

#include "SphericalHarmonicMapMesh.h"

namespace MeshLib {

/*! \brief CSphericalHarmonicMap class
 *
 *   Spherical harmonic map algorithm that maps a topological sphere to 
 *   the unit sphere with minimal harmonic energy.
 */
    class CSphericalHarmonicMap {
    public:
        /*!
         *  CSphericalHarmonicMap constructor
         */
        CSphericalHarmonicMap() : m_pMesh(NULL) {};

        /*!
         *  Set mesh and initialization
         *  \param pMesh input mesh
         */
        void set_mesh(CSHMMesh *pMesh);

        /*!
         *  Take next one step
         *  \param step step length
         *  \return harmonic energy
         */
        double step_one(int steps = 1, double step_length = 5e-1);

        /*!
         *  The spherical harmonic map
         *  \param step step length
         *  \param epsilon error threshold
         */
        void map(double step_length = 5e-1, double epsilon = 1e-3);

        /*
        *   normalize the mapping
        */
        void _normalize();

    protected:
        /*!
         *  Compute vertex normal
         */
        void _calculate_normal();

        /*!
         *  Compute edge weight
         */
        void _calculate_edge_weight();

        /*!
         *  Compute harmonic energy
         *  \return harmonic energy
         */
        double _calculate_harmonic_energy();

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
        CSHMMesh *m_pMesh;
    };
}
#endif // !_SPHERICAL_HARMONIC_MAP_H_
