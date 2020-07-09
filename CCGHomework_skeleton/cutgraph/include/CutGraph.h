#ifndef _CUT_GRAPH_H_
#define _CUT_GRAPH_H_

#include "CutGraphMesh.h"

namespace MeshLib
{
/*! \brief CCutGraph class
 *
 *   Compute the spanning tree of the dual mesh,
 *   the edges whose duals are not on the tree form the cut locus,
 *   label the cut locus as the sharp edges.
 */
class CCutGraph
{
  public:
    /*! 
     *  CCutGraph constructor
     *  \param pMesh input closed mesh
     */
    CCutGraph(CCutGraphMesh* pMesh) { m_pMesh = pMesh; };

    /*! 
     * Compute the cut graph.
     */
    void cut_graph();

  protected:
    /*! 
     *  Input closed mesh.
     */
    CCutGraphMesh* m_pMesh;

    /*! 
     *  Compute the spanning tree of the dual mesh.
     */
    void _dual_spanning_tree();

    /*!
     * Prune the branches which attached to valence-1 nodes.
     */
    void _prune();
};
} // namespace MeshLib
#endif // !_CUT_GRAPH_H_
