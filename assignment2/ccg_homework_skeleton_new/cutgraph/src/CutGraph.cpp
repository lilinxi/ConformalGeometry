#include <queue>
#include "CutGraph.h"

void MeshLib::CCutGraph::cut_graph()
{
    _dual_spanning_tree();

    // The cut graph contains all edges which their duals are not
    // in the spanning tree.
    for (CCutGraphMesh::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter)
    {
        CCutGraphEdge* pE = *eiter;
        pE->sharp() = !pE->sharp();
    }

    _prune();
}

/*----------------------------------------------------------------------------

Modify the method CCutGraph::_dual_spanning_tree()

------------------------------------------------------------------------------*/

void MeshLib::CCutGraph::_dual_spanning_tree()
{
    // Mark all sharp flags false
    for (CCutGraphMesh::MeshEdgeIterator eiter(m_pMesh); !eiter.end(); ++eiter)
    {
        CCutGraphEdge* pE = *eiter;
        //insert your code here
    }

    // Mark all touched flags false, and select a face
    CCutGraphFace* pHeadFace = NULL;
    for (CCutGraphMesh::MeshFaceIterator fiter(m_pMesh); !fiter.end(); ++fiter)
    {
        CCutGraphFace* pF = *fiter;
        //insert your code here
        pHeadFace = pF;
    }
    
    // 1. Construct the dual mesh conceptually.

    // 2. Generate a minimal spanning tree of the vertices in the dual mesh.
    std::queue<CCutGraphFace*> fQueue;
    fQueue.push(pHeadFace);
    pHeadFace->touched() = true;
    while (!fQueue.empty())
    {
        CCutGraphFace* pF = fQueue.front();
        fQueue.pop();

        for (CCutGraphMesh::FaceHalfedgeIterator fhiter(pF); !fhiter.end(); ++fhiter)
        {
            CCutGraphHalfEdge* pH = *fhiter;
            CCutGraphHalfEdge* pSymH = m_pMesh->halfedgeSym(pH);
            if (pSymH != NULL)
            {
                CCutGraphFace* pSymFace = m_pMesh->halfedgeFace(pSymH);
                if (!pSymFace->touched())
                {
                    //insert your code here
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------

Modify the method _CCutGraph::_prune()

------------------------------------------------------------------------------*/

void MeshLib::CCutGraph::_prune() 
{
    // A queue used to store valence-1 vertices
    std::queue<CCutGraphVertex*> vQueue;

    // 1. Compute the valence of each vertex, and record all valence-1 vertices.
    for (CCutGraphMesh::MeshVertexIterator viter(m_pMesh); !viter.end(); ++viter)
    {
        CCutGraphVertex* pV = *viter;
        pV->valence() = 0;

        for (CCutGraphMesh::VertexEdgeIterator veiter(pV); !veiter.end(); ++veiter)
        {
            CCutGraphEdge* pE = *veiter;
            //insert your code here
        }

        if (pV->valence() == 1)
            vQueue.push(pV);
    }

    // 2. Remove the segments which attached to valence-1 vertices.
    while (!vQueue.empty())
    {
        CCutGraphVertex* pV = vQueue.front();
        vQueue.pop();

        for (CCutGraphMesh::VertexEdgeIterator veiter(pV); !veiter.end(); ++veiter)
        {
            CCutGraphEdge* pE = *veiter;
            CCutGraphVertex* pW = m_pMesh->edgeVertex1(pE) == pV ? 
                m_pMesh->edgeVertex2(pE) : m_pMesh->edgeVertex1(pE);
            if (pE->sharp())
            {
                //insert your code here
                break;
            }
        }
    }
}
