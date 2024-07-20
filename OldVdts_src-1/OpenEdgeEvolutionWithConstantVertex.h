#if !defined(AFX_OpenEdgeEvolutionWithConstantVertex_H_INCLUDED_)
#define AFX_OpenEdgeEvolutionWithConstantVertex_H_INCLUDED_
#include "OpenEdgeEvolution.h"
#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "SimDef.h"
#include "MESH.h"
#include "RNG.h"
#include "Energy.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 
 This is a class to change the size of an open edge surface. This only creats triangles from edges that connected.
 
 The four main the the most important functions are
 
    1) links* CreateALink(vertex *);
    2) bool KillALink(links *);
    3) triangle* CloseATriangleHole(vertex *v1);
    4) bool KillATriangle(links *l1);
 
 
 */

class State;
class OpenEdgeEvolutionWithConstantVertex: public OpenEdgeEvolution { // to use for
public:
    OpenEdgeEvolutionWithConstantVertex();
    OpenEdgeEvolutionWithConstantVertex(int rate, State *pState);
    ~OpenEdgeEvolutionWithConstantVertex();
    inline int GetRate() {return m_Rate;}

public:
    void Initialize();
    void MC_Move(RNG* rng, double lmin, double lmax, double maxangle);
    
private:
    std::vector<triangle>       m_GhostT; // Some trinagles for initial storing
    std::vector<links>          m_GhostL;
    std::vector<triangle*>      m_pGhostT; // Some trinagles ....
    std::vector<links*>         m_pGhostL;  // some edges for  ...
    int m_Rate;
    MESH* m_pMESH;
    Vec3D * m_pBox;
    State *m_pState;
    Energy *m_pEnergyCalculator;
    double m_Beta;

    void RemoveFromLinkList(links* z, std::vector<links*> &vect);
    void RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect);
    void RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect);
    void AddtoLinkList(links* z, std::vector<links*> &vect);
    void AddtoVertexList(vertex* z, std::vector<vertex*> &vect);
    void AddtoTriangleList(triangle* z, std::vector<triangle*> &vect);
    bool Linkisvalid(vertex *, double lmin, double lmax, double maxangle);
    double  SystemEnergy();


    // the main hard part of the code. 4 interesting function
private:  // this functions could be in princeple public, but no need for now
    
    links* CreateALink(vertex *);    // creates a link between two connected edge at the edge
    bool KillALink(links *);         // kills an edge that is on the edge. not an edge on the surface
    triangle* CloseATriangleHole(vertex *v1);  //if an edge contains only three links, it will close it
    bool KillATriangle(links *l1);      // kills a triangle anywhere

    
};
#endif
