#if !defined(AFX_OpenEdgeEvolutionWithConstantVertex_B_H_INCLUDED_)
#define AFX_OpenEdgeEvolutionWithConstantVertex_B_H_INCLUDED_
#include "OpenEdgeEvolution.h"
#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "SimDef.h"
#include "MESH.h"
#include "RNG.h"
#include "Energy.h"

class State;
class OpenEdgeEvolutionWithConstantVertex_B: public OpenEdgeEvolution { // to use for
public:
    OpenEdgeEvolutionWithConstantVertex_B();
    OpenEdgeEvolutionWithConstantVertex_B(int rate, State *pState);
    ~OpenEdgeEvolutionWithConstantVertex_B();
    inline int GetRate() {return m_Rate;}

    
public:
    double m_WholeSize;
    
public:
    std::vector<triangle*>      m_pGhostT; // Some trinagles ....
    std::vector<links*>         m_pGhostL;  // some edges for  ...
    void Initialize();
    void MC_Move(RNG* rng, double lmin, double lmax, double maxangle);
    
private:
    int m_Rate;
    std::vector<triangle>       m_GhostT; // Some trinagles for initial storing
    std::vector<links>          m_GhostL;
    MESH* m_pMESH;
    Vec3D * m_pBox;
    State *m_pState;
    Energy *m_pEnergyCalculator;
    void RemoveFromLinkList(links* z, std::vector<links*> &vect);
    void RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect);
    void RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect);
    
    void AddtoLinkList(links* z, std::vector<links*> &vect);
    void AddtoVertexList(vertex* z, std::vector<vertex*> &vect);
    void AddtoTriangleList(triangle* z, std::vector<triangle*> &vect);
    double m_Beta;
    bool Anglevalid4Vhole(vertex *v1, double minangle);



    
private:
    links* CreateALink(vertex *);
    void KillALink(links *);
    void KillALinkOnSurf(links *);
    bool Linkisvalid(vertex *, double lmin, double lmax, double maxangle);
    double  SystemEnergy();



    
    
};
#endif
