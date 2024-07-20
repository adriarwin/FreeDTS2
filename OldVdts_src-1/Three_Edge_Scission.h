#if !defined(THREE_EDGE_SCISSION_H_INCLUDED)
#define THREE_EDGE_SCISSION_H_INCLUDED

#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "SimDef.h"
#include "MESH.h"
#include "RNG.h"
#include "Energy.h"
#include "DynamicTopology.h"
struct pot_triangle {    // data structure for a potential triangle
        int id;
        int cid; // connected id
        vertex *pv1;
        vertex *pv2;
        vertex *pv3;
        links *pl1;
        links *pl2;
        links *pl3;
    };
    struct pair_pot_triangle {    // data structure for a potential triangle
        bool state;
        int id;
        pot_triangle PT1;
        pot_triangle PT2;
        std::vector <links *> ConnectingLinks;       // this does not include the mirror links
        std::vector <triangle *> ConnectingTriangles;

    };
class State;
class Three_Edge_Scission : public DynamicTopology { // to use for polymorphism


    
public:
    Three_Edge_Scission();
    Three_Edge_Scission(int period, State *pState);
    ~Three_Edge_Scission();
    void initialize();
    bool MCMove(int step, double * TotalEnergy, RNG *rng, GenerateCNTCells *pGenCNT );

private:
    std::vector<triangle>       m_GhostT; // Some trinagles for initial storing
    std::vector<links>          m_GhostL;
    MESH* m_pMESH;
    Vec3D * m_pBox;
    State *m_pState;
    Energy *m_pEnergyCalculator;
    void RemoveFromLinkList(links* z, std::vector<links*> &vect);
    void RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect);
    void AddtoLinkList(links* z, std::vector<links*> &vect);
    void AddtoTriangleList(triangle* z, std::vector<triangle*> &vect);
    template<typename T>
    bool AddtoVectorCarefully(T* z, std::vector<T*> &vect);
    
    bool Anglevalid4Vhole(vertex *v1, double minangle);
    bool CorrectOrientation(pot_triangle &p1,pot_triangle &p2);
    pair_pot_triangle connected_2pot_triangles(pot_triangle potT1, pot_triangle potT2);
    std::vector<pair_pot_triangle> FindPotentialTriangles(MESH* mesh);
    std::vector <triangle *> DoAScission(pair_pot_triangle &pair);
    bool ReverseAScission(pair_pot_triangle &pair, triangle *t1, triangle *t2);   // this is the exact reverse action of DoAScission; different from DoAFussion

    bool DoAFussion(pair_pot_triangle pair);
    triangle * CreateATriangleFromAPotentialTriangle(pot_triangle &p1);
    bool MCScissionMove(int step, double * TotalEnergy, RNG *rng);
    bool MCFussionMove(int step, double * TotalEnergy, RNG *rng);
    double m_Beta;
    int m_Period;
    double  UpdateEnergy();

public:
    std::vector<triangle*>      m_pGhostT; // Some trinagles ....
    std::vector<links*>         m_pGhostL;  // some edges for  ...


template<typename T>
void KeepOneOccurrence(std::vector<T*> &vect);


    
    
};
#endif
