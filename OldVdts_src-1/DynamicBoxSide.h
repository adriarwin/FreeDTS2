#if !defined(AFX_DynamicBoxSide_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_)
#define AFX_DynamicBoxSide_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_


#include "SimDef.h"
#include "links.h"
#include "Vec3D.h"
#include "Energy.h"
#include "GenerateCNTCells.h"
#include "SpringPotentialBetweenTwoGroups.h"
#include "CNTCell.h"
#include "triangle.h"
#include "MESH.h"
#include "DynamicBox.h"


class State;
class DynamicBoxSide : public DynamicBox {
public:
    DynamicBoxSide();
    DynamicBoxSide(int tau, double f,State *st);
	~DynamicBoxSide();

    
    bool GetCNTCondition();
    int GetTau();
    void initialize();
    bool MCMoveBoxChange(double dx, double * TotalEnergy, double temp, int step, GenerateCNTCells *pGenCNT );

    
    
private:
    //=== old functions
    void CheckCNTSize();
    double DistanceSquardBetweenTwoVertices(vertex *,vertex *,Vec3D& );
    bool CheckFaceAngle();
    bool CheckMaxLinkLength();
    void PerformMove();
    void RejectMove();
    void AcceptMove();
    bool CheckFaceAngle(links * l);
private:

    int m_Tau;
    double m_F0;
    double m_dx;
    bool m_UpdateCNT;
    bool m_Move;
    double m_Beta; // 1/k_BT
    double m_oldLx;
    double m_newLx;
    double m_Lxon;
    int m_step;
    double *m_pLmin2;
    double *m_pLmax2;
    double *m_pminAngle;
    SpringPotentialBetweenTwoGroups *m_pSPBTG;



private:
    //=== update since 2023
    bool CheckMinDistance();
    
    //=== updates aug 2023
    private:
    Energy *m_pEnergyCalculator;
    std::vector<CNTCell *> m_pAllCNT;
    GenerateCNTCells *m_pGenCNT;
    State *m_pState;
    Vec3D *m_pBox;
    
    //=== copy containeir
    std::vector<triangle> m_ActiveT;
    std::vector<vertex > m_ActiveV;
    std::vector<links > m_SurfL;
    std::vector<links > m_MSurfL;
    std::vector<links > m_EdgeL;
    MESH* m_pMESH;



};


#endif
