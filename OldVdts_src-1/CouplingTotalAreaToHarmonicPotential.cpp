#include "CouplingTotalAreaToHarmonicPotential.h"
#include "Nfunction.h"


CouplingTotalAreaToHarmonicPotential::CouplingTotalAreaToHarmonicPotential(int eqsteps, double gamma,  double K0)
{
    double pi = acos(-1);
     m_6SQPI = 1.0/(6.0*sqrt(pi));   /// 1/6pi^1/2    
    m_TotalArea = 0;
    m_NoEQStep = eqsteps;
    m_K0 = K0;
    m_State = true;
    m_Gamma  = gamma;
    m_K0 = m_K0/2;
    m_NT = 1;
    
    if(m_Gamma<0 || m_Gamma>1)
    {
    std::cout<<"---> error in constant area; gamma is bad; make sure you know what are you doing \n";
    exit(0);
    }

    
}
CouplingTotalAreaToHarmonicPotential::~CouplingTotalAreaToHarmonicPotential()
{
    
}
//==========================================================
void CouplingTotalAreaToHarmonicPotential::Initialize(std::vector<triangle *> &pTriangle)
{
    double A=0.0;
	m_NT = double(pTriangle.size());
    for (std::vector<triangle *>::iterator it = pTriangle.begin() ; it != pTriangle.end(); ++it)
        A+=(*it)->GetArea();

    m_TotalArea = A;

   double A0 = m_NT*sqrt(3)/4.0; /// a_t = sqrt(3)/4.0*l^2 ;; l=sqrt(1-3)
   
   
   m_A0 = (1+2*m_Gamma)*A0;   // selecting between 1-3
    
    m_K0 = m_K0/m_NT/3*16;       // making the K0 t dependent   E=0.5*N*K*(A/A0-1)^2; K=0.5*N*K/A0^2==> E=0.5*N*K/A0^2(A-A0)^2

}
double CouplingTotalAreaToHarmonicPotential::CalculateEnergyChange(int step, double oldarea, double newarea)
{

    double alpha=1;
    if(step<m_NoEQStep)
        alpha= double(step)/double(m_NoEQStep);
        
        double da = newarea - oldarea;
        
        //DE = (A+DA-A0)^2-(A-A0)^2 = 2(A-A0+DA/2)DA


       double DE = 2*(m_TotalArea+da/2-m_A0)*da;
       
       DE = alpha*m_K0*DE;



    return DE;
}
void CouplingTotalAreaToHarmonicPotential::UpdateArea(double oldarea,  double newarea)
{

    m_TotalArea+=newarea-oldarea;
}
double CouplingTotalAreaToHarmonicPotential::AreaofTrianglesAroundVertex(vertex *pv)
{
    double A=0.0;
    
    
    std::vector <triangle *> pvT=pv->GetVTraingleList();
    for (std::vector<triangle *>::iterator it = pvT.begin() ; it != pvT.end(); ++it)
        A+=(*it)->GetArea();

    
    return A;
}

