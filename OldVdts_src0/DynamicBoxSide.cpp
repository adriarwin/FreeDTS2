

#include <stdio.h>
#include "DynamicBoxSide.h"
#include "Nfunction.h"
#include "vertex.h"
#include "Curvature.h"
#include "State.h"

/*
===============================================================================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class changes the box in x direction to minimize the energy. It is called in input file.
=================================================================================================================
*/
DynamicBoxSide::DynamicBoxSide()
{
    
    
}
DynamicBoxSide::DynamicBoxSide(int tau, double f,State *pstate)
{
    m_F0 = f;
    m_pState = pstate;
    m_Tau = tau;
}
DynamicBoxSide::~DynamicBoxSide()
{
    
}
void DynamicBoxSide::initialize(){
    std::cout<<"---> the algorithm for box change involves applying a constant force on one side. \n";
    m_pEnergyCalculator = m_pState->GetEnergyCalculator();
    m_pBox=(m_pState->m_pMesh)->m_pBox;
    m_pMESH = m_pState->m_pMesh;
    m_Beta = m_pState->m_Beta;
    m_pLmin2 = &(m_pState->m_MinVerticesDistanceSquare);
    m_pLmax2 = &(m_pState->m_MaxLinkLengthSquare);
    m_pminAngle = &(m_pState->m_MinFaceAngle);
    m_pCFGC = m_pState->GetGlobalCurvature();
    m_step = 0;

}
bool DynamicBoxSide::GetCNTCondition(){
    return m_UpdateCNT;
}
int DynamicBoxSide::GetTau() {
    return m_Tau;
}
bool DynamicBoxSide::MCMoveBoxChange(double dx, double * tot_Energy, double temp, int step, GenerateCNTCells *pGenCNT) {

    m_UpdateCNT=true;   // should be set to true. true means do not update it.
    m_pGenCNT = pGenCNT;
    m_oldLx=(*m_pBox)(0);
    m_newLx=(*m_pBox)(0)+dx;
    m_Lxon = m_newLx/m_oldLx;

    m_step = step;
    CheckCNTSize();      // CNT cell should not be smaller then 1.8 after the move; we use 1.8 so no need for frequent changes

    
    if(CheckMinDistance()==false)
        return false;
    
    if(CheckMaxLinkLength()==false)
        return false;
    

    
    // === when is coupled to Fixed global curvature; before performing the move we need to calculate the area and total curvature
    // === before any move should be done
    double DetaR = 0;
    double DeltaA = 0;
    if(m_pCFGC->GetState()==true)
    {
        for (std::vector<vertex *>::iterator it = (m_pMESH->m_pActiveV).begin() ; it != (m_pMESH->m_pActiveV).end(); ++it)
        {
            std::vector<double> C=(*it)->GetCurvature();
            double area = (*it)->GetArea();
            DeltaA+=area;
            DetaR+=(C.at(0)+C.at(1))*area;

        }
    }

    
    
     // we do this before performing the move so that the force is correct
     double dE_nematic_force = 0;
     if((m_pState->m_pConstant_NematicForce)->m_F0!=0)
     {
        for (std::vector<vertex *>::iterator it = (m_pMESH->m_pActiveV).begin() ; it != (m_pMESH->m_pActiveV).end(); ++it)
        {
            double x_o = (*it)->GetVXPos();
            double x_n = (m_Lxon-1)*x_o;
            Vec3D dx(x_n,0,0);
            dE_nematic_force+= (m_pState->m_pConstant_NematicForce)->Energy_of_Force(*it, dx);
        }
     }

    
    // by calling "PerformMove" now the move is performed:
    // box changed; vertices, triangles and links are coppied and the vertices positions are updated (rescaled in x and y)
    PerformMove();

    if(CheckFaceAngle()==false)
    {
        RejectMove();
        m_Move=false;
        return false;
    }

    // we updated the trinagle area and normal
    for (std::vector<triangle *>::iterator it = (m_pMESH->m_pActiveT).begin() ; it != (m_pMESH->m_pActiveT).end(); ++it)
            (*it)->UpdateNormal_Area(m_pBox);

    // we updated the surface links normal and shape operator
    for (std::vector<links *>::iterator it = (m_pMESH->m_pHL).begin() ; it != (m_pMESH->m_pHL).end(); ++it)
        {
            (*it)->UpdateNormal();
            (*it)->UpdateShapeOperator(m_pBox);
        }

    // updating the edge links; they do not have a shape operator
    for (std::vector<links *>::iterator it = (m_pMESH->m_pEdgeL).begin() ; it != (m_pMESH->m_pEdgeL).end(); ++it)
    {
            (*it)->UpdateEdgeVector(m_pBox);
    }

    //== updating the surface curvature
    for (std::vector<vertex *>::iterator it = (m_pMESH->m_pActiveV).begin() ; it != (m_pMESH->m_pActiveV).end(); ++it)
    {
        if((*it)->m_VertexType==0)
        (m_pState->CurvatureCalculator())->SurfVertexCurvature(*it);
        else if((*it)->m_VertexType==1)
        (m_pState->CurvatureCalculator())->EdgeVertexCurvature(*it);
        else
        std::cout<<" error202--> vertex type is unknown "<<std::endl;
    }
    
    //== this is for applying constant area which is not recomended and might be removed in future
    double DE_totA = 0;
    double newarea_FixedTotArea = 0;
    double oldarea_FixedTotArea = 0;
    if((m_pState->GetApply_Constant_Area())->GetState()==true)
    {
        oldarea_FixedTotArea = (m_pState->GetApply_Constant_Area())->GetTotalArea();
        for (std::vector<triangle *>::iterator it = (m_pMESH->m_pActiveT).begin() ; it != (m_pMESH->m_pActiveT).end(); ++it)
            newarea_FixedTotArea+=(*it)->GetArea();

        DE_totA =(m_pState->GetApply_Constant_Area())->GetEnergyChange(m_step,oldarea_FixedTotArea,newarea_FixedTotArea);
    }
    //=== previous might be removed in future
    
    // energy old and new
    double DE=0;        // total energy change
    double DE_DeltaLX=0;    // nergy change of  side x
    
    double oldEnergy = (*tot_Energy);
    double newEnergy = m_pEnergyCalculator->TotalEnergy((m_pMESH->m_pSurfV), (m_pMESH->m_pHL));
    newEnergy+= m_pEnergyCalculator->TotalEnergy((m_pMESH->m_pEdgeV), (m_pMESH->m_pEdgeL));
    
    DE= newEnergy-oldEnergy+DE_totA;

    //=== energy change of  side x
    DE_DeltaLX = -m_F0*(m_newLx-m_oldLx);
    
    //=== energy of changes in gloabal curvature due to CouplingtoFixedGlobalCurvature
    //=== the value of DeltaA and DetaR were computed before now we subtract new area
    double eG = 0;
    if(m_pCFGC->GetState()==true)
    {
        for (std::vector<vertex *>::iterator it = (m_pMESH->m_pActiveV).begin() ; it != (m_pMESH->m_pActiveV).end(); ++it)
        {
            std::vector<double> C=(*it)->GetCurvature();
            double area = (*it)->GetArea();
            DeltaA-=area;
            DetaR-=(C.at(0)+C.at(1))*area;
        }

        eG = m_pCFGC->CalculateEnergyChange(-DeltaA,-DetaR);
    }


    double diff_energy = (DE+DE_DeltaLX+eG);
    int nv = (m_pMESH->m_pActiveV).size();
    
        //if(nv*log(AreaRatio)-m_Beta*diff_energy>log(temp) )
        if(pow((m_newLx/m_oldLx),nv)*exp(-m_Beta*(diff_energy+dE_nematic_force))>temp )
       {
        if(m_pCFGC->GetState()==true)
        m_pCFGC->UpdateEnergyChange(-DeltaA,-DetaR);
         
         (*tot_Energy)=(*tot_Energy)+DE;  // we only add DE because this is only elastic energy
         (m_pState->m_pConstant_NematicForce)->m_ActiveEnergy+=dE_nematic_force;
         // Update area for the constant area system
         if((m_pState->GetApply_Constant_Area())->GetState()==true)
         (m_pState->GetApply_Constant_Area())->UpdateArea(oldarea_FixedTotArea,newarea_FixedTotArea);

         m_Move=true;
           


      }
      else
      {
          RejectMove();
          m_Move=false;

      }//    if(pow((m_newLx/m_oldLx),nv)*exp(-m_Beta*diff_energy)>temp )
    return m_Move;
}
    

/*
 ======================================================
 Checking for validity of the CNT Cells
 =========================================================
 */
void DynamicBoxSide::CheckCNTSize()
{

    m_pAllCNT = m_pGenCNT->GetAllCNTCells();
    Nfunction f;
    CNTCell*  TemCNT=m_pAllCNT.at(0);
    Vec3D cntlength=TemCNT->GetCNTSidemax()-TemCNT->GetCNTSidemin();

    if(cntlength(0)<1.8)
    {
        m_UpdateCNT=false;
        m_Move=false;
        std::string sms=" CNT cell is small. C_Y= "+f.Int_to_String(cntlength(1))+" and C_X= "+f.Int_to_String(cntlength(0))+".  it is re-generated at step "+f.Int_to_String(m_step);
        f.Write_One_LogMessage(sms);
        
    }
    else
    {
        m_UpdateCNT=true;
    }

}
double DynamicBoxSide::DistanceSquardBetweenTwoVertices(vertex * v1,vertex * v2,Vec3D &Box)
{
    double dx=v2->GetVXPos()-v1->GetVXPos();
    double dy=v2->GetVYPos()-v1->GetVYPos();
    double dz=v2->GetVZPos()-v1->GetVZPos();
    dx = dx*m_Lxon;
    
    if(fabs(dx)>Box(0)/2.0)
    {
        if(dx<0)
            dx=Box(0)+dx;
        else if(dx>0)
            dx=dx-Box(0);
    }
    if(fabs(dy)>Box(1)/2.0)
    {
        if(dy<0)
            dy=Box(1)+dy;
        else if(dy>0)
            dy=dy-Box(1);
    }
    if(fabs(dz)>Box(2)/2.0)
    {
        if(dz<0)
            dz=Box(2)+dz;
        else if(dz>0)
            dz=dz-Box(2);
    }
    double l2=dx*dx+dy*dy+dz*dz;
    return l2;
}
/*======================================================================
function to check the angle between two faces of a link
 ======================================================================*/
bool DynamicBoxSide::CheckFaceAngle()
{
   for (std::vector<links *>::iterator it = (m_pMESH->m_pHL).begin() ; it != (m_pMESH->m_pHL).end(); ++it)
    {
        if (CheckFaceAngle((*it))==false)
        {
            m_Move=false;
            return false;
        }
    }
    
    return true;
}

/*======================================================================
 This function makes a copy of vertices, links and trinagles incase of rejection. 
 ======================================================================*/
void DynamicBoxSide::PerformMove()
{
//=== copy the objects to use them in case the move got rejected
    m_ActiveV.clear();
    m_SurfL.clear();
    m_MSurfL.clear();
    m_ActiveT.clear();
    m_EdgeL.clear();
    
    for (std::vector<links *>::iterator it = (m_pMESH->m_pHL).begin() ; it != (m_pMESH->m_pHL).end(); ++it)
    {
        m_SurfL.push_back(*(*it));

    }
    for (std::vector<links *>::iterator it = (m_pMESH->m_pMHL).begin() ; it != (m_pMESH->m_pMHL).end(); ++it)
    {
        m_MSurfL.push_back(*(*it));
    }
    for (std::vector<links *>::iterator it = (m_pMESH->m_pEdgeL).begin() ; it != (m_pMESH->m_pEdgeL).end(); ++it)
    {
        m_EdgeL.push_back(*(*it));
    }

    for (std::vector<triangle *>::iterator it = (m_pMESH->m_pActiveT).begin() ; it != (m_pMESH->m_pActiveT).end(); ++it)
    {
        m_ActiveT.push_back(*(*it));
    }

    for (std::vector<vertex *>::iterator it = (m_pMESH->m_pActiveV).begin() ; it != (m_pMESH->m_pActiveV).end(); ++it)
    {
        m_ActiveV.push_back(*(*it));
    }

    //=== perform the move: first update the box and then update the positions; geometric varialbes are updated elsewhere
        (*m_pBox)(0)=m_newLx;

        //=== update the position of all vertices now: the move has been performed
        for (std::vector<vertex *>::iterator it = (m_pMESH->m_pActiveV).begin() ; it != (m_pMESH->m_pActiveV).end(); ++it)
        {
            double x = (*it)->GetVXPos();
            (*it)->UpdateVXPos(m_Lxon*x);
        }
}



/*
 ======================================================
Function to accept the move and update everything
 =========================================================
 */
void DynamicBoxSide::AcceptMove()
{
    // we need to resize the cnt cells;
    for (std::vector<CNTCell *>::iterator it = m_pAllCNT.begin() ; it != m_pAllCNT.end(); ++it)
    {
        Vec3D side1=(*it)->GetCNTSidemax();
        Vec3D side2=(*it)->GetCNTSidemin();
        side1(0)=side1(0)*m_Lxon;
        side2(0)=side2(0)*m_Lxon;
        (*it)->UpdateCNTSidemax(side1);
        (*it)->UpdateCNTSidemin(side2);
        
    }

}

void DynamicBoxSide::RejectMove()
{
    //=== return the box to its orginal position
    (*m_pBox)(0)=m_oldLx;
	
    //=== all the vertices also get updated to what they were
        int i=0;
        for (std::vector<vertex *>::iterator it = (m_pMESH->m_pActiveV).begin() ; it != (m_pMESH->m_pActiveV).end(); ++it)
        {
                *(*it)=m_ActiveV[i];
                i++;
        }
    //=== all the surf links and their mirros
        i=0;
        for (std::vector<links *>::iterator it = (m_pMESH->m_pHL).begin() ; it != (m_pMESH->m_pHL).end(); ++it)
        {
            *(*it)=m_SurfL[i];
            i++;
        }
        i=0;
        for (std::vector<links *>::iterator it = (m_pMESH->m_pMHL).begin() ; it != (m_pMESH->m_pMHL).end(); ++it)
        {
            *(*it)=m_MSurfL[i];
            i++;
        }
    //=== all the edge links

        i=0;
        for (std::vector<links *>::iterator it = (m_pMESH->m_pEdgeL).begin() ; it != (m_pMESH->m_pEdgeL).end(); ++it)
        {
            *(*it)=m_EdgeL[i];
            i++;
        }
    //=== all the triangles
        i=0;
        for (std::vector<triangle *>::iterator it = (m_pMESH->m_pActiveT).begin() ; it != (m_pMESH->m_pActiveT).end(); ++it)
        {
            *(*it)=m_ActiveT[i];
            i++;
        }
}

bool   DynamicBoxSide::CheckFaceAngle(links * l)
{

    Vec3D n=(l->GetTriangle())->GetNormalVector();
    Vec3D n3=((l->GetMirrorLink())->GetTriangle())->GetNormalVector();

    if(n.dot(n,n3)<(*m_pminAngle))
        return false;

    return true;
}

// since 2023
/*======================================================================
 This function goes through CNT cells  checks the distance between vertices after the box change.
 It is very ineffiecnt, should be changed as soon as possible
 ======================================================================*/
bool DynamicBoxSide::CheckMinDistance()
{
    // Note: Function "DistanceSquardBetweenTwoVertices" uses rescaled distance so this somehow reperasent after move; yet the move has not been done
    // This could also be more optimised, we are using CNT cells maybe just direct commparation works
    Vec3D Box(m_newLx,(*m_pBox)(1),(*m_pBox)(2));


    
    for (std::vector<CNTCell *>::iterator it = m_pAllCNT.begin() ; it != m_pAllCNT.end(); ++it)
    {
        std::vector <vertex *> pver=(*it)->GetVertexList();
        for (int i=0;i<pver.size();i++)
            for (int j=i+1;j<pver.size();j++)
            {
                double l2=DistanceSquardBetweenTwoVertices(pver[i],pver[j],Box);
                if(l2<(*m_pLmin2))
                    return false;
            }

        std::vector <CNTCell *> Nib=(*it)->GetVNeighbourCNTCell();
        for (std::vector<CNTCell *>::iterator it2 = Nib.begin() ; it2 != Nib.end(); ++it2)
        {
            std::vector <vertex *> pverN=(*it2)->GetVertexList();
            
            for (std::vector<vertex *>::iterator itv1 = pver.begin() ; itv1 != pver.end(); ++itv1)
            for (std::vector<vertex *>::iterator itv2 = pverN.begin() ; itv2 != pverN.end(); ++itv2)
            {
                    
                    double l2=DistanceSquardBetweenTwoVertices(*itv1,*itv2,Box);
                    if(l2<(*m_pLmin2))
                        return false;
            }
            
        }
    }

return true;
}
/*======================================================================
 This function goes through all links  checks the distance between vertices after the box change.
 ======================================================================*/
bool DynamicBoxSide::CheckMaxLinkLength()
{
    
    Vec3D Box(m_newLx,(*m_pBox)(1),(*m_pBox)(2));

    for (std::vector<links *>::iterator it = (m_pMESH->m_pHL).begin() ; it != (m_pMESH->m_pHL).end(); ++it)
    {
        vertex *v1= (*it)->GetV1();
        vertex *v2= (*it)->GetV2();
                    double l2=DistanceSquardBetweenTwoVertices(v1,v2,Box);
                    if(l2>(*m_pLmax2))
                    return false;
    }
    for (std::vector<links *>::iterator it = (m_pMESH->m_pEdgeL).begin() ; it != (m_pMESH->m_pEdgeL).end(); ++it)
    {

        vertex *v1= (*it)->GetV1();
        vertex *v2= (*it)->GetV2();

                    double l2=DistanceSquardBetweenTwoVertices(v1,v2,Box);
                    if(l2>(*m_pLmax2))
                    return false;
    }
    return true;
}










