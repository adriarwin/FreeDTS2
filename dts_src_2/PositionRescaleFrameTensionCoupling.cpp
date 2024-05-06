

#include <stdio.h>
#include "PositionRescaleFrameTensionCoupling.h"
#include "Nfunction.h"
#include "vertex.h"
#include "State.h"
#include "Voxelization.h"


/*
===============================================================================================================
 Last update Aug 2023 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class changes the box in x and y direction to minimize the energy. It is called in input file by Frame_Tension  = on 0 2, where the first flag should be on to start this command and the second argument is frame tension and the third argument is the priodic of the operating.

The way it works is based on the changing the box in x and y direction by a small size of dr by calling "MCMoveBoxChange" function.
=================================================================================================================
*/
PositionRescaleFrameTensionCoupling::PositionRescaleFrameTensionCoupling(int period, double sigma, State *pState){
    m_SigmaP = sigma;
    m_pState = pState;
    m_Period = period;

}
PositionRescaleFrameTensionCoupling::~PositionRescaleFrameTensionCoupling()
{
    
}
void PositionRescaleFrameTensionCoupling::Initialize()
{
    std::cout<<"---> the algorithm for box change involves applying PositionRescaleFrameTensionCoupling. \n";
    
    m_pBox=(m_pState->GetMesh())->GetBox();
   // m_pLmin2 = &(m_pState->m_MinVerticesDistanceSquare);
   // m_pLmax2 = &(m_pState->m_MaxLinkLengthSquare);
  //  m_pminAngle = &(m_pState->m_MinFaceAngle);
   // m_Beta = m_pState->m_Beta;
}
bool PositionRescaleFrameTensionCoupling::ChangeBoxSize(int step){
    
//---> if does not match the preiod, return false
    if(step%m_Period != 0)
        return false;
    
    double dx = 0;
    double dy = 0;

//---> find the size of box change; isotropic method
    dx=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
    dx = m_DR*dx;
    dy = dx*((*m_pBox)(1))/(*m_pBox)(0);


//---> check if we do the move, the distance will be normal
    if(!VertexMoveIsFine(dx, dy, *m_pLmin2, *m_pLmax2)){
        return false;
    }
    
//---> check the angles
    
    
    //const std::vector<inclusion *>& pAllInclusion = m_pState->GetMesh()->GetInclusion();

    
    
    double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);

    
    return true;
}
bool PositionRescaleFrameTensionCoupling::VertexMoveIsFine(double dx,double dy, double mindist2, double maxdist2){
    
    
    
    return true;
}



/*
bool PositionRescaleFrameTensionCoupling::MCMoveBoxChange(double dr, double * tot_Energy, double temp, int step, GenerateCNTCells *pGenCNT)
{
//==== some updates Aug. 2023
//==========================
    

    





    
    if(CheckMinDistance()==false)
        return false;
    
    if(CheckMaxLinkLength()==false)
        return false;
    

    
    // === when is coupled to Fixed global curvature; before performing the move we need to calculate the area and total curvature
    // === before any move should be done
    double DetaR = 0;
    double DeltaA = 0;
    if(m_pState->GetGlobalCurvature()->GetState()==true)
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
            double y_o = (*it)->GetVYPos();
            double x_n = m_Lnox*x_o;
            double y_n = m_Lnoy*y_o;
            Vec3D dx(x_n-x_o,y_n-y_o,0);
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
    if((m_pState->GetTotalAreaCoupling())->GetState()==true)
    {
        oldarea_FixedTotArea = (m_pState->GetTotalAreaCoupling())->GetTotalArea();
        for (std::vector<triangle *>::iterator it = (m_pMESH->m_pActiveT).begin() ; it != (m_pMESH->m_pActiveT).end(); ++it)
            newarea_FixedTotArea+=(*it)->GetArea();

        DE_totA =(m_pState->GetTotalAreaCoupling())->CalculateEnergyChange(m_step,oldarea_FixedTotArea,newarea_FixedTotArea);
    }
    //=== previous might be removed in future
    
    // energy old and new
    double DE=0;        // total energy change
    double DEArea=0;    // energy contribuion from projected area change
    
    double oldEnergy = (*tot_Energy);
    double newEnergy = m_pEnergyCalculator->TotalEnergy((m_pMESH->m_pSurfV), (m_pMESH->m_pHL));
    newEnergy+= m_pEnergyCalculator->TotalEnergy((m_pMESH->m_pEdgeV), (m_pMESH->m_pEdgeL));
    
    DE= newEnergy-oldEnergy+DE_totA;

    //=== energy change of constant projected area
    DEArea = -m_SigmaP*(m_newLx*m_newLy-m_oldLx*m_oldLy);
    //=== energy of changes in gloabal curvature due to CouplingtoFixedGlobalCurvature
    //=== the value of DeltaA and DetaR were computed before now we subtract new area
    double eG = 0;
    if(m_pState->GetGlobalCurvature()->GetState()==true)
    {
        for (std::vector<vertex *>::iterator it = (m_pMESH->m_pActiveV).begin() ; it != (m_pMESH->m_pActiveV).end(); ++it)
        {
            std::vector<double> C=(*it)->GetCurvature();
            double area = (*it)->GetArea();
            DeltaA-=area;
            DetaR-=(C.at(0)+C.at(1))*area;
        }

        eG = m_pState->GetGlobalCurvature()->CalculateEnergyChange(-DeltaA,-DetaR);
    }
    
    
     ///=== maybe in future
     // harmonic ponettailas should not be done in x and y direction many problems appear so we assume this
     //==============
    if(m_pSPBTG->GetState()==true)
    {
     double harmonicde=0;

        m_pSPBTG->CalculateEnergy(m_step);
        double e1 = m_pSPBTG->GetEnergy();
        for (std::vector<vertex *>::iterator it = m_ActiveV.begin() ; it != m_ActiveV.end(); ++it)
        {
            double x = (*it)->GetVXPos();
            double y = (*it)->GetVYPos();
            Vec3D dr(m_Lnox*x-x,m_Lnoy*y-y,0);
            m_pSPBTG->MovingVertex((*it), dr);
        }
        m_pSPBTG->CalculateEnergy(m_step);
        double e2 = m_pSPBTG->GetEnergy();
        harmonicde = e2-e1;
    }





    double diff_energy = (DE+DEArea+eG);
    int nv = (m_pMESH->m_pActiveV).size();
    
        //if(nv*log(AreaRatio)-m_Beta*diff_energy>log(temp) )
        if(pow((AreaRatio),nv)*exp(-m_Beta*(diff_energy+dE_nematic_force))>temp )
       {
        if(m_pState->GetGlobalCurvature()->GetState()==true)
            m_pState->GetGlobalCurvature()->UpdateEnergyChange(-DeltaA,-DetaR);
         
         (*tot_Energy)=(*tot_Energy)+DE;  // we only add DE because this is only elastic energy
         (m_pState->m_pConstant_NematicForce)->m_ActiveEnergy+=dE_nematic_force;
         // Update area for the constant area system
         if((m_pState->GetTotalAreaCoupling())->GetState()==true)
         (m_pState->GetTotalAreaCoupling())->UpdateArea(oldarea_FixedTotArea,newarea_FixedTotArea);

         m_Move=true;
           


      }
      else
      {
          RejectMove();
          m_Move=false;


      }//    if(pow((AreaRatio),nv)*exp(-m_Beta*diff_energy)>temp )
    return m_Move;
}
    


void PositionRescaleFrameTensionCoupling::CheckCNTSize()
{

    m_pAllCNT = m_pGenCNT->GetAllCNTCells();
    Nfunction f;
    CNTCell*  TemCNT=m_pAllCNT.at(0);
    Vec3D cntlength=TemCNT->GetCNTSidemax()-TemCNT->GetCNTSidemin();

    if(cntlength(1)<1.8 || cntlength(0)<1.8)
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
double PositionRescaleFrameTensionCoupling::DistanceSquardBetweenTwoVertices(vertex * v1,vertex * v2,Vec3D Box)
{
    double x2=v2->GetVXPos();
    double y2=v2->GetVYPos();
    double z2=v2->GetVZPos();
    double x1=v1->GetVXPos();
    double y1=v1->GetVYPos();
    double z1=v1->GetVZPos();
	x2=m_Lnox*x2;
	x1=m_Lnox*x1;
	y2=m_Lnoy*y2;
	y1=m_Lnoy*y1;


    double dx=x2-x1;
    double dy=y2-y1;
    double dz=z2-z1;
    
    
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

bool PositionRescaleFrameTensionCoupling::CheckFaceAngle()
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

void PositionRescaleFrameTensionCoupling::PerformMove()
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
        (*m_pBox)(1)=m_newLy;

        //=== update the position of all vertices now: the move has been performed
        for (std::vector<vertex *>::iterator it = (m_pMESH->m_pActiveV).begin() ; it != (m_pMESH->m_pActiveV).end(); ++it)
        {
            double x = (*it)->GetVXPos();
            (*it)->UpdateVXPos(m_Lnox*x);
            double y = (*it)->GetVYPos();
            (*it)->UpdateVYPos(m_Lnoy*y);
            
        }
}

void PositionRescaleFrameTensionCoupling::AcceptMove()
{
    // we need to resize the cnt cells;
    for (std::vector<CNTCell *>::iterator it = m_pAllCNT.begin() ; it != m_pAllCNT.end(); ++it)
    {
        Vec3D side1=(*it)->GetCNTSidemax();
        Vec3D side2=(*it)->GetCNTSidemin();
        side1(0)=side1(0)*m_Lnox;
        side1(1)=side1(1)*m_Lnoy;
        side2(0)=side2(0)*m_Lnox;
        side2(1)=side2(1)*m_Lnoy;
        (*it)->UpdateCNTSidemax(side1);
        (*it)->UpdateCNTSidemin(side2);
        
    }

}

void PositionRescaleFrameTensionCoupling::RejectMove()
{
    //=== return the box to its orginal position
    (*m_pBox)(0)=(*m_pBox)(0)-m_drx;
    (*m_pBox)(1)=(*m_pBox)(1)-m_dry;
	
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

bool   PositionRescaleFrameTensionCoupling::CheckFaceAngle(links * l)
{

    Vec3D n=(l->GetTriangle())->GetNormalVector();
    Vec3D n3=((l->GetMirrorLink())->GetTriangle())->GetNormalVector();

    if(n.dot(n,n3)<(*m_pminAngle))
        return false;

    return true;
}


bool PositionRescaleFrameTensionCoupling::CheckMinDistance()
{
    // Note: Function "DistanceSquardBetweenTwoVertices" uses rescaled distance so this somehow reperasent after move; yet the move has not been done
    // This could also be more optimised, we are using CNT cells maybe just direct commparation works
    Vec3D Box;
    Box(0)=m_newLx;
    Box(1)=m_newLy;
    Box(2)=(*m_pBox)(2);

    
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

bool PositionRescaleFrameTensionCoupling::CheckMaxLinkLength()
{
    Vec3D Box;
    Box(0)=m_newLx;
    Box(1)=m_newLy;
    Box(2)=(*m_pBox)(2);

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

*/








