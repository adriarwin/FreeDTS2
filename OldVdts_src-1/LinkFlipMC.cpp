

#include <stdio.h>
#include "LinkFlipMC.h"
#include "Curvature.h"
#include "State.h"
LinkFlipMC::LinkFlipMC()
{
    
}
LinkFlipMC::LinkFlipMC(State *pState)
{
    m_pState = pState;
    m_pBox = (pState->m_pMesh)->m_pBox;
    m_pminAngle = &(m_pState->m_MinFaceAngle);
    m_pLmin2    = &(m_pState->m_MinVerticesDistanceSquare);
    m_pLmax2    = &(m_pState->m_MaxLinkLengthSquare);
    m_pInc     = m_pState->m_pinc_ForceField;
    m_pTotEnergy=&(m_pState->m_TotEnergy);
    m_step = 0;
    m_Beta = m_pState->m_Beta;
}
LinkFlipMC::~LinkFlipMC()
{
    
}
void LinkFlipMC::MC_FlipALink(int step, links *plinks,  double temp)
{
    
    if(plinks->m_LinkType==1)
    {
        std::cout<<" error---> an edge link cannot be flipped "<<plinks->GetMirrorFlag()<<" \n";
        if(plinks->GetMirrorFlag()==true)
        std::cout<<" error---> source code error: mirror flag and edge type are not consistent \n";
        exit(0);
    }
    
    
    m_AVer.clear();
    m_AT.clear();
    m_AL.clear();
    m_NeighborLinks.clear();
    m_LIntEChange.clear();
    m_pLIntEChange.clear();
    
    
    Vec3D *Box = m_pBox;
m_step = step;
m_oldEnergy=0.0;
m_pLinks=plinks;
m_DetaR = 0;
m_DeltaA = 0;
m_Thermal=temp;
m_Mirror=m_pLinks->GetMirrorLink();
m_face=true;
m_MoveValidity=0;
m_L1=m_pLinks->GetNeighborLink1();
m_L2=m_pLinks->GetNeighborLink2();
m_L3=m_Mirror->GetNeighborLink1();
m_L4=m_Mirror->GetNeighborLink2();
m_T1=m_pLinks->GetTriangle();
m_T2=m_Mirror->GetTriangle();
    
m_V1=m_pLinks->GetV1();
m_V2=m_pLinks->GetV2();
m_V3=m_pLinks->GetV3();
m_V4=m_Mirror->GetV3();
    

    if(m_pState->GetGlobalCurvature()->GetState()==true)
    {
        {std::vector<double> C=m_V1->GetCurvature();
            double area = m_V1->GetArea();
            m_DeltaA+=area;
            m_DetaR+=(C.at(0)+C.at(1))*area;}
        {std::vector<double> C=m_V2->GetCurvature();
            double area = m_V2->GetArea();
            m_DeltaA+=area;
            m_DetaR+=(C.at(0)+C.at(1))*area;}
        {std::vector<double> C=m_V3->GetCurvature();
            double area = m_V3->GetArea();
            m_DeltaA+=area;
            m_DetaR+=(C.at(0)+C.at(1))*area;}
        {std::vector<double> C=m_V4->GetCurvature();
            double area = m_V4->GetArea();
            m_DeltaA+=area;
            m_DetaR+=(C.at(0)+C.at(1))*area;}
    }
    
    
    Vec3D R3(m_V3->GetVXPos(),m_V3->GetVYPos(),m_V3->GetVZPos());
    Vec3D R4(m_V4->GetVXPos(),m_V4->GetVYPos(),m_V4->GetVZPos());
    //==== checking if the R4-R3 is larger then box size, meaning that it has PBC problem
    
    Vec3D DR=(R4-R3);
    if(R3.dot(DR,DR)>3)
    {
        double dx=DR(0);
        double dy=DR(1);
        double dz=DR(2);
        if(fabs(dx)>(*Box)(0)/2.0)
        {
            if(dx<0)
                dx=(*Box)(0)+dx;
            else if(dx>0)
                dx=dx-(*Box)(0);
        }
        if(fabs(dy)>(*Box)(1)/2.0)
        {
 
            if(dy<0)
                dy=(*Box)(1)+dy;
            else if(dy>0)
                dy=dy-(*Box)(1);
 
        }
        if(fabs(dz)>(*Box)(2)/2.0)
        {
            if(dz<0)
                dz=(*Box)(2)+dz;
            else if(dz>0)
                dz=dz-(*Box)(2);
        }
        
         DR(0) = dx;
         DR(1) = dy;
         DR(2) = dz;
    }

       // m_DV=R3.dot(DR,m_T1->GetAreaVector())/6;

    //=========== This can be more optimised maybe
    m_simplexarea = 0;
    m_simplexvolume = 0;
    
    if((m_pState->GetVolumeCoupling())->GetState()==true || (m_pState->GetTotalAreaCoupling())->GetState()==true)
    {
            if((m_pState->GetVolumeCoupling())->GetState()==true)
            {
            m_simplexvolume+=(m_pState->GetVolumeCoupling())->SingleTriangleVolume(m_T1);
            m_simplexvolume+=(m_pState->GetVolumeCoupling())->SingleTriangleVolume(m_T2);
            }
            m_simplexarea+=(m_T1)->GetArea();
            m_simplexarea+=(m_T2)->GetArea();
    }

    EnergyDifference();

}
void LinkFlipMC::EnergyDifference()
{
    double DE=0.0;
    if(CheckFlipCondition()==false)
        return;

        PerformMove();
        m_T1->UpdateNormal_Area(m_pBox);
        m_T2->UpdateNormal_Area(m_pBox);
    
        m_face=CheckFaceAngle();

        if(m_face==false)
        {
            RejectMove();
            return;
        }


    
    // we should check if they have a mirror
        if(m_pLinks->m_LinkType==0)
        {
            m_pLinks->UpdateNormal();
        	m_pLinks->UpdateShapeOperator(m_pBox);
        }
        else
        {
            m_pLinks->UpdateEdgeVector(m_pBox);
        }
        if(m_L1->m_LinkType==0)
        {
            m_L1->UpdateNormal();
        	m_L1->UpdateShapeOperator(m_pBox);
        }
        else
        {
            m_L1->UpdateEdgeVector(m_pBox);
        }
        if(m_L2->m_LinkType==0)
        {
            m_L2->UpdateNormal();
            m_L2->UpdateShapeOperator(m_pBox);
        }
        else
        {
            m_L2->UpdateEdgeVector(m_pBox);
        }
        if(m_L3->m_LinkType==0)
        {
            m_L3->UpdateNormal();
            m_L3->UpdateShapeOperator(m_pBox);
        }
        else
        {
            m_L3->UpdateEdgeVector(m_pBox);
        }
        if(m_L4->m_LinkType==0)
        {
            m_L4->UpdateNormal();
            m_L4->UpdateShapeOperator(m_pBox);
        }
        else
        {
            m_L4->UpdateEdgeVector(m_pBox);
        }



        if(m_V1->m_VertexType==0)
            (m_pState->CurvatureCalculator())->SurfVertexCurvature(m_V1);
        else
            (m_pState->CurvatureCalculator())->EdgeVertexCurvature(m_V1);
        
        if(m_V2->m_VertexType==0)
            (m_pState->CurvatureCalculator())->SurfVertexCurvature(m_V2);
        else
            (m_pState->CurvatureCalculator())->EdgeVertexCurvature(m_V2);
    
        if(m_V3->m_VertexType==0)
            (m_pState->CurvatureCalculator())->SurfVertexCurvature(m_V3);
        else
            (m_pState->CurvatureCalculator())->EdgeVertexCurvature(m_V3);
    
        if(m_V4->m_VertexType==0)
            (m_pState->CurvatureCalculator())->SurfVertexCurvature(m_V4);
        else
            (m_pState->CurvatureCalculator())->EdgeVertexCurvature(m_V4);
    
	Energy EE(m_pInc);
	double NewEnergy=EE.Energy_OneLinkFlip(m_pLinks);
    
    for (std::vector<links *>::iterator it = m_pLIntEChange.begin() ; it != m_pLIntEChange.end(); ++it)
        NewEnergy+=EE.TwoInclusionsInteractionEnergy(*it);
    
    
	DE=NewEnergy-m_oldEnergy;
	m_EnergyDifference=DE;
    
    
    double eG = 0;
    if(m_pState->GetGlobalCurvature()->GetState()==true)
    {
        
        {std::vector<double> C=m_V1->GetCurvature();
            double area = m_V1->GetArea();
            m_DeltaA-=area;
            m_DetaR-=(C.at(0)+C.at(1))*area;}
        {std::vector<double> C=m_V2->GetCurvature();
            double area = m_V2->GetArea();
            m_DeltaA-=area;
            m_DetaR-=(C.at(0)+C.at(1))*area;}
        {std::vector<double> C=m_V3->GetCurvature();
            double area = m_V3->GetArea();
            m_DeltaA-=area;
            m_DetaR-=(C.at(0)+C.at(1))*area;}
        {std::vector<double> C=m_V4->GetCurvature();
            double area = m_V4->GetArea();
            m_DeltaA-=area;
            m_DetaR-=(C.at(0)+C.at(1))*area;}
        eG = m_pState->GetGlobalCurvature()->CalculateEnergyChange(-m_DeltaA,-m_DetaR);

    }
    
    //=========== This can be more optimised maybe
    double DEPV = 0;    //
    double DE_totA = 0; // energy asscoaited with change in total area
    double newsimplexarea = 0;
    double newsimplexvolume = 0;
    
    if((m_pState->GetVolumeCoupling())->GetState()==true || (m_pState->GetTotalAreaCoupling())->GetState()==true)
    {
        if((m_pState->GetVolumeCoupling())->GetState()==true)
        {
            newsimplexvolume+=(m_pState->GetVolumeCoupling())->SingleTriangleVolume(m_T1);
            newsimplexvolume+=(m_pState->GetVolumeCoupling())->SingleTriangleVolume(m_T2);
        }


            newsimplexarea+=m_T1->GetArea();
            newsimplexarea+=m_T2->GetArea();
            
            
	        if((m_pState->GetVolumeCoupling())->GetState()==true)
            DEPV = (m_pState->GetVolumeCoupling())->GetEnergyChange(m_step,m_simplexarea,m_simplexvolume,newsimplexarea,newsimplexvolume);
        
           if((m_pState->GetTotalAreaCoupling())->GetState()==true)
            DE_totA = (m_pState->GetTotalAreaCoupling())->CalculateEnergyChange(m_step,m_simplexarea,newsimplexarea);
    }
                double diff_energy = m_Beta*(DE+eG+DEPV+DE_totA);
                
                //std::cout<<DE<<"  "<<eG<<"  "<<DEPV<<"   "<<DE_OP<<"  \n";
                                
            	if(diff_energy<=0 )
            	{
                	AccpetMove();
                	(*m_pTotEnergy)=(*m_pTotEnergy)+DE;
                    m_MoveValidity=1;
                    (m_pState->GetVolumeCoupling())->UpdateArea_Volume(m_simplexarea,m_simplexvolume,newsimplexarea,newsimplexvolume);
                    (m_pState->GetTotalAreaCoupling())->UpdateArea(m_simplexarea,newsimplexarea);
            	}
            	else if(exp(-diff_energy)>m_Thermal )
             	{
                 	AccpetMove();
                 	(*m_pTotEnergy)=(*m_pTotEnergy)+DE;
                    m_MoveValidity=1;
                    (m_pState->GetVolumeCoupling())->UpdateArea_Volume(m_simplexarea,m_simplexvolume,newsimplexarea,newsimplexvolume);
                    (m_pState->GetTotalAreaCoupling())->UpdateArea(m_simplexarea,newsimplexarea);
             	}
            	else 
            	{

                	RejectMove();

            	}
}
void LinkFlipMC::AccpetMove()
{
    if(m_pState->GetGlobalCurvature()->GetState()==true)
    {
        m_pState->GetGlobalCurvature()->UpdateEnergyChange(-m_DeltaA,-m_DetaR);
    }

}
void LinkFlipMC::PerformMove()
{
    m_NeighborLinks.push_back(*m_L1);
    m_NeighborLinks.push_back(*m_L2);
    m_NeighborLinks.push_back(*m_L3);
    m_NeighborLinks.push_back(*m_L4);
    
    if(m_L1->m_LinkType==0)
    m_NeighborLinks.push_back(*(m_L1->GetMirrorLink()));
    if(m_L2->m_LinkType==0)
    m_NeighborLinks.push_back(*(m_L2->GetMirrorLink()));
    if(m_L3->m_LinkType==0)
    m_NeighborLinks.push_back(*(m_L3->GetMirrorLink()));
    if(m_L4->m_LinkType==0)
    m_NeighborLinks.push_back(*(m_L4->GetMirrorLink()));

    m_AVer.push_back(*m_V1);
    m_AVer.push_back(*m_V2);
    m_AVer.push_back(*m_V3);
    m_AVer.push_back(*m_V4);
    
    m_AT.push_back(*m_T1);
    m_AT.push_back(*m_T2);
    
    m_AL.push_back(*m_pLinks);
    m_AL.push_back(*m_Mirror);
    
    m_oldEnergy+=m_V1->GetEnergy();
    m_oldEnergy+=m_V2->GetEnergy();
    m_oldEnergy+=m_V3->GetEnergy();
    m_oldEnergy+=m_V4->GetEnergy();
    
    
    //// this part of the code can become more efficient
    
    std::vector <links*> temlinklist;
        if(m_V1->VertexOwnInclusion()==true)
        {
            std::vector<links *> ltem=m_V1->GetVLinkList();
            temlinklist.insert(temlinklist.end(), ltem.begin(), ltem.end());
            if(m_V1->m_VertexType==1)
            temlinklist.push_back(m_V1->m_pPrecedingEdgeLink);
        }
        if(m_V2->VertexOwnInclusion()==true)
        {
            std::vector<links *> ltem=m_V2->GetVLinkList();
            temlinklist.insert(temlinklist.end(), ltem.begin(), ltem.end());
            if(m_V2->m_VertexType==1)
            temlinklist.push_back(m_V2->m_pPrecedingEdgeLink);
        }
        if(m_V3->VertexOwnInclusion()==true)
        {
            std::vector<links *> ltem=m_V3->GetVLinkList();
            temlinklist.insert(temlinklist.end(), ltem.begin(), ltem.end());
            if(m_V3->m_VertexType==1)
            temlinklist.push_back(m_V3->m_pPrecedingEdgeLink);
        }
        if(m_V4->VertexOwnInclusion()==true)
        {
            std::vector<links *> ltem=m_V4->GetVLinkList();
            temlinklist.insert(temlinklist.end(), ltem.begin(), ltem.end());
            if(m_V4->m_VertexType==1)
            temlinklist.push_back(m_V4->m_pPrecedingEdgeLink);
        }

    
    for (std::vector<links *>::iterator it = temlinklist.begin() ; it != temlinklist.end(); ++it)
    {
        bool addit=true;
        for (std::vector<links *>::iterator it2 = m_pLIntEChange.begin() ; it2 != m_pLIntEChange.end(); ++it2)
        {
            if((*it2)->GetID()==(*it)->GetID())
                addit=false;
            else if((*it2)->GetMirrorFlag()==true)
            if(((*it2)->GetMirrorLink())->GetID()==(*it)->GetID())
                addit=false;


        }
        if(addit==true)
            m_pLIntEChange.push_back((*it));

    }
    
    for (std::vector<links *>::iterator it = m_pLIntEChange.begin() ; it != m_pLIntEChange.end(); ++it)
    {
        m_oldEnergy+=2*((*it)->GetIntEnergy());
        m_LIntEChange.push_back(*(*it));
    }

m_pLinks->Flip();
}


bool LinkFlipMC::CheckFlipCondition()
{
Vec3D Box=(*m_pBox);
//==================== check if the vertex has more then three link
	std::vector <vertex *> list1=m_V1->GetVNeighbourVertex();
    std::vector <vertex *> list2=m_V2->GetVNeighbourVertex();
    //std::vector <vertex *> list3=m_V3->GetVNeighbourVertex();
    //std::vector <vertex *> list4=m_V4->GetVNeighbourVertex();

	if(list1.size()<4 || list2.size()<4 )
	 return false;

//	if(list3.size()>9 || list4.size()>9 )
//	condition=false;

//==================== check if this link already exist
    std::vector <vertex *> mvn=m_V3->GetVNeighbourVertex();
    for (std::vector<vertex *>::iterator it = mvn.begin() ; it != mvn.end(); ++it)
    {
        if((*it)->GetVID()==m_V4->GetVID())
            return false;
    }
//==================== check the length of the new link

 		double x1=m_V3->GetVXPos();
		double y1=m_V3->GetVYPos();
		double z1=m_V3->GetVZPos();
    		double x2=m_V4->GetVXPos();
    		double y2=m_V4->GetVYPos();
    		double z2=m_V4->GetVZPos();
		double dx=x2-x1;
        if(fabs(dx)>Box(0)/2.0)
        {
            if(dx<0)
                dx=Box(0)+dx;
            else if(dx>0)
                dx=dx-Box(0);
        }
		double dy=y2-y1;
        if(fabs(dy)>Box(1)/2.0)
        {
            if(dy<0)
                dy=Box(1)+dy;
            else if(dy>0)
                dy=dy-Box(1);
        }
		double dz=z2-z1;
        if(fabs(dz)>Box(2)/2.0)
        {
            if(dz<0)
                dz=Box(2)+dz;
            else if(dz>0)
                dz=dz-Box(2);
        }
		double l2=dx*dx+dy*dy+dz*dz;

       if(l2>(*m_pLmax2) || l2<(*m_pLmin2))
	      return false;
    
//==================== check if the v3 of the links is same
return true;
}
void LinkFlipMC::RejectMove()
{
     *(m_V1)=m_AVer.at(0);
     *(m_V2)=m_AVer.at(1);
     *(m_V3)=m_AVer.at(2);
     *(m_V4)=m_AVer.at(3);
    
    *(m_T1)=m_AT.at(0);
    *(m_T2)=m_AT.at(1);

    *(m_pLinks)=m_AL.at(0);
    *(m_Mirror)=m_AL.at(1);
    
    
     *(m_L1)=m_NeighborLinks.at(0);
     *(m_L2)=m_NeighborLinks.at(1);
     *(m_L3)=m_NeighborLinks.at(2);
     *(m_L4)=m_NeighborLinks.at(3);
    
    int i=4;
    if(m_L1->m_LinkType==0)
    {
     *(m_L1->GetMirrorLink())=m_NeighborLinks.at(i);
        i++;
    }
    if(m_L2->m_LinkType==0)
    {
     *(m_L2->GetMirrorLink())=m_NeighborLinks.at(i);
        i++;
    }
    if(m_L3->m_LinkType==0)
    {
     *(m_L3->GetMirrorLink())=m_NeighborLinks.at(i);
        i++;
    }
    if(m_L4->m_LinkType==0)
    {
     *(m_L4->GetMirrorLink())=m_NeighborLinks.at(i);
        i++;
    }
    
    std::vector<links>::iterator itr=m_LIntEChange.begin();
    for (std::vector<links *>::iterator it = m_pLIntEChange.begin() ; it != m_pLIntEChange.end(); ++it)
    {
        double en=(*itr).GetIntEnergy();
        (*it)->UpdateIntEnergy(en);
        if((*it)->GetMirrorFlag()==true)
        ((*it)->GetMirrorLink())->UpdateIntEnergy(en);
#if TEST_MODE == Enabled
        
        if((*itr).GetID()!=(*it)->GetID())
            std::cout<<"  "<<(*itr).GetID()<<"  "<<(*it)->GetID()<<"  Link flip: These link ids should be same but it is not \n";
        
#endif
        ++itr;
    }

}
//====
bool   LinkFlipMC::CheckFaceAngle()
{
    Vec3D n1=(m_pLinks->GetTriangle())->GetNormalVector();
    Vec3D n2=((m_pLinks->GetMirrorLink())->GetTriangle())->GetNormalVector();
    
    if(n1.dot(n1,n2)<(*m_pminAngle))
        return false;

    if((m_pLinks->GetNeighborLink1())->m_LinkType==0)
    {
        Vec3D n11=(((m_pLinks->GetNeighborLink1())->GetMirrorLink())->GetTriangle())->GetNormalVector();
        if (n1.dot(n1,n11)<(*m_pminAngle))
            return false;
    }
    

    if((m_pLinks->GetNeighborLink2())->m_LinkType==0)
    {
        Vec3D n12=(((m_pLinks->GetNeighborLink2())->GetMirrorLink())->GetTriangle())->GetNormalVector();
        if (n1.dot(n1,n12)<(*m_pminAngle))
            return false;
    }

    
    links *mlink = m_pLinks->GetMirrorLink();
    
    if((mlink->GetNeighborLink1())->m_LinkType==0)
    {
        Vec3D n21=(((mlink->GetNeighborLink1())->GetMirrorLink())->GetTriangle())->GetNormalVector();
        if (n2.dot(n2,n21)<(*m_pminAngle))
            return false;
    }
    
    if((mlink->GetNeighborLink2())->m_LinkType==0)
    {
        Vec3D n22=(((mlink->GetNeighborLink2())->GetMirrorLink())->GetTriangle())->GetNormalVector();
        if (n2.dot(n2,n22)<(*m_pminAngle))
            return false;
    }

    return true;
}
