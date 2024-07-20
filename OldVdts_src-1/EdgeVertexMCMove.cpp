

#include <stdio.h>
#include "EdgeVertexMCMove.h"
#include "State.h"
EdgeVertexMCMove::EdgeVertexMCMove()
{
}
EdgeVertexMCMove::EdgeVertexMCMove(State *pState)
{
    m_pState = pState;
    m_pBox = (pState->m_pMesh)->m_pBox;
    m_pminAngle = &(m_pState->m_MinFaceAngle);
    m_pLmin2    = &(m_pState->m_MinVerticesDistanceSquare);
    m_pLmax2    = &(m_pState->m_MaxLinkLengthSquare);
    m_pInc     = m_pState->m_pinc_ForceField;
    m_pTotEnergy=&(m_pState->m_TotEnergy);
    m_pSPBTG = m_pState->Get2GroupHarmonic();
    m_Beta =   m_pState->m_Beta;
    m_step = 0;
}
EdgeVertexMCMove::~EdgeVertexMCMove()
{
    
}
void EdgeVertexMCMove::MC_MoveAVertex(int step, vertex *pvertex, double dx, double dy, double dz,double temp)
{
    
    // 1) first check the cnt and see after such a move what will happen
    // 2) Constant volume is meaning less
    // 3) constant global curvature ???
    // 4) pull force can be
    //
    
    
    if((m_pState->GetVolumeCoupling())->GetState()==true)
    {
        std::cout<<" error---> We are moving an edge vertex, therefore coupling the system to a volume-controlling algorithm is meaningless \n";
        exit(0);
    }
    if ((m_pState->GetApply_Constant_Area())->GetState()==true)
    {
        std::cout<<" error---> We are moving an edge vertex, therefore coupling the system to a fix global curvature is not possible (perhaps meaningless) \n";
        exit(0);

    }
    
    
    m_LIntEChange.clear();
    m_AVer.clear();
    m_RingLinks.clear();
    m_nLinks.clear();
    m_Triangle.clear();
    m_pLIntEChange.clear();
    m_pAVer.clear();
m_step = step;
m_pvertex=pvertex;
    
    m_DetaR = 0;
    m_DeltaA = 0;
m_dx=dx;
m_dy=dy;
m_dz=dz;
m_oldEnergy=0.0;
m_pAVer=m_pvertex->GetVNeighbourVertex();
    
m_ChangeCNT=false;
m_Thermal=temp;
m_MoveValidity=0;
    m_oldX=m_pvertex->GetVXPos();
    m_oldY=m_pvertex->GetVYPos();
    m_oldZ=m_pvertex->GetVZPos();
	m_X=m_oldX+m_dx;
	m_Y=m_oldY+m_dy;
	m_Z=m_oldZ+m_dz;
       
   EnergyDifference();
    m_face=true;

}
void EdgeVertexMCMove::EnergyDifference()
{
    
    if(CheckDistnace()==false)
        return;
    
    double DE=0.0;

    Move();   // old energy get updated here
    

    //==== Update the ring triangle normal and making a copy in case the move got rejected
    std::vector<triangle *> pT=m_pvertex->GetVTraingleList();
        for (std::vector<triangle *>::iterator it = pT.begin() ; it != pT.end(); ++it)
        {
             m_Triangle.push_back(*(*it));
             (*it)->UpdateNormal_Area(m_pBox);

        }
    std::vector<links *> TpL=m_pvertex->GetVLinkList();
    for (std::vector<links *>::iterator it = TpL.begin() ; it != TpL.end(); ++it)
        {
            if(CheckFaceAngle(*it)==false)
            {
                m_face=false;
                break;
            }
        }
    //==== Update the ring link normal and shape operator
    if(m_face==true)
    {

        std::vector<links *> pL=m_pvertex->GetVLinkList();
        for (std::vector<links *>::iterator it = pL.begin() ; it != pL.end(); ++it)
        {
            links *nl=(*it)->GetNeighborLink1();
           if((*it)->m_LinkType == 0)
           {
            (*it)->UpdateNormal();
            (*it)->UpdateShapeOperator(m_pBox);
           }
           else
            (*it)->UpdateEdgeVector(m_pBox);

           if(nl->m_LinkType == 0)
           {
            nl->UpdateNormal();
            nl->UpdateShapeOperator(m_pBox);
           }
           else
            nl->UpdateEdgeVector(m_pBox);
        }
        // the PrecedingEdge should be updated directly, it is not a part of vertex link and does not have a mirror
        (m_pvertex->m_pPrecedingEdgeLink)->UpdateEdgeVector(m_pBox);

    //==== Update vertexes

        (m_pState->CurvatureCalculator())->EdgeVertexCurvature(m_pvertex);
        for (std::vector<vertex *>::iterator it = m_pAVer.begin() ; it != m_pAVer.end(); ++it)
        {
            if((*it)->m_VertexType == 0)
            (m_pState->CurvatureCalculator())->SurfVertexCurvature(*it);
            else
            (m_pState->CurvatureCalculator())->EdgeVertexCurvature(*it);

        }
    }
if(m_face==true)
{
Energy EE (m_pInc);

double NewEnergy = EE.Energy_OneVertexMove(m_pvertex);
    for (std::vector<links *>::iterator it = m_pLIntEChange.begin() ; it != m_pLIntEChange.end(); ++it)
       NewEnergy+=EE.TwoInclusionsInteractionEnergy(*it);
// if there is spring force potential so we calculate the energy associated with this move
    double harmonicde=0;
    if(m_pSPBTG->GetState()==true)
    {
        m_pSPBTG->CalculateEnergy(m_step);
        double e1 = m_pSPBTG->GetEnergy();
        Vec3D dr(m_dx,m_dy,m_dz);
        m_pSPBTG->MovingVertex(m_pvertex, dr);
        m_pSPBTG->CalculateEnergy(m_step);
        double e2 = m_pSPBTG->GetEnergy();
        harmonicde = e2-e1;
    }
    double dE_nematic_force = 0;
    if((m_pState->m_pConstant_NematicForce)->m_F0!=0)
    {
        Vec3D dx(m_dx,m_dy,m_dz);
        dE_nematic_force = (m_pState->m_pConstant_NematicForce)->Energy_of_Force(m_pvertex, dx);
    }
    
                DE=NewEnergy-m_oldEnergy;
                double diff_energy = (DE+harmonicde);
            	if(diff_energy+dE_nematic_force<=0 )
            	{
                     AccpetMove();
                    (m_pState->m_pConstant_NematicForce)->m_ActiveEnergy+=dE_nematic_force;
                    (*m_pTotEnergy)=(*m_pTotEnergy)+DE;
                    m_MoveValidity=1;

            	}
            	else if(exp(-m_Beta*(diff_energy+dE_nematic_force))>m_Thermal )
             	{
                 	AccpetMove();
                    (m_pState->m_pConstant_NematicForce)->m_ActiveEnergy+=dE_nematic_force;
                 	(*m_pTotEnergy)=(*m_pTotEnergy)+DE;
                    m_MoveValidity=1;
             	}
            	else 
            	{
                	RejectMove();
                    // Reject the harmonic potential move
                    // we should reject this move here because it does not happen at thee beginning
                    if(m_pSPBTG->GetState()==true)
                    {
                        Vec3D dr(m_dx,m_dy,m_dz);
                        m_pSPBTG->RejectMovingVertex(m_pvertex, dr);
                    }
                    //
            	}

    m_EnergyDifference=DE;
    }
    else  /// meaning that face is false
    {
        RejectMove();
    }
}

void EdgeVertexMCMove::Move()
{

	m_oldEnergy=m_pvertex->GetEnergy();

	m_vertex=(*m_pvertex);
    // Making a copy of all orginal vertices: note, we cannot use this list to replace the orginal one, all pointing pointers will go wrong, just for later use to go back to what it was
        for (std::vector<vertex *>::iterator it = m_pAVer.begin() ; it != m_pAVer.end(); ++it)
        {
		    m_AVer.push_back(*(*it));
        	m_oldEnergy+=(*it)->GetEnergy();
        }


    m_pvertex->UpdateVXPos(m_X);
	m_pvertex->UpdateVYPos(m_Y);
	m_pvertex->UpdateVZPos(m_Z);

    std::vector<links *> pL=m_pvertex->GetVLinkList();
    for (std::vector<links *>::iterator it = pL.begin() ; it != pL.end(); ++it)
    {
        links *nl=(*it)->GetNeighborLink1();
        m_nLinks.push_back(*(*it));
        m_RingLinks.push_back(*nl);
    }
//// this part of the code can become more efficient

    std::vector <links*> temlinklist;
    if(m_pvertex->m_VertexType==1)
        temlinklist.push_back(m_pvertex->m_pPrecedingEdgeLink);
    for (std::vector<vertex *>::iterator it = m_pAVer.begin() ; it != m_pAVer.end(); ++it)
    {
            if((*it)->VertexOwnInclusion()==true)
            {
                std::vector<links *> ltem=(*it)->GetVLinkList();
                temlinklist.insert(temlinklist.end(), ltem.begin(), ltem.end());
                if((*it)->m_VertexType==1)
                    temlinklist.push_back((*it)->m_pPrecedingEdgeLink);
            }
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
}
void EdgeVertexMCMove::AccpetMove()
{

    if(m_ChangeCNT==true)
    {
    
        UppdateVertexCNTCell();
    }
}

void EdgeVertexMCMove::UppdateVertexCNTCell()
{

m_OldCNT->RemoveFromVertexList(m_pvertex);
m_NewCNT->AddtoVertexList(m_pvertex);
m_pvertex->UpdateVCNTCell(m_NewCNT);

//// If the move is accepted we should update cnt cell

}


//========
//== Check Length
bool EdgeVertexMCMove::CheckDistnace()
{
    
    Vec3D Pos(m_X,m_Y,m_Z);
    // check if the nighbouring vertices are good
    std::vector<vertex *> nv= m_pvertex->GetVNeighbourVertex();
        for (std::vector<vertex *>::iterator it = m_pAVer.begin() ; it != m_pAVer.end(); ++it)
        {
            double l2 = Length2WithAVertex(Pos, *it);
            if(l2> (*m_pLmax2) || l2<(*m_pLmin2))
                return false;
        }
    // check other vertices within cnt
    m_OldCNT=m_pvertex->GetVCNTCell();
    
    bool length=true;
	if(length==true)
	{		
		bool moveoutcnt=false;
		m_OldCNT=m_pvertex->GetVCNTCell();
		Vec3D Smin=m_OldCNT->GetCNTSidemin();
		Vec3D Smax=m_OldCNT->GetCNTSidemax();
		if(m_X>=Smin(0) && m_Y>=Smin(1) && m_Z>=Smin(2) && m_X<Smax(0) && m_Y<Smax(1) && m_Z<Smax(2) )
		{
            moveoutcnt=false;

		}
		else
		{
            moveoutcnt=true;m_ChangeCNT=true;
            
           // std::cout<<"goes out of cnt \n";
		}
		if(moveoutcnt==false)
		{
			std::vector <vertex *> ver=m_OldCNT->GetVertexList();
            
           
			    for (std::vector<vertex *>::iterator it = ver.begin() ; it != ver.end(); ++it)// check the move is ok in the cell
				{
					if((*it)->GetVID()!=m_pvertex->GetVID())
						if(CheckLengthBetweenTwoVertex((*it))==false)
						{
							length=false;
							break;
						}
					
				}
				if(length==true) // check if the move is ok with the neighboring cells
			 	{
					std::vector <CNTCell *> cnts=m_OldCNT->GetVNeighbourCNTCell();
					
                    
					for (std::vector<CNTCell *>::iterator it1 = cnts.begin() ; it1 != cnts.end(); ++it1)
					{
						std::vector <vertex *> ve=(*it1)->GetVertexList();
			    			for (std::vector<vertex *>::iterator it = ve.begin() ; it != ve.end(); ++it)
						{
							
                                if((*it)->GetVID()!=m_pvertex->GetVID()) /// if you have only one cnt (same size of the box) will cause problems therefore we need this condition 
                                {
                                    if(CheckLengthBetweenTwoVertex((*it))==false)
                                    {
                                        length=false;
                                        break;
                                    }
                                }

					
						}
					}
				}
		}
		else  
		{
			//===== We should first find out where the vertex has gone?
			int cx=1;
			int cy=1;
			int cz=1;
			std::vector <CNTCell *> vcnt=m_OldCNT->GetVNeighbourCNTCell();
			if(m_X<Smin(0))
			{
				cx=0;
			}
			else if(m_X>=Smax(0))
			{
				cx=2;
			}
			if(m_Y<Smin(1))
			{
				cy=0;
			}
			else if(m_Y>=Smax(1))
			{
				cy=2;
			}
			if(m_Z<Smin(2))
			{
				cz=0;
			}
			else if(m_Z>=Smax(2))
			{
				cz=2;
			}
			int n=0;
				// find n;
			n=cx+3*cy+9*cz;

			    if(n>13)
			    {
			        n=n-1;
			    }
			m_NewCNT=vcnt.at(n);
            
            


				std::vector <vertex *> ver=m_NewCNT->GetVertexList();
			    for (std::vector<vertex *>::iterator it = ver.begin() ; it != ver.end(); ++it)
				{
					
						if(CheckLengthBetweenTwoVertex((*it))==false)
						{
							length=false;
							break;
						}
					
				}

				if(length==true)
			 	{
					std::vector <CNTCell *> cnts=m_NewCNT->GetVNeighbourCNTCell();
					
					for (std::vector<CNTCell *>::iterator it1 = cnts.begin() ; it1 != cnts.end(); ++it1)
					{
						std::vector <vertex *> ve=(*it1)->GetVertexList();
			    			for (std::vector<vertex *>::iterator it = ve.begin() ; it != ve.end(); ++it)
						{
							
							if((*it)->GetVID()!=m_pvertex->GetVID())
							if(CheckLengthBetweenTwoVertex((*it))==false)
							{
							length=false;
							break;
							}
					
						}
					}
				}

			
		}


	}
//=========================
//====== check if the faces are ok



return length;
}

void EdgeVertexMCMove::RejectMove()
{


(*m_pvertex)=m_vertex;

int i=0;
        for (std::vector<vertex>::iterator it = m_AVer.begin() ; it != m_AVer.end(); ++it)
        {

		*(m_pAVer.at(i))=(*it);
		i++;
        }

    i=0;
    std::vector<links *> pL=m_pvertex->GetVLinkList();
    for (std::vector<links *>::iterator it = pL.begin() ; it != pL.end(); ++it)
    {
        links *nl=(*it)->GetNeighborLink1();
        
        
        *(*it)=m_nLinks.at(i);
        
        if ((*it)->m_LinkType == 0)
        {
        (*it)->GetMirrorLink()->PutNormal((*it)->GetNormal());
        (*it)->GetMirrorLink()->PutShapeOperator((*it)->GetBe(),(*it)->GetHe());
        }
        
        *nl=m_RingLinks.at(i);
        if(nl->GetMirrorFlag()==true)
        {
            nl->GetMirrorLink()->PutNormal(nl->GetNormal());
            nl->GetMirrorLink()->PutShapeOperator(nl->GetBe(),nl->GetHe());
        }
        i++;
    }
    i=0;
    std::vector<triangle *> pT=m_pvertex->GetVTraingleList();
    for (std::vector<triangle *>::iterator it = pT.begin() ; it != pT.end(); ++it)
    {
        *(*it)=m_Triangle.at(i);
        
        i++;
    }
    std::vector<links>::iterator itr=m_LIntEChange.begin();
    for (std::vector<links *>::iterator it = m_pLIntEChange.begin() ; it != m_pLIntEChange.end(); ++it)
    {
            double en=(*itr).GetIntEnergy();
            (*it)->UpdateIntEnergy(en);
        if ((*it)->m_LinkType == 0)
        ((*it)->GetMirrorLink())->UpdateIntEnergy(en);
#if TEST_MODE == Enabled

        if((*itr).GetID()!=(*it)->GetID())
            std::cout<<"  "<<(*itr).GetID()<<"  "<<(*it)->GetID()<<"  These link ids should be same but it is not \n";

#endif
        ++itr;
    }
    // we need this because it is not a part of vertex links and does not have a mirror
    (m_pvertex->m_pPrecedingEdgeLink)->UpdateEdgeVector(m_pBox);

}
double EdgeVertexMCMove::Length2WithAVertex(Vec3D v1, vertex* ver)
{
    double x2=ver->GetVXPos();
    double y2=ver->GetVYPos();
    double z2=ver->GetVZPos();
    Vec3D v2(x2,y2,z2);
    Vec3D dx = v2-v1;
    for (int i=0;i<3;i++)
    {
        if(fabs(dx(i))>(*m_pBox)(i)/2.0)
        {
            if(dx(i)<0)
            dx(i)=(*m_pBox)(i)+dx(i);
            else if(dx(i)>0)
            dx(i)=dx(i)-(*m_pBox)(i);
        }
    }
    double l2=dx.dot(dx,dx);
    return l2;
}
bool EdgeVertexMCMove::CheckLengthBetweenTwoVertex( vertex* v2)
{

bool flag=true;
Vec3D Box=(*m_pBox);

 		double x2=v2->GetVXPos();
		double y2=v2->GetVYPos();
		double z2=v2->GetVZPos();
 		double x1=m_X;
		double y1=m_Y;
		double z1=m_Z;
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
        
           if(l2<(*m_pLmin2))
		  {            
            
			flag=false;
		  }




return flag;
}
bool   EdgeVertexMCMove::CheckFaceAngle(links * l)
{
    double t=0;
    Vec3D n=(l->GetTriangle())->GetNormalVector();
    if(l->GetMirrorFlag()==true)
    {
        Vec3D n3=((l->GetMirrorLink())->GetTriangle())->GetNormalVector();
        if(n.dot(n,n3)<(*m_pminAngle))
            return false;
    }
    if((l->GetNeighborLink1())->GetMirrorFlag()==true)
    {
        Vec3D n1=(((l->GetNeighborLink1())->GetMirrorLink())->GetTriangle())->GetNormalVector();
                t=n.dot(n,n1);
            if (t<(*m_pminAngle))
                return false;

    }

    return true;
}



