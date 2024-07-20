#include "OpenEdgeEvolutionWithConstantVertex_B.h"
#include "State.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian

 Edge treatment, this is a new development since June 2023;
 What it does:
 
 // 1. We need to give error if there is a edge; we cannot have osmotic pressure ...
 // 2. We can create a hole closer function; instead of invoking the making links twice 
 
 */
OpenEdgeEvolutionWithConstantVertex_B::OpenEdgeEvolutionWithConstantVertex_B()
{

}
OpenEdgeEvolutionWithConstantVertex_B::OpenEdgeEvolutionWithConstantVertex_B(int rate, State *pState)
{
    m_Rate = rate;
    m_WholeSize = 0;
    m_pState = pState;
    m_Beta =   m_pState->m_Beta;
    m_pEnergyCalculator = pState->GetEnergyCalculator();
}
OpenEdgeEvolutionWithConstantVertex_B::~OpenEdgeEvolutionWithConstantVertex_B()
{
    
}
void OpenEdgeEvolutionWithConstantVertex_B::Initialize()
{
  // creating some trinagle and links to later create from them
    m_pMESH = m_pState->m_pMesh;
    int Nt = m_pMESH->m_pEdgeV.size();
    m_pBox = m_pMESH->m_pBox;
    // due to mirro links 2*(Nt-3)+Nt
    int lid = 2*((m_pMESH->m_pMHL).size())+(m_pMESH->m_pEdgeL).size();
    for (int i=0;i<3*Nt-6;i++)
    {
        links teml(lid);
        m_GhostL.push_back(teml);
        lid++;
    }
    int tid = (m_pMESH->m_pActiveT).size() ;
    for (int i=0;i<Nt-2;i++)
    {
        triangle temt(tid);
        m_GhostT.push_back(temt);
        tid++;
    }
    for (std::vector<links>::iterator it = m_GhostL.begin() ; it != m_GhostL.end(); ++it)
        m_pGhostL.push_back(&(*it));
    
    for (std::vector<triangle>::iterator it = m_GhostT.begin() ; it != m_GhostT.end(); ++it)
        m_pGhostT.push_back(&(*it));
    
}
void OpenEdgeEvolutionWithConstantVertex_B::MC_Move(RNG* rng, double lmin, double lmax, double minangle)
{

    double createorkill = rng->UniformRNG(1.0);
    double thermal = rng->UniformRNG(1.0);
     
#if DEBUG_MODE == Enabled
    std::cout<<" mess_d: MC_Move Function  for the open edge: Starts \n";
#endif
    int NL = (m_pMESH->m_pEdgeV).size();
    int NS = (m_pMESH->m_pSurfV).size();
if(createorkill<0.5)
{
        // an atempt to create a link                       //       v3
                                                            //      /   \
                                                          //     v1-----v2
    if((m_pMESH->m_pEdgeL).size()>4)
    {
        int n = rng->IntRNG(m_pMESH->m_pEdgeV.size());
        vertex *v1 =m_pMESH->m_pEdgeV[n];

        if(Linkisvalid(v1,lmin, lmax, minangle)==false)
            return;

        //=== edge and vertices who has changed energy
        links* l2 = v1->m_pEdgeLink;
        vertex *v3 = l2->GetV2();
        links* l1 = v3->m_pEdgeLink;
        vertex *v2 = l1->GetV2();



        double eold = v1->GetEnergy();
        eold+= v2->GetEnergy();
        eold+= v3->GetEnergy();

        //=== inclusion interaction energy;
        
        std::vector <links *> nvl1 = v1->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();

        std::vector <links *> nvl2 = v2->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();

        std::vector <links *> nvl3 = v3->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();
        


        // create a link (this also updates the gemotry)
         links *newlink = CreateALink(v1);

        
        // new energy
        double enew = m_pEnergyCalculator->SingleVertexEnergy(v1);
      //  std::cout<<"\n new "<<enew<<"  ";

        enew+= m_pEnergyCalculator->SingleVertexEnergy(v2);
      //  std::cout<<" "<<m_pEnergyCalculator->SingleVertexEnergy(v2)<<"  ";

        enew+=m_pEnergyCalculator->SingleVertexEnergy(v3);
      //  std::cout<<" "<<m_pEnergyCalculator->SingleVertexEnergy(v3)<<" length  "<<v3->m_VLength<<" ";

        //=== inclusion interaction energy
        
        for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
      //  std::cout<<" "<<enew<<"  ";

        for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
      //  std::cout<<" "<<enew<<"  ";

        for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
      //  std::cout<<" "<<enew<<"  ";

            // the new created link
           enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(newlink);

      //  std::cout<<" "<<enew<<" \n ";

        double *glo_energy=&(m_pState->m_TotEnergy);
        double DE = (enew-eold);
       // std::cout<<" DE "<<DE<<" \n ";
      //  std::cout<<" edge "<<(m_pMESH->m_pEdgeV).size()<<" \n ";

        //if(double(NL)/double(NS+1)*exp(-m_Beta*DE)>thermal )
        if(exp(-m_Beta*DE)>thermal )
        {
             (*glo_energy)=(*glo_energy)+DE;
        }
        else
        {
            KillALink(newlink);
            double e = m_pEnergyCalculator->SingleVertexEnergy(v1);
            e = m_pEnergyCalculator->SingleVertexEnergy(v2);
            e = m_pEnergyCalculator->SingleVertexEnergy(v3);
            {
            for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
            for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
            for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
            }
           // std::cout<<"created rejected \n";

            
        }
        // int energy;
        //double eint = TwoInclusionsInteractionEnergy(l0);

        // we have not included the edge energy yet, in particular if there are proteins on v1 and v2; new interaction should appear.
        // also many more geo vertices should be updated
  
        
    }
    else if((m_pMESH->m_pEdgeL).size()==4 )
    {
            int n = rng->IntRNG(m_pMESH->m_pEdgeV.size());
            vertex *v1 =m_pMESH->m_pEdgeV.at(n);
        if(Linkisvalid(v1,lmin, lmax, minangle)==false)
            return;
        
        if(Anglevalid4Vhole(v1, minangle)==false)        // we should check if the two new trinagules has valid angle
            return;

        
        

        links* l2 = v1->m_pEdgeLink;
        vertex *v3 = l2->GetV2();
        links* l1 = v3->m_pEdgeLink;
        vertex *v2 = l1->GetV2();
        links* l3 = v2->m_pEdgeLink;
        vertex *v4 = l3->GetV2();
        links* l4 = v4->m_pEdgeLink;

        double eold = v1->GetEnergy();
        eold+= v2->GetEnergy();
        eold+= v3->GetEnergy();
        eold+= v4->GetEnergy();
        
        //=== inclusion interaction energy;
        
        std::vector <links *> nvl1 = v1->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();
        
        std::vector <links *> nvl2 = v2->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();
        
        std::vector <links *> nvl3 = v3->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();
        
        std::vector <links *> nvl4 = v4->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl4.begin() ; it != nvl4.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();
        
        
        
        double *glo_energy=&(m_pState->m_TotEnergy);

        
       // double testE1 = SystemEnergy();
        //std::cout<<testE1<<"  energy before "<<*glo_energy<<"\n";
        
        links *newL1 = CreateALink(v1);
        links *mnewL1 = CreateALink(v2);
        
        
        // now we merge these two links
        newL1->UpdateMirrorFlag(true);
        mnewL1->UpdateMirrorFlag(true);
        newL1->UpdateMirrorLink(mnewL1);
        mnewL1->UpdateMirrorLink(newL1);
        newL1->m_LinkType = 0;
        mnewL1->m_LinkType = 0;
        RemoveFromLinkList(newL1,m_pMESH->m_pEdgeL);
        RemoveFromLinkList(mnewL1,m_pMESH->m_pEdgeL);
        AddtoLinkList(newL1,m_pMESH->m_pHL);
        AddtoLinkList(mnewL1,m_pMESH->m_pMHL);
        
        // v1 and v2 should also become a surface vertex
        RemoveFromVertexList(v1,m_pMESH->m_pEdgeV);
        RemoveFromVertexList(v2,m_pMESH->m_pEdgeV);
        AddtoVertexList(v1,m_pMESH->m_pSurfV);
        AddtoVertexList(v2,m_pMESH->m_pSurfV);
        v1->m_VertexType = 0;
        v2->m_VertexType = 0;
        // we need a sperate create function for this; it causes probelems like below which is ugly solutions 
        v2->RemoveFromNeighbourVertex(v1);  // since this has been done twice we removed it one
        v1->RemoveFromNeighbourVertex(v2);   // since this has been done twice we removed it one
        v1->AddtoNeighbourVertex(v2);       // because remove function removes both repeated and therefore should be added again
        v2->AddtoNeighbourVertex(v1);      // // because remove function removes both repeated and therefore should be added again

        // geometry update
        // the mirror will be updated normally
        newL1->UpdateNormal();   // normal of the links with mirror should be updated
        newL1->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
              
      (m_pState->CurvatureCalculator())->SurfVertexCurvature(v1);  // v1 is still an edge vertex
      (m_pState->CurvatureCalculator())->SurfVertexCurvature(v2);  // // v2 is still an edge vertex

        
        double enew = m_pEnergyCalculator->SingleVertexEnergy(v1);
        enew+= m_pEnergyCalculator->SingleVertexEnergy(v2);
        enew+=m_pEnergyCalculator->SingleVertexEnergy(v3);
        enew+=m_pEnergyCalculator->SingleVertexEnergy(v4);
        
        //=== inclusion interaction energy
        
        for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
        for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
        for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
        for (std::vector<links *>::iterator it = nvl4.begin() ; it != nvl4.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
        
        // the new created link
       enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(newL1);
        
        double DE = (enew-eold);


        //if(double(NL)/double(NS+1)*exp(-m_Beta*DE)>thermal )
        if(exp(-m_Beta*DE)>thermal )
        {
             (*glo_energy)=(*glo_energy)+DE;
        }
        else
        {
            KillALinkOnSurf(newL1);

            double e = m_pEnergyCalculator->SingleVertexEnergy(v1);
            e = m_pEnergyCalculator->SingleVertexEnergy(v2);
            e = m_pEnergyCalculator->SingleVertexEnergy(v3);
            e = m_pEnergyCalculator->SingleVertexEnergy(v4);
            
            {
            for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
            for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
            for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
            for (std::vector<links *>::iterator it = nvl4.begin() ; it != nvl4.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
            }
         //   (*glo_energy)=(*glo_energy)+DE;

        }
    }
    else if((m_pMESH->m_pEdgeL).size()==0)
    {
        
    }
    else
    {
        std::cout<<(m_pMESH->m_pEdgeL).size()<<"\n";
        std::cout<<"error---> it should not happen, this is not expected \n";
    }
}
else  // if(createorkill<0.5)
{
        // an atempt to kill a link
    if(m_pMESH->m_pEdgeL.size()!=0 )
    {
    
            int n = rng->IntRNG(m_pMESH->m_pEdgeL.size());
            links *plink =m_pMESH->m_pEdgeL.at(n);
            vertex *v1 = plink->GetV1();
            vertex *v2 = plink->GetV2();
            vertex *v3 = plink->GetV3();

        if(v1->GetVLinkList().size()<2 || v2->GetVLinkList().size()<2 || v3->m_VertexType == 1)
            return;
        
        // calculate the old energy
        double eold = v1->GetEnergy();
        eold+= v2->GetEnergy();
        eold+= v3->GetEnergy();
        
        
        //
        {
        std::vector <links *> nvl1 = v1->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();
        
        std::vector <links *> nvl2 = v2->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();
        
        std::vector <links *> nvl3 = v3->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();
            
            // because these two links were counted two times
            eold-=2*(plink->GetNeighborLink1())->GetIntEnergy();
            eold-=2*(plink->GetNeighborLink2())->GetIntEnergy();

        }
        

        
        
        
        // we kill a link and update the geomotry
        KillALink(plink);
 
        
        // new energy
        double enew = m_pEnergyCalculator->SingleVertexEnergy(v1);
        enew+= m_pEnergyCalculator->SingleVertexEnergy(v2);
        enew+=m_pEnergyCalculator->SingleVertexEnergy(v3);
        
        
        {
        std::vector <links *> nvl1 = v1->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
            
        std::vector <links *> nvl2 = v2->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
            
        std::vector <links *> nvl3 = v3->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
        }
        
        double *glo_energy=&(m_pState->m_TotEnergy);

        double DE = (enew-eold);
               
        //if(double(NS)/double(NL+1)*(exp(-m_Beta*DE)>thermal ))
        if((exp(-m_Beta*DE)>thermal ))
        {
             (*glo_energy)=(*glo_energy)+DE;
        }
        else
        {
            CreateALink(v1);

            
            
            {
            double e = m_pEnergyCalculator->SingleVertexEnergy(v1);
            e = m_pEnergyCalculator->SingleVertexEnergy(v2);
            e = m_pEnergyCalculator->SingleVertexEnergy(v3);
            
            std::vector <links *> nvl1 = v1->GetVLinkList();
            for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
                
            std::vector <links *> nvl2 = v2->GetVLinkList();
            for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
                
            std::vector <links *> nvl3 = v3->GetVLinkList();
            for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);

            }
           // std::cout<<"killed rejected \n";

        }
        // int energy;
        //double eint = TwoInclusionsInteractionEnergy(l0);

        // we have not included the edge energy yet, in particular if there are proteins on v1 and v2; new interaction should appear.
        // also many more geo vertices should be updated
        
        
    }
}
  
#if DEBUG_MODE == Enabled
    std::cout<<" mess_d: MC_Move Function  for the open edge: Ends \n";
#endif
}
void OpenEdgeEvolutionWithConstantVertex_B::RemoveFromLinkList(links* z, std::vector<links*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void OpenEdgeEvolutionWithConstantVertex_B::RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect)
{

    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void OpenEdgeEvolutionWithConstantVertex_B::RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect)
{

    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
// this can be written as template
void OpenEdgeEvolutionWithConstantVertex_B::AddtoLinkList(links* z, std::vector<links*> &vect)
{
    vect.push_back(z);
}
void OpenEdgeEvolutionWithConstantVertex_B::AddtoVertexList(vertex* z, std::vector<vertex*> &vect)
{
    vect.push_back(z);
}
void OpenEdgeEvolutionWithConstantVertex_B::AddtoTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.push_back(z);
}
links* OpenEdgeEvolutionWithConstantVertex_B::CreateALink(vertex *v1)
{
    // an atempt to create a link                       //       v3
    // between v1 and v2                                //      /   \
                                                      //     v1-----v2

    links* l2 = v1->m_pEdgeLink;
    vertex *v3 = l2->GetV2();
    links* l1 = v3->m_pEdgeLink;
    vertex *v2 = l1->GetV2();

    links *outlink;

        // we added the triangle
        if(m_pGhostT.size()<1)
        {
            std::cout<<" error 882---> we do not have enough storage for the new trinagles \n";
            exit(0);
        }
        triangle *t = m_pGhostT.at(0);  // one trinagule is created
        t->UpdateVertex(v1,v2,v3);
        RemoveFromTriangleList(t,m_pGhostT);
        AddtoTriangleList(t,m_pMESH->m_pActiveT);
        
        // add the trinagle to the vertcies list
        v1->AddtoTraingleList(t);
        v2->AddtoTraingleList(t);
        v3->AddtoTraingleList(t);
        
        // v3 is no longer an edge vertex
        RemoveFromVertexList(v3,m_pMESH->m_pEdgeV);
        AddtoVertexList(v3,m_pMESH->m_pSurfV);
        v3->m_VertexType = 0; // meaning it is not an edge vertex anymore
        
        // we create three links, two mirror and one edge
        if(m_pGhostL.size()<3)
        {
            std::cout<<" error 883---> we do not have enough storage for the new links \n";
            exit(0);
        }
        RemoveFromLinkList(l1,m_pMESH->m_pEdgeL);
        RemoveFromLinkList(l2,m_pMESH->m_pEdgeL);
        
        links *l0 = m_pGhostL[0];
        links *ml1 = m_pGhostL[1];
        links *ml2 = m_pGhostL[2];
        l0->m_LinkType = 1; // only this needs it as the rest will have 0 in constrcuter
        RemoveFromLinkList(ml1,m_pGhostL);
        RemoveFromLinkList(ml2,m_pGhostL);
        RemoveFromLinkList(l0,m_pGhostL);
        AddtoLinkList(l0,m_pMESH->m_pEdgeL);
        AddtoLinkList(ml1,m_pMESH->m_pActiveL);
        AddtoLinkList(ml2,m_pMESH->m_pActiveL);

        v1->m_pEdgeLink = l0;
        v2->m_pPrecedingEdgeLink = l0;
        l0->UpdateV(v1,v2,v3);
        l0->UpdateTriangle(t);

        
        ml1->UpdateV(v2,v3,v1);
        ml2->UpdateV(v3,v1,v2);
        ml1->UpdateTriangle(t);
        ml2->UpdateTriangle(t);
        
        ml1->UpdateMirrorLink(l1);
        ml2->UpdateMirrorLink(l2);
        l1->UpdateMirrorLink(ml1);
        l2->UpdateMirrorLink(ml2);
        l1->m_LinkType = 0;
        l2->m_LinkType = 0;
        l1->UpdateMirrorFlag(true);
        l2->UpdateMirrorFlag(true);
        ml1->UpdateMirrorFlag(true);
        ml2->UpdateMirrorFlag(true);
        ml1->UpdateNeighborLink1(ml2);
        ml1->UpdateNeighborLink2(l0);
        ml2->UpdateNeighborLink1(l0);
        ml2->UpdateNeighborLink2(ml1);
        l0->UpdateNeighborLink1(ml1);
        l0->UpdateNeighborLink2(ml2);
        
        v1->AddtoNeighbourVertex(v2);
        v2->AddtoNeighbourVertex(v1);
        v1->AddtoLinkList(l0);
        v3->AddtoLinkList(ml2);
        v2->AddtoLinkList(ml1);

        // adding ml1 and ml2 and l1 and l2 to m_pMHL and m_pHL
        AddtoLinkList(l1,m_pMESH->m_pMHL);
        AddtoLinkList(l2,m_pMESH->m_pMHL);
        AddtoLinkList(ml1,m_pMESH->m_pHL);
        AddtoLinkList(ml2,m_pMESH->m_pHL);
        outlink = l0;
    
      // now we should update the geomtry of the affected v, l, t
    t->UpdateNormal_Area(m_pBox);   //trinagule normal and area should be found
    l1->UpdateNormal();   // normal of the links with mirror should be updated
    l2->UpdateNormal();   // normal of the links with mirror should be updated

    // their mirror will be updated by the function within
    l1->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    l2->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    l0->UpdateEdgeVector(m_pBox);   // l0 is an edge link, we need only the edge vector and length
    
    (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v1);  // v1 is still an edge vertex
    (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v2);  // // v2 is still an edge vertex
    (m_pState->CurvatureCalculator())->SurfVertexCurvature(v3);  // v3 is now a surface vertex

    return outlink;
}
void OpenEdgeEvolutionWithConstantVertex_B::KillALinkOnSurf(links *plink)
{
   // std::cout<<" this function has not been impliminted yet \n";
    links *mplink = plink->GetMirrorLink();
    triangle *t1 = plink->GetTriangle();
    triangle *t2 = mplink->GetTriangle();
    vertex *v1 = plink->GetV1();
    vertex *v2 = plink->GetV2();
    vertex *v3 = plink->GetV3();
    vertex *v4 = mplink->GetV3();
    links *l1 = plink->GetNeighborLink1();
    links *l2 = plink->GetNeighborLink2();
    links *l3 = mplink->GetNeighborLink1();
    links *l4 = mplink->GetNeighborLink2();
    links *ml1 = l1->GetMirrorLink();
    links *ml2 = l2->GetMirrorLink();
    links *ml3 = l3->GetMirrorLink();
    links *ml4 = l4->GetMirrorLink();
    
    // make the plink and mplink a ghost
    RemoveFromLinkList(plink,m_pMESH->m_pHL);           // too expensive
    RemoveFromLinkList(mplink,m_pMESH->m_pHL);          // too expensive
    RemoveFromLinkList(plink,m_pMESH->m_pMHL);          // too expensive
    RemoveFromLinkList(mplink,m_pMESH->m_pMHL);      // too expensive
    RemoveFromLinkList(plink,m_pMESH->m_pActiveL);   // too expensive
    RemoveFromLinkList(mplink,m_pMESH->m_pActiveL);  // too expensive
    AddtoLinkList(mplink,m_pGhostL);
    AddtoLinkList(plink,m_pGhostL);

    
    // adding the 4 mirrors into the edge
    AddtoLinkList(ml1,m_pMESH->m_pEdgeL);
    AddtoLinkList(ml2,m_pMESH->m_pEdgeL);
    AddtoLinkList(ml3,m_pMESH->m_pEdgeL);
    AddtoLinkList(ml4,m_pMESH->m_pEdgeL);
    RemoveFromLinkList(ml1,m_pMESH->m_pHL);  // too expensive; I should find a better way
    RemoveFromLinkList(ml2,m_pMESH->m_pHL);  // too expensive
    RemoveFromLinkList(ml1,m_pMESH->m_pMHL);  // too expensive
    RemoveFromLinkList(ml2,m_pMESH->m_pMHL);  // too expensive
    RemoveFromLinkList(ml3,m_pMESH->m_pHL);  // too expensive; I should find a better way
    RemoveFromLinkList(ml4,m_pMESH->m_pHL);  // too expensive
    RemoveFromLinkList(ml3,m_pMESH->m_pMHL);  // too expensive
    RemoveFromLinkList(ml4,m_pMESH->m_pMHL);  // too expensive
    
    // the 4 edge lose mirror
    ml1->UpdateMirrorFlag(false);
    ml2->UpdateMirrorFlag(false);
    ml3->UpdateMirrorFlag(false);
    ml4->UpdateMirrorFlag(false);
    ml1->m_LinkType = 1;
    ml2->m_LinkType = 1;
    ml3->m_LinkType = 1;
    ml4->m_LinkType = 1;
    
    // all 4 v becomes  edge vertices
    RemoveFromVertexList(v1,m_pMESH->m_pSurfV);  // too expensive
    RemoveFromVertexList(v2,m_pMESH->m_pSurfV);  // too expensive
    RemoveFromVertexList(v3,m_pMESH->m_pSurfV);  // too expensive
    RemoveFromVertexList(v4,m_pMESH->m_pSurfV);  // too expensive
    AddtoVertexList(v1,m_pMESH->m_pEdgeV);
    AddtoVertexList(v2,m_pMESH->m_pEdgeV);
    AddtoVertexList(v3,m_pMESH->m_pEdgeV);
    AddtoVertexList(v4,m_pMESH->m_pEdgeV);
    v1->m_VertexType = 1; // v1 is now an edge vertex; it used to be surf vertex
    v2->m_VertexType = 1; // v2 is now an edge vertex; it used to be surf vertex
    v3->m_VertexType = 1; // v3 is now an edge vertex; it used to be surf vertex
    v4->m_VertexType = 1; // v4 is now an edge vertex; it used to be surf vertex

    v1->m_pEdgeLink = ml2;
    v2->m_pEdgeLink = ml4;
    v3->m_pEdgeLink = ml1;
    v4->m_pEdgeLink = ml3;

    v1->m_pPrecedingEdgeLink = ml3;
    v2->m_pPrecedingEdgeLink = ml1;
    v3->m_pPrecedingEdgeLink = ml2;
    v4->m_pPrecedingEdgeLink = ml4;
    
    // now 4  links need to be removed as well
    RemoveFromLinkList(l1,m_pMESH->m_pActiveL);  // too expensive
    RemoveFromLinkList(l2,m_pMESH->m_pActiveL);  // too expensive
    RemoveFromLinkList(l3,m_pMESH->m_pActiveL);  // too expensive
    RemoveFromLinkList(l4,m_pMESH->m_pActiveL);  // too expensive
    RemoveFromLinkList(l1,m_pMESH->m_pHL);  // too expensive
    RemoveFromLinkList(l2,m_pMESH->m_pHL);  // too expensive
    RemoveFromLinkList(l3,m_pMESH->m_pHL);  // too expensive
    RemoveFromLinkList(l4,m_pMESH->m_pHL);  // too expensive
    RemoveFromLinkList(l1,m_pMESH->m_pMHL);  // too expensive
    RemoveFromLinkList(l2,m_pMESH->m_pMHL);  // too expensive
    RemoveFromLinkList(l3,m_pMESH->m_pMHL);  // too expensive
    RemoveFromLinkList(l4,m_pMESH->m_pMHL);  // too expensive
    AddtoLinkList(l1,m_pGhostL);
    AddtoLinkList(l2,m_pGhostL);
    AddtoLinkList(l3,m_pGhostL);
    AddtoLinkList(l4,m_pGhostL);
    
    // removing the edge from vertex list
     v1->RemoveFromLinkList(plink);
     v2->RemoveFromLinkList(mplink);
     v1->RemoveFromLinkList(l3);
     v2->RemoveFromLinkList(l1);
     v4->RemoveFromLinkList(l4);
     v3->RemoveFromLinkList(l2);

    // v1 and v2 should also not be connected and they are nit nighbour anymore
    v1->RemoveFromNeighbourVertex(v2);
    v2->RemoveFromNeighbourVertex(v1);
    
    // the trinagle will be killed
    RemoveFromTriangleList(t1,m_pMESH->m_pActiveT); // too expensive
    AddtoTriangleList(t1,m_pGhostT);
    RemoveFromTriangleList(t2,m_pMESH->m_pActiveT); // too expensive
    AddtoTriangleList(t2,m_pGhostT);
    // remove the trinagle from all the vertcies list
    v1->RemoveFromTraingleList(t1);
    v2->RemoveFromTraingleList(t1);
    v3->RemoveFromTraingleList(t1);
    v1->RemoveFromTraingleList(t2);
    v2->RemoveFromTraingleList(t2);
    v3->RemoveFromTraingleList(t2);
    
    // now we should update the geomtry of the affected v, l
ml1->UpdateEdgeVector(m_pBox);   // edge vector should be updated
ml2->UpdateEdgeVector(m_pBox);   // edge vector should be updated
ml3->UpdateEdgeVector(m_pBox);   // edge vector should be updated
ml4->UpdateEdgeVector(m_pBox);   // edge vector should be updated

(m_pState->CurvatureCalculator())->EdgeVertexCurvature(v1);  // v1 is an edge vertex
(m_pState->CurvatureCalculator())->EdgeVertexCurvature(v2);  // // v2 is  an edge vertex
(m_pState->CurvatureCalculator())->EdgeVertexCurvature(v3);  // v3 is now an edge vertex
(m_pState->CurvatureCalculator())->EdgeVertexCurvature(v4);  // v3 is now an edge vertex

}
void OpenEdgeEvolutionWithConstantVertex_B::KillALink(links *plink)
{
        triangle *t = plink->GetTriangle();
        vertex *v1 = plink->GetV1();
        vertex *v2 = plink->GetV2();
        vertex *v3 = plink->GetV3();
        links *l1 = plink->GetNeighborLink1();
        links *l2 = plink->GetNeighborLink2();
        links *ml1 = l1->GetMirrorLink();
        links *ml2 = l2->GetMirrorLink();
        // make the plink a ghost
        RemoveFromLinkList(plink,m_pMESH->m_pEdgeL);
        AddtoLinkList(plink,m_pGhostL);
        // adding the two mirror into the edge
        AddtoLinkList(ml1,m_pMESH->m_pEdgeL);
        AddtoLinkList(ml2,m_pMESH->m_pEdgeL);
        RemoveFromLinkList(ml1,m_pMESH->m_pHL);  // too expensive; I should find a better way
        RemoveFromLinkList(ml2,m_pMESH->m_pHL);  // too expensive
        RemoveFromLinkList(ml1,m_pMESH->m_pMHL);  // too expensive
        RemoveFromLinkList(ml2,m_pMESH->m_pMHL);  // too expensive
        
        // the two edge lose mirror
        ml1->UpdateMirrorFlag(false);
        ml2->UpdateMirrorFlag(false);
        ml1->m_LinkType = 1;
        ml2->m_LinkType = 1;
    
        // v3 becomes an edge vertex
        RemoveFromVertexList(v3,m_pMESH->m_pSurfV);  // too expensive
        AddtoVertexList(v3,m_pMESH->m_pEdgeV);
        v3->m_VertexType = 1; // v3 is not an edge vertex; it used to be surf vertex
        v3->m_pEdgeLink = ml1;
        v3->m_pPrecedingEdgeLink = ml2;
    
        // edge link of v1 and PrecedingEdgeLink of v2 should be updated
        v1->m_pEdgeLink = ml2;
        v2->m_pPrecedingEdgeLink = ml1;
        
        // now two  links need to be removed as well
        RemoveFromLinkList(l1,m_pMESH->m_pActiveL);  // too expensive
        RemoveFromLinkList(l2,m_pMESH->m_pActiveL);  // too expensive
        RemoveFromLinkList(l1,m_pMESH->m_pHL);  // too expensive
        RemoveFromLinkList(l2,m_pMESH->m_pHL);  // too expensive
        RemoveFromLinkList(l1,m_pMESH->m_pMHL);  // too expensive
        RemoveFromLinkList(l2,m_pMESH->m_pMHL);  // too expensive
        AddtoLinkList(l1,m_pGhostL);
        AddtoLinkList(l2,m_pGhostL);
        // we need to also remove them from m_pHL and m_pMHL
        

       // removing the edge from vertex list
        v1->RemoveFromLinkList(plink);
        v3->RemoveFromLinkList(l2);
        v2->RemoveFromLinkList(l1);
        // v1 and v2 should also not be connected and they are nit nighbour anymore
        v1->RemoveFromNeighbourVertex(v2);
        v2->RemoveFromNeighbourVertex(v1);

        // the trinagle will be killed
        RemoveFromTriangleList(t,m_pMESH->m_pActiveT); // too expensive
        AddtoTriangleList(t,m_pGhostT);
        // remove the trinagle from all the vertcies list
        v1->RemoveFromTraingleList(t);
        v2->RemoveFromTraingleList(t);
        v3->RemoveFromTraingleList(t);
    
    
            // now we should update the geomtry of the affected v, l
        ml1->UpdateEdgeVector(m_pBox);   // edge vector should be updated
        ml2->UpdateEdgeVector(m_pBox);   // edge vector should be updated
  
        (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v1);  // v1 is still an edge vertex
        (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v2);  // // v2 is still an edge vertex
        (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v3);  // v3 is now an edge vertex
    
}
bool OpenEdgeEvolutionWithConstantVertex_B::Linkisvalid(vertex *v1, double lmin, double lmax, double minangle)
{
    

    // check if the new link length is within the allowed range and also if the angle of the new trinagule is fine with respect to the two other trinagules
    bool isvalid = true;
    
    links* l2 = v1->m_pEdgeLink;
    vertex *v3 = l2->GetV2();
    links* l1 = v3->m_pEdgeLink;
    vertex *v2 = l1->GetV2();
    
    triangle t(0,v1,v2,v3);
    Vec3D *Box = v1->GetBox();
    t.UpdateNormal_Area(Box);
    Vec3D n1 = t.GetNormalVector();
    Vec3D n2 = (l2->GetTriangle())->GetNormalVector();
    Vec3D n3 = (l1->GetTriangle())->GetNormalVector();

    Vec3D P1(v1->GetVXPos(),v1->GetVYPos(),v1->GetVZPos());
    Vec3D P2(v2->GetVXPos(),v2->GetVYPos(),v2->GetVZPos());
    Vec3D DP = P2-P1;
    
    for (int i=0;i<3;i++)
    if(fabs(DP(i))>(*Box)(i)/2.0)
    {
        if(DP(i)<0)
            DP(i)=(*Box)(i)+DP(i);
        else if(DP(i)>0)
            DP(i)=DP(i)-(*Box)(i);
    }

    double dist2 = P1.dot(DP,DP);
    
    
    if(dist2<lmin || dist2>lmax)
        return false;
    if( n1.dot(n1,n2)<minangle)
        return false;
    if( n1.dot(n1,n3)<minangle)
        return false;
    
    return isvalid;
    
}
double  OpenEdgeEvolutionWithConstantVertex_B::SystemEnergy()
{
    std::vector<vertex *> ActiveV = m_pMESH->m_pSurfV;
    std::vector<triangle *> pActiveT = m_pMESH->m_pActiveT;
    std::vector<links *> mLink = m_pMESH->m_pHL;
    std::vector<links *>  pEdgeL = m_pMESH->m_pEdgeL;
    std::vector<vertex *> EdgeV  = m_pMESH->m_pEdgeV;
    double en = 0;
    

    for (std::vector<triangle *>::iterator it = pActiveT.begin() ; it != pActiveT.end(); ++it)
        (*it)->UpdateNormal_Area(m_pBox);
    
    
    for (std::vector<links *>::iterator it = mLink.begin() ; it != mLink.end(); ++it)
    {
        (*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
    }
    
    
    for (std::vector<vertex *>::iterator it = ActiveV.begin() ; it != ActiveV.end(); ++it)
        (m_pState->CurvatureCalculator())->SurfVertexCurvature(*it);

    //====== edge links should be updated
    for (std::vector<links *>::iterator it = pEdgeL.begin() ; it != pEdgeL.end(); ++it)
            (*it)->UpdateEdgeVector(m_pBox);

    for (std::vector<vertex *>::iterator it = EdgeV.begin() ; it != EdgeV.end(); ++it)
            (m_pState->CurvatureCalculator())->EdgeVertexCurvature(*it);
    
    en=m_pEnergyCalculator->TotalEnergy(ActiveV,mLink);
    en=en+ m_pEnergyCalculator->TotalEnergy(EdgeV,pEdgeL);
    
    
    return en;
}
bool OpenEdgeEvolutionWithConstantVertex_B::Anglevalid4Vhole(vertex *v1, double minangle)
{
    

    // Only to check when only 4 vertices in the hole to see the angle is valid
    bool isvalid = true;
    
    links* l2 = v1->m_pEdgeLink;
    vertex *v3 = l2->GetV2();
    links* l1 = v3->m_pEdgeLink;
    vertex *v2 = l1->GetV2();
    links* l3 = v2->m_pEdgeLink;
    vertex *v4 = l3->GetV2();;

    triangle t1(0,v1,v2,v3);
    triangle t2(0,v1,v4,v2);

    Vec3D *Box = v1->GetBox();
    t1.UpdateNormal_Area(Box);
    Vec3D n1 = t1.GetNormalVector();
    t2.UpdateNormal_Area(Box);
    Vec3D n2 = t2.GetNormalVector();
    

    if( n1.dot(n1,n2)<minangle)
        return false;

    
    return isvalid;
    
}
