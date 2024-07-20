#include "OpenEdgeEvolutionWithConstantVertex.h"
#include "State.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian

 Edge treatment, this is a new development since June 2023;
 
 What it does:
 
 // 1. We need to give error if there is a edge; we cannot have osmotic pressure ...
 // 2. We can create a hole closer function; instead of invoking the making links twice 
 
 */
OpenEdgeEvolutionWithConstantVertex::OpenEdgeEvolutionWithConstantVertex()
{

}
OpenEdgeEvolutionWithConstantVertex::OpenEdgeEvolutionWithConstantVertex(int rate, State *pState)
{
    m_Rate = rate;
    m_pState = pState;
    m_Beta =   m_pState->m_Beta;
    m_pEnergyCalculator = pState->GetEnergyCalculator();
}
OpenEdgeEvolutionWithConstantVertex::~OpenEdgeEvolutionWithConstantVertex(){
    
}
void OpenEdgeEvolutionWithConstantVertex::Initialize(){
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
void OpenEdgeEvolutionWithConstantVertex::MC_Move(RNG* rng, double lmin, double lmax, double minangle)
{
    double createorkill = rng->UniformRNG(1.0);
    double thermal = rng->UniformRNG(1.0);
    int NL = (m_pMESH->m_pEdgeV).size();
    int NS = (m_pMESH->m_pSurfV).size();


if(createorkill<0.5){
        // an atempt to create a link                       //       v3
                                                            //      /   \
                                                          //     v1-----v2
    if((m_pMESH->m_pEdgeL).size()>3)
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
 //----   for constant global direction type of moves
        if(v1->VertexOwnInclusion()){
            if(!(v1->GetInclusion())->UpdateGlobalDirectionFromLocal())
                return;
        }
        if(v2->VertexOwnInclusion()){
            if(!(v2->GetInclusion())->UpdateGlobalDirectionFromLocal())
                return;
        }
        if(v3->VertexOwnInclusion()){
            if(!(v3->GetInclusion())->UpdateGlobalDirectionFromLocal())
                return;
        }

        double eold = v1->GetEnergy();
        eold+= v2->GetEnergy();
        eold+= v3->GetEnergy();

        //=== inclusion interaction energy; it is multipled by two because only half is assinged to each edge
        std::vector <links *> nvl1 = v1->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();

        std::vector <links *> nvl2 = v2->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();

        std::vector <links *> nvl3 = v3->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
            eold+=2*(*it)->GetIntEnergy();

        eold+=2*(v1->m_pPrecedingEdgeLink)->GetIntEnergy();

        // create a link (this also updates the gemotry)
         links *newlink = CreateALink(v1);

        //----   for constant global direction type of moves
        if (v1->VertexOwnInclusion() || v2->VertexOwnInclusion() || v3->VertexOwnInclusion()) {
            Vec3D LD1, LD2, LD3;
            if(v1->VertexOwnInclusion()){
                LD1 = (v1->GetInclusion())->GetGDirection();
                LD1 = (v1->GetG2LTransferMatrix())*LD1;
                if(LD1.isbad()){
                    KillALink(newlink);
                    return;
                }
                LD1(2) = 0;
                LD1 = LD1*(1/LD1.norm());
            }
            if(v2->VertexOwnInclusion()){
                LD2 = (v2->GetInclusion())->GetGDirection();
                LD2 = (v2->GetG2LTransferMatrix())*LD2;
                if(LD2.isbad()){
                    //== reject the move
                    KillALink(newlink);
                    return;
                }
                LD2(2) = 0;
                LD2 = LD2*(1/LD2.norm());
            }
            if(v3->VertexOwnInclusion()){
                LD3 = (v3->GetInclusion())->GetGDirection();
                LD3 = (v3->GetG2LTransferMatrix())*LD3;
                if(LD3.isbad()){
                    //== reject the move
                    KillALink(newlink);
                    return;
                }
                LD3(2) = 0;
                LD3 = LD3*(1/LD3.norm());
            }
//-- this should happen at the end, otherwise this early links get bad number
            if(v1->VertexOwnInclusion())
            (v1->GetInclusion())->UpdateLocalDirection(LD1);
            if(v2->VertexOwnInclusion())
            (v2->GetInclusion())->UpdateLocalDirection(LD2);
            if(v3->VertexOwnInclusion())
            (v3->GetInclusion())->UpdateLocalDirection(LD3);
        }
        
        // new energy
        double enew = m_pEnergyCalculator->SingleVertexEnergy(v1);
        enew+= m_pEnergyCalculator->SingleVertexEnergy(v2);
        enew+=m_pEnergyCalculator->SingleVertexEnergy(v3);

        //=== inclusion interaction energy
        for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
        for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
        for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);

            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(v1->m_pPrecedingEdgeLink);

            // the new created link
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(newlink);

        double *glo_energy=&(m_pState->m_TotEnergy);
        double DE = (enew-eold);

        //if(double(NL)/double(NS+1)*exp(-m_Beta*DE)>thermal )
        if(exp(-m_Beta*DE)>thermal )
        {
             (*glo_energy)=(*glo_energy)+DE;
        }
        else
        {
            KillALink(newlink);
                if(v1->VertexOwnInclusion()){
                    Tensor2  G2L = v1->GetG2LTransferMatrix();
                    Vec3D LD = G2L*((v1->GetInclusion())->GetGDirection());
                    if(fabs(LD(2))>0.00001){
                        std::cout<<" something is wrong here, this should not happen. vector should be on the plane \n";
                    }
                    (v1->GetInclusion())->UpdateLocalDirection(LD);
                }
                if(v2->VertexOwnInclusion()){
                    Tensor2  G2L = v2->GetG2LTransferMatrix();
                    Vec3D LD = G2L*((v2->GetInclusion())->GetGDirection());
                    if(fabs(LD(2))>0.00001){
                        std::cout<<" something is wrong here, this should not happen. vector should be on the plane \n";
                    }
                    (v2->GetInclusion())->UpdateLocalDirection(LD);
                }
                if(v3->VertexOwnInclusion()){
                    Tensor2  G2L = v3->GetG2LTransferMatrix();
                    Vec3D LD = G2L*((v3->GetInclusion())->GetGDirection());
                    if(fabs(LD(2))>0.00001){
                        std::cout<<" something is wrong here, this should not happen. vector should be on the plane \n";
                    }
                    (v3->GetInclusion())->UpdateLocalDirection(LD);
                }
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
            e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(v1->m_pPrecedingEdgeLink);

        }
    }
    else if((m_pMESH->m_pEdgeL).size()==3){
        vertex *v1 =m_pMESH->m_pEdgeV[0];
        links* l1 = v1->m_pEdgeLink;
        vertex *v2 = l1->GetV2();
        links* l2 = v2->m_pEdgeLink;
        vertex *v3 = l2->GetV2();
        links* l3 = v3->m_pEdgeLink;

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

        double *glo_energy=&(m_pState->m_TotEnergy);

        //--- to have the inclsuion direction constant in 3d
            if(v1->VertexOwnInclusion()){
                if(!(v1->GetInclusion())->UpdateGlobalDirectionFromLocal())
                    return;
            }
            if(v2->VertexOwnInclusion()){
                if(!(v2->GetInclusion())->UpdateGlobalDirectionFromLocal())
                    return;
            }
            if(v3->VertexOwnInclusion()){
                if(!(v3->GetInclusion())->UpdateGlobalDirectionFromLocal())
                    return;
            }
  
        //==== creating the triangles
        triangle *T = CloseATriangleHole(v1);

        if (v1->VertexOwnInclusion() || v2->VertexOwnInclusion() || v3->VertexOwnInclusion()) {
            Vec3D LD1, LD2, LD3;
            if(v1->VertexOwnInclusion()){
                LD1 = (v1->GetInclusion())->GetGDirection();
                LD1 = (v1->GetG2LTransferMatrix())*LD1;
                if(LD1.isbad()){
                    KillATriangle(l1);
                    return;
                }
                LD1(2) = 0;
                LD1 = LD1*(1/LD1.norm());
            }
            if(v2->VertexOwnInclusion()){
                LD2 = (v2->GetInclusion())->GetGDirection();
                LD2 = (v2->GetG2LTransferMatrix())*LD2;
                if(LD2.isbad()){
                    KillATriangle(l1);
                    return;
                }
                LD2(2) = 0;
                LD2 = LD2*(1/LD2.norm());
            }
            if(v3->VertexOwnInclusion()){
                LD3 = (v3->GetInclusion())->GetGDirection();
                LD3 = (v3->GetG2LTransferMatrix())*LD3;
                if(LD3.isbad()){
                    KillATriangle(l1);
                    return;
                }
                LD3(2) = 0;
                LD3 = LD3*(1/LD3.norm());
            }
//-- this should happen at the end, otherwise this early links get bad number
            if(v1->VertexOwnInclusion())
            (v1->GetInclusion())->UpdateLocalDirection(LD1);
            if(v2->VertexOwnInclusion())
            (v2->GetInclusion())->UpdateLocalDirection(LD2);
            if(v3->VertexOwnInclusion())
            (v3->GetInclusion())->UpdateLocalDirection(LD3);
        }

        
        double enew = m_pEnergyCalculator->SingleVertexEnergy(v1);
        enew+= m_pEnergyCalculator->SingleVertexEnergy(v2);
        enew+=m_pEnergyCalculator->SingleVertexEnergy(v3);
        
        //=== inclusion interaction energy
        for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
        for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
        for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
            enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);

        double DE = (enew-eold);

        if(exp(-m_Beta*DE)>thermal )
        {
            (*glo_energy)=(*glo_energy)+DE;
        }
        else{
            KillATriangle(l1);
            if(v1->VertexOwnInclusion()){
                Tensor2  G2L = v1->GetG2LTransferMatrix();
                Vec3D LD = G2L*((v1->GetInclusion())->GetGDirection());
                if(fabs(LD(2))>0.00001){
                    std::cout<<"--> error: something is wrong here, this should not happen. vector should be on the plane \n";
                }
                (v1->GetInclusion())->UpdateLocalDirection(LD);
            }
            if(v2->VertexOwnInclusion()){
                Tensor2  G2L = v2->GetG2LTransferMatrix();
                Vec3D LD = G2L*((v2->GetInclusion())->GetGDirection());
                if(fabs(LD(2))>0.00001){
                    std::cout<<"--> error: something is wrong here, this should not happen. vector should be on the plane \n";
                }
                (v2->GetInclusion())->UpdateLocalDirection(LD);
            }
            if(v3->VertexOwnInclusion()){
                Tensor2  G2L = v3->GetG2LTransferMatrix();
                Vec3D LD = G2L*((v3->GetInclusion())->GetGDirection());
                if(fabs(LD(2))>0.00001){
                    std::cout<<"--> error: something is wrong here, this should not happen. vector should be on the plane \n";
                }
                (v3->GetInclusion())->UpdateLocalDirection(LD);
            }
            
            double e = m_pEnergyCalculator->SingleVertexEnergy(v1);
            e = m_pEnergyCalculator->SingleVertexEnergy(v2);
            e = m_pEnergyCalculator->SingleVertexEnergy(v3);
            
        
            for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
            for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
            for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
                e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(*it);
        }
    }
    else if((m_pMESH->m_pEdgeL).size()==0){
        return;
    }
    else if((m_pMESH->m_pEdgeL).size()!=3){
        std::cout<<(m_pMESH->m_pEdgeL).size()<<"\n";
        std::cout<<"error---> it should not happen, this is not expected \n";
    }
    return;
}
else// if (createorkill>0.5)
{
    // an atempt to kill a link
    if(m_pMESH->m_pEdgeL.size()!=0 ){
    
            int n = rng->IntRNG(m_pMESH->m_pEdgeL.size());
            links *plink =m_pMESH->m_pEdgeL[n];
    
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
            
            //--- this link does not exist in the nb of any vertices
            eold+=2*(v1->m_pPrecedingEdgeLink)->GetIntEnergy();

            // because these two links were counted two times
            eold-=2*(plink->GetNeighborLink1())->GetIntEnergy();
            eold-=2*(plink->GetNeighborLink2())->GetIntEnergy();
        }
        //== this is a new piece to use global inclusion direction as a constant variable during the move: note global direction should be projected on the new local plane
        if(v1->VertexOwnInclusion()){
            if(!(v1->GetInclusion())->UpdateGlobalDirectionFromLocal())
                return;
        }
        if(v2->VertexOwnInclusion()){
            if(!(v2->GetInclusion())->UpdateGlobalDirectionFromLocal())
                return;
        }
        if(v3->VertexOwnInclusion()){
            if(!(v3->GetInclusion())->UpdateGlobalDirectionFromLocal())
                return;
        }

        // we kill a link and update the geomotry
        KillALink(plink);

        // updating the local inc direction from the global one
        if (v1->VertexOwnInclusion() || v2->VertexOwnInclusion() || v3->VertexOwnInclusion()) {
            Vec3D LD1, LD2, LD3;
            if(v1->VertexOwnInclusion()){
                LD1 = (v1->GetInclusion())->GetGDirection();
                LD1 = (v1->GetG2LTransferMatrix())*LD1;
                if(LD1.isbad()){
                    CreateALink(v1);
                    return;
                }
                LD1(2) = 0;
                LD1 = LD1*(1/LD1.norm());
            }
            if(v2->VertexOwnInclusion()){
                LD2 = (v2->GetInclusion())->GetGDirection();
                LD2 = (v2->GetG2LTransferMatrix())*LD2;
                if(LD2.isbad()){
                    CreateALink(v1);
                    return;
                }
                LD2(2) = 0;
                LD2 = LD2*(1/LD2.norm());
            }
            if(v3->VertexOwnInclusion()){
                LD3 = (v3->GetInclusion())->GetGDirection();
                LD3 = (v3->GetG2LTransferMatrix())*LD3;
                if(LD3.isbad()){
                    CreateALink(v1);
                    return;
                }
                LD3(2) = 0;
                LD3 = LD3*(1/LD3.norm());
            }
//-- this should happen at the end, otherwise this early links get bad number
            if(v1->VertexOwnInclusion())
            (v1->GetInclusion())->UpdateLocalDirection(LD1);
            if(v2->VertexOwnInclusion())
            (v2->GetInclusion())->UpdateLocalDirection(LD2);
            if(v3->VertexOwnInclusion())
            (v3->GetInclusion())->UpdateLocalDirection(LD3);
        }
 
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
        enew+=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(v1->m_pPrecedingEdgeLink);

        
        double *glo_energy=&(m_pState->m_TotEnergy);
        double DE = (enew-eold);
        //if(double(NS)/double(NL+1)*(exp(-m_Beta*DE)>thermal ))
        if(exp(-m_Beta*DE)>thermal) {
             (*glo_energy)=(*glo_energy)+DE;
        }
        else{       // reject the move
                    CreateALink(v1);
                if(v1->VertexOwnInclusion())
                {
                    Tensor2  G2L = v1->GetG2LTransferMatrix();
                    Vec3D LD = G2L*((v1->GetInclusion())->GetGDirection());
                    if(fabs(LD(2))>0.00001){
                        std::cout<<" something is wrong here, this should not happen. vector should be on the plane \n";
                    }
                    (v1->GetInclusion())->UpdateLocalDirection(LD);
                }
                if(v2->VertexOwnInclusion())
                {
                    Tensor2  G2L = v2->GetG2LTransferMatrix();
                    Vec3D LD = G2L*((v2->GetInclusion())->GetGDirection());
                    if(fabs(LD(2))>0.00001){
                        std::cout<<" something is wrong here, this should not happen. vector should be on the plane \n";
                    }
                    (v2->GetInclusion())->UpdateLocalDirection(LD);
                }
                if(v3->VertexOwnInclusion())
                {
                    Tensor2  G2L = v3->GetG2LTransferMatrix();
                    Vec3D LD = G2L*((v3->GetInclusion())->GetGDirection());
                    if(fabs(LD(2))>0.00001){
                        std::cout<<" something is wrong here, this should not happen. vector should be on the plane \n";
                    }
                    (v3->GetInclusion())->UpdateLocalDirection(LD);
                }
            
            
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
            
            e=m_pEnergyCalculator->TwoInclusionsInteractionEnergy(v1->m_pPrecedingEdgeLink);
            }
           // std::cout<<"killed rejected \n";

        }
        // int energy;
        //double eint = TwoInclusionsInteractionEnergy(l0);
    }
}
  
#if DEBUG_MODE == Enabled
    std::cout<<" mess_d: MC_Move Function  for the open edge: Ends \n";
#endif
}
void OpenEdgeEvolutionWithConstantVertex::RemoveFromLinkList(links* z, std::vector<links*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void OpenEdgeEvolutionWithConstantVertex::RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect)
{

    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void OpenEdgeEvolutionWithConstantVertex::RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect)
{

    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
// this can be written as template
void OpenEdgeEvolutionWithConstantVertex::AddtoLinkList(links* z, std::vector<links*> &vect)
{
    vect.push_back(z);
}
void OpenEdgeEvolutionWithConstantVertex::AddtoVertexList(vertex* z, std::vector<vertex*> &vect)
{
    vect.push_back(z);
}
void OpenEdgeEvolutionWithConstantVertex::AddtoTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.push_back(z);
}
// an atempt to create a link                       //        v3                      v3
// between v1 and v2                                // ml3  /   \  ml2---> ml3,l3  // T \\  l2, ml2   ml3(v1,v3)  ml2(v3,v2)
                                                  //     v1      v2               v1---->v2
links* OpenEdgeEvolutionWithConstantVertex::CreateALink(vertex *v1){           //     l1
// l1 will be created and
    if(m_pGhostT.size()<1 && m_pGhostL.size()<3){
        std::cout<<" error 882---> we do not have enough storage for the new trinagle and links \n";
        exit(0);
    }
    
    links* ml3 = v1->m_pEdgeLink;
    vertex *v3 = ml3->GetV2();
    links* ml2 = v3->m_pEdgeLink;
    vertex *v2 = ml2->GetV2();
    links *l1 = m_pGhostL[0];
    links *l2 = m_pGhostL[1];
    links *l3 = m_pGhostL[2];
    triangle *tre = m_pGhostT[0];  // one trinagule is created

//---- create the triangle
    tre->UpdateVertex(v1,v2,v3);
    m_pGhostT.erase(m_pGhostT.begin());
    AddtoTriangleList(tre,m_pMESH->m_pActiveT);
    l1->UpdateTriangle(tre);
    l2->UpdateTriangle(tre);
    l3->UpdateTriangle(tre);
    v1->AddtoTraingleList(tre);
    v2->AddtoTraingleList(tre);
    v3->AddtoTraingleList(tre);
    
//----- update the links
      //-- update l1 and l2 that have became a surface edge
    RemoveFromLinkList(ml2,m_pMESH->m_pEdgeL);
    RemoveFromLinkList(ml3,m_pMESH->m_pEdgeL);
    AddtoLinkList(ml2,m_pMESH->m_pHL);
    AddtoLinkList(ml3,m_pMESH->m_pHL);
    ml2->m_LinkType = 0;
    ml3->m_LinkType = 0;
    l2->m_LinkType = 0;
    l3->m_LinkType = 0;
    l1->m_LinkType = 1;
    ml2->UpdateMirrorFlag(true);
    ml3->UpdateMirrorFlag(true);
    l2->UpdateMirrorFlag(true);
    l3->UpdateMirrorFlag(true);
    l1->UpdateMirrorFlag(false);
    ml2->UpdateMirrorLink(l2);
    ml3->UpdateMirrorLink(l3);
    l2->UpdateMirrorLink(ml2);
    l3->UpdateMirrorLink(ml3);
    AddtoLinkList(l2,m_pMESH->m_pMHL);   // since we just added ml1 to m_pHL
    AddtoLinkList(l3,m_pMESH->m_pMHL);
    m_pGhostL.erase(m_pGhostL.begin());   // three times for l1,l2,l1
    m_pGhostL.erase(m_pGhostL.begin());   // three times for l1,l2,l1
    m_pGhostL.erase(m_pGhostL.begin());    // three times for l1,l2,l1
    AddtoLinkList(l1,m_pMESH->m_pEdgeL);
    AddtoLinkList(l1,m_pMESH->m_pActiveL);
    AddtoLinkList(l2,m_pMESH->m_pActiveL);
    AddtoLinkList(l3,m_pMESH->m_pActiveL);

    l1->UpdateNeighborLink1(l2);
    l1->UpdateNeighborLink2(l3);
    l2->UpdateNeighborLink1(l3);
    l2->UpdateNeighborLink2(l1);
    l3->UpdateNeighborLink1(l1);
    l3->UpdateNeighborLink2(l2);
    l1->UpdateV(v1,v2,v3);
    l2->UpdateV(v2,v3,v1);
    l3->UpdateV(v3,v1,v2);
//----update vertices
    //----update v3 that has become a surf
    RemoveFromVertexList(v3,m_pMESH->m_pEdgeV);   // only this vertex becames a surf
    AddtoVertexList(v3,m_pMESH->m_pSurfV);
    v3->m_VertexType = 0;                          // meaning it is not an edge vertex anymore
    //----update other vertices
    v1->m_pEdgeLink = l1;
    v2->m_pPrecedingEdgeLink = l1;
    v1->AddtoNeighbourVertex(v2);
    v2->AddtoNeighbourVertex(v1);
    v1->AddtoLinkList(l1);
    v2->AddtoLinkList(l2);
    v3->AddtoLinkList(l3);
    
 //----- now we should update the geomtry of the affected v, l, t
    tre->UpdateNormal_Area(m_pBox);   //trinagule normal and area should be found
    l2->UpdateNormal();   // normal of the links with mirror should be updated
    l3->UpdateNormal();   // normal of the links with mirror should be updated

    // their mirror will be updated by the function within
    l1->UpdateEdgeVector(m_pBox);   // l1 is an edge link, we need only the edge vector and length
    l2->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    l3->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    
    (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v1);  // v1 is still an edge vertex
    (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v2);  // // v2 is still an edge vertex
    (m_pState->CurvatureCalculator())->SurfVertexCurvature(v3);  // v3 is now a surface vertex

    return l1;
}

//         v3                    v2
//  l3  // T1 \\ l2  -->  ml3  /   \  ml2
//      v1-l1- v2             v1    v3
bool OpenEdgeEvolutionWithConstantVertex::KillALink(links *l1)
{
        triangle *tri = l1->GetTriangle();
        vertex *v1 = l1->GetV1();
        vertex *v2 = l1->GetV2();
        vertex *v3 = l1->GetV3();
        links *l2 = l1->GetNeighborLink1();
        links *l3 = l2->GetNeighborLink1();
        links *ml2 = l2->GetMirrorLink();
        links *ml3 = l3->GetMirrorLink();
    
//--- remove the triangle
    RemoveFromTriangleList(tri,m_pMESH->m_pActiveT);
    AddtoTriangleList(tri,m_pGhostT);
          //--- remove the trinagle from all the vertcies list
    v1->RemoveFromTraingleList(tri);
    v2->RemoveFromTraingleList(tri);
    v3->RemoveFromTraingleList(tri);
    
//---- remove the links l1, l2 and l3 ; note the mirors remain allive but will become an edge link
    RemoveFromLinkList(l1,m_pMESH->m_pActiveL);  // too expensive
    RemoveFromLinkList(l2,m_pMESH->m_pActiveL);  // too expensive
    RemoveFromLinkList(l3,m_pMESH->m_pActiveL);  // too expensive
    RemoveFromLinkList(l2,m_pMESH->m_pHL);  // too expensive
    RemoveFromLinkList(l3,m_pMESH->m_pHL);  // too expensive
    RemoveFromLinkList(l2,m_pMESH->m_pMHL);  // too expensive
    RemoveFromLinkList(l3,m_pMESH->m_pMHL);  // too expensive
    RemoveFromLinkList(l1,m_pMESH->m_pEdgeL);    // this link does not exist in m_pMHL and m_pHL
    AddtoLinkList(l1,m_pGhostL);
    AddtoLinkList(l2,m_pGhostL);
    AddtoLinkList(l3,m_pGhostL);
    v1->RemoveFromLinkList(l1);
    v2->RemoveFromLinkList(l2);
    v3->RemoveFromLinkList(l3);
    
//--- convert the links into edge links ml2 and ml3
    RemoveFromLinkList(ml2,m_pMESH->m_pHL);  // too expensive; I should find a better way
    RemoveFromLinkList(ml3,m_pMESH->m_pHL);  // too expensive
    RemoveFromLinkList(ml2,m_pMESH->m_pMHL);  // too expensive
    RemoveFromLinkList(ml3,m_pMESH->m_pMHL);  // too expensive
    ml2->UpdateMirrorFlag(false);
    ml3->UpdateMirrorFlag(false);
    ml2->m_LinkType = 1;
    ml3->m_LinkType = 1;
    AddtoLinkList(ml2,m_pMESH->m_pEdgeL);    // adding the two mirror into the edge
    AddtoLinkList(ml3,m_pMESH->m_pEdgeL);   // adding the two mirror into the edge

//--- convert the three vertices into edge vertex,
        //--- only v2 needs an update: v1 and v3 are already an edge vertex.
    RemoveFromVertexList(v3,m_pMESH->m_pSurfV);  // too expensive
    AddtoVertexList(v3,m_pMESH->m_pEdgeV);
    v3->m_VertexType = 1; //
       //-- v1 and v2 are no longer connected
    v1->RemoveFromNeighbourVertex(v2);
    v2->RemoveFromNeighbourVertex(v1);
       // now we updated the vertex edge links and preceding link; note: v3 edge links will not be changed and also v1  preceding link will not be changed
        v1->m_pEdgeLink = ml3;
        v3->m_pEdgeLink = ml2;
        v3->m_pPrecedingEdgeLink = ml3;
        v2->m_pPrecedingEdgeLink = ml2 ;
//--- now we should update the geomtry of the affected v, l
        ml2->UpdateEdgeVector(m_pBox);   // edge vector should be updated
        ml3->UpdateEdgeVector(m_pBox);   // edge vector should be updated
  
        (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v1);  // v1 is still an edge vertex
        (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v2);  // // v2 is still an edge vertex
        (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v3);  // v3 is now an edge vertex
    
    return true;
    
}

// an atempt to create a triangle                    //       v3              v3
// using v1                                          // l3   /   \ l2  -->  // T \\
                                                  //       v1---v2        v1 ====v2
                                                    //        l1
triangle* OpenEdgeEvolutionWithConstantVertex::CloseATriangleHole(vertex *v1)
{


    links* l1 = v1->m_pEdgeLink;
    vertex *v2 = l1->GetV2();
    links* l2 = v2->m_pEdgeLink;
    vertex *v3 = l2->GetV2();
    links* l3 = v3->m_pEdgeLink;

        // we added the triangle
    if(m_pGhostT.size()<1 || m_pGhostL.size()<3){
        std::cout<<" error 882---> we do not have enough storage for the new trinagles or links \n";
        exit(0);
    }

//-- new triangle and links that need to be created.
    links *ml1 = m_pGhostL[0];
    links *ml2 = m_pGhostL[1];
    links *ml3 = m_pGhostL[2];
    triangle *outtriangle = m_pGhostT[0];
//-- build the new triangle
    m_pGhostT.erase(m_pGhostT.begin());
    AddtoTriangleList(outtriangle,m_pMESH->m_pActiveT);
    // this trinagle is made of v1, v2 and v3 but aniwise
    outtriangle->UpdateVertex(v2,v1,v3);
    // add the trinagle to the vertcies list
    v1->AddtoTraingleList(outtriangle);
    v2->AddtoTraingleList(outtriangle);
    v3->AddtoTraingleList(outtriangle);
    // note, the exsisting links belong to different trinagles. ml1, ml2, ml3 making this trinagle
    ml1->UpdateTriangle(outtriangle);
    ml2->UpdateTriangle(outtriangle);
    ml3->UpdateTriangle(outtriangle);
//-- build the ml links
    //-- removing from ghost containers
    m_pGhostL.erase(m_pGhostL.begin());  // rm ml1
    m_pGhostL.erase(m_pGhostL.begin());  // rm ml2
    m_pGhostL.erase(m_pGhostL.begin());  // rm ml3
    //--add ml links to active
    AddtoLinkList(ml1,m_pMESH->m_pActiveL);
    AddtoLinkList(ml2,m_pMESH->m_pActiveL);
    AddtoLinkList(ml3,m_pMESH->m_pActiveL);
    //--- to mhl
    AddtoLinkList(ml1,m_pMESH->m_pMHL);  // note: this is true because l1 is added to hl
    AddtoLinkList(ml2,m_pMESH->m_pMHL);
    AddtoLinkList(ml3,m_pMESH->m_pMHL);
    //-- next link
    ml1->UpdateNeighborLink1(ml3);
    ml1->UpdateNeighborLink2(ml2);
    ml2->UpdateNeighborLink1(ml1);
    ml2->UpdateNeighborLink2(ml3);
    ml3->UpdateNeighborLink1(ml2);
    ml3->UpdateNeighborLink2(ml1);
    //-- update their vertex
    ml1->UpdateV(v2,v1,v3);
    ml2->UpdateV(v3,v2,v1);
    ml3->UpdateV(v1,v3,v2);
    //== adding the m links to vertices list
    v1->AddtoLinkList(ml3);
    v3->AddtoLinkList(ml2);
    v2->AddtoLinkList(ml1);
    ml1->m_LinkType = 0; //
    ml2->m_LinkType = 0; //
    ml3->m_LinkType = 0; //
    ml1->UpdateMirrorFlag(true);
    ml2->UpdateMirrorFlag(true);
    ml3->UpdateMirrorFlag(true);
    ml1->UpdateMirrorLink(l1);
    ml2->UpdateMirrorLink(l2);
    ml3->UpdateMirrorLink(l3);
//-- update l1,l2,l3
    l1->m_LinkType = 0; //
    l2->m_LinkType = 0; //
    l3->m_LinkType = 0; //
    l1->UpdateMirrorFlag(true);
    l2->UpdateMirrorFlag(true);
    l3->UpdateMirrorFlag(true);
    l1->UpdateMirrorLink(ml1);
    l2->UpdateMirrorLink(ml2);
    l3->UpdateMirrorLink(ml3);
    //  l1 and l2 to  m_pHL
    AddtoLinkList(l1,m_pMESH->m_pHL);  //note  ml1  ml2 and ml3 are in mhl
    AddtoLinkList(l2,m_pMESH->m_pHL);
    AddtoLinkList(l3,m_pMESH->m_pHL);
    RemoveFromLinkList(l1,m_pMESH->m_pEdgeL);     // they are part of normal links now
    RemoveFromLinkList(l2,m_pMESH->m_pEdgeL);
    RemoveFromLinkList(l3,m_pMESH->m_pEdgeL);
//== uodate vertices
    v1->m_VertexType = 0; // meaning it is not an edge vertex anymore
    v2->m_VertexType = 0; // meaning it is not an edge vertex anymore
    v3->m_VertexType = 0; // meaning it is not an edge vertex anymore
    //-- removing from vedge containers
    RemoveFromVertexList(v1,m_pMESH->m_pEdgeV);   // v1, v2 , v3 is no longer an edge vertex
    RemoveFromVertexList(v2,m_pMESH->m_pEdgeV);
    RemoveFromVertexList(v3,m_pMESH->m_pEdgeV);
    //--- adding to the containers
    AddtoVertexList(v1,m_pMESH->m_pSurfV);
    AddtoVertexList(v2,m_pMESH->m_pSurfV);
    AddtoVertexList(v3,m_pMESH->m_pSurfV);

// updating the geomtry
      // triangle area and normal
    outtriangle->UpdateNormal_Area(m_pBox);   //trinagule normal and area should be found
    ml1->UpdateNormal();   // normal of the links with mirror should be updated
    ml2->UpdateNormal();   // normal of the links with mirror should be updated
    ml3->UpdateNormal();   // normal of the links with mirror should be updated

    // their mirror will be updated by the function within
    ml1->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    ml2->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    ml3->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    
    (m_pState->CurvatureCalculator())->SurfVertexCurvature(v1);  // v1 is still an edge vertex
    (m_pState->CurvatureCalculator())->SurfVertexCurvature(v2);  // // v2 is still an edge vertex
    (m_pState->CurvatureCalculator())->SurfVertexCurvature(v3);  // v3 is now a surface vertex


    return outtriangle;
}
// the triangle should belong to l1 (l1 get killed at the end)
// ---------------------------                      //               v3               v3
// between v1 and v2                                //     l3     //  \\ l2 --->     /   \
                                                  //             v1 ====v2          v1---v2
                                                         //          l1

// this function kills target_triangle and the assosciated links l1,l2,l3. It send them to the ghost area
bool OpenEdgeEvolutionWithConstantVertex::KillATriangle(links *l1){
//=== note: the mirror will be killed not l1,
    // this are the objects that will be killed
    links *ml1 = l1->GetMirrorLink();
    links *ml3 = ml1->GetNeighborLink1();   //    in the reverse function ml1->UpdateNeighborLink1(ml3);
    links *ml2 = ml3->GetNeighborLink1();   //    in the reverse function ml3->UpdateNeighborLink1(ml2);
    triangle *target_triangle = ml1->GetTriangle();

    //-- other objects
    links *l2 = ml2->GetMirrorLink();
    links *l3 = ml3->GetMirrorLink();
    vertex *v1 = ml1->GetV2();
    vertex *v2 = ml1->GetV1();      //  in the reverse function ml2->UpdateV(v3,v2,v1);
    vertex *v3 = ml1->GetV3();

//-- remove the triangle, only from active t and verices. Note the links that have this triangle will die anyway
    AddtoTriangleList(target_triangle,m_pGhostT);
    RemoveFromTriangleList(target_triangle,m_pMESH->m_pActiveT);
    v1->RemoveFromTraingleList(target_triangle);
    v2->RemoveFromTraingleList(target_triangle);
    v3->RemoveFromTraingleList(target_triangle);
    
//--- remove the links ml1, ml2, ml3 links
    RemoveFromLinkList(ml1,m_pMESH->m_pActiveL);
    RemoveFromLinkList(ml2,m_pMESH->m_pActiveL);
    RemoveFromLinkList(ml3,m_pMESH->m_pActiveL);
    RemoveFromLinkList(ml1,m_pMESH->m_pMHL);
    RemoveFromLinkList(ml2,m_pMESH->m_pMHL);
    RemoveFromLinkList(ml3,m_pMESH->m_pMHL);
    RemoveFromLinkList(ml1,m_pMESH->m_pHL);  // we know that they were added to mhl but for more generality
    RemoveFromLinkList(ml2,m_pMESH->m_pHL);
    RemoveFromLinkList(ml3,m_pMESH->m_pHL);
    AddtoLinkList(ml1,m_pGhostL);
    AddtoLinkList(ml2,m_pGhostL);
    AddtoLinkList(ml3,m_pGhostL);
// update l1,l2 and l3 links
    RemoveFromLinkList(l1,m_pMESH->m_pHL);  // this link also has became an edge links, should be removed form the list
    RemoveFromLinkList(l2,m_pMESH->m_pHL); // this link also has became an edge links, should be removed form the list
    RemoveFromLinkList(l3,m_pMESH->m_pHL); // this link also has became an edge links, should be removed form the list
    RemoveFromLinkList(l1,m_pMESH->m_pMHL); // we know they were added to phl but for generality
    RemoveFromLinkList(l2,m_pMESH->m_pMHL); // this link also has became an edge links, should be removed form the list
    RemoveFromLinkList(l3,m_pMESH->m_pMHL); // this link also has became an edge links, should be removed form the list
    AddtoLinkList(l1,m_pMESH->m_pEdgeL);    // so they became an edge link
    AddtoLinkList(l2,m_pMESH->m_pEdgeL);
    AddtoLinkList(l3,m_pMESH->m_pEdgeL);
    l1->UpdateMirrorFlag(false);
    l2->UpdateMirrorFlag(false);
    l3->UpdateMirrorFlag(false);
    l1->m_LinkType = 1;
    l2->m_LinkType = 1;
    l3->m_LinkType = 1;
//-- convert the vertices into edge
    RemoveFromVertexList(v1,m_pMESH->m_pSurfV);
    RemoveFromVertexList(v2,m_pMESH->m_pSurfV);
    RemoveFromVertexList(v3,m_pMESH->m_pSurfV);
    AddtoVertexList(v1,m_pMESH->m_pEdgeV);
    AddtoVertexList(v2,m_pMESH->m_pEdgeV);
    AddtoVertexList(v3,m_pMESH->m_pEdgeV);
    v1->m_VertexType = 1;
    v2->m_VertexType = 1;
    v3->m_VertexType = 1;
    
    v1->RemoveFromLinkList(ml3);    // in the counter function    v1->AddtoLinkList(ml3);
    v2->RemoveFromLinkList(ml1);    // in the counter function    v2->AddtoLinkList(ml1);
    v3->RemoveFromLinkList(ml2);    // in the counter function    v3->AddtoLinkList(ml2);

    v1->m_pEdgeLink = l1;    // for here is not needed (since v1 was an edge vector). But in a general case, v1 may not have an edge vector
    v2->m_pEdgeLink = l2;
    v3->m_pEdgeLink = l3;
    v1->m_pPrecedingEdgeLink = l3;
    v2->m_pPrecedingEdgeLink = l1;
    v3->m_pPrecedingEdgeLink = l2;
      // now we should update the geomtry of the affected v1,v2,v3, ml1,ml2,ml3,
    l1->UpdateEdgeVector(m_pBox);   // this is an edge link, we need only the edge vector and length
    l2->UpdateEdgeVector(m_pBox);   // this is an edge link, we need only the edge vector and length
    l3->UpdateEdgeVector(m_pBox);   // this is an edge link, we need only the edge vector and length
    
    (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v1);  // v1 is still an edge vertex
    (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v2);  // // v2 is still an edge vertex
    (m_pState->CurvatureCalculator())->EdgeVertexCurvature(v3);  // v3 is now a surface vertex

    return true;
}
double  OpenEdgeEvolutionWithConstantVertex::SystemEnergy()
{
    std::vector<vertex *> ActiveV = m_pMESH->m_pSurfV;
    std::vector<triangle *> pActiveT = m_pMESH->m_pActiveT;
    std::vector<links *> mLink = m_pMESH->m_pHL;
    std::vector<links *>  pEdgeL = m_pMESH->m_pEdgeL;
    std::vector<vertex *> EdgeV  = m_pMESH->m_pEdgeV;
    double en = 0;
    

    for (std::vector<triangle *>::iterator it = pActiveT.begin() ; it != pActiveT.end(); ++it)
        (*it)->UpdateNormal_Area(m_pBox);
    
    
    for (std::vector<links *>::iterator it = (m_pMESH->m_pHL).begin() ; it != (m_pMESH->m_pHL).end(); ++it){
        (*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
    }
    
    for (std::vector<vertex *>::iterator it = (m_pMESH->m_pSurfV).begin() ; it != (m_pMESH->m_pSurfV).end(); ++it)
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
//           va -- v3----vb
//            \n2 /n1\n3/
//             v1--- v2
//                DP
bool OpenEdgeEvolutionWithConstantVertex::Linkisvalid(vertex *v1, double lmin, double lmax, double minangle)
{
    // check if the new link length is within the allowed range and also if the angle of the new trinagule is fine with respect to the two other trinagules
    links* l2 = v1->m_pEdgeLink;
    vertex *v3 = l2->GetV2();
    links* l1 = v3->m_pEdgeLink;
    vertex *v2 = l1->GetV2();
    Vec3D *Box = v1->GetBox();

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
    if(dist2<lmin || dist2>lmax){
        return false;
    }
    
    triangle t(0,v1,v2,v3);
    t.UpdateNormal_Area(Box);
    Vec3D n1 = t.GetNormalVector();
    Vec3D n2 = (l2->GetTriangle())->GetNormalVector();
    Vec3D n3 = (l1->GetTriangle())->GetNormalVector();
    
    if( n1.dot(n1,n2)<minangle){
        return false;
    }
    if( n1.dot(n1,n3)<minangle){
        return false;
    }
    
    return true;
    
}
