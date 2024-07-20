

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <fstream>
#include "Convert.h"
#include "TSIFile.h"
#include "QFile.h"
#include "WritevtuFiles.h"
#include "VMDOutput.h"

Convert::Convert()
{
}
Convert::Convert(std::vector <std::string> argument)
{

    m_Argument = argument;
    m_Healthy =true;
    m_Seed =36723;
    m_MinFaceAngle = -0.5;
    m_MinVerticesDistanceSquare = 1.0;
    m_MaxLinkLengthSquare = 3.0;
    m_OutputFilename = "out.q";
    m_InputFilename = "in.tsi";
    m_Box(0) = 10;
    m_Box(1) = 10;
    m_Box(2) = 10;
    m_Zoom(0) = 1;
    m_Zoom(1) = 1;
    m_Zoom(2) = 1;
    m_center = false;
    m_CutType = "direction";
    m_Translate(0)=0; m_Translate(1)=0; m_Translate(2)=0;
    m_NoVF = 0;
    Nfunction f;
    ExploreArguments();     // read the input data

    std::string ext = m_InputFilename.substr(m_InputFilename.find_last_of(".") + 1);
    QFile Q;
    TSIFile tsi;
    WritevtuFiles vtu;
    VMDOutput  gro;
    if(ext==TSExt)
    {
        m_BluePrint = Q.Read(m_InputFilename);
    }
    else if(ext==TSIExt)
    {
        m_BluePrint = tsi.ReadTSI(m_InputFilename);
    }
    else
    {
        std::cout<<"---> Error: input file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
    }
    //==============================================================
    
    // make the changes into the mesh
    
        std::vector<Vertex_Map> bvertex = m_BluePrint.bvertex;
        m_Box = m_BluePrint.simbox;
        ExploreArguments();     // read the input data
    
    // first we do the translation
    if(m_Translate.norm() != 0 )
    for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it)
    {
         double x = it->x + m_Translate(0);
         double y = it->y + m_Translate(1);
         double z = it->z + m_Translate(2);
            x=x*(m_Zoom(0));
            y=y*(m_Zoom(1));
            z=z*(m_Zoom(2));
        m_Box(0) = m_Box(0)*fabs(m_Zoom(0));
        m_Box(1) = m_Box(1)*fabs(m_Zoom(1));
        m_Box(2) = m_Box(2)*fabs(m_Zoom(2));

         if(x>m_Box(0))
             x = x-(int(x/m_Box(0)))*m_Box(0);
         else if(x<0)
             x= x+(int(fabs(x)/m_Box(0))+1)*m_Box(0);
        
        if(y>m_Box(1))
            y = y-(int(y/m_Box(1)))*m_Box(1);
        else if(y<0)
            y= y+(int(fabs(y)/m_Box(1))+1)*m_Box(1);
        
        if(z>m_Box(2))
            z = z-(int(z/m_Box(2)))*m_Box(2);
        else if(z<0)
            z= z+(int(fabs(z)/m_Box(2))+1)*m_Box(2);
        
        it->x = x;
        it->y = y;
        it->z = z;
        

    }
    // we do box change
    m_BluePrint.simbox = m_Box;
    
    // we do centering;
    double xcm = 0;    double ycm = 0; double zcm = 0;

    for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it)
    {
        xcm+= it->x;
        ycm+= it->y;
        zcm+= it->z;
    }
    xcm = xcm/double(bvertex.size());
    ycm = ycm/double(bvertex.size());
    zcm = zcm/double(bvertex.size());

    if(m_center == true)
    for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it)
    {
        it->x = it->x - xcm + m_Box(0)/2;
        it->y = it->y - ycm + m_Box(1)/2;
        it->z = it->z - zcm + m_Box(2)/2;
        
    }
    m_BluePrint.bvertex = bvertex;

    

    std::vector<Vertex_Map> newVlist;
    std::vector<Vertex_Map> RemVlist;
//======================== make the cut
if(m_CutType == "direction"){
    double cutsize = 0.2;
    Vec3D Dir(1,0,0);
    std::cout<<" enter the direction vector: \n";
    std::cin>>Dir(0)>>Dir(1)>>Dir(2);
    Dir = Dir*(1.0/Dir.norm());
    std::cout<<" enter the cut percentage%: \n";
    std::cin>>cutsize;
    
    double dmin = 10000;
    double dmax = -100000;

for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it){
    double x = it->x -m_Box(0)/2;
    double y = it->y -m_Box(1)/2;
    double z = it->z -m_Box(2)/2;
    Vec3D X(x,y,z);
    
    double cos = X.dot(X,Dir);
    if(cos>dmax)
        dmax = cos;
    
    if(cos<dmin)
        dmin = cos;
}
    
    double size = dmax-dmin;
    size = size*(1-cutsize);
    
for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it){
        double x = it->x -m_Box(0)/2;
        double y = it->y -m_Box(1)/2;
        double z = it->z -m_Box(2)/2;
        Vec3D X(x,y,z);
        double cos = X.dot(X,Dir);
        double f = cos - dmin;
        if(size<f)
        {
        RemVlist.push_back(*it);
        }
        else
        {
        newVlist.push_back(*it);
        }
        
    }
}
else if(m_CutType == "vertex"){

    int vid;
    std::cout<<" enter id of the vertex \n";
    std::cin>>vid;
    double x0 = bvertex[0].x;// -m_Box(0)/2;
    double y0 = bvertex[0].y; //-m_Box(1)/2;
    double z0 = bvertex[0].z;// -m_Box(2)/2;
    
    Vec3D A(1,1,1);
    std::cout<<" enter the a b c: \n";
    std::cin>>A(0)>>A(1)>>A(2);
    A = A*(1.0/A.norm());

    double R0= 8;
    std::cout<<" enter the R \n";
    std::cin>>R0;
    

    


    
    for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it){
        double x = it->x ;//-m_Box(0)/2;
        double y = it->y ;//-m_Box(1)/2;
        double z = it->z ;//-m_Box(2)/2;
        
        
        double d = A(0)*(x-x0)*(x-x0)+A(1)*(y-y0)*(y-y0)+A(2)*(z-z0)*(z-z0);
        
        if(d<R0*R0)
        {
        RemVlist.push_back(*it);
        }
        else
        {
        newVlist.push_back(*it);
        }
    }
    
}
else if(m_CutType == "center"){

    double x0 = m_Box(0)/2;
    double y0 = m_Box(1)/2;
    double z0 = m_Box(2)/2;
    
    Vec3D A(1,1,1);
    std::cout<<" enter the a b c: \n";
    std::cin>>A(0)>>A(1)>>A(2);
    A = A*(1.0/A.norm());

    double R0= 8;
    std::cout<<" enter the R \n";
    std::cin>>R0;
    

    


    
    for (std::vector<Vertex_Map>::iterator it = bvertex.begin() ; it != bvertex.end(); ++it){
        double x = it->x ;//-m_Box(0)/2;
        double y = it->y ;//-m_Box(1)/2;
        double z = it->z ;//-m_Box(2)/2;
        
        
        double d = A(0)*(x-x0)*(x-x0)+A(1)*(y-y0)*(y-y0)+A(2)*(z-z0)*(z-z0);
        
        if(d<R0*R0)
        {
        RemVlist.push_back(*it);
        }
        else
        {
        newVlist.push_back(*it);
        }
    }
    
}
//================
    std::map<int, int> newVmap;  //index to vid
    int i=0;
     for (std::vector<Vertex_Map>::iterator it = newVlist.begin() ; it != newVlist.end(); ++it)
    {
        int vid = it->id;
        newVmap.insert(std::make_pair(vid, i));
        it->id = i;
        i++;
    }


     std::vector<Triangle_Map> newTr;
     std::vector<Inclusion_Map> newInc;
     
        for (std::vector<Triangle_Map>::iterator it1 = (m_BluePrint.btriangle).begin() ; it1 != (m_BluePrint.btriangle).end(); ++it1)
    	{
    		int idv1 = it1->v1;
    		int idv2 = it1->v2;
    		int idv3 = it1->v3;
    		bool good = true;
   		for (std::vector<Vertex_Map>::iterator it = RemVlist.begin() ; it != RemVlist.end(); ++it)
   		 {
  		  	int idv = it->id;
    			if(idv==idv1 || idv==idv2||  idv==idv3 )
	    		good = false;
        	}
        	if(good==true)
        	{
        		newTr.push_back(*it1);
        	}
    	}   
    	
    	 for (std::vector<Triangle_Map>::iterator it1 = newTr.begin() ; it1 != newTr.end(); ++it1)
    	{

        		it1->v1 = (newVmap.find(it1->v1))->second;
        		it1->v2 = (newVmap.find(it1->v2))->second;
        		it1->v3 = (newVmap.find(it1->v3))->second;
    	} 

    	
     
       	for (std::vector<Inclusion_Map>::iterator it2 = (m_BluePrint.binclusion).begin() ; it2 != (m_BluePrint.binclusion).end(); ++it2)
    	{  
    	    		int idi = it2->vid;
    	    		bool good = true;
    		for (std::vector<Vertex_Map>::iterator it = RemVlist.begin() ; it != RemVlist.end(); ++it)
    		{
    			int idv = it->id;
    			if(idv==idi)
    	    		 good = false;
        
    		}     
    		if(good==true)
    		{
    		it2->vid = (newVmap.find(it2->vid))->second;
        	newInc.push_back(*it2);
        	}
    }
    
        m_BluePrint.bvertex = newVlist;
        m_BluePrint.btriangle = newTr;
        m_BluePrint.binclusion = newInc;
        m_BluePrint.number_vector_field = m_NoVF;
    //--- we generate vector fields
    std::vector<VectorField_Map> bVF;
    srand(3847);

    for (std::vector<Vertex_Map>::iterator it = newVlist.begin() ; it != newVlist.end(); ++it) {
        std::string s_vf;
        for  (int i = 0; i<m_NoVF; i++){
            double x = (double(rand() % 2000000) / 2000000.0);
            double y = (double(rand() % 2000000) / 2000000.0);
            double size = sqrt(x*x+y*y);
            x = x/size;
            y = y/size;
            s_vf += " "+ f.Int_to_String(i+1) +" "+ f.Int_to_String(x) +" "+ f.Int_to_String(y) ;
        }

        VectorField_Map tvf;
        tvf.data_line = s_vf;
        bVF.push_back(tvf);
   }
    m_BluePrint.bvectorfields  = bVF;
    
    
    
    
    
    //----
    
    
    
    std::cout<<" eveything is good here \n";
    //=======================================================
    std::string outext = m_OutputFilename.substr(m_OutputFilename.find_last_of(".") + 1);
    if(outext==TSExt)
    {
        Q.Write(m_OutputFilename, m_BluePrint);
    }
    else if(outext==TSIExt)
    {
        tsi.WriteTSI(m_OutputFilename, m_BluePrint);
    }
    else if(outext=="vtu")
    {
        vtu.Write(m_BluePrint, m_OutputFilename);
    }
    else if(outext=="gro")
    {
        gro.Write(m_OutputFilename, m_BluePrint);
    }
    else
    {
        std::cout<<"---> Error: input file with "<<ext<<" extension is not recognized. "<<std::endl;
    }

    
}
Convert::~Convert()
{
    
}
void Convert::ExploreArguments()
{
    Nfunction f;
    for (long i=1;i<m_Argument.size();i=i+2)
    {
        std::string Arg1 = m_Argument.at(i);
        if(Arg1=="-h")
        {
            HelpMessage();
            exit(0);
            break;
        }
        else if(Arg1=="-o")
        {
            m_OutputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-in")
        {
            m_InputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-zoom")
        {
            m_Zoom(0) = f.String_to_Double(m_Argument.at(i+1));
            m_Zoom(1) = f.String_to_Double(m_Argument.at(i+2));
            m_Zoom(2) = f.String_to_Double(m_Argument.at(i+3));
            i++;i++;
        }
        else if(Arg1=="-defout")
        {
            m_GeneralOutputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-angle")
        {
            m_MinFaceAngle = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-maxDist")
        {
            m_MaxLinkLengthSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-type")
        {
            m_CutType = m_Argument.at(i+1);
        }
        else if(Arg1=="-minDist")
        {
            m_MinVerticesDistanceSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-translate")
        {
            m_Translate(0) = f.String_to_Double(m_Argument.at(i+1));
            m_Translate(1) = f.String_to_Double(m_Argument.at(i+2));
            m_Translate(2) = f.String_to_Double(m_Argument.at(i+3));
            i++;i++;

        }
        else if(Arg1=="-Box")
        {
            m_Box(0) = f.String_to_Double(m_Argument.at(i+1));
            m_Box(1) = f.String_to_Double(m_Argument.at(i+2));
            m_Box(2) = f.String_to_Double(m_Argument.at(i+3));
            i++;i++;
            
        }
        else if(Arg1=="-NoVF")
        {
            m_NoVF = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-center")
        {
            m_center = true;
            i=i-1;
            
        }
        else if(Arg1=="-minDist")
        {
            m_MinVerticesDistanceSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else
        {
            std::cout << "---> Error: "<<Arg1<<" ... bad argument ";
            std::cout<<"\n"<<"For more information and tips run "<< m_Argument.at(0) <<" -h"<<"\n";
            m_Healthy =false;
            exit(0);
            break;
        }
    }

}
void Convert::HelpMessage()
{
    std::cout<<"--------------------copyright: Weria Pezeshkian------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"--------------version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"------------simple example for exacuting  -------------------"<<"\n";
    std::cout<<" ./CNV -in out.q  -o out.tsi"<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  option    type        default            description "<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  -in         string       out.q               input file name, could be tsi or q file formats "<<"\n";
    std::cout<<"  -o          string       out.tsi             output file name, could be tsi/vtu/gro or q file formats "<<"\n";
    std::cout<<"  -type       string       direction           vertex/direction "<<"\n";

    
    
    std::cout<<"=========================================================================="<<"\n";
    std::cout<<"=========================================================================="<<"\n";
    std::cout<<"------------------ version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
}
