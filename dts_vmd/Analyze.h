#if !defined(AFX_Analyze_H_7F4B21B8_D13C_9321_QF23_124095086231__INCLUDED_)
#define AFX_Analyze_H_7F4B21B8_D13C_9321_QF23_124095086231__INCLUDED_

#include "SimDef.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "State.h"
#include "xdrfile.h"
#include "xdrfile_xtc.h"

class Analyze
{
public:
    
	Analyze(State *state);
	 ~Analyze();

public:


private:
    State *m_pState;


    void Write_SURF_XTC_Frame(std::vector<vertex*> pV, XDRFILE * fout,int &step,float &time,matrix &box,float & prec);
    
    void WriteGro(std::vector<vertex*> pV,std::string filename, Vec3D Box);
    void  UpdateGeometry(MESH *pmesh, Vec3D *pBox);
    void  GreatePro(std::vector<vertex*> &pV);
    
    std::vector<vertex> m_MoreV;


};


#endif
