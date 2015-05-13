#ifndef UTILITIES_H
#define UTILITIES_H
#include <functional>
#include <boost/serialization/access.hpp>
// include input and output archivers
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::exit(EXIT_FAILURE); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif


class JumpEvent{
public:
	int x, y;
	int phi;
	double prob;
	JumpEvent(){}
	JumpEvent(int x0, int y0, int phi0, double prob0){
		x = x0;
		y = y0;
		phi = phi0;	
		prob = prob0;
	}
	
};

class CoorPair{
public:
	int x;
	int y;

   friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & x;
      ar & y;    
    }


	CoorPair(){};

	CoorPair(int x0,int y0){x=x0;y=y0;}

};

typedef struct
{
	std::size_t operator() (const CoorPair & CP) const {
		std::size_t h1=std::hash<int>()(CP.x);
		std::size_t h2 = std::hash<int>()(CP.y);
		return h1^(h2<<1);
	}
}CoorPairHash;

typedef struct
{
	bool operator() (const CoorPair & CP1,const CoorPair & CP2) const {
		return (CP1.x==CP2.x)&&(CP1.y==CP2.y);
	}
}CoorPairEqual;

class Config{
	public:
	int x;
	int y;
	int phi;
	
	   friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & x;
      ar & y;
      ar & phi;    
    }
	
	Config(){};
	Config(int x0, int y0, int phi0){ x = x0; y = y0; phi = phi0; }
};

typedef struct
{
	std::size_t operator() (const Config & CP) const {
		std::size_t h1 = std::hash<int>()(CP.x);
		std::size_t h2 = std::hash<int>()(CP.y);
		std::size_t h3 = std::hash<int>()(CP.phi);
		return h1 ^ (h2 << 1) ^(h3 <<2);
	}
}ConfigHash;

typedef struct
{
	bool operator() (const Config & CP1, const Config & CP2) const {
		return (CP1.x == CP2.x) && (CP1.y == CP2.y)&&(CP1.phi == CP2.phi);
	}
}ConfigEqual;

class Edge{
public:
	int from;
	int to;
	double cost;
	int controlOpt;
	double transProb;
	
	   friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & from;
      ar & to;
      ar & cost;
      ar & controlOpt;
      ar & transProb;    
    }
	
	Edge(){};

	Edge(int v,int w,double c, int u, double trans){
	from=v;
	to=w;
	cost=c;
	controlOpt=u;
	transProb = trans;
	}
};

class GNode{
public:
	int phi;
	int index;
	CoorPair pos;
	double JFunc, JFunc_old;
	int OptControl, OptControl_old;
        int isolation;
	
	   friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & phi;
      ar & index;
      ar & pos;
      ar & JFunc;
      ar & JFunc_old;
      ar & OptControl;
      ar & OptControl_old;
      ar & isolation;
    }
	

	GNode(){};

	GNode(int x0,int y0,int phi0){
		pos.x = x0;
		pos.y = y0;
		phi=phi0;
                
	}
};

#endif
