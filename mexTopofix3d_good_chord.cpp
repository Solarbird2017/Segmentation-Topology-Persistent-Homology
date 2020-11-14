//mex CC=gcc-4.3 LD=g++-4.3 -I/usr/local/Cellar/boost/1.66.0/include/ -O mexTopofix3d.cpp
#define MAX_PERS_PTS 1000	//maximum of numbers of persistence pts
#define BIG_INT 0xFFFFFFF

#define drawinmerge_comp
//#define Draw_criPoint_On_perturbM
//#define Draw_criPoint_On_critM





//========== Functions needs  =============
#define MATLAB
#ifdef MATLAB
#include "mex.h"
#endif

#include <exception>
#include <math.h>
#include <queue>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
//add the BGL ****************************************
#include <utility>                   // for std::pair
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/depth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
using namespace boost;
//add the BGL ****************************************

using namespace std;

#ifdef MATLAB
extern void _main();

// Function declarations.
// -----------------------------------------------------------------
double getMatlabScalar (const mxArray* ptr) {

  // Make sure the input argument is a scalar in double-precision.
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != 1)
    mexErrMsgTxt("The input argument must be a double-precision scalar");

  return *mxGetPr(ptr);
}
double& createMatlabScalar (mxArray*& ptr) { 
  ptr = mxCreateDoubleMatrix(1,1,mxREAL);
  return *mxGetPr(ptr);
}

#define LOG_FILE "pers_MRF_2D_log.txt"
void myMessage (const string msg,bool showtime){
//	mexWarnMsgTxt(msg.c_str());
//	mexPrintf("%s\n",msg.c_str());
	
	 time_t now;
	 time(&now);

	fstream filestr;
	filestr.open (LOG_FILE, fstream::in | fstream::out | fstream::ate);
	
	if (showtime){
	  filestr << ctime(&now) << "-------" << msg.c_str() << endl;
	}else{
	  filestr << msg.c_str() << endl;
	}
	filestr.close();

}
#else
void myMessage (const string msg,bool showtime){
	time_t now;
	 time(&now);
	if (showtime){
	  cout << ctime(&now) << "-------" << msg.c_str() << endl;
	}else{
	  cout<<msg.c_str()<<endl;
	}
}
#endif

#define  OUTPUT_MSG(MSG)  {stringstream * ss=new stringstream(stringstream::in | stringstream::out); *ss << MSG << ends; myMessage(ss->str(),true); delete ss; }

#define  OUTPUT_NOTIME_MSG(MSG)  {stringstream * ss=new stringstream(stringstream::in | stringstream::out); *ss << MSG << ends; myMessage(ss->str(),false); delete ss; }

// Function definitions.
// -----------------------------------------------------------------
// myassert will print message if bool a is false.
void myassert(bool a, int source=0){
	if (!a){
	       	int b=5;
		if(source==0){
			OUTPUT_MSG("ASSERTION FAILED!!!!!!!");
		}else if(source==2){
			b=6;
			OUTPUT_MSG("ASSERTION FAILED!!---- Due to hole--killing!!!!!!!!-----captured");
		}else{
			OUTPUT_MSG("ASSERTION FAILED!!---- Due to hole--killing!!!!!!!!");
		}
	}
//	assert(a);
	return;
}


// -----------------------  xudong: need to modify the 2D myDoubleMatrix to 3D ------------------------------------//
//class myDoubleMatrix {
//	public:
//		int nrow;
//		int ncol;
//		vector< vector< double > > data;    // data is 2d vector.
//    
//  myDoubleMatrix(int m,int n,double v=0.0) { 
//	  nrow=m;
//	  ncol=n;
//	  int i,j;
//	  for (i=0;i<nrow;i++)
//		  data.push_back(vector< double >(ncol,v));
//	  return;
//  } 
//  double get(int r, int c) { 
//	 	myassert((0<=r)&&(r<nrow));
//	 	myassert((0<=c)&&(c<ncol));
//	     return data[r][c];
//  }
//  void put(int r, int c, double v) { 
//	 	myassert((0<=r)&&(r<nrow));
//	 	myassert((0<=c)&&(c<ncol));
//	        data[r][c]=v;
//  }
//    
//    // give input *ptr to data.
//  void input1D(double *ptr){
//	  int i,j;
//	  for (i=0;i<nrow;i++)
//		  for (j=0;j<ncol;j++)
//			  data[i][j]=ptr[j*nrow+i]; // input *ptr is counted along column, so have such code.
//  }
//  void output1D(double *ptr){
//	  int i,j;
//	  for (i=0;i<nrow;i++)
//		  for (j=0;j<ncol;j++)
//			  ptr[j*nrow+i]=data[i][j];
//  }
//};
void print3dmatrix(vector< vector< vector< double > > > data){
//    cout<<"2. x dimension: "<<data.size()<<", y dimension: "<<data[0].size()<<", z dimension: "<<data[0][0].size()<<endl;
    for(int x=0; x<data.size();x++){
        for(int y=0; y<data[0].size();y++){
            for(int z=0; z<data[0][0].size();z++){
                cout<<data[x][y][z]<<"  ";
            }
            cout<<endl;
        }
        cout<<"-------------- the above is slice: "<<x<<" ---------------------"<<endl;
    }
}

//xudong add the following code.
class my3dMatrix {
public:
    int nrow;
    int ncol;
    int nslice;
    vector< vector< vector< double > > > data;    // data is 2d vector.
    
    my3dMatrix(int m,int n, int l, double v=0.0) {
        nrow=m;
        ncol=n;
        nslice=l;
        for(int i=0;i<nslice;i++){
            data.push_back(vector < vector < double > >(nrow,vector<double>(ncol,v)));   //xudong need to push_back 2d vectors.
//            return;   //return function here lead to dimension of data is 1X2X3 which causes the program crash.
        }
//        cout<<"the 3d matrix after executing the constructor:"<<endl;
//        print3dmatrix(data);
        return;
    }
//    double get(int r, int c, int l) {
//        myassert((0<=r)&&(r<nrow));
//        myassert((0<=c)&&(c<ncol));
//        return data[r][c][l];
//    }
//    void put(int r, int c, int l, double v) {
//        myassert((0<=r)&&(r<nrow));
//        myassert((0<=c)&&(c<ncol));
//        data[r][c][l]=v;
//    }
    void clearTozero(){
        int i,j,k;
        for(k=0;k<nslice;k++){
            for (i=0;i<nrow;i++)
                for (j=0;j<ncol;j++){
                    data[k][i][j]=0;
                }
        }
    }
    
    // give input *ptr to data.
    void input1D(double *ptr){
//        cout<<"1. x dimension: "<<data.size()<<", y dimension: "<<data[0].size()<<", z dimension: "<<data[0][0].size()<<endl;
        int i,j,k;
        for(k=0;k<nslice;k++){
            for (i=0;i<nrow;i++)
                for (j=0;j<ncol;j++){
                    data[k][i][j]=ptr[k*ncol*nrow + j*nrow + i]; // input *ptr is counted along column, so have such code.
//                    print3dmatrix(data);
                }
        }
//        cout<<"------------ 1. The following are data after input ------------"<<endl;  //ok here.
//        print3dmatrix(data);         // xudong: ok till here.
//        cout<<"------------ 1. The above are data after input ------------"<<endl<<endl;  //ok here.
    }
    void output1D(double *ptr){
//        cout<<"------------ 2. The following are data after entering output1D function ------------"<<endl;  //???
//        print3dmatrix(data);
//        cout<<"------------ 2. The above are data after entering output1D function  ------------"<<endl<<endl;  //???
        int i,j,k;
        for(k=0;k<nslice;k++){
            for (i=0;i<nrow;i++)
                for (j=0;j<ncol;j++)
                    ptr[k*ncol*nrow + j*nrow + i]=data[k][i][j];
        }

    
    }
    void copy3d(my3dMatrix *ptr){
//        cout<<"ptr->data.size(): "<<ptr->data.size()<<endl;
//        cout<<"ptr->data[0].size(): "<<ptr->data[0].size()<<endl;
//        cout<<"ptr->data[0][0].size(): "<<ptr->data[0][0].size()<<endl;
        if(nslice != ptr->data.size())
            cout<<"z dimension are different"<<endl;
        if(nrow != ptr->data[0].size())
            cout<<"x dimension are different"<<endl;
        if(ncol != ptr->data[0][0].size())
            cout<<"y dimension are different"<<endl;
        
        int i,j,k;
        for(k=0;k<nslice;k++){
            for (i=0;i<nrow;i++)
                for (j=0;j<ncol;j++){
                    this->data[k][i][j]=ptr->data[k][i][j];
                }
        }
        
        
    }
    void copy3d_doNotCopyPointWithIntensity0(my3dMatrix *ptr){
        //        cout<<"ptr->data.size(): "<<ptr->data.size()<<endl;
        //        cout<<"ptr->data[0].size(): "<<ptr->data[0].size()<<endl;
        //        cout<<"ptr->data[0][0].size(): "<<ptr->data[0][0].size()<<endl;
        if(nslice != ptr->data.size())
            cout<<"z dimension are different"<<endl;
        if(nrow != ptr->data[0].size())
            cout<<"x dimension are different"<<endl;
        if(ncol != ptr->data[0][0].size())
            cout<<"y dimension are different"<<endl;
        
        int i,j,k;
        for(k=0;k<nslice;k++){
            for (i=0;i<nrow;i++)
                for (j=0;j<ncol;j++){
                    if(0 == ptr->data[k][i][j])
                        continue;
                    this->data[k][i][j] = ptr->data[k][i][j];
                }
        }
        
        
    }
};




class my4dMatrix {
public:
    int nrow;
    int ncol;
    int nslice;
    int nstack;
    vector< vector< vector< vector< double > > > > data4d;    // data is 2d vector.
    
//    my4dMatrix(int m=0,int n=0, int l=0, int z=0, double v=0.0) {
//        nrow=m;
//        ncol=n;
//        nslice=l;
//        nstack = z;
//        
//        ////    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ xudong need to modify the following code ~~~~~~~~~~~~~~~~~~
//        for(int i=0;i<nstack;i++){
//            data.push_back(vector< vector < vector < double > > >(nslice, vector < vector < double > >(nrow,vector<double>(ncol,v))));   //xudong need to push_back 2d vectors.
//            //            return;   //return function here lead to dimension of data is 1X2X3 which causes the program crash.
//        }
//        //        cout<<"the 3d matrix after executing the constructor:"<<endl;
//        //        print3dmatrix(data);
//        return;
//    }
    void output1D(double *ptr){
        nrow = data4d[0][0].size();
        ncol = data4d[0][0][0].size();
        nslice = data4d[0].size();
        nstack = data4d.size();
        int i,j,k,h;
        for(h=0;h<nstack;h++){
            for(k=0;k<nslice;k++){
                for (i=0;i<nrow;i++){
                    for (j=0;j<ncol;j++){
                        ptr[h*nslice*ncol*nrow + k*ncol*nrow + j*nrow + i]=data4d[h][k][i][j];
//                        if(data4d[h][k][i][j] != 0){cout<<"Find element value isn't 0 !!!"<<endl;}
//                        cout<<"data4d[h][k][i][j]: "<<data4d[h][k][i][j]<<endl;
                    }
                }
            }
        }
    }

};



class mVpM_2dMatrix {
public:
    int nrow;
    int ncol;
    vector< vector< int > > coord;
    vector< int > range = vector< int >(6,0);   //the x,y,z range.
    vector<double> centroid = vector< double >(3,0);
    
    mVpM_2dMatrix(int m,int n, int v = -1) {
        nrow=m;
        ncol=n;
        int i,j;
        for (i=0;i<nrow;i++)
            coord.push_back(vector< int >(ncol,v));
        return;
    }
    
    double get(int r, int c) {
        myassert((0<=r)&&(r<nrow));
        myassert((0<=c)&&(c<ncol));
        return coord[r][c];
    }
    void put(int r, int c, double v) {
        myassert((0<=r)&&(r<nrow));
        myassert((0<=c)&&(c<ncol));
        coord[r][c]=v;
    }
    
    // give input *ptr to data.
    void input1D(double *ptr){
        int i,j;
        double temp_x=0, temp_y=0, temp_z=0;
        for (i=0;i<nrow;i++){
            for (j=0;j<ncol;j++){
                coord[i][j]=ptr[j*nrow+i]; // input *ptr is counted along column, so have such code.
            }
        }

        range[0] = coord[0][0];
        range[1] = coord[0][0];
        range[2] = coord[0][1];
        range[3] = coord[0][1];
        range[4] = coord[0][2];
        range[5] = coord[0][2];
        for(i=1;i<coord.size();i++){
            if(range[0] > coord[i][0]){range[0]=coord[i][0];}
            if(range[1] < coord[i][0]){range[1]=coord[i][0];}
            if(range[2] > coord[i][1]){range[2]=coord[i][1];}
            if(range[3] < coord[i][1]){range[3]=coord[i][1];}
            if(range[4] > coord[i][2]){range[4]=coord[i][2];}
            if(range[5] < coord[i][2]){range[5]=coord[i][2];}
        }
        
//        cout<<"The range of chordRange (x_min, x_max, y_min, y_max, z_min, z_max): "<<endl;
//        for (i=0;i<range.size();i++){
//            cout<<range[i]<<", ";
//        }cout<<endl;
        
        
        //calculate the centroid as followings.
        for(i=0;i<coord.size();i++){
            temp_x += coord[i][0];
            temp_y += coord[i][1];
            temp_z += coord[i][2];
        }
//        cout<<"temp_x: "<<temp_x<<", temp_y: "<<temp_y<<", temp_z: "<<temp_z<<endl;
        centroid[0] = (double)temp_x/coord.size();
        centroid[1] = (double)temp_y/coord.size();
        centroid[2] = (double)temp_z/coord.size();
//        cout<<"The coord.size() is: "<<coord.size()<<endl;
//        cout<<"The centroid is: "<<centroid[0]<<", "<<centroid[1]<<", "<<centroid[2]<<endl;
    }
    void output1D(double *ptr){
        int i,j;
        for (i=0;i<nrow;i++)
            for (j=0;j<ncol;j++)
                ptr[j*nrow+i]=coord[i][j];
    }
};



class Vertex{   //xudong modify the Vertex to 3D case.
	public:
    int zidx;   //xudong add this line.
	int xidx;
	int yidx;
    int mapZCoord;  //xudong add this line.
	int mapXCoord;  //mapXCoord is the X coordinate in cell order which is equal to 2*xi.
	int mapYCoord;  ////mapXCoord is the Y coordinate in cell order which is equal to 2*yi.
    
	Vertex(int zi, int xi,int yi) : zidx(zi),xidx(xi),yidx(yi),mapZCoord(zi*2),mapXCoord(xi*2),mapYCoord(yi*2){}    //xudong need to check the mapZCoord?????!??.
	Vertex() : zidx(-1),xidx(-1),yidx(-1),mapZCoord(-1),mapXCoord(-1),mapYCoord(-1){}
};

class Edge{
	public:
	int v1_order;
	int v2_order;
	int mapXCoord;
	int mapYCoord;
    int mapZCoord;  //xudong add this line.
	Edge(int v1o,int v2o, int mz, int mx, int my) : v1_order(min(v1o,v2o)),v2_order(max(v1o,v2o)),mapZCoord(mz),mapXCoord(mx),mapYCoord(my){}
	Edge() : v1_order(-1),v2_order(-1),mapZCoord(-1),mapXCoord(-1),mapYCoord(-1){}
};

class Triangle{
	public:
	int v1_order;
	int v2_order;
	int v3_order;
	int v4_order;
	int e1_order;
	int e2_order;
	int e3_order;
	int e4_order;
	int mapXCoord;
	int mapYCoord;
    int mapZCoord;  //xudong add this line.
    
	Triangle(int v1o,int v2o,int v3o,int v4o,int e1o,int e2o,int e3o,int e4o,int mz,int mx,int my):
	mapZCoord(mz),mapXCoord(mx),mapYCoord(my){
		vector< int >tmpvec (4,0);
		tmpvec[0]=v1o;
		tmpvec[1]=v2o;
		tmpvec[2]=v3o;
		tmpvec[3]=v4o;
		sort(tmpvec.begin(),tmpvec.end());
		v1_order=tmpvec[0];
		v2_order=tmpvec[1];
		v3_order=tmpvec[2];
		v4_order=tmpvec[3];

		vector< int >tmpEdgevec (4,0);
		tmpEdgevec[0]=e1o;
		tmpEdgevec[1]=e2o;
		tmpEdgevec[2]=e3o;
		tmpEdgevec[3]=e4o;
		sort(tmpEdgevec.begin(),tmpEdgevec.end());
		e1_order=tmpEdgevec[0];
		e2_order=tmpEdgevec[1];
		e3_order=tmpEdgevec[2];
		e4_order=tmpEdgevec[3];
	}
	Triangle() : v1_order(-1),v2_order(-1),v3_order(-1),v4_order(-1),mapZCoord(-1),mapXCoord(-1),mapYCoord(-1){}
};

// -----------------------  xudong add class cube{}  ------------------------------------//
class Cube{ //Cube has 8 vertices, 12 edges, 6 faces.
    public:
        int v1_order;
        int v2_order;
        int v3_order;
        int v4_order;
        int v5_order;
        int v6_order;
        int v7_order;
        int v8_order;
 
        int e1_order;
        int e2_order;
        int e3_order;
        int e4_order;
        int e5_order;
        int e6_order;
        int e7_order;
        int e8_order;
        int e9_order;
        int e10_order;
        int e11_order;
        int e12_order;
 
        int t1_order;
        int t2_order;
        int t3_order;
        int t4_order;
        int t5_order;
        int t6_order;
 
        int mapXCoord;
        int mapYCoord;
        int mapZCoord;  //xudong add this line.
    
        Cube(int v1o,int v2o,int v3o,int v4o,int v5o,int v6o,int v7o,int v8o,int e1o,int e2o,int e3o,int e4o,int e5o,int e6o,int e7o,int e8o,int e9o,int e10o,int e11o,int e12o, int t1o, int t2o, int t3o, int t4o, int t5o, int t6o, int mz,int mx,int my): mapZCoord(mz),mapXCoord(mx),mapYCoord(my){
        vector< int >tmpvec (8,0);
        tmpvec[0]=v1o;
        tmpvec[1]=v2o;
        tmpvec[2]=v3o;
        tmpvec[3]=v4o;
        tmpvec[4]=v5o;
        tmpvec[5]=v6o;
        tmpvec[6]=v7o;
        tmpvec[7]=v8o;
        sort(tmpvec.begin(),tmpvec.end());
        v1_order=tmpvec[0];
        v2_order=tmpvec[1];
        v3_order=tmpvec[2];
        v4_order=tmpvec[3];
        v5_order=tmpvec[4];
        v6_order=tmpvec[5];
        v7_order=tmpvec[6];
        v8_order=tmpvec[7];
 
        vector< int >tmpEdgevec (12,0);
        tmpEdgevec[0]=e1o;
        tmpEdgevec[1]=e2o;
        tmpEdgevec[2]=e3o;
        tmpEdgevec[3]=e4o;
        tmpEdgevec[4]=e5o;
        tmpEdgevec[5]=e6o;
        tmpEdgevec[6]=e7o;
        tmpEdgevec[7]=e8o;
        tmpEdgevec[8]=e9o;
        tmpEdgevec[9]=e10o;
        tmpEdgevec[10]=e11o;
        tmpEdgevec[11]=e12o;
        sort(tmpEdgevec.begin(),tmpEdgevec.end());
        e1_order=tmpEdgevec[0];
        e2_order=tmpEdgevec[1];
        e3_order=tmpEdgevec[2];
        e4_order=tmpEdgevec[3];
        e5_order=tmpEdgevec[4];
        e6_order=tmpEdgevec[5];
        e7_order=tmpEdgevec[6];
        e8_order=tmpEdgevec[7];
        e9_order=tmpEdgevec[8];
        e10_order=tmpEdgevec[9];
        e11_order=tmpEdgevec[10];
        e12_order=tmpEdgevec[11];
            
        vector< int >tmpTrivec (6,0);
        tmpTrivec[0]=t1o;
        tmpTrivec[1]=t2o;
        tmpTrivec[2]=t3o;
        tmpTrivec[3]=t4o;
        tmpTrivec[4]=t5o;
        tmpTrivec[5]=t6o;
        sort(tmpTrivec.begin(), tmpTrivec.end());
        t1_order=tmpTrivec[0];
        t2_order=tmpTrivec[1];
        t3_order=tmpTrivec[2];
        t4_order=tmpTrivec[3];
        t5_order=tmpTrivec[4];
        t6_order=tmpTrivec[5];
        }
    
        Cube() : v1_order(-1),v2_order(-1),v3_order(-1),v4_order(-1),v5_order(-1),v6_order(-1),v7_order(-1),v8_order(-1),mapZCoord(-1),mapXCoord(-1),mapYCoord(-1){}
 
};

//myDoubleMatrix * phi;
my3dMatrix * phi3d; //xudong add this line.

bool vCompVal(Vertex a, Vertex b){  // sort vertices accoding to intensity of each point.
    return phi3d->data[a.zidx][a.xidx][a.yidx] < phi3d->data[b.zidx][b.xidx][b.yidx];   // ascending order. xudong: needs to change phi to phi3d.
}

bool eComp(Edge a, Edge b){ // sort edge according to the order ot v2, then v1.
	if(a.v2_order!=b.v2_order)
		return a.v2_order<b.v2_order;
	if(a.v1_order!=b.v1_order)
		return a.v1_order<b.v1_order;
}

bool trigComp(Triangle a, Triangle b){ 
	if(a.e4_order!=b.e4_order)
		return a.e4_order<b.e4_order;
	if(a.e3_order!=b.e3_order)
		return a.e3_order<b.e3_order;
	if(a.e2_order!=b.e2_order)
		return a.e2_order<b.e2_order;
	if(a.e1_order!=b.e1_order)
		return a.e1_order<b.e1_order;
}
//xudong add the following code.
bool CubeComp(Cube a, Cube b){
    if(a.t6_order!=b.t6_order)
        return a.t6_order<b.t6_order;
    if(a.t5_order!=b.t5_order)
        return a.t5_order<b.t5_order;
    if(a.t4_order!=b.t4_order)
        return a.t4_order<b.t4_order;
    if(a.t3_order!=b.t3_order)
        return a.t3_order<b.t3_order;
    if(a.t2_order!=b.t2_order)
        return a.t2_order<b.t2_order;
    if(a.t1_order!=b.t1_order)
        return a.t1_order<b.t1_order;
}
//xudong add the code above.


template <class Container>
struct Counter : public std::iterator <std::output_iterator_tag, void, void, void, void>{
	size_t &cnt;

    Counter(size_t &x) : cnt(x) {}	
 
	template<typename t>
    Counter& operator=(t)
	{        
        return *this;
    }
    
    Counter& operator* () 
	{
        return *this;
    }
    
    Counter& operator++(int) 
	{
		++cnt;
        return *this;
    }    

	Counter& operator++() 
	{
		++cnt;
        return *this;
    }    
};

// We avoid excessive allocations by calculating the size of the resulting list.
// Then we resize the result and populate it with the actual values.
vector< int > list_sym_diff(vector< int > &sa, vector< int > &sb){  //return a vector.
//assume inputs are both sorted increasingly
	size_t count = 0;
//    cout<<count<<"  --> 1"<<endl;
	Counter< vector< int > > counter(count);    //counter(0).
//    cout<<count<<"  --> 2"<<endl;
	set_symmetric_difference(sa.begin(), sa.end(), sb.begin(), sb.end(), counter);  ///求集合A，B的对称差（即(A-B)并(B-A))
    
//    cout<<count<<"  --> 3"<<endl;
	vector< int > out;
	out.reserve(count); //pre-allocate count memory space.
	set_symmetric_difference(sa.begin(), sa.end(), sb.begin(), sb.end(), back_inserter(out));   ///求集合A，B的对称差（即(A-B)并(B-A))

	return out;
}


//-----------------------------------------------------
//vertex-edge pair and persistence
//edge-trig pair and persistence
//Trig-cube pair and persistence    //xudong add this line.
//-----------------------------------------------------
class VertEdgePair{
  public:
  int vbidx;
  int edidx;
  double robustness;
  double birth;
  double death;

  //initialize coordinates using the input vertices and persistence
  VertEdgePair(int vbi, int edi, double rob, double b, double d) : vbidx(vbi),edidx(edi), robustness(rob),birth(b),death(d){}
  
  bool operator<(const VertEdgePair &rhs) const{
    return (this->robustness >= rhs.robustness);
  }
};

class EdgeTrigPair{
  public:
  int ebidx;
  int tdidx;
  double robustness;
  double birth;
  double death;
  
  //initialize coordinates using the input vertices and persistence
  EdgeTrigPair( int ebi, int tdi, double rob,double b,double d) : ebidx(ebi),tdidx(tdi),robustness(rob),birth(b),death(d){}

  bool operator<(const EdgeTrigPair &rhs) const{    // the EdgeTrigPair pair will be ranked fisrtly if the robustness ranked
    return (this->robustness >= rhs.robustness);
  }
};

//xudong add the following code.
class TrigCubePair{
public:
    int tbidx;
    int cdidx;
    double robustness;
    double birth;
    double death;
    
    //initialize coordinates using the input vertices and persistence
    TrigCubePair( int tbi, int cdi, double rob,double b,double d) : tbidx(tbi),cdidx(cdi),robustness(rob),birth(b),death(d){}
    
    bool operator<(const TrigCubePair &rhs) const{
        return (this->robustness >= rhs.robustness);
    }
};

//xudong add the class path.
class path{
public:
    int zcoord;
    int xcoord;
    int ycoord;
    double max_intensity;
    double min_intensity;
    path(int z, int x, int y, double max, double min) : zcoord(z), xcoord(x), ycoord(y), max_intensity(max), min_intensity(min) {}
    path() : zcoord(-1), xcoord(-1), ycoord(-1), max_intensity(0), min_intensity(0) {}
    
    //for sorting.
/*
    bool operator<(const path & other) const{
        // two coordinated are equal
        if((this->zcoord==other.zcoord)&&(this->xcoord==other.xcoord)){
//        return ((this->zcoord<=other.zcoord)&&(this->xcoord<=other.xcoord)&&(this->ycoord<=other.ycoord));
            return (this->ycoord<=other.ycoord);
        }
        else if((this->zcoord==other.zcoord)&&(this->ycoord==other.ycoord)){
            return (this->xcoord<other.xcoord);
        }
        else if((this->xcoord==other.xcoord)&&(this->ycoord==other.ycoord)){
            return (this->zcoord<other.zcoord);
        }
        
        //only one coordinate is equal.
        else if((this->zcoord==other.zcoord)&&(this->xcoord!=other.xcoord)&&(this->ycoord!=other.ycoord)){
            //        return ((this->zcoord<=other.zcoord)&&(this->xcoord<=other.xcoord)&&(this->ycoord<=other.ycoord));
            return (this->xcoord<=other.xcoord);
        }
        else if((this->zcoord!=other.zcoord)&&(this->xcoord==other.xcoord)&&(this->ycoord!=other.ycoord)){
            return (this->zcoord<other.zcoord);
        }
        else if((this->zcoord!=other.zcoord)&&(this->xcoord!=other.xcoord)&&(this->ycoord==other.ycoord)){
            return (this->zcoord<other.zcoord);
        }
        
    }

    bool operator<(const path & other) const{
        // two coordinated are equal
        if(this->zcoord==other.zcoord){
            //        return ((this->zcoord<=other.zcoord)&&(this->xcoord<=other.xcoord)&&(this->ycoord<=other.ycoord));
            if(this->xcoord==other.xcoord){
                if(this->ycoord==other.ycoord){
                    // z=z, x=x, y=y -> do nothing.
                }
                else{
                    return (this->ycoord<=other.ycoord);    // z=z, x=x, y!=y
                }
            }
            else{
                if(this->ycoord==other.ycoord){
                    return (this->xcoord<=other.xcoord);// z=z, x!=x, y=y
                }
                else{
                    return ((this->xcoord<=other.xcoord)&&(this->ycoord<=other.ycoord));// z=z, x!=x, y!=y
                }
                
            }
        }
        else {
            if(this->xcoord==other.xcoord){
                if(this->ycoord==other.ycoord){
                    return (this->zcoord<=other.zcoord);// z!=z, x=x, y=y
                }
                else{
                   return ((this->zcoord<=other.zcoord)&&(this->ycoord<=other.ycoord));// z!=z, x=x, y!=y
                }
            }
            else{
                if(this->ycoord==other.ycoord){
                    return ((this->zcoord<=other.zcoord)&&(this->xcoord<=other.xcoord));// z!=z, x!=x, y=y
                }
                else{
                    return ((this->zcoord<=other.zcoord)&&(this->xcoord<=other.xcoord)&&(this->ycoord<=other.ycoord));// z!=z, x!=x, y!=y
                }
            }
            
        }
        
    }
*/
    // for function erase().
    bool operator==(const path & other) const {
        return ((this->zcoord==other.zcoord)&&(this->xcoord==other.xcoord)&&(this->ycoord==other.ycoord));
    }
 
//    int x() const { return x_; }
//    int y() const { return y_; }

};



//xudong: show informaiton by using the following functions.
void show_ve(VertEdgePair ve){
    cout<<"vbidx: "<<ve.vbidx<<", edidx: "<<ve.edidx<<", robustness: "<<ve.robustness<<", birth: "<<ve.birth<<", death: "<<ve.death<<endl;
}
void show_et(EdgeTrigPair et){
    cout<<"ebidx: "<<et.ebidx<<", tdidx: "<<et.tdidx<<", robustness: "<<et.robustness<<", birth: "<<et.birth<<", death: "<<et.death<<endl;
}
void show_tc(TrigCubePair tc){
    cout<<"tbidx: "<<tc.tbidx<<", cdidx: "<<tc.cdidx<<", robustness: "<<tc.robustness<<", birth: "<<tc.birth<<", death: "<<tc.death<<endl;
}
//show the xidx and yidx in the 1D Vertex vector.
void showvList(vector<Vertex> * vL){
    for(int i=0; i<vL->size(); i++){
        Vertex temp11 = (*vL)[i];
        cout<<"xidx: "<<temp11.xidx<<", yidx: "<<temp11.yidx<<", order="<<i<<endl;
        //            Vertex v1 = (*vList)[cellOrder[mapxcoord-1][mapycoord-1]];
    }
    
}
void showeList(vector<Edge> * eL){
    for(int i=0; i<eL->size(); i++){
        Edge temp11 = (*eL)[i];
        cout<<"Order= "<<i<<" -> v1_order: "<<temp11.v1_order<<", v2_order: "<<temp11.v2_order<<", mapZCoord: "<<temp11.mapZCoord<<", mapXCoord: "<<temp11.mapXCoord<<", mapYCoord: "<<temp11.mapYCoord<<i<<endl;
        //            Vertex v1 = (*vList)[cellOrder[mapxcoord-1][mapycoord-1]];
    }
    
}
void showtList(vector<Triangle> * tL){
    for(int i=0; i<tL->size(); i++){
        Triangle temp11 = (*tL)[i];
        cout<<"Order= "<<i<<" -> v1_order: "<<temp11.v1_order<<", v2_order: "<<temp11.v2_order<<", v3_order: "<<temp11.v3_order<<", v4_order: "<<temp11.v4_order<<", mapZCoord: "<<temp11.mapZCoord<<", mapXCoord: "<<temp11.mapXCoord<<", mapYCoord: "<<temp11.mapYCoord<<endl;
    }
    
}
void showt(Triangle tL){
        Triangle temp11 = tL;
        cout<<" -> v1_order: "<<temp11.v1_order<<", v2_order: "<<temp11.v2_order<<", v3_order: "<<temp11.v3_order<<", v4_order: "<<temp11.v4_order<<", mapZCoord: "<<temp11.mapZCoord<<", mapXCoord: "<<temp11.mapXCoord<<", mapYCoord: "<<temp11.mapYCoord<<endl;
    
}


// compare Vertex 1D vector, if idfferent?
void comList (vector<Vertex> * vL1, vector<Vertex> * vL2){
    vector<Vertex> * dif;
    for(int j=0; j<vL1->size(); j++){
        Vertex temp_1 = (*vL1)[j], temp_2 = (*vL2)[j];
        if(temp_1.xidx!=temp_2.xidx || temp_1.yidx!=temp_2.yidx){
            dif->push_back(temp_1);
            dif->push_back(temp_2);
        }
    }
    if(dif->size()!=0){
        showvList(dif);}
    else
        cout<<"Those two vList are same."<<endl;
}
//for( vector< int >::iterator myiiter=circle_elist.begin(); myiiter!=circle_elist.end(); myiiter++, counter_vertex++ ){
// show 2D vector:
void show2dvector(vector< vector < int > > TwoDVector){
    cout<<"size: "<<TwoDVector.size()<<", capacity: "<<TwoDVector.capacity()<<endl;
    vector< vector < int > >::iterator i;
    vector< int >::iterator j;
    int cot_i;
    for(i=TwoDVector.begin(); i!=TwoDVector.end(); i++){
        cot_i = i - TwoDVector.begin();
        cout<<"--> "<<cot_i<<": ";
        if((*i).empty()){continue;}
        for(j=i->begin(); j!=i->end(); j++){
            cout<<*j<<" ";
        }
        
        cout<<endl;
    }
    
}
void show2dvector(vector< vector < int > > *TwoDVector){
    vector< vector < int > >::iterator i;
    vector< int >::iterator j;
    
    for(i=(*TwoDVector).begin(); i!=(*TwoDVector).end(); i++){
        if((*i).empty()){continue;}
        for(j=i->begin(); j!=i->end(); j++){
            cout<<*j<<" ";
        }
        
        cout<<endl;
    }
    
}
void show1dchords(vector < path > temp){
    for(int j=0; j<temp.size(); j++){
        cout<<j<<" -> z: "<<temp[j].zcoord+1<<", x: "<<temp[j].xcoord+1<<", y:"<<temp[j].ycoord+1<<endl;
    }

}
void show1dvectorPath(vector < path > temp, vector< vector< vector< double > > > data){
        int z,x,y;
        if(0==temp.size()) cout<<"size of the path is 0 !!!!"<<endl;
        else{
            cout<<"------------------------------------- 1dvectorPath size: "<<temp.size()<<"-------------------------------------"<<endl;
            for(int j=0; j<temp.size(); j++){
                z = temp[j].zcoord;
                x = temp[j].xcoord;
                y = temp[j].ycoord;
                cout<<j<<" -> z: "<<z<<", x: "<<x<<", y:"<<y<<", intensity: "<<data[z][x][y]<<endl;
            }
        }
}


void show2dvectorPath(vector< vector < path > > temp, vector< vector< vector< double > > > data){
    int z,x,y;
    for(int i=0; i<temp.size(); i++){   //row# != 0, but not element in first row.
        if(0==temp[i].size()) continue;
        else{
            cout<<"-------------------------------------"<<i<<", temp["<<i<<"].size(): "<<temp[i].size()<<"-------------------------------------"<<endl;
            for(int j=0; j<temp[i].size(); j++){
                z = temp[i][j].zcoord;
                x = temp[i][j].xcoord;
                y = temp[i][j].ycoord;
                cout<<j<<" -> z: "<<z<<", x: "<<x<<", y:"<<y<<", intensity: "<<data[z][x][y]<<endl;
            }
        }
        cout<<endl;
    }
}



void show2dchords(vector< vector < path > > temp){
    for(int i=0; i<temp.size(); i++){   //row# != 0, but not element in first row.
        if(0==temp[i].size()) continue;
        else{
            cout<<"-------------------------------------"<<i<<", temp["<<i<<"].size(): "<<temp[i].size()<<"-------------------------------------"<<endl;
            for(int j=0; j<temp[i].size(); j++){
                cout<<j<<" -> z: "<<temp[i][j].zcoord+1<<", x: "<<temp[i][j].xcoord+1<<", y:"<<temp[i][j].ycoord+1<<endl;
            }
        }
        cout<<endl;
    }
    
}
void show2dchords2(vector< vector < path > > temp){
    for(int i=0; i<temp.size(); i++){   //row# != 0, but not element in first row.
        if(0==temp[i].size()) continue;
        else{
            cout<<"-------------------------------------"<<i<<", temp["<<i<<"].size(): "<<temp[i].size()<<"-------------------------------------"<<endl;
            for(int j=0; j<temp[i].size(); j++){
                cout<<j<<"-> x: "<<temp[i][j].xcoord+1<<", y:"<<temp[i][j].ycoord+1<<", z: "<<temp[i][j].zcoord+1<<endl;
            }
        }
        cout<<endl;
    }
    
}
void show1dvector(vector<double> OneDVector){
    //    cout<<"enter 1"<<endl;
    for(int i=0; i<OneDVector.size(); i++){
        cout<<i<<"-> "<<OneDVector[i]<<endl;
    }
}
int show1dvector(vector<int> OneDVector, int showmode=1){
    //    cout<<"enter 2"<<endl;
    int counter_temp = 0;
    for(int i=0; i<OneDVector.size(); i++){
        if(1==showmode){
            cout<<OneDVector[i]<<"  ";
        }
        else if(2==showmode){
            if(i%100==0) cout<<endl;
            cout<<OneDVector[i]<<" ";
        }
        else if(3==showmode){
            if(-1!=OneDVector[i]) counter_temp++;
            
        }
        else if(4==showmode){
            cout<<OneDVector[i]<<" ";
            
        }
        else{
            cout<<OneDVector[i]<<endl;
        }
    }
    if(4!=showmode){
        cout<<endl;
    }
    return counter_temp;
}
void intenOneChord(vector< vector< path > > chords){
    for(int j=0; j<chords.size(); j++){
        double intensity = 0;
        int i;
        for(i=0; i<chords[j].size(); i++){
            intensity += phi3d->data[chords[j][i].zcoord][chords[j][i].xcoord][chords[j][i].ycoord];
        }
        cout<<j<<" -> Total Intensity is: "<<intensity<<", from "<<i<<" vertices"<<endl;
    }
    
    
};


void showedge(Edge temp){
    cout<<"v1_order: "<<temp.v1_order<<", v2_order: "<<temp.v2_order<<", k: "<<temp.mapXCoord<<", i: "<<temp.mapXCoord<<", j: "<<temp.mapXCoord<<endl;
}
void showtri(Triangle temp){
    cout<<"v1_order: "<<temp.v1_order<<", v2_order: "<<temp.v2_order<<", v3_order: "<<temp.v3_order<<", v4_order: "<<temp.v4_order<<endl;
    cout<<"e1_order: "<<temp.e1_order<<", e2_order: "<<temp.e2_order<<", e3_order: "<<temp.e3_order<<", e4_order: "<<temp.e4_order<<endl;
    cout<<", k: "<<temp.mapXCoord<<", i: "<<temp.mapXCoord<<", j: "<<temp.mapXCoord<<endl;
}
void showcube(Cube temp){
    cout<<"v1_order: "<<temp.v1_order<<", v2_order: "<<temp.v2_order<<", v3_order: "<<temp.v3_order<<", v4_order: "<<temp.v4_order<<"v, 5_order: "<<temp.v5_order<<", v6_order: "<<temp.v6_order<<"v7_order: "<<temp.v7_order<<", v8_order: "<<temp.v8_order<<endl;
    cout<<"e1_order: "<<temp.e1_order<<", e2_order: "<<temp.e2_order<<", e3_order: "<<temp.e3_order<<", e4_order: "<<temp.e4_order<<endl;
    cout<<"e5_order: "<<temp.e5_order<<", e6_order: "<<temp.e6_order<<", e7_order: "<<temp.e7_order<<", e8_order: "<<temp.e8_order<<endl;
    cout<<"e9_order: "<<temp.e9_order<<", e10_order: "<<temp.e10_order<<", e11_order: "<<temp.e11_order<<", e12_order: "<<temp.e12_order<<endl;
    cout<<"t1_order: "<<temp.t1_order<<", t2_order: "<<temp.t2_order<<", t3_order: "<<temp.t3_order<<endl;
    cout<<"t4_order: "<<temp.t4_order<<", t5_order: "<<temp.t5_order<<", t6_order: "<<temp.t6_order<<endl;
    cout<<", k: "<<temp.mapXCoord<<", i: "<<temp.mapXCoord<<", j: "<<temp.mapXCoord<<endl;
}
//xudong: show informaiton by using the above functions.

//-----------------------------------------------------
//compute 2D persistence
// m,n: size of the two dimensions
// pers_thd: threshold of persistence (only bigger persistence would be recorded
// rob_thd: threshold of robustness
// levelset_val: the image value of the levelset (0 in image segmentation), xudong: which means that phi(x)<0.
// persistenceM: persistence flow, +pers to creator and -pers to destroyer
// robustnessM: robustness flow, +pers to creator or -pers to destroyer, depending on which is closer to the levelset_val 
// veList: vertex-edge pair, together with corresponding persistence
// etrigList: edge-triangle pair, together with corresponding persistence
//
// assume the global variable phi is already available (which stores the height function)
//-----------------------------------------------------


enum CellTypeEnum {CT_UNDEFINED, VERTEX, EDGEVERT, EDGEHORI, ZEDGE, TRIG, TRIGHORI, TRIGVERT, CUBE};   //xudong add CUBE.
//enum EdgePersTypeEnum {EP_UNDEFINED, DESTROYER, CREATOR};
enum ETrigPersTypeEnum {EP_UNDEFINED, DESTROYER, CREATOR};
enum CellFixTypeEnum {CF_UNDEFINED, MOVEDOWN, MOVEUP};

void print3dmatrix(vector< vector< vector < CellTypeEnum > > > data){
    //    cout<<"2. x dimension: "<<data.size()<<", y dimension: "<<data[0].size()<<", z dimension: "<<data[0][0].size()<<endl;
    for(int x=0; x<data.size();x++){
        for(int y=0; y<data[0].size();y++){
            for(int z=0; z<data[0][0].size();z++){
                cout<<data[x][y][z]<<"  ";
            }
            cout<<endl;
        }
        cout<<"-------------- the above is slice: "<<x<<" ---------------------"<<endl;
    }
}
void print3dmatrix(vector< vector< vector < int > > > data){
    //    cout<<"2. x dimension: "<<data.size()<<", y dimension: "<<data[0].size()<<", z dimension: "<<data[0][0].size()<<endl;
    for(int x=0; x<data.size();x++){
        for(int y=0; y<data[0].size();y++){
            for(int z=0; z<data[0][0].size();z++){
                cout<<data[x][y][z]<<"  ";
            }
            cout<<endl;
        }
        cout<<"-------------- the above is slice: "<<x<<" ---------------------"<<endl;
    }
}


//calculate the euclide distance between (zcoord,xcoord,ycoord) and everypoint in marker_coord. Return true if there exist on distance less than the n_pix.
bool vertexRangeChecker(const int n_pix,int tolerance,const int zcoord, const int xcoord, const int ycoord, vector< vector< int > > marker_coord){
    int z, x, y, temp_c, temp_a, temp_b;
    int bias = n_pix + tolerance;
    bool temp=false;
    for(int iter=0; iter<marker_coord.size(); iter++){
        z = marker_coord[iter][2];
        x = marker_coord[iter][0];
        y = marker_coord[iter][1];
        temp_c = (zcoord-z)*(zcoord-z);
        temp_a = (xcoord-x)*(xcoord-x);
        temp_b = (ycoord-y)*(ycoord-y);
        //        cout<<"zcoord: "<<zcoord<<", xcoord: "<<ycoord<<", ycoord: "<<xcoord<<endl;
        //        cout<<iter<<": marker-> z-bias: "<<z-bias<<", z+bias: "<<z+bias<<endl;
        //        cout<<iter<<": marker-> x-bias: "<<y-bias<<", x+bias: "<<y+bias<<endl;
        //        cout<<iter<<": marker-> y-bias: "<<x-bias<<", y+bias: "<<x+bias<<endl;
//        if((z-bias <= zcoord)&&(zcoord <= z+bias)&&(x-bias <= xcoord)&&(xcoord <= x+bias)&&(y-bias <= ycoord)&&(ycoord <= y+bias)){
        if(temp_c + temp_a + temp_b <= n_pix*n_pix){
            temp=true;
            break;
        }
        //        else{cout<<"return false"<<endl; return false;}
    }
    return temp;
}

//just check one piece of position (x,y,z).
//calculate the euclide distance between (zcoord,xcoord,ycoord) and only one point in marker_coord. Return true if there exist on distance less than the n_pix.
bool vertexRangeChecker(const int n_pix,int tolerance,const int zcoord, const int xcoord, const int ycoord, vector< int > marker_coord){
    int z, x, y, temp_c, temp_a, temp_b;
    int bias = n_pix + tolerance;
    bool temp=false;
    z = marker_coord[2];
    x = marker_coord[0];
    y = marker_coord[1];
    temp_c = (zcoord-z)*(zcoord-z);
    temp_a = (xcoord-x)*(xcoord-x);
    temp_b = (ycoord-y)*(ycoord-y);

//    if((z-bias <= zcoord)&&(zcoord <= z+bias)&&(x-bias <= xcoord)&&(xcoord <= x+bias)&&(y-bias <= ycoord)&&(ycoord <= y+bias)){
    if(temp_c + temp_a + temp_b <= n_pix*n_pix){
        temp=true;
    }
    //        else{cout<<"return false"<<endl; return false;}
    return temp;
}

//checking whether the chords were inside the range
//(zcoord,xcoord,ycoord) within the chordRange(z,x,y)+bias(z_bias,x_bias,y_bias)..
bool insideMarkerBox(const int z_bias,const int x_bias,const int y_bias,const int zcoord, const int xcoord, const int ycoord, const vector< int > chordRange){
    if(chordRange.size() != 6){cout<<"Wrong!!! -> chordRange.size() != 6"<<endl; return false;}
    int z_min = chordRange[0];
    int z_max = chordRange[1];
    int x_min = chordRange[2];
    int x_max = chordRange[3];
    int y_min = chordRange[4];
    int y_max = chordRange[5];
    
    if((z_min-z_bias <= zcoord)&&(zcoord <= z_max+z_bias)&&
       (x_min-x_bias <= xcoord)&&(xcoord <= x_max+x_bias)&&
       (y_min-y_bias <= ycoord)&&(ycoord <= y_max+y_bias)){
        return true;
    }
    else{
        return false;
    }
}

//change the intensity on point (zcoord, xcoord, ycoord) and its neighbors (with radius=n_pix) to max_total.
void changeInteCore(const int n_pix, vector< vector< vector< double > > > & data, const int zcoord, const int xcoord, const int ycoord, double max_total){
    const int NSlice = data.size();
    const int NRow = data[0].size();
    const int NCol = data[0][0].size();
    const int range = n_pix*2+1;
    
    if(0 == n_pix){
        data[zcoord][xcoord][ycoord]=max_total;   // only change the intensity of point at [zcoord][xcoord][ycoord].
    }
    else{
        for (int k=0;k<range;k++){ //make sure: 0 <= zcoord-n_pix+k <= mapNSlice-1 <=> 0 <= zcoord-n_pix+k < mapNSlice
            if ((zcoord-n_pix+k < 0)||(zcoord-n_pix+k >= NSlice)){ continue;}
            for (int i=0;i<range;i++){
                if ((xcoord-n_pix+i < 0)||(xcoord-n_pix+i >= NRow)){ continue;}
                for(int j=0;j<range;j++){
                    if ((ycoord-n_pix+j < 0)||(ycoord-n_pix+j >= NCol)){ continue;}
                    data[zcoord-n_pix+k][xcoord-n_pix+i][ycoord-n_pix+j]=max_total;
                }
            }
        }
    }
    
}
void drawBallOnTips(const int r_pix, vector< vector< vector< double > > > & data, const int z_coord, const int x_coord, const int y_coord, double toIntensity){
    int temp_c, temp_a, temp_b;
    const int NSlice = data.size();
    const int NRow = data[0].size();
    const int NCol = data[0][0].size();
    if(0 == r_pix){
        data[z_coord][x_coord][y_coord]=toIntensity;   // only change the intensity of point at [zcoord][xcoord][ycoord].
    }
    else{
        for(int k=z_coord-r_pix; k<=z_coord+r_pix; k++){
            if ((k < 0)||(k >= NSlice)){ continue;}
            for(int i=x_coord-r_pix; i<=x_coord+r_pix; i++){
                if ((i < 0)||(i >= NRow)){ continue;}
                for(int j=y_coord-r_pix; j<=y_coord+r_pix; j++){
                    if ((j < 0)||(j >= NCol)){ continue;}
                    temp_c = (z_coord-k)*(z_coord-k);
                    temp_a = (x_coord-i)*(x_coord-i);
                    temp_b = (y_coord-j)*(y_coord-j);
                    if(temp_c + temp_a + temp_b <= r_pix*r_pix){
                        data[k][i][j]=toIntensity;
                    }
                }
            }
        }
    }
}
void oneSortedChordsToTwoChords(vector< vector< path > > sortedChords,
                                vector< vector< int > > pM_marker,
                                vector< vector< int > >mV_marker,
                                vector< vector< path > > & halfSortedChords1,
                                vector< vector< path > > & halfSortedChords2
                                ){
    
    int pM_marker_z, pM_marker_x, pM_marker_y, mV_marker_z, mV_marker_x, mV_marker_y;
    
    int z, x, y, iter, i, j;
    for(i=0; i<sortedChords.size(); i++){
        int pM_dis = 0xFFFFFFF, mV_dis = 0xFFFFFFF; //int mV_dis, pM_dis1, pM_dis2;
        int pM_dis_temp, mV_dis_temp;
//        cout<<"pM_dis: "<<pM_dis<<endl; cout<<"mV_dis: "<<mV_dis<<endl;
        path end1, end2;
        int end1_id, end2_id;
//        cout<<"-------------------------------------------------------"<<endl;
        for(j=0; j<sortedChords[i].size(); j++){
            z = sortedChords[i][j].zcoord;
            x = sortedChords[i][j].xcoord;
            y = sortedChords[i][j].ycoord;
            
            //calculat the shorted distance from the chord to pM.
            for(iter=0; iter<pM_marker.size(); iter++){
                pM_marker_z = pM_marker[iter][2];
                pM_marker_x = pM_marker[iter][0];
                pM_marker_y = pM_marker[iter][1];
                
                pM_dis_temp = (z-pM_marker_z)*(z-pM_marker_z) + (x-pM_marker_x)*(x-pM_marker_x) + (y-pM_marker_y)*(y-pM_marker_y);
                if(pM_dis_temp < pM_dis){
                    pM_dis = pM_dis_temp;
                    end1.zcoord = z;
                    end1.xcoord = x;
                    end1.ycoord = y;
                    end1_id = j;
                }
            }
            
            for(iter=0; iter<mV_marker.size(); iter++){
                mV_marker_z = mV_marker[iter][2];
                mV_marker_x = mV_marker[iter][0];
                mV_marker_y = mV_marker[iter][1];
                
                mV_dis_temp = (z-mV_marker_z)*(z-mV_marker_z) + (x-mV_marker_x)*(x-mV_marker_x) + (y-mV_marker_y)*(y-mV_marker_y);
                if(mV_dis_temp < mV_dis){
                    mV_dis = mV_dis_temp;
                    end2.zcoord = z;
                    end2.xcoord = x;
                    end2.ycoord = y;
                    end2_id = j;
                }
            }
            

        }
//        cout<<"Pass finding end1 and end2 "<<endl;
//        cout<<"end1_id: "<<end1_id<<endl; cout<<"end2_id: "<<end2_id<<endl; //end1_id: 3, end2_id: 89
//        cout<<"----------------------------"<<endl;
        int breakiter_1 = 0;
        for(j=end1_id; j<sortedChords[i].size(); j=(j+1)%(sortedChords[i].size()), breakiter_1++){
            halfSortedChords1[i].push_back(sortedChords[i][j]);
            if(j == end2_id){
                break;
            }
            if(breakiter_1>=sortedChords[i].size()){
                break;
            }
        }
//        cout<<"Pass halfSortedChords1 in oneSortedChordsToTwoChords()"<<endl;
        
        int breakiter_2 =0;
        for(j=end2_id; j<sortedChords[i].size(); j=(j+1)%(sortedChords[i].size()), breakiter_2++){
            halfSortedChords2[i].push_back(sortedChords[i][j]);
            if(j == end1_id){
                break;
            }
            if(breakiter_2>=sortedChords[i].size()){
                break;
            }
        }
//        cout<<"Pass halfSortedChords2 in oneSortedChordsToTwoChords()"<<endl;
//        cout<<"-------------------------------------------------------"<<endl;
        
    }
    
}
void greddyAlo(vector< vector< vector< double > > > data, vector<vector< path > >halfSortedChords, vector< vector< path > > & newhalfSortedChords, const int times){
    const int NSlice = data.size();
    const int NRow = data[0].size();
    const int NCol = data[0][0].size();
    //            const int range = n_pix*2+1;
    int row, col, z, x, y;
    srand( (unsigned)time(NULL) );
    
    //copy halfSortedChords to newhalfSortedChords:
    for(int temp_row = 0; temp_row < halfSortedChords.size(); temp_row++){
        for(int temp_col = 0; temp_col < halfSortedChords[temp_row].size(); temp_col++){
            newhalfSortedChords[temp_row].push_back(halfSortedChords[temp_row][temp_col]);
        }
//        copy ( myints, myints+7, myvector.begin() );
    }
//    show1dvectorPath(halfSortedChords[0], data);
//    show2dvectorPath(newhalfSortedChords, data);
//    cout<<"================================ *************************** ================================"<<endl;
//    cout<<"----- 1 ------"<<endl;
//    cout<<"newhalfSortedChords.size(): "<<newhalfSortedChords.size()<<endl; //=0?
    //greddy algorithm:
    for(row=0; row<newhalfSortedChords.size(); row++){
//    for(row=0; row<1; row++){   //only check the first chord.
        // repeat for multiple times.
//        cout<<"----- 2 ------"<<endl;
        for(int times_temp=0; times_temp<times; times_temp++){
            // two points randomly selected from 0 to size()-1 on the halfSortedChords.
            int pointIndex_1 = newhalfSortedChords[row].size()*(double)rand()/RAND_MAX-1;
            int pointIndex_2 = newhalfSortedChords[row].size()*(double)rand()/RAND_MAX-1;  //max pointIndex_2 <= randNum - 1.
            while(pointIndex_1 == pointIndex_2){
                pointIndex_2 = newhalfSortedChords[row].size()*(double)rand()/RAND_MAX;
            }
//            cout<<"----- 3 ------"<<endl;
            if(pointIndex_1 > pointIndex_2){
                // pointIndex_2 = halfSortedChords.size()*(double)rand()/RAND_MAX;
                int temp = pointIndex_1;
                pointIndex_1 = pointIndex_2;
                pointIndex_2 = temp;
            }
//            cout<<"-------------------------------------------------------"<<endl;
//            cout<<"halfSortedChords[row].size(): "<<halfSortedChords[row].size()<<", pointIndex_1: "<<pointIndex_1<<", pointIndex_2: "<<pointIndex_2<<endl;
            
            //        path firstPoint = halfSortedChords[row][pointIndex_1];
            //        newhalfSortedChords[row].push_back(prePoint); // push the first point(path) into newhalfSortedChords.
//            cout<<"----- 4 ------"<<endl;
            for(col=pointIndex_1; col<pointIndex_2; col++){
//                z = halfSortedChords[row][col].zcoord;
//                x = halfSortedChords[row][col].xcoord;
//                y = halfSortedChords[row][col].ycoord;
                
                z = newhalfSortedChords[row][col].zcoord;
                x = newhalfSortedChords[row][col].xcoord;
                y = newhalfSortedChords[row][col].ycoord;
                double lowestIntensity = 0xFFFFFFF*1.0;
                path betterpoint = path();
                
                //checking neighboring 26 points of (z,x,y).
                for (int k=z-1;k<z+2;k++){ //make sure: 0 <= zcoord-n_pix+k <= mapNSlice-1 <=> 0 <= zcoord-n_pix+k < mapNSlice
                    if ((k < 0)||(k >= NSlice)){ continue;}
                    for (int i=x-1;i<x+2;i++){
                        if ((i < 0)||(i >= NRow)){ continue;}
                        for(int j=y-1;j<y+2;j++){
                            if ((j < 0)||(j >= NCol)){ continue;}
                            //                        if(k==prePoint.zcoord || i==prePoint.xcoord || j==prePoint.ycoord){continue;}   // discard the starting point.
                            if(lowestIntensity > data[k][i][j]){
                                lowestIntensity = data[k][i][j];
                                betterpoint.zcoord = k;
                                betterpoint.xcoord = i;
                                betterpoint.ycoord = j;
                            }
                        }
                    }
                }
                newhalfSortedChords[row][col] = betterpoint;
            }
//            cout<<"----- 5 ------"<<endl;
//            show1dvectorPath(newhalfSortedChords[row], data);
//            cout<<"=========================================================================================="<<endl;
        }
        
    }
//    cout<<"================================ *************************** ================================"<<endl;
//    show1dvectorPath(halfSortedChords[0], data);
}

// find the pM_coord according to the ball around the tips.
vector< vector < int > > reconstruct_pM(vector< vector< int > > coord, const int r_pix,const int NSlice,const int NRow,const int NCol){
    if(!coord.empty()){
        int k, i, j;
        int z_coord, x_coord,  y_coord;
        int temp_a, temp_b, temp_c;
        vector<int> temp_coord_1d(3,0);
        vector< vector< int > > temp_coord_2d;
//        for(i=0; i<coord.size(); i++){
//            for(j=0; j<coord[0].size(); j++){
//                temp[i][j] = coord[i][j];
//            }
//        }
        
//        coord.clear();  // Clear 2d vector coord, then recreate it as followings.
        
        for(i=0; i<coord.size(); i++){
            z_coord = coord[i][2];
            x_coord = coord[i][0];
            y_coord = coord[i][1];
            for(int k=z_coord-r_pix; k<=z_coord+r_pix; k++){
                if ((k < 0)||(k >= NSlice)){ continue;}
                for(int i=x_coord-r_pix; i<=x_coord+r_pix; i++){
                    if ((i < 0)||(i >= NRow)){ continue;}
                    for(int j=y_coord-r_pix; j<=y_coord+r_pix; j++){
                        if ((j < 0)||(j >= NCol)){ continue;}
                        temp_c = (z_coord-k)*(z_coord-k);
                        temp_a = (x_coord-i)*(x_coord-i);
                        temp_b = (y_coord-j)*(y_coord-j);
                        if(temp_c + temp_a + temp_b <= r_pix*r_pix){
                            temp_coord_1d[0] = i; temp_coord_1d[1] = j; temp_coord_1d[2] = k;
                            temp_coord_2d.push_back(temp_coord_1d);
                        }
                    }
                }
            }
            
        }
        return temp_coord_2d;
  
    }
    else{
        cout<<"Wrong!!! -> vector< vector< int > > coord is empty "<<endl;
    }
}

// find the mV_coord according to the intensity.
vector< vector < int > > reconstruct_mV(vector< vector< vector< double > > > data, vector< vector< int > > coord, double intensity_max, double intensity_min){
    int nslice = data.size();
    int nrow = data[0].size();
    int ncol = data[0][0].size();
    int k, j, i;
    vector<int> temp_coord_1d(3,0);
    vector< vector < int > > temp_coord_2d;
//    cout<<"in Program -> size of phi2d_mV: "<<coord.size()<<endl;
    
    for(k=0;k<nslice;k++){
        for (i=0;i<nrow;i++){
            for (j=0;j<ncol;j++){
//                if(data[k][i][j] <= intensity_max && data[k][i][j] >= intensity_min){ // voxel within intensity range will add into the temp_coord_2d (more elements).
                if(data[k][i][j] == intensity_min){ // only the voxel with intensity_min will add into the temp_coord_2d (less elements).
                    temp_coord_1d[0] = i; temp_coord_1d[1] = j; temp_coord_1d[2] = k;
                    temp_coord_2d.push_back(temp_coord_1d);
//                    cout<<"Happened!!!"<<endl;
                }
            }
        }
    }
    return temp_coord_2d;
}

vector< vector< path > > inside_mVpM_box(const int z_bias,const int x_bias,const int y_bias, const vector< vector< path > > SortedPath, const vector< int > chordRange){
    if(chordRange.size() != 6){cout<<"Wrong!!! -> chordRange.size() != 6"<<endl;}
    vector< vector< path > > final_SortedPath;
    int i, j, z_coord, x_coord, y_coord;
    int z_min = chordRange[0];
    int z_max = chordRange[1];
    int x_min = chordRange[2];
    int x_max = chordRange[3];
    int y_min = chordRange[4];
    int y_max = chordRange[5];
    bool temp;
    
    for(i=0; i<SortedPath.size(); i++){
        temp = true;
        for(j=0; j<SortedPath[i].size(); j++){
            z_coord = SortedPath[i][j].zcoord;
            x_coord = SortedPath[i][j].xcoord;
            y_coord = SortedPath[i][j].ycoord;
            
            if((z_min-z_bias <= z_coord)&&(z_coord <= z_max+z_bias)&&
               (x_min-x_bias <= x_coord)&&(x_coord <= x_max+x_bias)&&
               (y_min-y_bias <= y_coord)&&(y_coord <= y_max+y_bias)){ }
//            else{ temp = false; cout<<"Find Bad Chord by inside_mVpM_box() "<<endl; break;}
            else{ temp = false; break;}
            
        }
        if(temp == true){final_SortedPath.push_back(SortedPath[i]);}
    }
    return final_SortedPath;
}




//-----------------************-------------------------------------------//
//------------- xudong: The following is the class CellMap ---------------//
//-----------------************-------------------------------------------//


bool checkcellorder(vector< vector< vector < int > > > NeiOrder, vector< vector< vector < int > > > neighborOrder){
    int k, i, j;
    for (k = 0; k < 5; k++){
        for (i = 0; i < 5; i++){
            for (j = 0; j < 5; j++){
                if (NeiOrder[k][i][j] != neighborOrder[k][i][j]){
                    return true;
                }
                
            }
        }
    }
    return false;
}
void datacomp(vector< vector< vector< double > > > a, vector< vector< vector< double > > > b){
    int i, j, k;
    bool temp=true;
    if(a.size() != b.size()){cout<<" Dim of a.size() and b.size() is different "<<endl;}
    if(a[0].size() != b[0].size()){cout<<" Dim of a[0].size() and b[0].size() is different "<<endl;}
    if(a[0][0].size() != b[0][0].size()){cout<<" Dim of a[0][0].size() and b[0][0].size() is different "<<endl;}
    for(i=0; i<a.size(); i++){
        for(j=0; j<a[i].size(); j++){
            for(k=0; k<a[i][j].size(); k++){
                if(a[i][j][k] != b[i][j][k]){
                    temp=false;
                    break;
                }
            }
        }
    }
    if(temp == true){cout<<"a and b are same!!!"<<endl;}
    else{cout<<" a and b are different!!!"<<endl;}
}


class CellMap{
public:
	int vertNum,edgeNum,trigNum, cubeNum;    // xudong add cubeNum.
    int mapNRow, mapNCol, mapNSlice;   // xudong add mapNSlice.
	int currEOrder, currTOrder, currCOrder;
    int currEOrder1, currTOrder1, currCOrder1;
    
    //create three 3D vectors.
	vector< vector< vector < int > > > cellOrder; // xudong modeify the cellOrder to 3D.
	vector< vector< vector < CellTypeEnum > > > cellType; // xudong modeify the cellType to 3D.
    
    // assign the dimension in constructor, and will be updated from function void setEPersType( vector< int > * low_2D_e2t ).
	vector< vector< vector < ETrigPersTypeEnum > > > edgePersType; // xudong modeify the edgePersType to 3D.
    vector< vector< vector < ETrigPersTypeEnum > > > trigPersType;   // xudong add the trigPersType for 3D image.
    
	vector< Vertex > * vList;   // give the value in the constructor.
	vector< Edge > * eList; // give the value in the constructor.
	vector< Triangle > * trigList;  // give the value in the constructor.
    vector< Cube > * cubeList;  //xudong add the cube List.

	vector< vector< vector < CellFixTypeEnum > > > cellFixType;   // xudong modeify the cellFixType to 3D.

	void setVertOrder(int sid, int rid, int cid, int vorder){
		
		myassert( (sid >=0)&&(sid < mapNSlice) && (rid >=0)&&(rid < mapNRow) && (cid >= 0)&&(cid < mapNCol) );
		myassert( cellType[sid][rid][cid] == VERTEX );  // how to check it's VERTEX, cellType 3d matrix hasn't been writen celltypes???????????????????
		myassert(cellOrder[sid][rid][cid] == -1); //this vert has not been specified, assign value by cellOrder.assign(mapNRow, vector< int >(mapNCol, -1));
        
        bool temp_b1 = (sid >=0)&&(sid < mapNSlice) && (rid >=0)&&(rid < mapNRow) && (cid >= 0)&&(cid < mapNCol);
        bool temp_b2 = (cellType[sid][rid][cid] == VERTEX);
        bool temp_b3 = (cellOrder[sid][rid][cid] == -1);
        
        if(!temp_b1){cout<<"error 1"<<endl;}
        if(!temp_b2){cout<<"error 1"<<endl;}
        if(!temp_b3){cout<<"error 1"<<endl;}
        
        
//        cout<<"rid: "<<rid<<", cid: "<<cid<<", vorder="<<vorder<<endl;
		cellOrder[sid][rid][cid] = vorder;  //write the oder into cellOrder 3D matrix, vorder is according to the intensity of that vertex.
        vector< vector< vector < int > > > neighborOrder(5, vector< vector< int > >(5,vector< int >(5,-1)));  // create 5X5X5 neighbor space (different types).
        vector< vector< vector < int > > > NeiOrder(5, vector< vector< int > >(5,vector< int >(5,-1)));
        
        if (checkcellorder(NeiOrder, neighborOrder)) cout<<"1. NeiOrder and neighborOrder are not the same !!!"<<endl;
        
		int k,i,j;
		// get the neighboring orders
        for (k=sid-2;k<=sid+2;k++){ // those two vertices in front and two vertices in back.
            if ((k<0)||(k>=mapNSlice)) continue;
            for (i=rid-2;i<=rid+2;i++){ // those two vertices above and two vertices below.
                if ((i<0)||(i>=mapNRow)) continue;
                for(j=cid-2;j<=cid+2;j++){  // those two vertices on the left and two vertices on the right.
                    if ((j<0)||(j>=mapNCol)) continue;
                    // copy the order in cellOrder to neighborOrder, but indices in neighborOrder was changed to 0~5.
                    neighborOrder[k+2-sid][i+2-rid][j+2-cid]=cellOrder[k][i][j];  //[0~4][0~4][0~4] = [sid-2 ~ sid+2][rid-2 ~ rid+2][cid-2 ~ cid+2].
                    NeiOrder[k+2-sid][i+2-rid][j+2-cid]=cellOrder[k][i][j];
                }
            }
        }
		
//        cout<<"------------------------ neighborOrder -----------------------"<<endl;
//		show2dvector(neighborOrder);   //--------->
		// update the neighboring orders
// 		int ue,de,le,re,ult,urt,dlt,drt;
// 		int uv,dv,lv,rv,ulv,urv,dlv,drv;
        
        //condition-> 26 neighboring vertices:
		int & LayerSelf_uv = neighborOrder[2][0][2];  int & LayerFront_uv = neighborOrder[0][0][2]; int & LayerBack_uv = neighborOrder[4][0][2];
		int & LayerSelf_dv = neighborOrder[2][4][2];  int & LayerFront_dv = neighborOrder[0][4][2]; int & LayerBack_dv = neighborOrder[4][4][2];
		int & LayerSelf_lv = neighborOrder[2][2][0];  int & LayerFront_lv = neighborOrder[0][2][0]; int & LayerBack_lv = neighborOrder[4][2][0];
		int & LayerSelf_rv = neighborOrder[2][2][4];  int & LayerFront_rv = neighborOrder[0][2][4]; int & LayerBack_rv = neighborOrder[4][2][4];
		int & LayerSelf_ulv = neighborOrder[2][0][0];int & LayerFront_ulv = neighborOrder[0][0][0];int & LayerBack_ulv = neighborOrder[4][0][0];
		int & LayerSelf_urv = neighborOrder[2][0][4];int & LayerFront_urv = neighborOrder[0][0][4];int & LayerBack_urv = neighborOrder[4][0][4];
		int & LayerSelf_dlv = neighborOrder[2][4][0];int & LayerFront_dlv = neighborOrder[0][4][0];int & LayerBack_dlv = neighborOrder[4][4][0];
		int & LayerSelf_drv = neighborOrder[2][4][4];int & LayerFront_drv = neighborOrder[0][4][4];int & LayerBack_drv = neighborOrder[4][4][4];
                                                 int & LayerFront_front_v = neighborOrder[0][2][2];int & LayerBack_back_v = neighborOrder[4][2][2];
        
        //results-> 6 neighboring edges:
		int & ue = neighborOrder[2][1][2];
		int & de = neighborOrder[2][3][2];
		int & le = neighborOrder[2][2][1];
		int & re = neighborOrder[2][2][3];
        int & front_e = neighborOrder[1][2][2];
        int & back_e = neighborOrder[3][2][2];
        
        //results-> 12 neighboring triangles:
		int & LayerSelf_ult = neighborOrder[2][1][1];int & LayerFront_ut = neighborOrder[1][1][2];int & LayerBack_ut = neighborOrder[3][1][2];
		int & LayerSelf_urt = neighborOrder[2][1][3];int & LayerFront_dt = neighborOrder[1][3][2];int & LayerBack_dt = neighborOrder[3][3][2];//correct here!!!
		int & LayerSelf_dlt = neighborOrder[2][3][1];int & LayerFront_lt = neighborOrder[1][2][1];int & LayerBack_lt = neighborOrder[3][2][1];
		int & LayerSelf_drt = neighborOrder[2][3][3];int & LayerFront_rt = neighborOrder[1][2][3];int & LayerBack_rt = neighborOrder[3][2][3];
        
        //results-> 8 neighboring cubes:
        int & LayerFront_ulc = neighborOrder[1][1][1];int & LayerBack_ulc = neighborOrder[3][1][1];
        int & LayerFront_dlc = neighborOrder[1][3][1];int & LayerBack_dlc = neighborOrder[3][3][1];
        int & LayerFront_urc = neighborOrder[1][1][3];int & LayerBack_urc = neighborOrder[3][1][3];
        int & LayerFront_drc = neighborOrder[1][3][3];int & LayerBack_drc = neighborOrder[3][3][3];

        //============ update 6 neighboring edges order =============
		if (LayerSelf_uv>=0){   //int & LayerSelf_uv = neighborOrder[2][0][2]
			ue = currEOrder;    //int & ue = neighborOrder[2][1][2];
			currEOrder++;
		}
		if (LayerSelf_dv>=0){   //int & LayerSelf_dv = neighborOrder[2][4][2];
			de = currEOrder;    //int & de = neighborOrder[2][3][2]
			currEOrder++;
		}
		if (LayerSelf_lv>=0){   //int & LayerSelf_lv = neighborOrder[2][2][0]
			le = currEOrder;    //int & le = neighborOrder[2][2][1]
			currEOrder++;
		}
		if (LayerSelf_rv>=0){   //int & LayerSelf_rv = neighborOrder[2][2][4]
			re = currEOrder;    //int & re = neighborOrder[2][2][3]
			currEOrder++;
		}
        if (LayerFront_front_v>=0){ //int & LayerFront_front_v = neighborOrder[0][2][2]
            front_e = currEOrder;   //int & front_e = neighborOrder[1][2][2]
            currEOrder++;
        }
        if (LayerBack_back_v>=0){   //int & LayerBack_back_v = neighborOrder[4][2][2]
            back_e = currEOrder;    //int & back_e = neighborOrder[3][2][2]
            currEOrder++;
        }
        
        
        
        if (NeiOrder[2][0][2]>=0){
            NeiOrder[2][1][2] = currEOrder1;
            currEOrder1++;
        }
        if (NeiOrder[2][4][2]>=0){
            NeiOrder[2][3][2] = currEOrder1;
            currEOrder1++;
        }
        if (NeiOrder[2][2][0]>=0){
            NeiOrder[2][2][1] = currEOrder1;
            currEOrder1++;
        }
        if (NeiOrder[2][2][4]>=0){
            NeiOrder[2][2][3] = currEOrder1;
            currEOrder1++;
        }
        if (NeiOrder[0][2][2]>=0){
            NeiOrder[1][2][2] = currEOrder1;
            currEOrder1++;
        }
        if (NeiOrder[4][2][2]>=0){
            NeiOrder[3][2][2] = currEOrder1;
            currEOrder1++;
        }
        
        if (checkcellorder(NeiOrder, neighborOrder)) cout<<"2. NeiOrder and neighborOrder are not the same !!!"<<endl;
        
        
        
        
        if ((LayerSelf_uv>=0)&&(LayerSelf_lv>=0)&&(LayerSelf_ulv>=0)){  //202, 220, 200 -> 211
			LayerSelf_ult = currTOrder;
			currTOrder++;
		}
		if ((LayerSelf_uv>=0)&&(LayerSelf_rv>=0)&&(LayerSelf_urv>=0)){  //202, 224, 204 -> 213
			LayerSelf_urt = currTOrder;
			currTOrder++;
		}
        if ((LayerSelf_uv>=0)&&(LayerFront_uv>=0)&&(LayerFront_front_v>=0)){    //202, 002, 022 -> 112
            LayerFront_ut = currTOrder;
            currTOrder++;
        }
        if ((LayerSelf_uv>=0)&&(LayerBack_uv>=0)&&(LayerBack_back_v>=0)){   //202, 402, 422 -> 312
            LayerBack_ut = currTOrder;
            currTOrder++;
        }
        if ((NeiOrder[2][0][2]>=0)&&(NeiOrder[2][2][0]>=0)&&(NeiOrder[2][0][0]>=0)){
            NeiOrder[2][1][1] = currTOrder1;
            currTOrder1++;
        }
        if ((NeiOrder[2][2][4]>=0)&&(NeiOrder[2][0][4]>=0)&&(NeiOrder[2][0][2]>=0)){
            NeiOrder[2][1][3] = currTOrder1;
            currTOrder1++;
        }
        if ((NeiOrder[0][0][2]>=0)&&(NeiOrder[2][0][2]>=0)&&(NeiOrder[0][2][2]>=0)){
            NeiOrder[1][1][2] = currTOrder1;
            currTOrder1++;
        }
        if ((NeiOrder[4][0][2]>=0)&&(NeiOrder[2][0][2]>=0)&&(NeiOrder[4][2][2]>=0)){
            NeiOrder[3][1][2] = currTOrder1;
            currTOrder1++;
        }
        //-------------------------------------------------------------
        
        if ((LayerSelf_dv>=0)&&(LayerSelf_lv>=0)&&(LayerSelf_dlv>=0)){  //242,220, 240 -> 231
            LayerSelf_dlt = currTOrder;
            currTOrder++;
        }
        if ((LayerSelf_dv>=0)&&(LayerSelf_rv>=0)&&(LayerSelf_drv>=0)){  //242, 224, 244 -> 233
            LayerSelf_drt = currTOrder;
            currTOrder++;
        }
        if ((LayerSelf_dv>=0)&&(LayerFront_dv>=0)&&(LayerFront_front_v>=0)){    //242, 042, 022 -> 132
            LayerFront_dt = currTOrder;
            currTOrder++;
        }
        if ((LayerSelf_dv>=0)&&(LayerBack_dv>=0)&&(LayerBack_back_v>=0)){   //242, 442, 422 -> 332
            LayerBack_dt = currTOrder;
            currTOrder++;
        }
        if ((NeiOrder[2][2][0]>=0)&&(NeiOrder[2][4][0]>=0)&&(NeiOrder[2][4][2]>=0)){
            NeiOrder[2][3][1] = currTOrder1;
            currTOrder1++;
        }
        if ((NeiOrder[2][2][4]>=0)&&(NeiOrder[2][4][4]>=0)&&(NeiOrder[2][4][2]>=0)){
            NeiOrder[2][3][3] = currTOrder1;
            currTOrder1++;
        }
        if ((NeiOrder[0][4][2]>=0)&&(NeiOrder[2][4][2]>=0)&&(NeiOrder[0][2][2]>=0)){
            NeiOrder[1][3][2] = currTOrder1;
            currTOrder1++;
        }
        if ((NeiOrder[4][4][2]>=0)&&(NeiOrder[2][4][2]>=0)&&(NeiOrder[4][2][2]>=0)){
            NeiOrder[3][3][2] = currTOrder1;
            currTOrder1++;
        }
        //-------------------------------------------------------------
        if ((LayerSelf_lv>=0)&&(LayerFront_lv>=0)&&(LayerFront_front_v>=0)){    //220, 020, 022 -> 121
            LayerFront_lt = currTOrder;
            currTOrder++;
        }
        if ((LayerSelf_lv>=0)&&(LayerBack_lv>=0)&&(LayerBack_back_v>=0)){   //220, 420, 422 -> 321
            LayerBack_lt = currTOrder;
            currTOrder++;
        }
        if ((LayerSelf_rv>=0)&&(LayerFront_rv>=0)&&(LayerFront_front_v>=0)){    //224, 024, 022 -> 123
            LayerFront_rt = currTOrder;
            currTOrder++;
        }
        if ((LayerSelf_rv>=0)&&(LayerBack_rv>=0)&&(LayerBack_back_v>=0)){   //224, 424, 422 -> 323
            LayerBack_rt = currTOrder;
            currTOrder++;
        }
        if ((NeiOrder[2][2][0]>=0)&&(NeiOrder[0][2][0]>=0)&&(NeiOrder[0][2][2]>=0)){
            NeiOrder[1][2][1] = currTOrder1;
            currTOrder1++;
        }
        if ((NeiOrder[2][2][0]>=0)&&(NeiOrder[4][2][0]>=0)&&(NeiOrder[4][2][2]>=0)){
            NeiOrder[3][2][1] = currTOrder1;
            currTOrder1++;
        }
        if ((NeiOrder[0][2][2]>=0)&&(NeiOrder[0][2][4]>=0)&&(NeiOrder[2][2][4]>=0)){
            NeiOrder[1][2][3] = currTOrder1;
            currTOrder1++;
        }
        if ((NeiOrder[2][2][4]>=0)&&(NeiOrder[4][2][4]>=0)&&(NeiOrder[4][2][2]>=0)){
            NeiOrder[3][2][3] = currTOrder1;
            currTOrder1++;
        }
        
        if (checkcellorder(NeiOrder, neighborOrder)) cout<<"3. NeiOrder and neighborOrder are not the same !!!"<<endl;
        
        

        if ((LayerFront_uv>=0)&&(LayerFront_ulv>=0)&&(LayerFront_lv>=0)&&(LayerFront_front_v>=0) && (LayerSelf_uv>=0)&&(LayerSelf_ulv>=0)&&(LayerSelf_lv>=0)){
            //  002, 000, 020, 022, 202, 200, 220 -> 111
            LayerFront_ulc = currCOrder;
            currCOrder++;
        }
        if ((NeiOrder[0][0][2]>=0) && (NeiOrder[0][0][0]>=0)&&(NeiOrder[0][2][0]>=0)&&(NeiOrder[0][2][2]>=0)&&(NeiOrder[2][0][0]>=0)&&(NeiOrder[2][0][2]>=0)&& (NeiOrder[2][2][0]>=0)){
            NeiOrder[1][1][1] = currCOrder1;
            currCOrder1++;
        }
        if (checkcellorder(NeiOrder, neighborOrder)) cout<<"4.1 NeiOrder and neighborOrder are not the same !!!"<<endl;
        
        if ((LayerBack_uv>=0)&&(LayerBack_ulv>=0)&&(LayerBack_lv>=0)&&(LayerBack_back_v>=0) && (LayerSelf_uv>=0)&&(LayerSelf_ulv>=0)&&(LayerSelf_lv>=0)){
            // 402, 400, 420, 422, 202, 200, 220 -> 311
            LayerBack_ulc = currCOrder;
            currCOrder++;
        }
        if ((NeiOrder[2][0][0]>=0)&&(NeiOrder[4][0][0]>=0)&&(NeiOrder[4][0][2]>=0)&&(NeiOrder[2][0][2]>=0) && (NeiOrder[2][2][0]>=0)&&(NeiOrder[4][2][0]>=0)&&(NeiOrder[4][2][2]>=0)){
            NeiOrder[3][1][1] = currCOrder1;
            currCOrder1++;
        }
        if (checkcellorder(NeiOrder, neighborOrder)) cout<<"4.2 NeiOrder and neighborOrder are not the same !!!"<<endl;
        
        if ((LayerFront_uv>=0)&&(LayerFront_urv>=0)&&(LayerFront_rv>=0)&&(LayerFront_front_v>=0) && (LayerSelf_uv>=0)&&(LayerSelf_urv>=0)&&(LayerSelf_rv>=0)){
            // 002, 004, 024, 022, 202, 204, 224 -> 113
            LayerFront_urc = currCOrder;
            currCOrder++;
        }
        if ((NeiOrder[0][0][2]>=0)&&(NeiOrder[0][0][4]>=0)&&(NeiOrder[0][2][4]>=0)&&(NeiOrder[0][2][2]>=0)&&(NeiOrder[2][0][2]>=0)&&(NeiOrder[2][0][4]>=0) && (NeiOrder[2][2][4]>=0)){
            NeiOrder[1][1][3] = currCOrder1;
            currCOrder1++;
        }
        if (checkcellorder(NeiOrder, neighborOrder)) cout<<"4.3 NeiOrder and neighborOrder are not the same !!!"<<endl;
        
        if ((LayerBack_uv>=0)&&(LayerBack_urv>=0)&&(LayerBack_rv>=0)&&(LayerBack_back_v>=0) && (LayerSelf_uv>=0)&&(LayerSelf_urv>=0)&&(LayerSelf_rv>=0)){
            // 402, 404, 424, 422, 202, 204, 224 -> 313
            LayerBack_urc = currCOrder;
            currCOrder++;
        }
        if ((NeiOrder[4][0][2]>=0)&&(NeiOrder[4][0][4]>=0)&&(NeiOrder[4][2][4]>=0)&&(NeiOrder[4][2][2]>=0)&&(NeiOrder[2][0][2]>=0)&&(NeiOrder[2][0][4]>=0) && (NeiOrder[2][2][4]>=0)){
            NeiOrder[3][1][3] = currCOrder1;
            currCOrder1++;
        }
        if (checkcellorder(NeiOrder, neighborOrder)) cout<<"4.4 NeiOrder and neighborOrder are not the same !!!"<<endl;
        
        
        if ((LayerFront_dv>=0)&&(LayerFront_lv>=0)&&(LayerFront_dlv>=0)&&(LayerFront_front_v>=0) && (LayerSelf_dv>=0)&&(LayerSelf_lv>=0)&&(LayerSelf_dlv>=0)){
            // 042, 020, 040, 022, 242, 220, 240 -> 131
            LayerFront_dlc = currCOrder;
            currCOrder++;
        }
        if ((NeiOrder[0][4][2]>=0)&&(NeiOrder[0][2][0]>=0)&&(NeiOrder[0][4][0]>=0)&&(NeiOrder[0][2][2]>=0)&&(NeiOrder[2][4][0]>=0)&&(NeiOrder[2][4][2]>=0)&&(NeiOrder[2][2][0]>=0)){
            NeiOrder[1][3][1] = currCOrder1;
            currCOrder1++;
        }
        if (checkcellorder(NeiOrder, neighborOrder)) cout<<"4.5 NeiOrder and neighborOrder are not the same !!!"<<endl;
        
        if ((LayerBack_dv>=0)&&(LayerBack_lv>=0)&&(LayerBack_dlv>=0)&&(LayerBack_back_v>=0) && (LayerSelf_dv>=0)&&(LayerSelf_lv>=0)&&(LayerSelf_dlv>=0)){
            //  442, 420, 440, 422, 242, 220, 240 -> 331
            LayerBack_dlc = currCOrder;
            currCOrder++;
        }
        if ((NeiOrder[4][4][2]>=0)&&(NeiOrder[4][2][0]>=0)&&(NeiOrder[4][2][2]>=0)&&(NeiOrder[4][4][0]>=0)&&(NeiOrder[2][4][2]>=0) &&(NeiOrder[2][2][0]>=0)&& (NeiOrder[2][4][0]>=0)){
            NeiOrder[3][3][1] = currCOrder1;
            currCOrder1++;
        }
        if (checkcellorder(NeiOrder, neighborOrder)) cout<<"4.6 NeiOrder and neighborOrder are not the same !!!"<<endl;
        
        if ((LayerFront_dv>=0)&&(LayerFront_rv>=0)&&(LayerFront_drv>=0)&&(LayerFront_front_v>=0) && (LayerSelf_dv>=0)&&(LayerSelf_rv>=0)&&(LayerSelf_drv>=0)){
            //  042, 024, 044, 022, 242, 224, 244 -> 133
            LayerFront_drc = currCOrder;
            currCOrder++;
        }
        if ((NeiOrder[0][4][2]>=0)&&(NeiOrder[0][2][4]>=0)&&(NeiOrder[0][4][4]>=0)&&(NeiOrder[0][2][2]>=0)&&(NeiOrder[2][4][2]>=0)&&(NeiOrder[2][4][4]>=0) && (NeiOrder[2][2][4]>=0)){
            NeiOrder[1][3][3] = currCOrder1;
            currCOrder1++;
        }
        if (checkcellorder(NeiOrder, neighborOrder)) cout<<"4.7 NeiOrder and neighborOrder are not the same !!!"<<endl;
        
        if ((LayerBack_dv>=0)&&(LayerBack_rv>=0)&&(LayerBack_drv>=0)&&(LayerBack_back_v>=0) && (LayerSelf_dv>=0)&&(LayerSelf_rv>=0)&&(LayerSelf_drv>=0)){
            //  442, 424, 444, 422, 242, 224, 244 -> 333
            LayerBack_drc = currCOrder;
            currCOrder++;
        }
        if ((NeiOrder[4][4][2]>=0)&&(NeiOrder[4][2][4]>=0)&&(NeiOrder[4][4][4]>=0)&&(NeiOrder[4][2][2]>=0)&&(NeiOrder[2][4][2]>=0)&&(NeiOrder[2][2][4]>=0) && (NeiOrder[2][4][4]>=0)){
            NeiOrder[3][3][3] = currCOrder1;
            currCOrder1++;
        }
        if (checkcellorder(NeiOrder, neighborOrder)) cout<<"4.8 NeiOrder and neighborOrder are not the same !!!"<<endl;
        
        
        
        //===================================================================================================================================
        
        
        

        

		// copy the neighboring orders back
//		for (i=rid-2;i<=rid+2;i++){
//			if ((i<0)||(i>=mapNRow)) continue;
//			for(j=cid-2;j<=cid+2;j++){
//				if ((j<0)||(j>=mapNCol)) continue;
//				if (cellOrder[i][j]!=neighborOrder[i+2-rid][j+2-cid]){  //if =, then skip the following operation to imporve the efficiency.
//					myassert( cellOrder[i][j]==-1 );
//					myassert( cellType[i][j]!=VERTEX );
//					cellOrder[i][j]=neighborOrder[i+2-rid][j+2-cid];
//				}
//			}
//		}
        
        for (k=sid-2;k<=sid+2;k++){ // those two vertices in front and two vertices in back.
            if ((k<0)||(k>=mapNSlice)) continue;
            for (i=rid-2;i<=rid+2;i++){ // those two vertices above and two vertices below.
                if ((i<0)||(i>=mapNRow)) continue;
                for(j=cid-2;j<=cid+2;j++){  // those two vertices on the left and two vertices on the right.
                    if ((j<0)||(j>=mapNCol)) continue;
                    // copy the order in cellOrder to neighborOrder, but indices in neighborOrder was changed to 0~5.
//                    if(neighborOrder[k+2-sid][i+2-rid][j+2-cid]==-1){
//                        cout<<"Find -1 in neighborHorder !!!!!!!!!! with sid: "<<sid<<", rid: "<<rid<<", cid: "<<cid<<endl;
//                        
//                    }
                    cellOrder[k][i][j]=neighborOrder[k+2-sid][i+2-rid][j+2-cid];  //[0~4][0~4][0~4] = [sid-2 ~ sid+2][rid-2 ~ rid+2][cid-2 ~ cid+2].
                }
            }
        }
        for (k=sid-2;k<=sid+2;k++){ // those two vertices in front and two vertices in back.
            if ((k<0)||(k>=mapNSlice)) continue;
            for (i=rid-2;i<=rid+2;i++){ // those two vertices above and two vertices below.
                if ((i<0)||(i>=mapNRow)) continue;
                for(j=cid-2;j<=cid+2;j++){  // those two vertices on the left and two vertices on the right.
                    if ((j<0)||(j>=mapNCol)) continue;
                    // copy the order in cellOrder to neighborOrder, but indices in neighborOrder was changed to 0~5.
                    //                    if(neighborOrder[k+2-sid][i+2-rid][j+2-cid]==-1){
                    //                        cout<<"Find -1 in neighborHorder !!!!!!!!!! with sid: "<<sid<<", rid: "<<rid<<", cid: "<<cid<<endl;
                    //
                    //                    }
                    cellOrder_xudong[k][i][j]=NeiOrder[k+2-sid][i+2-rid][j+2-cid];  //[0~4][0~4][0~4] = [sid-2 ~ sid+2][rid-2 ~ rid+2][cid-2 ~ cid+2].
                }
            }
        }
        
//        cout<<"------------------------ 2nd. cellOrder -----------------------"<<endl;
//		show2dvector(cellOrder);   //------->

	}
    
    vector< vector< vector < int > > > cellOrder_xudong;
    
    void xudongcheckeleincellOrder(vector< vector< vector < int > > > cellOrder, int vertNum, int edgeNum, int trigNum, int cubeNum){
        cout<<"Number of Slides: "<<cellOrder.size()<<endl;
        cout<<"Number of rows: "<<cellOrder[0].size()<<endl;
        cout<<"Number of cols: "<<cellOrder[0][0].size()<<endl;
        int totalelements = cellOrder.size() * cellOrder[0].size() * cellOrder[0][0].size();
        vector<int> temp(totalelements, 0);
        for(int k=0;k<cellOrder.size();k++){
            for (int i=0;i<cellOrder[0].size();i++){
                for (int j=0;j<cellOrder[0][0].size();j++){
                    temp[cellOrder[k][i][j]]++;
                }
            }
        }
        cout<<"# of vertex = "<<vertNum<<", # of edge = "<<edgeNum<<", # of trig = "<<trigNum<<", # of cube = "<<cubeNum<<endl;
        for(int k=0;k<totalelements-1;k++){
            if (temp[k]==4 && temp[k+1]!=4){
                cout<<"1st step is at: "<<k+1<<endl;
            }
            else if (temp[k]==3 && temp[k+1]!=3){
                cout<<"2nd step is at: "<<k+1<<endl;
            }
            else if (temp[k]==2 && temp[k+1]!=2){
                cout<<"3rd step is at: "<<k+1<<endl;
            }
            else if (temp[k]==1 && temp[k+1]!=1){
                cout<<"4th step is at: "<<k+1<<endl;
            }
        }
//        show1dvector(temp,2);
        
    }
    //---------------- xudong checked the function setVertOrder(), and it's good. ------------------>

	CellMap(vector< Vertex > * vL, vector< Edge > * eL, vector< Triangle > * tL, vector< Cube > * cL){
		int k,i,j,idx;    //idx was defined here.

		vList=vL;
		eList=eL;
		trigList=tL;
        cubeList=cL;    //xudong add this line.
        
		myassert(eList->size()==0);
		myassert(trigList->size()==0);
		
//		vertNum = vList->size();    //2D case.
//		edgeNum = phi3d->nrow *(phi3d->ncol - 1) + phi3d->ncol*(phi3d->nrow - 1);   //2D case.
//		trigNum = (phi3d->ncol - 1) * (phi3d->nrow - 1);    //2D case.
        
        int mm = phi3d->nrow;
        int nn = phi3d->ncol;
        int ll = phi3d->nslice;
        vertNum = vList->size();
//        cout<<"The vertNum from vListSize: "<<vertNum<<", vertNum from m*n*l: "<<mm*nn*ll<<endl; //3D case.
        edgeNum = (mm*(nn-1)+nn*(mm-1))*ll + mm*nn*(ll-1);
//        cout<<"The edgeNum: "<<edgeNum<<endl;  //3D case.
        trigNum = (mm-1)*(nn-1)*ll + (ll-1)*(mm*(nn-1)+nn*(mm-1));
//        cout<<"The triNum: "<<trigNum<<endl; //3D case.
        cubeNum = (mm-1)*(nn-1)*(ll-1);
//        cout<<"The cubeNum: "<<cubeNum<<endl; //3D case.
//        cout<<"Total number of elements from sum is: "<<vertNum+edgeNum+trigNum+cubeNum<<endl;
//        cout<<"Total number of elements from calculating is: "<<(2*mm-1)*(2*nn-1)*(2*ll-1)<<endl;

		mapNRow = 2*mm-1;  // mapNRow = 2m-1, vertex|edgehori|vertex...
		mapNCol = 2*nn-1;  // mapNCol = 2n-1, edgevert|triangle|edgevert...
        mapNSlice = 2*ll-1;
		currEOrder=0;
		currTOrder=0;
        currCOrder=0;
        
        currEOrder1=0;
        currTOrder1=0;
        currCOrder1=0;

//		int i, j;
        // assign dimension and element types for the following three 3D matrices.
		cellOrder.assign(mapNSlice, vector< vector<int> >(mapNRow, vector< int >(mapNCol, -1)));    //each entry of cellOrder was filled will -1.
        cellOrder_xudong.assign(mapNSlice, vector< vector<int> >(mapNRow, vector< int >(mapNCol, -1)));
		edgePersType.assign(mapNSlice, vector< vector<ETrigPersTypeEnum> >(mapNRow, vector< ETrigPersTypeEnum >(mapNCol, EP_UNDEFINED)));
        trigPersType.assign(mapNSlice, vector< vector<ETrigPersTypeEnum> >(mapNRow, vector< ETrigPersTypeEnum >(mapNCol, EP_UNDEFINED)));
		cellFixType.assign(mapNSlice, vector< vector<CellFixTypeEnum> >(mapNRow, vector< CellFixTypeEnum >(mapNCol, CF_UNDEFINED)));

		// create cellType row by row for different cases. --> vector< vector< vector < CellTypeEnum > > > cellType;
        //--> CellTypeEnum {0.CT_UNDEFINED, 1.VERTEX, 2.EDGEVERT, 3.EDGEHORI, 4.ZEDGE, 5.TRIG, 6.TRIGHORI, 7.TRIGVERT, 8.CUBE};
        vector< CellTypeEnum > 	zeroRow_slice_1( mapNCol, VERTEX );
        vector< CellTypeEnum > 	firstRow_slice_1( mapNCol, EDGEVERT );
        vector< CellTypeEnum > 	zeroRow_slice_2( mapNCol, ZEDGE );
        vector< CellTypeEnum > 	firstRow_slice_2( mapNCol, TRIGVERT );
        for (i = 0; i < mapNCol; i++){
            if ( i % 2 == 1 ){
                zeroRow_slice_1[i] = EDGEHORI;    // odd entry will be assigned to EDGEHORI.
                firstRow_slice_1[i] = TRIG;   // odd entry will be assigned to TRIG.
                zeroRow_slice_2[i] = TRIGHORI;   // odd entry will be assigned to TRIG.
                firstRow_slice_2[i] = CUBE;   // odd entry will be assigned to CUBE.
            }
        }
        //write different types into 2D matrix slice_1 and slice_2.
        vector< vector <CellTypeEnum> > slice_1;
        vector< vector <CellTypeEnum> > slice_2;
        for (i = 0; i < mapNRow; i++){
            if (i % 2 == 0){
                slice_1.push_back(zeroRow_slice_1);
                slice_2.push_back(zeroRow_slice_2);
            }
            else{
                slice_1.push_back(firstRow_slice_1);
                slice_2.push_back(firstRow_slice_2);
            }
                
        }
        //write different types into cellType.
        for (k = 0; k < mapNSlice; k++)
            if (k % 2 == 0){
                //xudong: write types for even slice.
                cellType.push_back(slice_1);
            }
            else{
                //xudong: write types for odd slice.
                cellType.push_back(slice_2);
            }
//        cout<<"-------------------------------------- print cellType below -------------------------------------"<<endl;
//        print3dmatrix(cellType);        // xudong: ok till here., xudong: check the logic of code and it's good.
//        cout<<"-------------------------------------- print cellType above -------------------------------------"<<endl;
		

//        showvList(vList);
//        vector<Vertex> *temp1 = vList;
//        show2dvector(cellOrder);   // all the all the (2m-1)*(2n-1) were filled with -1.
        
        for( i = 0; i<vertNum; i++){
//        for( i = 0; i<1000; i++){
			setVertOrder((* vList)[i].zidx*2, (* vList)[i].xidx*2 ,(* vList)[i].yidx*2, i);  //write vertexorder,edgeorder,trigorder,cubeorder into cellOrder 3D matrix.
//            cout<<"------------------------------------- "<<i<<" ----------------------------------------"<<endl;
        }
        

        if (checkcellorder(cellOrder_xudong, cellOrder)) cout<<" cellOrder_xudong and cellOrder are different "<<endl;
      
        
//        cout<<"currEOrder: "<<currEOrder<<endl;cout<<"edgeNum: "<<edgeNum<<endl;    //full 3D image ->currEOrder:251068416, edgeNum: 251068416 -> cellType[][][] is good
//        cout<<"currTOrder: "<<currTOrder<<endl;cout<<"trigNum: "<<trigNum<<endl;    //full 3D image ->currTOrder:250479936, trigNum: 250479936 -> cellType[][][] is good.
//        cout<<"currCOrder: "<<currCOrder<<endl;cout<<"cubeNum: "<<cubeNum<<endl;
        if(currEOrder != edgeNum){cout<<"BAD! currEOrder != edgeNum"<<endl;}
        if(currTOrder != trigNum){cout<<"BAD! currTOrder != trigNum"<<endl;}
        if(currCOrder != cubeNum){cout<<"BAD! currCOrder != cubeNum"<<endl;}
        if(currEOrder1 != edgeNum){cout<<"BAD! currEOrder1 != edgeNum"<<endl;}
        if(currTOrder1 != trigNum){cout<<"BAD! currTOrder1 != trigNum"<<endl;}
        if(currCOrder1 != cubeNum){cout<<"BAD! currCOrder1 != cubeNum"<<endl;}
        
		myassert(currEOrder == edgeNum);
		myassert(currTOrder == trigNum);
//        cout<<endl<<"-------------------------------------- print cellOrder below-------------------------------------"<<endl;
//        print3dmatrix(cellOrder);            // xudong: ok till here.
//        cout<<"-------------------------------------- print cellOrder above -------------------------------------"<<endl<<endl;
        
//        vector<Vertex> *temp2 = vList;    // vList does not change.
//        comList(temp1, temp2);
//        showvList(vList);
//        cout<<"Good"<<endl;
//        showedge(e1); //checking the constructor of class Edge, it's good.
//        showtri(t1);  //checking the constructor of class triangle, it's good.
//        showcube(c1); //checking the constructor of class cube, it's good.
//        xudongcheckeleincellOrder(cellOrder, vertNum, edgeNum, trigNum, cubeNum);   // check those 4 sequences of element order,it's good.
cout<<"======================== 2. pass ===================="<<endl;
		eList->assign(edgeNum,Edge());
		trigList->assign(trigNum,Triangle());   // assign dimension and element types. all the values are -1, cannot be used, since mapZCoord = -1, which is wrong !!!
        cubeList->assign(cubeNum,Cube());
 
//		int idx;
//		Line 590~609: build eList and trigList,
//      CellTypeEnum {CT_UNDEFINED, 1.VERTEX, 2.EDGEVERT, 3.EDGEHORI, 4.ZEDGE, 5.TRIG, 6.TRIGHORI, 7.TRIGVERT, 8.CUBE};
//      Vertex(int zi,int xi,int yi)
        int edgeNum_from_cellOrder =0;
        int triNum_from_cellOrder =0;
        int cubeNum_from_cellOrder =0;
        
        for(k=0;k<mapNSlice;k++){
            for (i=0;i<mapNRow;i++){
                for (j=0;j<mapNCol;j++){
                    idx = cellOrder[k][i][j];  //idx should be the order of one of vertex, (edgehori & edgevert), trig, or cube.
                    myassert(idx>=0); if(idx < 0) cout<<" Error !!! -> \"idx < 0\" Happened."<<endl;
//                    Vertex(int zi, int xi,int yi) : zidx(zi),xidx(xi),yidx(yi),mapZCoord(zi*2),mapXCoord(xi*2),mapYCoord(yi*2){}
//                    Edge(int v1o,int v2o, int mz, int mx, int my) : v1_order(min(v1o,v2o)),v2_order(max(v1o,v2o)),mapZCoord(mz),mapXCoord(mx),mapYCoord(my){}
                    if (cellType[k][i][j]==EDGEHORI){
                        myassert(idx<edgeNum);  if(idx >= edgeNum) cout<<"Error idx >= edgeNum"<<endl;
                        edgeNum_from_cellOrder++;   // counting edge number.
                        // store right and left vertices as (min, max), and the position in cellorder matrix.
//                        (* eList)[idx]=Edge(left,right,k,i,j);
                        (* eList)[idx]=Edge(cellOrder[k][i][j-1],cellOrder[k][i][j+1],k,i,j);
                    }
                    else if(cellType[k][i][j]==EDGEVERT){
                        myassert(idx<edgeNum);  if(idx >= edgeNum) cout<<"Error idx >= edgeNum"<<endl;
                        edgeNum_from_cellOrder++;   // counting edge number.
                        // store up and down vertices as (min, max), and the position in cellorder matrix.
//                        (* eList)[idx]=Edge(up,down,k,i,j);
                        (* eList)[idx]=Edge(cellOrder[k][i-1][j],cellOrder[k][i+1][j],k,i,j);
                    }
                    else if(cellType[k][i][j]==ZEDGE){
                        myassert(idx<edgeNum);  if(idx >= edgeNum) cout<<"Error idx >= edgeNum"<<endl;
                        edgeNum_from_cellOrder++;   // counting edge number.
                        // store front and back vertices as (min, max), and the position in cellorder matrix.
//                        (* eList)[idx]=Edge(front,back,k,i,j);
                        (* eList)[idx]=Edge(cellOrder[k-1][i][j],cellOrder[k+1][i][j],k,i,j);
                    }
                }
            }
            
        }
        if(edgeNum > edgeNum_from_cellOrder){cout<<"Missing edge when extracted from cellOrder 3D matrix : "<<endl;}
		sort(eList->begin(),eList->end(),eComp);
        for(i=0;i<edgeNum;i++){ // sort(eList) -> edge order=0~edgeNum.
			cellOrder[(*eList)[i].mapZCoord][(*eList)[i].mapXCoord][(*eList)[i].mapYCoord]=i;  // update cellOrder according to the sorted eList.
        }
//        cout<<"======================== 3 xudong: pass ===================="<<endl;

        
        
        
        
        
        for(k=0;k<mapNSlice;k++){
            for (i=0;i<mapNRow;i++){
                for (j=0;j<mapNCol;j++){
                    idx = cellOrder[k][i][j];  //idx should be the order of one of vertex, (edgehori & edgevert), trig, or cube.
                    myassert(idx>=0); if(idx < 0) cout<<" Error !!! -> \"idx < 0\" Happened."<<endl;
                    
                    //                    Triangle(int v1o,int v2o,int v3o,int v4o,int e1o,int e2o,int e3o,int e4o,int mz,int mx,int my)
                    else if(cellType[k][i][j]==TRIG){
                        myassert(idx<trigNum);  if(idx >= trigNum) cout<<"Error idx >= trigNum"<<endl;
                        //                        cout<<"k="<<k<<", i="<<i<<", j="<<j<<endl;  //check the output, this line ok.
                        //                        cout<<up_left<<endl;cout<<down_left<<endl;cout<<up_right<<endl;cout<<down_right<<endl;  //check the output, this line ok.
                        //                        cout<<"k="<<k<<", i-1="<<i-1<<", j="<<j<<", up="<<up<<endl; // this line wrong.     int up = cellOrder[k][i-1][j];
                        //                        cout<<"k="<<k<<", i-1="<<i-1<<", j="<<j<<", cell type is (3): "<<cellType[k][i-1][j]<<", cellOrder: "<<cellOrder[k][i-1][j]<<", up: "<<up<<endl;
                        //                        cout<<down<<endl;   // this line wrong.     int down = cellOrder[k][i+1][j];
                        //                        cout<<left<<endl;   //this line ok.
                        //                        cout<<right<<endl;  //this line ok.
                        //                        (* trigList)[idx]=Triangle(cellOrder[i-1][j-1],cellOrder[i+1][j-1],cellOrder[i-1][j+1],cellOrder[i+1][j+1],
                        //                                                   cellOrder[i-1][j],cellOrder[i+1][j],cellOrder[i][j-1],cellOrder[i][j+1],
                        //                                                   i,j);
                        
                        //                        (* trigList)[idx]=Triangle(up_left,down_left,up_right,down_right,up,down,left,right,k,i,j);
                        triNum_from_cellOrder++;  // counting triangle number.
                        (* trigList)[idx]=Triangle(cellOrder[k][i-1][j-1],cellOrder[k][i+1][j-1],cellOrder[k][i-1][j+1],cellOrder[k][i+1][j+1],cellOrder[k][i-1][j],cellOrder[k][i+1][j],cellOrder[k][i][j-1],cellOrder[k][i][j+1],k,i,j);
                        //                        cout<<triNum_from_cellOrder<<", Enter TRIG  idx: "<<idx; showt((* trigList)[idx]);
                    }
                    else if(cellType[k][i][j]==TRIGHORI){
                        myassert(idx<trigNum);  if(idx >= trigNum) cout<<"Error idx >= trigNum"<<endl;
                        //                        (* trigList)[idx]=Triangle(front_left,front_right,back_left,back_right,front,back,left,right,k,i,j);
                        triNum_from_cellOrder++;  // counting triangle number.
                        (* trigList)[idx]=Triangle(cellOrder[k-1][i][j-1],cellOrder[k-1][i][j+1],cellOrder[k+1][i][j-1],cellOrder[k+1][i][j+1],cellOrder[k-1][i][j],cellOrder[k+1][i][j],cellOrder[k][i][j-1],cellOrder[k][i][j+1],k,i,j);
                        //                        cout<<triNum_from_cellOrder<<", Enter TRIGHORI  idx: "<<idx; showt((* trigList)[idx]);
                    }
                    else if(cellType[k][i][j]==TRIGVERT){
                        myassert(idx<trigNum);  if(idx >= trigNum) cout<<"Error idx >= trigNum"<<endl;
                        //                        (* trigList)[idx]=Triangle(front_up,front_down,back_up,back_down, up,down,front,back,k,i,j);
                        triNum_from_cellOrder++;  // counting triangle number.
                        (* trigList)[idx]=Triangle(cellOrder[k-1][i-1][j],cellOrder[k-1][i+1][j],cellOrder[k+1][i-1][j],cellOrder[k+1][i+1][j], cellOrder[k][i-1][j],cellOrder[k][i+1][j],cellOrder[k-1][i][j],cellOrder[k+1][i][j],k,i,j);
                        //                        cout<<triNum_from_cellOrder<<", Enter TRIGVERT  idx: "<<idx; showt((* trigList)[idx]);
                    }
                }
            }
            
        }
        

        if(trigNum > triNum_from_cellOrder) {cout<<"Missing trig when extracted from cellOrder 3D matrix"<<endl;}

        sort(trigList->begin(),trigList->end(),trigComp);
        for(i=0;i<trigNum;i++){
            //            cout<<"Ok when i = "<<i<<endl;  // when i=255 -> mapZCoord = -1, which is wrong !!!
            cellOrder[(*trigList)[i].mapZCoord][(*trigList)[i].mapXCoord][(*trigList)[i].mapYCoord]=i;    // update cellOrder according to the sorted trigList.
        }
        //cout<<"======================== 4 xudong: pass ===================="<<endl;
        
        
        
        
        for(k=0;k<mapNSlice;k++){
            for (i=0;i<mapNRow;i++){
                for (j=0;j<mapNCol;j++){
                    idx = cellOrder[k][i][j];  //idx should be the order of one of vertex, (edgehori & edgevert), trig, or cube.
                    myassert(idx>=0); if(idx < 0) cout<<" Error !!! -> \"idx < 0\" Happened."<<endl;
                    
//                    Format of cube is: Cube(int v1o,int v2o,int v3o,int v4o,int v5o,int v6o,int v7o,int v8o,
//                    int e1o,int e2o,int e3o,int e4o,int e5o,int e6o,int e7o,int e8o,int e9o,int e10o,int e11o,int e12o,
//                    int t1o, int t2o, int t3o, int t4o, int t5o, int t6o, int mz,int mx,int my)
                    else if(cellType[k][i][j]==CUBE){
                        myassert(idx<cubeNum);  if(idx >= cubeNum) cout<<"Error idx >= cubeNum "<<endl;
                        cubeNum_from_cellOrder++;   // counting cube number.
                        
                        (* cubeList)[idx]=Cube(cellOrder[k-1][i-1][j-1],cellOrder[k-1][i-1][j+1],cellOrder[k-1][i+1][j-1],cellOrder[k-1][i+1][j+1],   cellOrder[k+1][i-1][j-1],cellOrder[k+1][i-1][j+1],cellOrder[k+1][i+1][j-1],cellOrder[k+1][i+1][j+1], cellOrder[k][i-1][j-1],cellOrder[k][i+1][j-1],cellOrder[k][i-1][j+1],cellOrder[k][i+1][j+1],   cellOrder[k-1][i-1][j],cellOrder[k-1][i+1][j],cellOrder[k-1][i][j-1],cellOrder[k-1][i][j+1], cellOrder[k+1][i-1][j],cellOrder[k+1][i+1][j],cellOrder[k+1][i][j-1],cellOrder[k+1][i][j+1],cellOrder[k][i-1][j],cellOrder[k][i+1][j],cellOrder[k][i][j-1],cellOrder[k][i][j+1],cellOrder[k-1][i][j],cellOrder[k+1][i][j],k,i,j);
                    }
                }
            }
            
        }
        
        if(cubeNum > cubeNum_from_cellOrder){cout<<"Missing cube when extracted from cellOrder 3D matrix"<<endl;}
        sort(cubeList->begin(),cubeList->end(),CubeComp);
        
        for(i=0;i<cubeNum;i++){
            cellOrder[(*cubeList)[i].mapZCoord][(*cubeList)[i].mapXCoord][(*cubeList)[i].mapYCoord]=i;    // update cellOrder according to the sorted trigList.
        }
//        cout<<"======================== 5 xudong: pass ===================="<<endl;
        

        
        
        
        
        
        
        

	}

    
    
    vector<int> birthPoint1D(vector<int> chordRange, mVpM_2dMatrix * const phi2d_mV, const int NofPoints, mVpM_2dMatrix * const phi2d_pM){
        vector<int> birPoi_1D;
        int mV_index, temp_z, temp_x, temp_y, temp;
//         vector< vector< int > > coord;
        set<int> index_list;
        int mV_list_size = phi2d_mV->coord.size();
        while(index_list.size()<NofPoints){
            mV_index = rand() % mV_list_size;
            temp_x = (phi2d_mV->coord)[mV_index][0];
            temp_y = (phi2d_mV->coord)[mV_index][1];
            temp_z = (phi2d_mV->coord)[mV_index][2];
            
            if(chordRange[2]<=temp_x && temp_x<=chordRange[3] && chordRange[4]<=temp_y && temp_y<=chordRange[5] && chordRange[0]<=temp_z && temp_z<=chordRange[1])
                index_list.insert(mV_index);
        }
        
        for(set<int>::iterator ii=index_list.begin(); ii!=index_list.end(); ii++){
            cout<<"index_list: "<<*ii<<endl;
        }
        
        
        for(set<int>::iterator ii=index_list.begin(); ii!=index_list.end(); ii++){
            temp_z = (phi2d_mV->coord)[*ii][2];
            temp_x = (phi2d_mV->coord)[*ii][0];
            temp_y = (phi2d_mV->coord)[*ii][1];
            cout<<"mV -> "<<" z: "<<temp_z+1<<", x: "<<temp_x+1<<", y: "<<temp_y+1<<endl;
            cout<<"The intensity of mV is: "<<phi3d->data[temp_z][temp_x][temp_y]<<endl;
            temp = cellOrder[temp_z*2][temp_x*2][temp_y*2];
            birPoi_1D.push_back(temp);
            
        }
        
        
        
        for(int i=0;i<phi2d_pM->coord.size();i++){
            temp_z = (phi2d_pM->coord)[i][2];
            temp_x = (phi2d_pM->coord)[i][0];
            temp_y = (phi2d_pM->coord)[i][1];
            cout<<"pM_tips -> "<<" z: "<<temp_z+1<<", x: "<<temp_x+1<<", y: "<<temp_y+1<<endl;
            cout<<"The intensity of tips is: "<<phi3d->data[temp_z][temp_x][temp_y]<<endl;
            temp = cellOrder[temp_z*2][temp_x*2][temp_y*2];
            birPoi_1D.push_back(temp);
            
        }
        return birPoi_1D;
        
    }
    
    
    void buildBoundary3D(vector<vector < int > > * boundary_3D){    // square the boundary, which is edge.
        int i, j, idx;
        for (i=0; i<cubeNum; i++){
            (* boundary_3D)[i].push_back( (* cubeList)[i].t1_order );
            (* boundary_3D)[i].push_back( (* cubeList)[i].t2_order );
            (* boundary_3D)[i].push_back( (* cubeList)[i].t3_order );
            (* boundary_3D)[i].push_back( (* cubeList)[i].t4_order );
            (* boundary_3D)[i].push_back( (* cubeList)[i].t5_order );
            (* boundary_3D)[i].push_back( (* cubeList)[i].t6_order );

        }
    }
    //*************************** xudong add this code ***********************
    void showCoord(int e){
        int v1 = (* eList)[e].v1_order;
        int v2 = (* eList)[e].v2_order;
        cout<<"v1-> x: "<<(*vList)[v1].xidx<<", y: "<<(*vList)[v1].yidx<<", z: "<<(*vList)[v1].zidx<<endl;
        cout<<"v2-> x: "<<(*vList)[v2].xidx<<", y: "<<(*vList)[v2].yidx<<", z: "<<(*vList)[v2].zidx<<endl;
    }
    //*************************** xudong add this code ***********************
	void buildBoundary2D(vector<vector < int > > * boundary_2D){    // square the boundary, which is edge.
		int i, j, idx;
        bool status=false;
		for (i=0; i<trigNum; i++){
            (* boundary_2D)[i].push_back( (* trigList)[i].e1_order );
			(* boundary_2D)[i].push_back( (* trigList)[i].e2_order );
			(* boundary_2D)[i].push_back( (* trigList)[i].e3_order );
			(* boundary_2D)[i].push_back( (* trigList)[i].e4_order );
//            if((* trigList)[i].e1_order == 1370037){cout<<"The i of e1_order is: "<<i<<endl;status=true;}   //xudong add.
//            else if((* trigList)[i].e2_order == 1370037){cout<<"The i of e2_order is: "<<i<<endl;status=true;}  //xudong add.
//            else if((* trigList)[i].e3_order == 1370037){cout<<"The i of e3_order is: "<<i<<endl;status=true;}  //xudong add.
//            else if((* trigList)[i].e4_order == 1370037){cout<<"The i of e4_order is: "<<i<<endl;status=true;}  //xudong add.
//            if(true == status){
//                showCoord((* trigList)[i].e1_order );   //xudong add.
//                showCoord((* trigList)[i].e2_order );   //xudong add.
//                showCoord((* trigList)[i].e3_order );   //xudong add.
//                showCoord((* trigList)[i].e4_order );   //xudong add.
//                status = false;
//            }
		}
	}
    
	void buildBoundary1D(vector<vector < int > > * boundary_1D){
		int i, j, idx;
		for (i=0; i<edgeNum; i++){
			(* boundary_1D)[i].push_back( (* eList)[i].v1_order );
			(* boundary_1D)[i].push_back( (* eList)[i].v2_order );
		}
	}

    //identify which edge-tri pair was creator, which as destroyer.
//	void setEPersType( vector< int > * low_2D_e2t ){    //xudong: original code for 2D image was deleted.
    void setEPersType( vector< int > * low_2D_e2t ){
        myassert(low_2D_e2t->size() == edgeNum);
        int k,i,j,idx;
        for(k=0;k<mapNSlice;k++){
            for (i=0; i< mapNRow; i++){
                for (j=0; j< mapNCol; j++){
                    if ((cellType[k][i][j]==EDGEVERT)||(cellType[k][i][j]==EDGEHORI)||(cellType[k][i][j]==ZEDGE)){
                        idx=cellOrder[k][i][j];    //obtain the index of edge.
                        if ((* low_2D_e2t)[idx]!=-1)
                            edgePersType[k][i][j] = CREATOR;
                        else
                            edgePersType[k][i][j] = DESTROYER;
                    }
                }
            }
        }
    
    }
    
    
    void setTPersType( vector< int > * low_3D_t2c ){    //xudong: code for 3D image.
		myassert(low_3D_t2c->size() == edgeNum);
		int k,i,j,idx;
        for(k=0;k<mapNSlice;k++){
            for (i=0; i< mapNRow; i++){
                for (j=0; j< mapNCol; j++){
//                    enum CellTypeEnum {CT_UNDEFINED, VERTEX, EDGEVERT, EDGEHORI, ZEDGE, TRIG, TRIGHORI, TRIGVERT, CUBE};
                    if ((cellType[k][i][j]==TRIG)||(cellType[k][i][j]==TRIGHORI)||(cellType[k][i][j]==TRIGVERT)){
                        idx=cellOrder[k][i][j];    //obtain the index of edge.
                        if ((* low_3D_t2c)[idx]!=-1)
                            trigPersType[k][i][j] = CREATOR;    //edgePersType[k][i][j] will be used in merge-comp, flood_comp, merge_hole, flood_hole.
                        else
                            trigPersType[k][i][j] = DESTROYER;  //edgePersType[k][i][j] will be used in merge-comp, flood_comp, merge_hole, flood_hole.
                    }
                }
            }
        }
        
	}
 /*
    void setVertFixType(int rid, int cid, CellFixTypeEnum ftype){
        myassert(ftype != CF_UNDEFINED);
        
        myassert( (rid >=0)&&(rid < mapNRow)&&(cid >= 0)&&(cid < mapNCol) );
        
        myassert( cellType[rid][cid] == VERTEX );
        
        if(cellFixType[rid][cid] == ftype) return;	//alreay fixed
        
        myassert(cellFixType[rid][cid] == CF_UNDEFINED);	//cell was either unknown, or the same as ftype
        
        cellFixType[rid][cid] = ftype;	//fix the vertex first
        
        vector< vector < CellFixTypeEnum > > neighborFixType(5, vector< CellFixTypeEnum >(5,CF_UNDEFINED));
        
        int i,j;
        
        // get the neighboring orders
        for (i=rid-2;i<=rid+2;i++){
            if ((i<0)||(i>=mapNRow)) continue;
            for(j=cid-2;j<=cid+2;j++){
                if ((j<0)||(j>=mapNCol)) continue;
                neighborFixType[i+2-rid][j+2-cid]=cellFixType[i][j];
            }
        }
        
        // 		CellFixTypeEnum ue,de,le,re,ult,urt,dlt,drt;
        // 		CellFixTypeEnum uv,dv,lv,rv,ulv,urv,dlv,drv;
        CellFixTypeEnum & uv = neighborFixType[0][2];
        CellFixTypeEnum & dv = neighborFixType[4][2];
        CellFixTypeEnum & lv = neighborFixType[2][0];
        CellFixTypeEnum & rv = neighborFixType[2][4];
        CellFixTypeEnum & ulv = neighborFixType[0][0];
        CellFixTypeEnum & urv = neighborFixType[0][4];
        CellFixTypeEnum & dlv = neighborFixType[4][0];
        CellFixTypeEnum & drv = neighborFixType[4][4];
        
        CellFixTypeEnum & ue = neighborFixType[1][2];
        CellFixTypeEnum & de = neighborFixType[3][2];
        CellFixTypeEnum & le = neighborFixType[2][1];
        CellFixTypeEnum & re = neighborFixType[2][3];
        CellFixTypeEnum & ult = neighborFixType[1][1];
        CellFixTypeEnum & urt = neighborFixType[1][3];
        CellFixTypeEnum & dlt = neighborFixType[3][1];
        CellFixTypeEnum & drt = neighborFixType[3][3];
        
        if (uv==ftype){
            myassert(ue==CF_UNDEFINED);
            ue = ftype;
        }
        if (dv==ftype){
            myassert(de==CF_UNDEFINED);
            de = ftype;
        }
        if (lv==ftype){
            myassert(le==CF_UNDEFINED);
            le = ftype;
        }
        if (rv==ftype){
            myassert(re==CF_UNDEFINED);
            re = ftype;
        }
        if ((uv==ftype)&&(ulv==ftype)&&(lv==ftype)){
            myassert(ult == CF_UNDEFINED);
            ult = ftype;
        }
        if ((uv==ftype)&&(urv==ftype)&&(rv==ftype)){
            myassert(urt == CF_UNDEFINED);
            urt = ftype;
        }
        if ((dv==ftype)&&(dlv==ftype)&&(lv==ftype)){
            myassert(dlt == CF_UNDEFINED);
            dlt = ftype;
        }
        if ((dv==ftype)&&(drv==ftype)&&(rv==ftype)){
            myassert(drt == CF_UNDEFINED);
            drt = ftype;
        }
        
        // copy the neighboring orders back
        for (i=rid-2;i<=rid+2;i++){
            if ((i<0)||(i>=mapNRow)) continue;
            for(j=cid-2;j<=cid+2;j++){
                if ((j<0)||(j>=mapNCol)) continue;
                if (cellFixType[i][j]!=neighborFixType[i+2-rid][j+2-cid]){
                    myassert( cellFixType[i][j]==CF_UNDEFINED );
                    myassert( cellType[i][j]!=VERTEX );
                    cellFixType[i][j]=neighborFixType[i+2-rid][j+2-cid];
                }
            }
        }
        
    }
 */

    //xudong: getFixType() won't be used in finding the chords.
	CellFixTypeEnum getFixType( int dim, int idx ){
		int mapZ, mapX, mapY;
		if(dim == 0){
            mapZ = (* vList)[idx].mapZCoord;
			mapX = (* vList)[idx].mapXCoord;
			mapY = (* vList)[idx].mapYCoord;
			return cellFixType[mapZ][mapX][mapY];
		};
		if(dim == 1){
            mapZ = (* eList)[idx].mapZCoord;
			mapX = (* eList)[idx].mapXCoord;
			mapY = (* eList)[idx].mapYCoord;
			return cellFixType[mapZ][mapX][mapY];
		};
		if(dim == 2){
            mapZ = (* trigList)[idx].mapZCoord;
			mapX = (* trigList)[idx].mapXCoord;
			mapY = (* trigList)[idx].mapYCoord;
			return cellFixType[mapZ][mapX][mapY];
		};
        if(dim == 3){
            mapZ = (* cubeList)[idx].mapZCoord;
            mapX = (* cubeList)[idx].mapXCoord;
            mapY = (* cubeList)[idx].mapYCoord;
            return cellFixType[mapZ][mapX][mapY];
        };
		myassert(false);
	}
    
    
    
    //-----------------************-------------------------------------------------------------------------------//
    // xudong has modified the code above for suiting 3d case.-------------------------~~~------------------------//
    //-----------------************-------------------------------------------------------------------------------//
    
//    //change the intensity on point (zcoord, xcoord, ycoord) and its neighbors (with radius=n_pix) to max_total.
//    void changeIntensity(const int n_pix,vector< vector< vector< double > > > &data, int zcoord, int xcoord, int ycoord, double max_total){
////        double tempMax = max(max1, max2);
//        if(0 == n_pix){
//            data[zcoord][xcoord][ycoord]=max_total;   // only change the intensity of point at [zcoord][xcoord][ycoord].
//        }
//        else {
//            int range = n_pix*2+1;
//            for (int k=0;k<range;k++){ //make sure: 0 <= zcoord-n_pix+k <= mapNSlice-1 <=> 0 <= zcoord-n_pix+k < mapNSlice
//                if ((zcoord-n_pix+k < 0)||(zcoord-n_pix+k >= mapNSlice)){ continue;}
//                for (int i=0;i<range;i++){
//                    if ((xcoord-n_pix+i < 0)||(xcoord-n_pix+i >= mapNRow)){ continue;}
//                    for(int j=0;j<range;j++){
//                        if ((ycoord-n_pix+j < 0)||(ycoord-n_pix+j >= mapNCol)){ continue;}
//                        data[zcoord-n_pix+k][xcoord-n_pix+i][ycoord-n_pix+j]=max_total;  //
//                    }
//                }
//            }
//        }
//        
//    }
    
    
    
    //========================== From merge_comp() ============================================================
    //  I used three method to terminate the depth-first search.
    //  1. if the id of tmpv is less than the id of vBirth, then stop.
    //  2. if the intensity of tmpv is less than the intensity of vBirth, then stop.
    //  3. if the intensity of tmpv is less than the min(changeToIntensity_mV, changeToIntensity_pM), then stop.
    //  Method 1 and 2 could find very similar chords, but 2 will search less points than 1, since
    //  there are many voxels have the same intensity., hence method 1 will search more voxels.
    //   But, m1 and m2 needs to assign at lease two point, one from mV and another was from pM.
    //  Since, m2 use intensity as the searching stop criterion, I use m3 (do not need to use two point, but)
    //  intensity to be the searching stop criterion.
    //  m1: if (tmpv<=vBirth) { break;}
    //  m2: if (tmpv_intensity <= vBirth_intensity) { break;}
    //  m3 if (tmpv_intensity <= min_change) { break;}
    //==========================================================================================================

//#define showinfo
//#define vBirth_1
////    #define vBirth_intensity_2
////    #define min_change_3   // if set two birthPoints, one from mV and another one from pM. Then can find two identical chords.

	vector<path> merge_comp(
                    int vBirth, //vBirth is a vertex selected from mV or pM.
                    double min_change,    //
                    int eDeath,

                   my3dMatrix * const perturbM,

                   double changeToIntensity,
                   int & counter_xu
//                   bool &chords_temp1_write,
//                   int Len_chords_ub,
//                   int Len_chords_lb,
//                   vector< vector< int > > pM_marker,
//                   vector< vector< int > > mV_marker,
//                   vector< path > &critchords1d,
//                   bool & critchords1d_write,
//                   vector< int > chordRange,
//                   const int tol
                   )
    {
        vector<path> chord;
        vector< path > temp_chords_1, temp_chords_2;    //store the path in two directions into temp_chords_1 and temp_chords_2.
//#ifdef showinfo
//        cout<<"0. eDeath: "<<eDeath<<endl;
//        cout<<"1. vBirth: "<<vBirth<<", with intensity : "<<phi3d->data[(* vList)[vBirth].zidx][(* vList)[vBirth].xidx][(* vList)[vBirth].yidx]<<endl;
//        cout<<"2. intensity of min_change : "<<min_change<<endl;
//
//#endif

        
        
        int pointCounter = 0;
//        chords_temp1_write = false;   //prevent writing short chord in next loop.
//        critchords1d_write = false; //prevent writing wrong chord in next loop.
        
		//starting from edeath, go through a path of destroyer edges,  connecting vBirth and an component born earlier, and only go with edge with lower order
		// for example (38, 104355, 10.0561, -19.9998, 10.0561) --> vbidx: 38, edidx: 104355, robustness: 10.0561, birth: -19.9998, death: 10.0561
		vector< bool > vert_visited(vertNum,false);
		stack< int > vert_stack;

        int tmpv1 = (* eList)[eDeath].v1_order;
//        cout<<"tmpv1 id:"<<tmpv1<<endl;
        int tmpv2 = (* eList)[eDeath].v2_order;
//        cout<<"tmpv2 id:"<<tmpv2<<endl;
        int zcoord,xcoord,ycoord, mapzcoord,mapxcoord,mapycoord,tmpv;
        vector< int > parent_v(vertNum,-1);
        double tmpv_intensity, vBirth_intensity;
//        double tmpv_intensity;
        int childVertex;
        int xudongcounter1 = 0;
        double max1, max2, max_total, min_total;
        

        
        
        //xudong: add the following code.
//        double vbirth_intensity=phi->data[(* vList)[vBirth].xidx][(* vList)[vBirth].yidx];
//        cout<<"vBirth: "<<vBirth<<", xcoord: "<<(* vList)[vBirth].xidx<<", ycoord: "<<(* vList)[vBirth].yidx<<", vbirth_intensity: "<<vbirth_intensity<<endl; //-20
//        double tempv1_intensity=phi->data[(* vList)[tmpv1].xidx][(* vList)[tmpv1].yidx];
//        double tempv2_intensity=phi->data[(* vList)[tmpv2].xidx][(* vList)[tmpv2].yidx];
//        cout<<"tmpv1: "<<tmpv1<<", xcoord: "<<(* vList)[tmpv1].xidx<<", ycoord: "<<(* vList)[tmpv1].yidx<<", tempv1_intensity: "<<tempv1_intensity<<endl; //-20<<endl;
//        cout<<"tmpv2: "<<tmpv2<<", xcoord: "<<(* vList)[tmpv2].xidx<<", ycoord: "<<(* vList)[tmpv2].yidx<<", tempv2_intensity: "<<tempv2_intensity<<endl;
        //xudong: add the above code.
  
        
        
        
//        //            tmpv_intensity = phi3d->data[(* vList)[tmpv].zidx][(* vList)[tmpv].xidx][(* vList)[tmpv].yidx];
//        //            vBirth_intensity = phi3d->data[(* vList)[vBirth].zidx][(* vList)[vBirth].xidx][(* vList)[vBirth].yidx];
//        //            cout<<"tmpv_intensity: "<<tmpv_intensity<<", vBirth_intensity: "<<vBirth_intensity<<endl;
//        //            if (tmpv<=vBirth) { //1.
//        if (tmpv_intensity <= -900) { //2.

//            //                cout<<"======================== 18: break!!! happened in finding tmpv1. ===================="<<endl;
//            //                cout<<"************ tmpv<=vBirth happpen, then break!!! ***************"<<endl;
//            //                cout<<"tmpv-> "<<tmpv<<": "<<phi3d->data[(* vList)[tmpv].zidx][(* vList)[tmpv].xidx][(* vList)[tmpv].yidx]<<endl;
//            //                cout<<"vBirth-> "<<vBirth<<": "<<phi3d->data[(* vList)[vBirth].zidx][(* vList)[vBirth].xidx][(* vList)[vBirth].yidx]<<endl;
//            counter_xu++;
//            break;  //found the first end, then break the while loop. xudong: the first end was not pop fromt he stack!!!!
//        }
        
        
        
        
        
        
        
        
        
        // 1.find the path from tmpv1.
//        cout<<"======================== 17: start find the path from tmpv1. ===================="<<endl;
        vert_stack.push(tmpv1);
        while(! vert_stack.empty()){
            xudongcounter1++;
            tmpv = vert_stack.top();
//            if(tmpv2==tmpv){ cout<<"tmpv2==tmpv"<<endl;}
            if (vert_visited[tmpv]){    // if this vertex has been visited, pop it then visit the vertex stored below.
                vert_stack.pop();
                continue;
            };
//            if(xudongcounter1%1000==0){cout<<"father id:"<<tmpv<<endl;}
            if(tmpv2==tmpv){break;}   //1. Good
//            if(tmpv2==tmpv){ continue;} // BAD->non-stop loop, since always pick up tmpv2 from stack top, then loop again at this command line.
//            if(tmpv2==tmpv){ vert_stack.pop();continue;}   //1. BAD->even this could finish, but visit too much point, then tmpv2 does not have voxel to visit in next step.
            
            vert_visited[tmpv]=true;    // mark this vertex to be true if haven't been visited.
//            if(vert_visited[tmpv2]) {cout<<"tmpv2 was set to be true on last command line."<<endl; break;} //2.
            zcoord=(* vList)[tmpv].zidx;
            xcoord=(* vList)[tmpv].xidx;
            ycoord=(* vList)[tmpv].yidx;
            
            mapzcoord=(* vList)[tmpv].mapZCoord;
            mapxcoord=(* vList)[tmpv].mapXCoord;
            mapycoord=(* vList)[tmpv].mapYCoord;
            myassert(cellType[mapzcoord][mapxcoord][mapycoord]==VERTEX);
            myassert((cellFixType[mapzcoord][mapxcoord][mapycoord]==CF_UNDEFINED)||(cellFixType[mapzcoord][mapxcoord][mapycoord]==MOVEDOWN));
            
            // 			if ((tmpv==vBirth)||(tmpv==vHigh)) break;	//found the first end
            tmpv_intensity = phi3d->data[(* vList)[tmpv].zidx][(* vList)[tmpv].xidx][(* vList)[tmpv].yidx];
            vBirth_intensity = phi3d->data[(* vList)[vBirth].zidx][(* vList)[vBirth].xidx][(* vList)[vBirth].yidx];
//#ifdef vBirth_1
//            if (tmpv<=vBirth) {
//                cout<<"1. define vBirth"<<endl;
//#endif
//#ifdef vBirth_intensity_2
//            if (tmpv_intensity <= vBirth_intensity) {
//                cout<<"2. define vBirth_intensity: "<<vBirth_intensity<<endl;
//#endif
//#ifdef min_change_3
            if (tmpv_intensity <= min_change) {
//                cout<<"3. define value_900"<<endl;
//#endif
                
                //                cout<<"======================== 18: break!!! happened in finding tmpv1. ===================="<<endl;
                //                cout<<"************ tmpv<=vBirth happpen, then break!!! ***************"<<endl;
                //                cout<<"tmpv: "<<phi->data[(* vList)[tmpv].xidx][(* vList)[tmpv].yidx]<<endl;
                //                cout<<"vBirth: "<<phi->data[(* vList)[vBirth].xidx][(* vList)[vBirth].yidx]<<endl;;
                counter_xu++;
                break;  //found the first end, then break the while loop. xudong: the first end was not pop fromt he stack!!!!
            }
            

            
            vert_stack.pop();
            //            cout<<"check the vertex: "<<tmpv<<endl; //53101, 52444,
            //            cout<<"---> eDeath:"<<eDeath<<endl; //==104355
            //            cout<<"---> up edge id:"<<cellOrder[mapxcoord-1][mapycoord]<<", edgePersType:"<<edgePersType[mapxcoord-1][mapycoord]<<endl;
            //            cout<<"---> down edge id:"<<cellOrder[mapxcoord+1][mapycoord]<<", edgePersType:"<<edgePersType[mapxcoord+1][mapycoord]<<endl;
            //            cout<<"---> left edge id:"<<cellOrder[mapxcoord][mapycoord-1]<<", edgePersType:"<<edgePersType[mapxcoord][mapycoord-1]<<endl;
            //            cout<<"---> right edge id:"<<cellOrder[mapxcoord][mapycoord+1]<<", edgePersType:"<<edgePersType[mapxcoord][mapycoord+1]<<endl;
            //check all neighbor edges, push relevant vertices into the queue;
            
            //up vertex
            //            if ((mapxcoord>0)&&(edgePersType[mapzcoord][mapxcoord-1][mapycoord]==DESTROYER)&&(cellOrder[mapzcoord][mapxcoord-1][mapycoord]<eDeath)){
            if ((mapxcoord-2>=0)&&(edgePersType[mapzcoord][mapxcoord-1][mapycoord]==DESTROYER)&&(cellOrder[mapzcoord][mapxcoord-1][mapycoord]<eDeath)){
                //                cout<<"enter up"<<endl;
                childVertex = cellOrder[mapzcoord][mapxcoord-2][mapycoord];
                if (! vert_visited[childVertex]){ //check whether the previous vertex has already visited according to vert_visited[][].
                    //                    int childVertex = cellOrder[mapzcoord][mapxcoord-2][mapycoord];
                    myassert(parent_v[childVertex]<0);
                    parent_v[childVertex]=tmpv;
                    vert_stack.push(childVertex);
                    //                    cout<<"---> up Happen, add vertex: "<<childVertex<<", with intensity: "<<phi->data[(* vList)[childVertex].xidx][(* vList)[childVertex].yidx]<<endl;
                }
            }
            //down vertex
            if ((mapxcoord+2<mapNRow)&&(edgePersType[mapzcoord][mapxcoord+1][mapycoord]==DESTROYER)&&(cellOrder[mapzcoord][mapxcoord+1][mapycoord]<eDeath)){
                //                cout<<"enter down"<<endl;
                childVertex = cellOrder[mapzcoord][mapxcoord+2][mapycoord];
                if (! vert_visited[childVertex]){ //check whether the previous vertex has already visited according to vert_visited[][].
                    //                    int childVertex = cellOrder[mapzcoord][mapxcoord+2][mapycoord];
                    myassert(parent_v[childVertex]<0);
                    parent_v[childVertex]=tmpv;
                    vert_stack.push(childVertex);
                    //                    cout<<"---> down Happen, add vertex: "<<childVertex<<", with intensity: "<<phi->data[(* vList)[childVertex].xidx][(* vList)[childVertex].yidx]<<endl;
                }
            }
            //left vertex
            if ((mapycoord-2>=0)&&(edgePersType[mapzcoord][mapxcoord][mapycoord-1]==DESTROYER)&&(cellOrder[mapzcoord][mapxcoord][mapycoord-1]<eDeath)){
                //               cout<<"enter left"<<endl;
                childVertex = cellOrder[mapzcoord][mapxcoord][mapycoord-2];
                if (! vert_visited[childVertex]){ //check whether the previous vertex has already visited according to vert_visited[][].
                    //                    int childVertex = cellOrder[mapzcoord][mapxcoord][mapycoord-2];
                    myassert(parent_v[childVertex]<0);
                    parent_v[childVertex]=tmpv;
                    vert_stack.push(childVertex);
                    //                    cout<<"---> left Happen, add vertex: "<<childVertex<<", with intensity: "<<phi->data[(* vList)[childVertex].xidx][(* vList)[childVertex].yidx]<<endl;
                }
            }
            //right vertex
            if ((mapycoord+2<mapNCol)&&(edgePersType[mapzcoord][mapxcoord][mapycoord+1]==DESTROYER)&&(cellOrder[mapzcoord][mapxcoord][mapycoord+1]<eDeath)){
                //                cout<<"enter right"<<endl;
                childVertex = cellOrder[mapzcoord][mapxcoord][mapycoord+2];
                if (! vert_visited[childVertex]){ //check whether the previous vertex has already visited according to vert_visited[][].
                    //                    int childVertex = cellOrder[mapzcoord][mapxcoord][mapycoord+2];
                    myassert(parent_v[childVertex]<0);
                    parent_v[childVertex]=tmpv;
                    vert_stack.push(childVertex);
                    //                    cout<<"---> right Happen, add vertex: "<<childVertex<<", with intensity: "<<phi->data[(* vList)[childVertex].xidx][(* vList)[childVertex].yidx]<<endl;
                }
            }
            //front vertex
            if ((mapzcoord-2>=0)&&(edgePersType[mapzcoord-1][mapxcoord][mapycoord]==DESTROYER)&&(cellOrder[mapzcoord-1][mapxcoord][mapycoord]<eDeath)){
                //                cout<<"enter right"<<endl;
                childVertex = cellOrder[mapzcoord-2][mapxcoord][mapycoord];
                if (! vert_visited[childVertex]){ //check whether the previous vertex has already visited according to vert_visited[][].
                    //                    int childVertex = cellOrder[mapzcoord-2][mapxcoord][mapycoord];
                    myassert(parent_v[childVertex]<0);
                    parent_v[childVertex]=tmpv;
                    vert_stack.push(childVertex);
                    //                    cout<<"---> right Happen, add vertex: "<<childVertex<<", with intensity: "<<phi->data[(* vList)[childVertex].xidx][(* vList)[childVertex].yidx]<<endl;
                }
            }
            //back vertex.
            if ((mapzcoord+2<mapNSlice)&&(edgePersType[mapzcoord+1][mapxcoord][mapycoord]==DESTROYER)&&(cellOrder[mapzcoord+1][mapxcoord][mapycoord]<eDeath)){
                //                cout<<"enter right"<<endl;
                childVertex = cellOrder[mapzcoord+2][mapxcoord][mapycoord];
                if (! vert_visited[childVertex]){ //check whether the previous vertex has already visited according to vert_visited[][].
                    //                    int childVertex = cellOrder[mapzcoord+2][mapxcoord][mapycoord];
                    myassert(parent_v[childVertex]<0);
                    parent_v[childVertex]=tmpv;
                    vert_stack.push(childVertex);
                    //                    cout<<"---> right Happen, add vertex: "<<childVertex<<", with intensity: "<<phi->data[(* vList)[childVertex].xidx][(* vList)[childVertex].yidx]<<endl;
                }
            }
                
            //            if(xudongcounter1 == 30)break;
            
            //            cout<<"/////////////////////////////////-> vert_stack.empty(): "<<!vert_stack.empty()<<endl;
        };	//end of while(! vert_stack.empty())
        //        cout<<"======================== 19: finish finding the path from tmpv1. ===================="<<endl;
        //        show1dvector(parent_v,2);
        //        cout<<"/////////////////////////////////-----------------------/////////////////////////"<<endl;
        //        cout<<"myassert: false:"<<endl;myassert(false); //no output.
        //        cout<<"myassert: ture:"<<endl;myassert(true);   //no output.
        myassert(! vert_stack.empty());		//vert_stack stores the first path.
        vector< bool > pathv_visited(vertNum,false);
        //		cout<<"1 ---> vert_stack is empty:  "<<vert_stack.empty()<<endl;
        int endV1 = -1;
        if(! vert_stack.empty()){
            endV1 = vert_stack.top();	//first end v still in the stack, connected to tmpv1
        }
        while(! vert_stack.empty()){ // clear vert_stack.
            int temptop = vert_stack.top();
            //            cout<<"The top element in vert_stack: "<<temptop<<", with intensity: "<<phi->data[(* vList)[temptop].xidx][(* vList)[temptop].yidx]<<endl;
            vert_stack.pop();
        }
        tmpv = endV1;
        int oldtmpv=tmpv;
        while(tmpv >= 0){
            //move down everything on the path starting from the first end point for all the vertices stroed in parent_v .
            myassert(! pathv_visited[tmpv]);
            pathv_visited[tmpv]=true;
            
            zcoord=(* vList)[tmpv].zidx;
            xcoord=(* vList)[tmpv].xidx;
            ycoord=(* vList)[tmpv].yidx;
            mapzcoord=(* vList)[tmpv].mapZCoord;
            mapxcoord=(* vList)[tmpv].mapXCoord;
            mapycoord=(* vList)[tmpv].mapYCoord;
            myassert(cellType[mapzcoord][mapxcoord][mapycoord]==VERTEX);
            myassert((cellFixType[mapzcoord][mapxcoord][mapycoord]==CF_UNDEFINED)||(cellFixType[mapzcoord][mapxcoord][mapycoord]==MOVEDOWN));
            

//            max1 = perturbM->data[tmpv1.zidx][tmpv1.xidx][tmpv1.yidx];
//            max2 = perturbM->data[zcoord][xcoord][ycoord];
//            max_total = max(max1,max2);
//            min_total = min(max1,max2);
#ifdef drawinmerge_comp
            changeInteCore(0,perturbM->data, zcoord, xcoord, ycoord, changeToIntensity);    // draw from endV1 to tmpv1 (eDeath_vertex1).
            pointCounter++;
#endif
            //            changeIntensity(n_pix,perturbM->data, zcoord, xcoord, ycoord, max1, max2);  //xudong: add this the change the intensity on 3 directions.
            temp_chords_1.push_back(path(zcoord, xcoord, ycoord, 0, 0)); //xudong add this line to store the chord path.--------------->>>>>
            
            oldtmpv=tmpv;
            tmpv = parent_v[tmpv];
        };
        myassert( oldtmpv == tmpv1 );
        
        
        
        // 2.find the path from tmpv2.
//                cout<<"======================== 20: start find the path from tmpv2. ===================="<<endl;
        //        cout<<"--------->>>>>> mapNSlice: "<<mapNSlice<<"--------->>>>>> mapNRow: "<<mapNRow<<"--------->>>>>> mapNCol: "<<mapNCol<<endl;
        myassert(! vert_visited[tmpv2]); if(vert_visited[tmpv2]!=false) cout<<"Wrong!!! tmpv2 already visited in previous calculation."<<endl;
        vert_stack.push(tmpv2);	//search path connecting tmpv2 to someone
        while(! vert_stack.empty()){
//                        cout<<"======================== 20.1 ===================="<<endl;
            tmpv = vert_stack.top();
            if (vert_visited[tmpv]){
//                                cout<<"======================== 20.2 ===================="<<endl;
                vert_stack.pop();
//                                cout<<"======================== 20.3 ===================="<<endl;
                continue;
            };
            vert_visited[tmpv]=true;
//                        cout<<"======================== 20.4 ===================="<<endl;
            zcoord=(* vList)[tmpv].zidx;
            xcoord=(* vList)[tmpv].xidx;
            ycoord=(* vList)[tmpv].yidx;
            
            mapzcoord=(* vList)[tmpv].mapZCoord;
            mapxcoord=(* vList)[tmpv].mapXCoord;
            mapycoord=(* vList)[tmpv].mapYCoord;
            myassert(cellType[mapzcoord][mapxcoord][mapycoord]==VERTEX);
            myassert((cellFixType[mapzcoord][mapxcoord][mapycoord]==CF_UNDEFINED)||(cellFixType[mapzcoord][mapxcoord][mapycoord]==MOVEDOWN));
            
            tmpv_intensity = phi3d->data[(* vList)[tmpv].zidx][(* vList)[tmpv].xidx][(* vList)[tmpv].yidx];
            vBirth_intensity = phi3d->data[(* vList)[vBirth].zidx][(* vList)[vBirth].xidx][(* vList)[vBirth].yidx];
//#ifdef vBirth_1
//            if (tmpv<=vBirth) {
////                cout<<"1. define vBirth"<<endl;
//#endif
//#ifdef vBirth_intensity_2
//            if (tmpv_intensity <= vBirth_intensity) {
////                cout<<"2. define vBirth_intensity"<<endl;
//#endif
//#ifdef min_change_3
            if (tmpv_intensity <= min_change) {
//                cout<<"3. define value_900"<<endl;
//#endif

                counter_xu++;
                break;	//found the first end
            }
            vert_stack.pop();
            
            //up vertex
            if ((mapxcoord-2 >= 0)&&(edgePersType[mapzcoord][mapxcoord-1][mapycoord]==DESTROYER)&&(cellOrder[mapzcoord][mapxcoord-1][mapycoord]<eDeath)){
                //                cout<<"edgePersType[mapzcoord][mapxcoord-1][mapycoord]: "<<mapzcoord<<", "<<mapxcoord-1<<", "<<mapycoord<<", "<<edgePersType[mapzcoord][mapxcoord-1][mapycoord]<<endl;
                //                cout<<"cellOrder[mapzcoord][mapxcoord-1][mapycoord]: "<<cellOrder[mapzcoord][mapxcoord-1][mapycoord]<<endl;
                //                cout<<"childVertex = cellOrder[mapzcoord][mapxcoord-2][mapycoord]: "<<cellOrder[mapzcoord][mapxcoord-2][mapycoord]<<endl;
                //                cout<<"enter up 1"<<endl;
                childVertex = cellOrder[mapzcoord][mapxcoord-2][mapycoord];
//                cout<<"enter up 2"<<endl;
                if (! vert_visited[childVertex]){ //check whether the previous vertex has already visited according to vert_visited[][].
                    //                    int childVertex = cellOrder[mapzcoord][mapxcoord-2][mapycoord];
                    myassert(parent_v[childVertex]<0);
                    parent_v[childVertex]=tmpv;
                    //                    cout<<"======================== 20.7 ===================="<<endl;
                    vert_stack.push(childVertex);
//                                        cout<<"======================== 20.8 ===================="<<endl;
                    //                    cout<<"---> up Happen, add vertex: "<<childVertex<<", with intensity: "<<phi->data[(* vList)[childVertex].xidx][(* vList)[childVertex].yidx]<<endl;
                }}
            
            
            
            //down vertexs
            //            cout<<"enter down 0"<<endl;
            //            if ((mapxcoord<mapNRow-1)&&(edgePersType[mapzcoord][mapxcoord+1][mapycoord]==DESTROYER)&&(cellOrder[mapzcoord][mapxcoord+1][mapycoord]<eDeath)){
            if ((mapxcoord+2 < mapNRow)&&(edgePersType[mapzcoord][mapxcoord+1][mapycoord]==DESTROYER)&&(cellOrder[mapzcoord][mapxcoord+1][mapycoord]<eDeath)){
                //                cout<<"mapNRow-1: "<<mapNRow-1<<endl;
                //                cout<<"edgePersType[mapzcoord][mapxcoord+1][mapycoord]: "<<mapzcoord<<", "<<mapxcoord+1<<", "<<mapycoord<<", "<<edgePersType[mapzcoord][mapxcoord+1][mapycoord]<<endl;
                //                cout<<"cellOrder[mapzcoord][mapxcoord+1][mapycoord]: "<<cellOrder[mapzcoord][mapxcoord+1][mapycoord]<<endl;
                //                cout<<"childVertex = cellOrder[mapzcoord][mapxcoord+2][mapycoord]: "<<cellOrder[mapzcoord][mapxcoord+2][mapycoord]<<endl;
                
//                cout<<"enter down 1"<<endl;
                childVertex = cellOrder[mapzcoord][mapxcoord+2][mapycoord];
                //                cout<<"enter down 2"<<endl;
                if (! vert_visited[childVertex]){ //check whether the previous vertex has already visited according to vert_visited[][].
                    //                    int childVertex = cellOrder[mapzcoord][mapxcoord+2][mapycoord];
                    myassert(parent_v[childVertex]<0);
                    parent_v[childVertex]=tmpv;
                    //                    cout<<"======================== 20.9 ===================="<<endl;
                    vert_stack.push(childVertex);
//                                        cout<<"======================== 20.10 ===================="<<endl;
                    //                    cout<<"---> down Happen, add vertex: "<<childVertex<<", with intensity: "<<phi->data[(* vList)[childVertex].xidx][(* vList)[childVertex].yidx]<<endl;
                }}
            
            
            
            //left vertex
            //            cout<<"enter left 0"<<endl;
            //            if ((mapycoord>0)&&(edgePersType[mapzcoord][mapxcoord][mapycoord-1]==DESTROYER)&&(cellOrder[mapzcoord][mapxcoord][mapycoord-1]<eDeath)){
            if ((mapycoord-2 >= 0)&&(edgePersType[mapzcoord][mapxcoord][mapycoord-1]==DESTROYER)&&(cellOrder[mapzcoord][mapxcoord][mapycoord-1]<eDeath)){
                //                cout<<"edgePersType[mapzcoord][mapxcoord][mapycoord-1]: "<<mapzcoord<<", "<<mapxcoord<<", "<<mapycoord-1<<", "<<edgePersType[mapzcoord][mapxcoord][mapycoord-1]<<endl;
                //                cout<<"cellOrder[mapzcoord][mapxcoord][mapycoord-1]: "<<cellOrder[mapzcoord][mapxcoord][mapycoord-1]<<endl;
                //                cout<<"childVertex = cellOrder[mapzcoord][mapxcoord][mapycoord-2]: "<<cellOrder[mapzcoord][mapxcoord][mapycoord-2]<<endl;
                //
//                                cout<<"enter left 1"<<endl;
                childVertex = cellOrder[mapzcoord][mapxcoord][mapycoord-2];
                //                cout<<"enter left 2"<<endl;
                if (! vert_visited[childVertex]){ //check whether the previous vertex has already visited according to vert_visited[][].
                    //                    int childVertex = cellOrder[mapzcoord][mapxcoord][mapycoord-2];
                    myassert(parent_v[childVertex]<0);
                    parent_v[childVertex]=tmpv;
                    //                    cout<<"======================== 20.11 ===================="<<endl;
                    vert_stack.push(childVertex);
//                                        cout<<"======================== 20.12 ===================="<<endl;
                    //                    cout<<"---> left Happen, add vertex: "<<childVertex<<", with intensity: "<<phi->data[(* vList)[childVertex].xidx][(* vList)[childVertex].yidx]<<endl;
                }}
            
            
            
            //right vertex
            //            cout<<"enter right 0"<<endl;
            //            if ((mapycoord<mapNCol-1)&&(edgePersType[mapzcoord][mapxcoord][mapycoord+1]==DESTROYER)&&(cellOrder[mapzcoord][mapxcoord][mapycoord+1]<eDeath)){
            if ((mapycoord+2 < mapNCol)&&(edgePersType[mapzcoord][mapxcoord][mapycoord+1]==DESTROYER)&&(cellOrder[mapzcoord][mapxcoord][mapycoord+1]<eDeath)){
                //                cout<<"mapNCol-1: "<<mapNCol-1<<endl;
                //                cout<<"edgePersType[mapzcoord][mapxcoord][mapycoord+1]: "<<mapzcoord<<", "<<mapxcoord<<", "<<mapycoord+1<<", "<<edgePersType[mapzcoord][mapxcoord][mapycoord+1]<<endl;
                //                cout<<"cellOrder[mapzcoord][mapxcoord][mapycoord+1]: "<<cellOrder[mapzcoord][mapxcoord][mapycoord+1]<<endl;
                //                cout<<"childVertex = cellOrder[mapzcoord][mapxcoord][mapycoord+2]: "<<cellOrder[mapzcoord][mapxcoord][mapycoord+2]<<endl;
                
//                                cout<<"enter right 1"<<endl;
                childVertex = cellOrder[mapzcoord][mapxcoord][mapycoord+2];
                //                cout<<"enter right 2"<<endl;
                if (! vert_visited[childVertex]){ //check whether the previous vertex has already visited according to vert_visited[][].
                    //                    int childVertex = cellOrder[mapzcoord][mapxcoord][mapycoord+2];
                    myassert(parent_v[childVertex]<0);
                    parent_v[childVertex]=tmpv;
                    //                    cout<<"======================== 20.13 ===================="<<endl;
                    vert_stack.push(childVertex);
//                                        cout<<"======================== 20.14 ===================="<<endl;
                    //                    cout<<"---> right Happen, add vertex: "<<childVertex<<", with intensity: "<<phi->data[(* vList)[childVertex].xidx][(* vList)[childVertex].yidx]<<endl;
                }}
            
            
            
            //front vertex
            //            cout<<"enter front 0"<<endl;
            //            if ((mapzcoord>0)&&(edgePersType[mapzcoord-1][mapxcoord][mapycoord]==DESTROYER)&&(cellOrder[mapzcoord-1][mapxcoord][mapycoord]<eDeath)){
            if ((mapzcoord-2 >= 0)&&(edgePersType[mapzcoord-1][mapxcoord][mapycoord]==DESTROYER)&&(cellOrder[mapzcoord-1][mapxcoord][mapycoord]<eDeath)){
                //                cout<<"edgePersType[mapzcoord-1][mapxcoord][mapycoord]: "<<mapzcoord-1<<", "<<mapxcoord<<", "<<mapycoord<<", "<<edgePersType[mapzcoord-1][mapxcoord][mapycoord]<<endl;
                //                cout<<"cellOrder[mapzcoord-1][mapxcoord][mapycoord]: "<<cellOrder[mapzcoord-1][mapxcoord][mapycoord]<<endl;
                //                cout<<"childVertex = cellOrder[mapzcoord-2][mapxcoord][mapycoord]: "<<cellOrder[mapzcoord-2][mapxcoord][mapycoord]<<endl;
                //
//                                cout<<"enter front 1"<<endl;
                childVertex = cellOrder[mapzcoord-2][mapxcoord][mapycoord];
                //                cout<<"enter front 2"<<endl;
                if (! vert_visited[childVertex]){ //check whether the previous vertex has already visited according to vert_visited[][].
                    //                    int childVertex = cellOrder[mapzcoord-2][mapxcoord][mapycoord];
                    myassert(parent_v[childVertex]<0);
                    parent_v[childVertex]=tmpv;
                    //                    cout<<"======================== 20.15 ===================="<<endl;
                    vert_stack.push(childVertex);
//                                        cout<<"======================== 20.16 ===================="<<endl;
                    //                    cout<<"---> right Happen, add vertex: "<<childVertex<<", with intensity: "<<phi->data[(* vList)[childVertex].xidx][(* vList)[childVertex].yidx]<<endl;
                }}
            //back vertex.
            //            cout<<"enter back 0"<<endl;
            //            if ((mapzcoord<mapNSlice-1)&&(edgePersType[mapzcoord+1][mapxcoord][mapycoord]==DESTROYER)&&(cellOrder[mapzcoord+1][mapxcoord][mapycoord]<eDeath)){
            if ((mapzcoord+2 < mapNSlice)&&(edgePersType[mapzcoord+1][mapxcoord][mapycoord]==DESTROYER)&&(cellOrder[mapzcoord+1][mapxcoord][mapycoord]<eDeath)){
                //                cout<<"mapNSlice-1: "<<mapNSlice-1<<endl;   //mapNSlice-1=190---> z=0~190, mapzcoord+1<<", "<<mapxcoord<<", "<<mapycoord = 191, 66, 222.
                //                cout<<"edgePersType[mapzcoord+1][mapxcoord][mapycoord]: "<<mapzcoord+1<<", "<<mapxcoord<<", "<<mapycoord<<", "<<edgePersType[mapzcoord+1][mapxcoord][mapycoord]<<endl;
                //                cout<<"cellOrder[mapzcoord+1][mapxcoord][mapycoord]: "<<cellOrder[mapzcoord+1][mapxcoord][mapycoord]<<endl;
                //                cout<<"childVertex = cellOrder[mapzcoord+2][mapxcoord][mapycoord]: "<<cellOrder[mapzcoord+2][mapxcoord][mapycoord]<<endl;
                
//                                cout<<"enter back 1"<<endl;
                childVertex = cellOrder[mapzcoord+2][mapxcoord][mapycoord];
                //                cout<<"enter back 2"<<endl;
                if (! vert_visited[childVertex]){ //check whether the previous vertex has already visited according to vert_visited[][].
                    //                    int childVertex = cellOrder[mapzcoord+2][mapxcoord][mapycoord];
                    myassert(parent_v[childVertex]<0);
                    parent_v[childVertex]=tmpv;
                    //                    cout<<"======================== 20.17 ===================="<<endl;
                    vert_stack.push(childVertex);
//                                        cout<<"======================== 20.18 ===================="<<endl;
                    //                    cout<<"---> right Happen, add vertex: "<<childVertex<<", with intensity: "<<phi->data[(* vList)[childVertex].xidx][(* vList)[childVertex].yidx]<<endl;
                }}
            //            if(xudongcounter1 == 30)break;
            
        };	//end of while(! vert_stack.empty())
        //        cout<<"======================== 22: finish finding the path from tmpv2. ===================="<<endl;
        
        
        myassert(! vert_stack.empty());		//vert_stack stores the second path.
        //        cout<<"1---> vert_stack is empty:  "<<vert_stack.empty()<<endl;
        int endV2 = -1;
        if(! vert_stack.empty()){
            endV2 = vert_stack.top();	//second end v, will connect to tmpv2
        }
        
        myassert(endV2!=endV1);
        myassert((endV2==vBirth)||(endV1==vBirth));
        // 		myassert((endV2==vHigh)||(endV1==vHigh));
        while(! vert_stack.empty()){
            vert_stack.pop();
        }
        tmpv = endV2;
        oldtmpv=tmpv;
        while(tmpv >= 0){
            //move down everything on the path
            myassert(! pathv_visited[tmpv]);
            pathv_visited[tmpv]=true;
            zcoord=(* vList)[tmpv].zidx;
            xcoord=(* vList)[tmpv].xidx;
            ycoord=(* vList)[tmpv].yidx;
            
            mapzcoord=(* vList)[tmpv].mapZCoord;
            mapxcoord=(* vList)[tmpv].mapXCoord;
            mapycoord=(* vList)[tmpv].mapYCoord;
            myassert(cellType[mapzcoord][mapxcoord][mapycoord]==VERTEX);
            myassert((cellFixType[mapzcoord][mapxcoord][mapycoord]==CF_UNDEFINED)||(cellFixType[mapzcoord][mapxcoord][mapycoord]==MOVEDOWN));
            
//            max1 = perturbM->data[tmpv2.zidx][tmpv2.xidx][tmpv2.yidx];;
//            max2 = perturbM->data[zcoord][xcoord][ycoord];
//            max_total = max(max1,max2);
//            min_total = min(max1,max2);
#ifdef drawinmerge_comp
            changeInteCore(0,perturbM->data, zcoord, xcoord, ycoord, changeToIntensity);    // draw from endV2 to tmpv2 (eDeath_vertex2).
            pointCounter++;
#endif
            ////            changeIntensity(n_pix,perturbM->data, zcoord, xcoord, ycoord, max1, max2);  //xudong: add this the change the intensity on 3 directions.
            temp_chords_2.push_back(path(zcoord, xcoord, ycoord, 0, 0)); //xudong add this line to store the 2nd half chord path.--------->
            ////            chords.push_back(path(zcoord, xcoord, ycoord, max_total)); //xudong add this line to store the chord path.--------------->>>>>
            

            
            ////			setVertFixType(mapxcoord,mapycoord,MOVEDOWN);   //xudong: delete this line, it might not help in merge-comp.
            
            oldtmpv=tmpv;
            tmpv=parent_v[tmpv];
        };	//end of while(! vert_stack.empty())
        myassert(oldtmpv==tmpv2);
        
        
        
        for(int i=0; i<temp_chords_1.size(); i++){
            chord.push_back(temp_chords_1[i]);
        }
        for(int i=temp_chords_2.size()-1; i>=0; i--){
            chord.push_back(temp_chords_2[i]);
        }
        //        cout<<"temp_chords_1: "<<temp_chords_1.size()<<", temp_chords_2.size(): "<<temp_chords_2.size()<<", chord.size(): "<<chord.size()<<endl;
        if(temp_chords_1.size()+temp_chords_2.size() != chord.size())
            cout<<"In merge_comp() size of temp_chords_1 + temp_chords_2 != size of chord"<<endl;
        
        return chord;
        
    }
    
    
 

    


    

};	// end of class CellMap

bool showedgeFrom2db(vector<int> vec, int target){
    for(int k=0; k<vec.size(); k++) {
        if(vec[k]==target){
//            cout<<"Find the target edge! "<<endl;
//            for(int d=0; d<vec.size(); d++) {
//                cout<<vec[d]<<", ";
//            }
            return true;
        }
    }
    return false;
}


// xudong: calcPers() is the function for calculating Persistent Homology .original:pair< int, int > calcPers
void calcPers(const int m,
              const int n,
              const int l,
              my3dMatrix * const img,
              const int n_pix,
//              const double rob_thd,
//              const double levelset_val,
//              double perturb_thd,
              my3dMatrix * const dijk3dimg,
              my3dMatrix * const perturbM,
              my3dMatrix * const critM,
              my3dMatrix * const dijM,
              my4dMatrix * const verifychord,
//              my3dMatrix * const verifychord_1,
//              my3dMatrix * const verifychord_2,
//              my3dMatrix * const verifychord_3,
//              my3dMatrix * const verifychord_4,
//              my3dMatrix * const verifychord_5,
//              my3dMatrix * const verifychord_6,
//              my3dMatrix * const verifychord_7,
//              my3dMatrix * const verifychord_8,
//              my3dMatrix * const verifychord_9,
//              my3dMatrix * const verifychord_10,
//       
              int ncomp_ub,
              const int chord_ub,
              const int chord_lb,
              double changeToIntensity,
              mVpM_2dMatrix * const phi2d_pM,
              const vector< vector < int > > pM_coord,
              mVpM_2dMatrix * const phi2d_mV,
              const vector< vector < int > > mV_coord,
              const int continue_reduction,
              const int death_1d,
              const double changeToIntensity_mV,
              const double changeToIntensity_pM){
    //double dist_from_path_thd, int remove_only, int kill_holes, int big_holes_remaining, double big_holes_thd){
    srand( (unsigned)time(NULL) );
	OUTPUT_MSG("Begin computing persistence");
    int k,i,j;
    //Creating the chordRange vector to show -> (z_min, z_max, x_min, x_max, y_min, y_max).
    //phi2d_pM->range and phi2d_mV->range are ranged follow x, y, z.
    vector<int> chordRange = vector<int>(6);
    chordRange[0] = max(min((phi2d_pM->range)[4],(phi2d_mV->range)[4]), 0);
    chordRange[1] = min(max((phi2d_pM->range)[5],(phi2d_mV->range)[5]), l);
    chordRange[2] = max(min((phi2d_pM->range)[0],(phi2d_mV->range)[0]), 0);
    chordRange[3] = min(max((phi2d_pM->range)[1],(phi2d_mV->range)[1]), m);
    chordRange[4] = max(min((phi2d_pM->range)[2],(phi2d_mV->range)[2]), 0);
    chordRange[5] = min(max((phi2d_pM->range)[3],(phi2d_mV->range)[3]), n);

    cout<<"The range of chordRange (z_min, z_max, x_min, x_max, y_min, y_max): "<<endl;
    for (i=0;i<chordRange.size();i++){
        cout<<chordRange[i]<<", ";
    }cout<<endl;
    
    

	//constructing vList
    vector<int> veQueueList;
	vector< Vertex > * vList=new vector< Vertex >;
	vector< Edge > * eList=new vector< Edge >;
	vector< Triangle > * trigList=new vector< Triangle >;
    vector< Cube > * cubeList = new vector <Cube>;  //xudong add this line.
	for (k=0;k<l;k++)   //xudong add the slice index k.
        for (i=0;i<m;i++)
            for (j=0;j<n;j++)
                vList->push_back(Vertex(k,i,j));    // just push_back Vertex with dimension as k by i by j.
    
	sort(vList->begin(), vList->end(), vCompVal);   // sort the vList according to the intensity of 3d points in small to large value (-1503 -->  1537).
    
	OUTPUT_MSG("--vList constructed and sorted");

	CellMap myCM(vList,eList,trigList,cubeList); // vList -> cellOrder -> eList -> trigList.

    cout<<"======================== 11: pass CellMap myCM() ===================="<<endl;    //--> Feb.10,2016 check the code till here. everything is good.
    
    
    

    
    
    
    
    
    
    
    
    //-***********************************************************************************************************************************************
#define dijkstra
//#define w1edge
#define w2edge
#define w3edge
    cout<<"======================== enter dijkstra===================="<<endl;
    typedef adjacency_list < listS, vecS, undirectedS,no_property, property < edge_weight_t, double > > graph_t;
    typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
    //    typedef std::pair<int, int> Edge;
    //    typedef pair<Vertex, Vertex> Edge3d;  //Wrong!!!, the BGL need the pair type has int type.
    typedef pair<int, int> Edge3d;
    
    int z,x,y,TotalAddNum=0;
    double w1=1, w2=1.414, w3=1.732;
//    double alf=1, beta=0.00001;
    double alf=0, beta=1;
    double intensityBias=0;
    double intensityBias_temp;
    
    for(k=0;k<l;k++){
        for (i=0;i<m;i++){
            for (j=0;j<n;j++){
                intensityBias_temp = dijk3dimg->data[k][i][j];
                if(intensityBias_temp < intensityBias){
                    intensityBias = intensityBias_temp;
                }  //find the min value of intensity among the 3d img.
            }
        }
    }
    if(intensityBias == 0){intensityBias= -1;}
    cout<<"intensityBias = "<<intensityBias<<endl;  // = -1077
    
    int inten_1,inten_2,inten_3,inten_4,inten_5,inten_6,inten_7,inten_8;
    //   vector<Edge3d> *edge_array3d = new vector<Edge3d>(l*m*n, Edge3d());
    //   vector<double> *weight3d = new vector<double>(l*m*n, 0.0);
    vector<Edge3d> *edge_array3d = new vector<Edge3d>;
    vector<double> *weight3d = new vector<double>;
    double w3_8_edge_ave_intensity=0;
    double divider = 1.0;
    
    for (z=0;z<l-1;z++){   //xudong add the slice index k.
        // if(k+1 >= l) continue;
        for (x=0;x<m-1;x++){
            // if(i+1 >= m) continue;
            for (y=0;y<n-1;y++){
                // if(j+1>dim_j-1) continue;
                
                // the following intensities need to be changed to positive. vertex id = m*n*z + n*x + y;
                inten_1 = (dijk3dimg->data[z][x][y] - intensityBias)/divider;     const int &v1 = 1000*1000*z +     1000*x + y;
                inten_2 = (dijk3dimg->data[z][x][y+1] - intensityBias)/divider;   const int &v2 = 1000*1000*z +     1000*x + (y+1);
                inten_3 = (dijk3dimg->data[z+1][x][y+1] - intensityBias)/divider; const int &v3 = 1000*1000*(z+1) + 1000*x + (y+1);
                inten_4 = (dijk3dimg->data[z+1][x][y] - intensityBias)/divider;   const int &v4 = 1000*1000*(z+1) + 1000*x + y;
                inten_5 = (dijk3dimg->data[z][x+1][y] - intensityBias)/divider;   const int &v5 = 1000*1000*z +     1000*(x+1) + y;
                inten_6 = (dijk3dimg->data[z][x+1][y+1] - intensityBias)/divider; const int &v6 = 1000*1000*z +     1000*(x+1) + (y+1);
                inten_7 = (dijk3dimg->data[z+1][x+1][y+1] - intensityBias)/divider;   const int &v7 = 1000*1000*(z+1) + 1000*(x+1) + (y+1);
                inten_8 = (dijk3dimg->data[z+1][x+1][y] - intensityBias)/divider; const int &v8 = 1000*1000*(z+1) + 1000*(x+1) + y;
                if(inten_1<0 || inten_2<0 || inten_3<0 || inten_4<0 || inten_5<0 || inten_6<0 || inten_7<0 || inten_8<0){
//                    cout<<"Negative intensity found !"<<endl;
                }
                //-******The vertex id for a single cube **********
                //                4-------3
                //              / |     / |
                //            1---|---2   |
                //            |   8---|---7
                //            | /    |  /
                //            5-------6
                //  inten_1 = phi3d->data[z][x][y]; const Vertex &v1 = Vertex(z,x,y);
                //  inten_2 = phi3d->data[z][x][y+1];   const Vertex &v2 = Vertex(z,x,y+1);
                //  inten_3 = phi3d->data[z+1][x][y+1]; const Vertex &v3 = Vertex(z+1,x,y+1);
                //  inten_4 = phi3d->data[z+1][x][y];   const Vertex &v4 = Vertex(z+1,x,y);
                //  inten_5 = phi3d->data[z][x+1][y];   const Vertex &v5 = Vertex(z,x+1,y);
                //  inten_6 = phi3d->data[z][x+1][y+1]; const Vertex &v6 = Vertex(z,x+1,y+1);
                //  inten_7 = phi3d->data[z+1][x+1][y+1];   const Vertex &v7 = Vertex(z+1,x+1,y+1);
                //  inten_8 = phi3d->data[z+1][x+1][y]; const Vertex &v8 = Vertex(z+1,x+1,y);
                //-************************************************
#ifdef w1edge
                // 12 w1 edges.
                edge_array3d->push_back(Edge3d(v1,v4)); weight3d->push_back(alf*w1 + beta*(inten_1 + inten_4)/2);  //back
                edge_array3d->push_back(Edge3d(v1,v5)); weight3d->push_back(alf*w1 + beta*(inten_1 + inten_5)/2);  //down
                edge_array3d->push_back(Edge3d(v1,v2)); weight3d->push_back(alf*w1 + beta*(inten_1 + inten_2)/2);  //right
                TotalAddNum += 3;
                if(z==l-2){
                    edge_array3d->push_back(Edge3d(v4,v3)); weight3d->push_back(alf*w1 + beta*(inten_4 + inten_3)/2);  //back-right
                    edge_array3d->push_back(Edge3d(v4,v8)); weight3d->push_back(alf*w1 + beta*(inten_4 + inten_8)/2);   //back-down
                    TotalAddNum += 2;
                }
                
                if(x==m-2){
                    edge_array3d->push_back(Edge3d(v5,v8)); weight3d->push_back(alf*w1 + beta*(inten_5 + inten_8)/2);//down-back
                    edge_array3d->push_back(Edge3d(v5,v6)); weight3d->push_back(alf*w1 + beta*(inten_5 + inten_6)/2);//down-right
                    TotalAddNum += 2;
                }
                
                if(y==n-2){
                    edge_array3d->push_back(Edge3d(v2,v3)); weight3d->push_back(alf*w1 + beta*(inten_2 + inten_3)/2);//right-back
                    edge_array3d->push_back(Edge3d(v2,v6)); weight3d->push_back(alf*w1 + beta*(inten_2 + inten_6)/2);//right-down
                    TotalAddNum += 2;
                }
                
                if(z==l-2 && y==n-2){TotalAddNum += 1;
                    edge_array3d->push_back(Edge3d(v7,v3)); weight3d->push_back(alf*w1 + beta*(inten_7 + inten_3)/2);}//up
                
                if(x==m-2 && y==n-2){TotalAddNum += 1;
                    edge_array3d->push_back(Edge3d(v7,v6)); weight3d->push_back(alf*w1 + beta*(inten_7 + inten_6)/2);}//front
                
                if(z==l-2 && x==m-2){TotalAddNum += 1;
                    edge_array3d->push_back(Edge3d(v7,v8)); weight3d->push_back(alf*w1 + beta*(inten_7 + inten_8)/2);}//left
#endif
                
#ifdef w2edge
                // 12 w2 edges.
                edge_array3d->push_back(Edge3d(v1,v6));
                weight3d->push_back(alf*w2 + beta*(inten_1 + inten_6)/2);
//                weight3d->push_back(alf*w2 + beta*(inten_1 + inten_6 + inten_2 + inten_5)/4);
                
                edge_array3d->push_back(Edge3d(v2,v5));
                weight3d->push_back(alf*w2 + beta*(inten_2 + inten_5)/2);
//                weight3d->push_back(alf*w2 + beta*(inten_2 + inten_5 + inten_1 + inten_6)/4);
                
                edge_array3d->push_back(Edge3d(v4,v5));
                weight3d->push_back(alf*w2 + beta*(inten_4 + inten_5)/2);
//                weight3d->push_back(alf*w2 + beta*(inten_4 + inten_5 + inten_1 + inten_8)/4);
                
                edge_array3d->push_back(Edge3d(v1,v8));
                weight3d->push_back(alf*w2 + beta*(inten_1 + inten_8)/2);
//                weight3d->push_back(alf*w2 + beta*(inten_1 + inten_8 + inten_4 + inten_5)/4);
                
                edge_array3d->push_back(Edge3d(v4,v2));
                weight3d->push_back(alf*w2 + beta*(inten_4 + inten_2)/2);
//                weight3d->push_back(alf*w2 + beta*(inten_4 + inten_2 + inten_1 + inten_3)/4);
                
                edge_array3d->push_back(Edge3d(v1,v3));
                weight3d->push_back(alf*w2 + beta*(inten_1 + inten_3)/2);
//                weight3d->push_back(alf*w2 + beta*(inten_1 + inten_3 + inten_4 + inten_2)/4);
                
                TotalAddNum += 6;
                
                if(z==l-2){
                    edge_array3d->push_back(Edge3d(v3,v8));
                    weight3d->push_back(alf*w2 + beta*(inten_3 + inten_8)/2);
//                    weight3d->push_back(alf*w2 + beta*(inten_3 + inten_8 + inten_4 + inten_7)/4);
                    
                    edge_array3d->push_back(Edge3d(v4,v7));
                    weight3d->push_back(alf*w2 + beta*(inten_4 + inten_7)/2);
//                    weight3d->push_back(alf*w2 + beta*(inten_4 + inten_7 + inten_3 + inten_8 )/4);
                    TotalAddNum += 2;
                }
                
                if(x==m-2){
                    edge_array3d->push_back(Edge3d(v5,v7));
                    weight3d->push_back(alf*w2 + beta*(inten_5 + inten_7)/2);
//                    weight3d->push_back(alf*w2 + beta*(inten_5 + inten_7 + inten_6 + inten_8)/4);
                    
                    edge_array3d->push_back(Edge3d(v6,v8));
                    weight3d->push_back(alf*w2 + beta*(inten_6 + inten_8)/2);
//                    weight3d->push_back(alf*w2 + beta*(inten_6 + inten_8 + inten_5 + inten_7)/4);
                    TotalAddNum += 2;
                }
                
                if(y==n-2){
                    edge_array3d->push_back(Edge3d(v2,v7));
                    weight3d->push_back(alf*w2 + beta*(inten_2 + inten_7)/2);
//                    weight3d->push_back(alf*w2 + beta*(inten_2 + inten_7 + inten_6 + inten_3)/4);
                    
                    edge_array3d->push_back(Edge3d(v6,v3));
                    weight3d->push_back(alf*w2 + beta*(inten_6 + inten_3)/2);
//                    weight3d->push_back(alf*w2 + beta*(inten_6 + inten_3 + inten_2 + inten_7)/4);
                    TotalAddNum += 2;
                }
#endif
#ifdef w3edge
                // 4 w3 edges.
                w3_8_edge_ave_intensity = (inten_1 + inten_2 + inten_3 + inten_4 + inten_5 + inten_6 + inten_7 + inten_8)/8;
                edge_array3d->push_back(Edge3d(v1,v7));
//                weight3d->push_back(alf*w3 + beta*w3_8_edge_ave_intensity);
                weight3d->push_back(alf*w3 + beta*(inten_1 + inten_7)/2);
                
                edge_array3d->push_back(Edge3d(v3,v5));
//                weight3d->push_back(alf*w3 + beta*w3_8_edge_ave_intensity);
                weight3d->push_back(alf*w3 + beta*(inten_3 + inten_5)/2);
                
                edge_array3d->push_back(Edge3d(v4,v6));
//                weight3d->push_back(alf*w3 + beta*w3_8_edge_ave_intensity);
                weight3d->push_back(alf*w3 + beta*(inten_4 + inten_6)/2);
                
                edge_array3d->push_back(Edge3d(v2,v8));
//                weight3d->push_back(alf*w3 + beta*w3_8_edge_ave_intensity);
                weight3d->push_back(alf*w3 + beta*(inten_2 + inten_8)/2);
                TotalAddNum += 4;
#endif
            }
        }
    }
    //        cout<<"TotalAddNum = "<<TotalAddNum<<endl;
    int edgeNum_temp = (m*(n-1)+n*(m-1))*l + m*n*(l-1);
    //        cout<<"calculate the edgeNum_temp: "<<edgeNum_temp<<endl;
    int edgeFaceNum_temp = ((m-1)*(n-1)*l + (l-1)*(m*(n-1)+n*(m-1)))*2;
    //        cout<<"calculate the edgeFaceNum_temp: "<<edgeFaceNum_temp<<endl;
    int edgeCubeNum_temp = (m-1)*(n-1)*(l-1)*4;
    //        cout<<"calculate the edgeCubeNum_temp: "<<edgeCubeNum_temp<<endl;
    if(TotalAddNum == (edgeNum_temp + edgeFaceNum_temp + edgeCubeNum_temp)) cout<<" Dijkstra: Find all edges for Graph. GOOD!!! "<<endl;
    else cout<<" Dijkstra: Lost edge!!! "<<endl;
    
    
    //--------------- destination Points for Dijkstra Algo ---------------
    // define two destination voxels.
    vector<int> mV_Des, sV;
    for (i=(phi2d_pM->coord).size()-2;i<(phi2d_pM->coord).size();i++){
        z = (phi2d_pM->coord)[i][2];
        x = (phi2d_pM->coord)[i][0];
        y = (phi2d_pM->coord)[i][1];
        //        cout<<"pM_Des-> "<<" z: "<<temp_z<<", x: "<<temp_x<<", y: "<<temp_y<<endl;
//        temp_des = 1000*1000*z + 1000*x + y;
        mV_Des.push_back(1000*1000*z + 1000*x + y);
        //        cout<<"pM_Des: "<<pM_Des[i]<<endl;
    }

//    for (i=0;i<(phi2d_pM->coord).size()-2;i++){
    for (i=3;i<6;i++){  // start from point_3, point_4, point_5.
        z = (phi2d_pM->coord)[i][2];
        x = (phi2d_pM->coord)[i][0];
        y = (phi2d_pM->coord)[i][1];
        //        cout<<"pM_Des-> "<<" z: "<<temp_z<<", x: "<<temp_x<<", y: "<<temp_y<<endl;
        //            temp_des = 1000*1000*temp_z + 1000*temp_x + temp_y;
        sV.push_back(1000*1000*z + 1000*x + y);
    }
    //--------------- destination Points for Dijkstra Algo ---------------


    int num_edge = edge_array3d->size();
    vector<double>::iterator weights_iterator = weight3d->begin();
    graph_t g(edge_array3d->begin(), edge_array3d->end(), weights_iterator, num_edge);
    vector< vector<path> > shortestPath(sV.size(),vector<path>());

    //    property_map<graph_t, edge_weight_t>::type weightmap = get(weight3d, g);
    for(j=0;j<mV_Des.size();j++){
        for(i=0;i<sV.size();i++){
            vector<vertex_descriptor> p(num_vertices(g));   // store the parent of *vi .
            vector<double> d(num_vertices(g)); // store the distance from sV to *vi.
            
            vertex_descriptor s = vertex(sV[i],g);
            //         dijkstra_shortest_paths(g, s,predecessor_map(boost::make_iterator_property_map(p.begin(),get(boost::vertex_index, g))).distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, g))));
            dijkstra_shortest_paths(g,s,predecessor_map(make_iterator_property_map(p.begin(),get(vertex_index, g))).distance_map(make_iterator_property_map(d.begin(), get(vertex_index, g))));
            
            for (int t = mV_Des[j]; t != sV[i]; t = p[t]){
                z = t/1000000;
                x = (t%1000000)/1000;
                y = t%1000;
                shortestPath[i].push_back(path(z,x,y,0,0));
            }
            z = sV[i]/1000000;
            x = (sV[i]%1000000)/1000;
            y = sV[i]%1000;
            shortestPath[i].push_back(path(z,x,y,0,0));
            
        }
    }
    
 
//            show2dchords(shortestPath);

    //        int temp_1, temp_2;
    //        Vertex temp_v1, temp_v2;

    for(i=0; i<shortestPath.size(); i++){
        for(j=0; j<shortestPath[i].size(); j++){
            
            changeInteCore(0,dijM->data, shortestPath[i][j].zcoord, shortestPath[i][j].xcoord, shortestPath[i][j].ycoord, changeToIntensity);
            
        }
    }

    
    delete edge_array3d;
    delete weight3d;
    
    cout<<"======================== pass dijkstra===================="<<endl;

    //-***********************************************************************************************************************************************
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
#ifdef rest
    int num_t_creator=0;  //for 3D case.
    int num_c_destroyer=0;    //for 3D case.
    int num_e_creator=0;
    int num_t_destroyer=0;
    int num_v_creator=0;// number of vertices creating non-essential class
    int num_e_destroyer=0;// number of edge destroyer (non-essential).
    
    //output edge-trig pairs whose persistence is bigger than pers_thd
    int vBirth,vDeath;
    int low;
    //    double tmp_pers,tmp_rob,tmp_double,tmp_death,tmp_birth;
    double tmp_rob,tmp_death,tmp_birth; //delete tmp_pers and tmp_double,
    //    int tmp_int;    // since it's useless, this tem_int was deleted.
    int tmptmp=0,tmptmptmp=0;;   // xudong: this line is useless.
    double eDeath_1D;//xudong add this parameters for setting the upper bound of triDeath.

//-----------------===
//Start building boundary_3D
//-----------------===


//    vector< Vertex >:: iterator iter_pM;
//    int cou = 0;
//    for(iter_pM = vList->begin(); iter_pM != vList->end();iter_pM++, cou++){
//
//        if(phi3d->data[(*iter_pM).zidx][(*iter_pM).xidx][(*iter_pM).yidx] <= changeToIntensity_pM){
//            cout<<"The order of vertex with intensity of changeToIntensity_pM is: "<<cou<<endl;
//            cout<<"vv->z: "<<(*iter_pM).zidx+1<<", vv->x: "<<(*iter_pM).xidx+1<<", vv->y: "<<(*iter_pM).yidx+1<<endl;
//        }
//    }
//    cout<<"------------------- The code above and below has the same output ---------------------------"<<endl;
    for(int iter_pM_1=0; iter_pM_1 < vList->size();iter_pM_1++){
        
        if(phi3d->data[(*vList)[iter_pM_1].zidx][(*vList)[iter_pM_1].xidx][(*vList)[iter_pM_1].yidx] <= changeToIntensity_pM){
//            cout<<"The order of vertex with intensity of changeToIntensity_pM is: "<<iter_pM_1<<endl;
//            cout<<"vv->z: "<<(*vList)[iter_pM_1].zidx+1<<", vv->x: "<<(*vList)[iter_pM_1].xidx+1<<", vv->y: "<<(*vList)[iter_pM_1].yidx+1<<endl;
        }
    }
    // xudong: need to complete the following code !!!!!!!!!!!
    
    
    //xudong: should obtain the edge order from cellmap, since I know the coordinate of those two tips.
    

    cout<<"======================== 12: Start building boundary_2D  ===================="<<endl;
//-----------------------------------====== Start building boundary_2D ======--------------------------------------
	//construct and reduce 2d boundary matrix
	vector< vector< int > > * boundary_2D=new vector< vector< int > >(myCM.trigNum, vector< int >());
	myCM.buildBoundary2D(boundary_2D);  // def of buildBoundary2D is at Line 577. push_back boundary to boudary_2D.
	vector< int > * low_2D_e2t=new vector< int >(myCM.edgeNum,-1);
    multiset< VertEdgePair > veQueue;
    multiset< VertEdgePair > veQueue_1;
	multiset< EdgeTrigPair > etQueue;
    int total_reduction_effect_2D = 0;

////    vector<reduction> records;
////    int iter_1 =0;  //xudong: xudong counting.
    ////------------------- Find the birthPoints -------------------
//    vector<int> birthPoint_1D = myCM.birthPoint1D(chordRange, phi2d_mV, 1, phi2d_pM);   // store the birthPoint in vector.
//    cout<<"The birthPoint_1D from the tips and mV are:"<<endl;
//    if (! birthPoint_1D.empty()){
//        for(int i=0;i<birthPoint_1D.size();i++){
//            cout<<"birthPoint_1D[i]: "<<birthPoint_1D[i]<<endl;
//            
//        }
//    }
    
    
    
    
    
	for (i=0;i<myCM.trigNum;i++){

		low = * (* boundary_2D)[i].rbegin();    // obtain the order of last edge for each triangle. boundary_2D was created from trigList, which contain 4 edge orders.
        
        
//
        while ( ( ! (* boundary_2D)[i].empty() ) && ( (* low_2D_e2t)[low]!=-1) ){

            (* boundary_2D)[i] = list_sym_diff((* boundary_2D)[i],(* boundary_2D)[(* low_2D_e2t)[low]]);  //case1: 282 336 355 356 358 359, case2: 312 327 358 359 361 362.

//            cout<<"---------------------------"<<endl<<endl;
            if(! (* boundary_2D)[i].empty()){
                low = * (* boundary_2D)[i].rbegin();  // change low value.
//                if( true == showedgeFrom2db((* boundary_2D)[i], 1370037)){cout<<"Enter while():  low= "<<low<<endl;}
            }
        }
        
        
//---============== separate the continuous reduction ==============
        
//        if (! (* boundary_2D)[i].empty()){
//            myassert(low>=0);
//            myassert((* low_2D_e2t)[low]==-1);
//            (* low_2D_e2t)[low]=i;
//        }
//    }
//    
//    
//    for (i=0;i<myCM.trigNum;i++){
//--=============== separate the continuous reduction =============
        
//        cout<<"======================== 12.2 ===================="<<endl;
        /************************************ continuous reduction *******************************************/
        //check the edges in (* boundary_2D)[i], if the edge has been the low before,
        //then combine the (* boundary_2D)[i] and (* boundary_2D)[(* low_2D_e2t)[low]], if
        //the number of edge can be decrease.
    
        int before_reduction_size = (* boundary_2D)[i].size();
        
        vector<int> tmp_bdry_v; //MatrixListType tmp_bdry_v;
        int new_pivot = low;
        vector<int>::iterator new_pivot_iter = lower_bound((* boundary_2D)[i].begin(),(* boundary_2D)[i].end(), new_pivot);
//        MatrixListType::iterator new_pivot_iter = lower_bound(boundary_upper[i].begin(),boundary_upper[i].end(), new_pivot);
        //   pos=lower_bound(myset.begin(),myset.end(),65);  //return  the element >= 65..
        assert(new_pivot_iter == (* boundary_2D)[i].end()-1);
        //        cout << "reducing column " << i << " size = " << boundary_upper[i].size() << endl;
//        cout<<"======================== 12.3 ===================="<<endl;
        while(new_pivot_iter != (* boundary_2D)[i].begin()){
            new_pivot_iter --;
            //            cout << new_pivot_iter - boundary_upper[i].begin() << " ";
            new_pivot = * new_pivot_iter;
//            cout<<"======================== 12.31 ===================="<<endl;
            if((* low_2D_e2t)[new_pivot] < BIG_INT){
                // the new_pivot is paired with some column before i, use it for reduction
                assert((* low_2D_e2t)[new_pivot] < i);
                
                if( continue_reduction == 1){
//                    (* boundary_2D)[i] = list_sym_diff((* boundary_2D)[i],(* boundary_2D)[(* low_2D_e2t)[low]]);  //->original code.
//                    cout<<"======================== 12.32 ===================="<<endl;
                    tmp_bdry_v = list_sym_diff((* boundary_2D)[i], (* boundary_2D)[(* low_2D_e2t)[new_pivot]]);
//                    cout<<"======================== 12.33 ===================="<<endl;
                    if(tmp_bdry_v.size() < (* boundary_2D)[i].size()){
                        (* boundary_2D)[i] =  tmp_bdry_v;                     //boundary_upper[i] = tmp_bdry_v;
                        // update the reduction list as well
//                        reduction_list[i]=list_sym_diff(reduction_list[i], reduction_list[(* low_2D_e2t)[new_pivot]]);
                    }
                }else if (continue_reduction == 2){
                    (* boundary_2D)[i] = list_sym_diff((* boundary_2D)[i], (* boundary_2D)[(* low_2D_e2t)[new_pivot]]);
//                    reduction_list[i]=list_sym_diff(reduction_list[i], reduction_list[(* low_2D_e2t)[new_pivot]]);
                }else{
                    assert(false);
                }
//                cout<<"======================== 12.34 ===================="<<endl;
            }
            
            new_pivot_iter = lower_bound((* boundary_2D)[i].begin(),(* boundary_2D)[i].end(), new_pivot);
        }
//        cout<<"======================== 12.4 ===================="<<endl;
        
        int after_reduction_size = (* boundary_2D)[i].size();
//        //        cout << "reduction effect, size = " << before_reduction_size << " --> " << after_reduction_size << endl;
//        //        cout << after_reduction_size - before_reduction_size << " ";
        total_reduction_effect_2D += after_reduction_size - before_reduction_size;
        
        /*************************************** continuous reduction ****************************************/
//        cout<<"======================== 12.5 ===================="<<endl;
        

        
       
        // <----the code above can merge two or more cycles to one cycle.---->
        //find one edge-triangle pair, which means find one more num_e_creator and one more num_t_destroyer.
        // But, need to decide whether the edge-triangle pair we found satisfy the conditions to be a component.
        if (! (* boundary_2D)[i].empty()){
            myassert(low>=0);
            myassert((* low_2D_e2t)[low]==-1);
            (* low_2D_e2t)[low]=i;

//            num_e_creator++;
//            num_t_destroyer++;
//            //record pair
            Edge edgeCreator=(* eList)[low];
            Triangle trigDestroyer=(* trigList)[i];
            vBirth= edgeCreator.v2_order;   // the intensity on the max vertex of low edge .
            vDeath= trigDestroyer.v4_order; // the intensity on the max vertex of this .
			tmp_death=phi3d->data[(* vList)[vDeath].zidx][(* vList)[vDeath].xidx][(* vList)[vDeath].yidx];    // get the intensity.
			tmp_birth=phi3d->data[(* vList)[vBirth].zidx][(* vList)[vBirth].xidx][(* vList)[vBirth].yidx];    // get the intensity.
//			tmp_rob=min(fabs(tmp_birth-levelset_val),fabs(tmp_death-levelset_val)); // calculate the robustness.
            int loopLen = (* boundary_2D)[i].size();
//            cout<<"In boundary-2D, et-pair: "<<i<<": vBirth-> "<<tmp_birth<<", vDeath-> "<<tmp_death<<endl;
//            if( (tmp_birth<levelset_val) && (tmp_death>levelset_val) && (tmp_rob>rob_thd) ){  //don't work for Xudong case.
//            if( ( tmp_birth < tmp_death ) && ( tmp_rob > rob_thd ) ){
//            if((-200<tmp_birth) && (tmp_birth < 200) && (tmp_birth < tmp_death) && (tmp_death > 0) && (tmp_death - tmp_birth > 60) && (chord_ub>=loopLen) && (loopLen>=chord_lb)){    //// 6/1/2018 -> this is good.
            
            int v1_index=edgeCreator.v1_order;
            int v2_index=edgeCreator.v2_order;
            Vertex v1 = (*vList)[v1_index];
            Vertex v2 = (*vList)[v2_index];
            //if(v2.zidx+1==28 && v2.xidx+1==25 && v2.yidx+1==23){
            
            // shouldn't have (chord_ub>=loopLen), since the chord found may be the multiple true chord connected together.
            if((-300<tmp_birth) && (tmp_birth < 300) && (tmp_birth < tmp_death) && (tmp_death > 0) && (tmp_death - tmp_birth > 60)  && (loopLen>=chord_lb)){
                etQueue.insert(EdgeTrigPair(low,i,tmp_rob, tmp_birth, tmp_death));
//                veQueue.insert(VertEdgePair(4,low,tmp_rob, tmp_birth, tmp_death));   // birth vbidx: 3902. 2429-> one tip
//                veQueue.insert(VertEdgePair(170,low,tmp_rob, tmp_birth, tmp_death));   // birth vbidx: 2572. 1001->one tip
//                veQueue_1.insert(VertEdgePair(1,low,tmp_rob, tmp_birth, tmp_death));   // birth vbidx: 3902. 2429-> one tip
//                veQueue_1.insert(VertEdgePair(199,low,tmp_rob, tmp_birth, tmp_death));   // birth vbidx: 2572. 1001->one tip
                
                
                // add multiple (at least two) points into the veQueue. In fact, we just need low (DeathTime).
//                if (! birthPoint_1D.empty()){
//                    for(int i=0;i<birthPoint_1D.size();i++){
//                        veQueue.insert(VertEdgePair(birthPoint_1D[i],low,tmp_rob, tmp_birth, tmp_death));
//                        veQueue_1.insert(VertEdgePair(birthPoint_1D[i],low,tmp_rob, tmp_birth, tmp_death));
//                        
//                    }
//                }

                veQueue.insert(VertEdgePair(0,low,tmp_rob, tmp_birth, tmp_death));
                veQueue_1.insert(VertEdgePair(0,low,tmp_rob, tmp_birth, tmp_death));
                
                /////************************************* Marker the critic points ******************************
#ifdef Draw_criPoint_On_perturbM
                    drawBallOnTips(2,perturbM->data, (* vList)[vBirth].zidx, (* vList)[vBirth].xidx, (* vList)[vBirth].yidx, changeToIntensity);
#endif
#ifdef Draw_criPoint_On_critM

                    drawBallOnTips(2,critM->data, (* vList)[vBirth].zidx, (* vList)[vBirth].xidx, (* vList)[vBirth].yidx, -800);
#endif
                /////************************************* Marker the critic points ******************************
                int v1_index=edgeCreator.v1_order;
                int v2_index=edgeCreator.v2_order;
                Vertex v1 = (*vList)[v1_index];
                Vertex v2 = (*vList)[v2_index];
//                cout<<"------------------------------------------------"<<endl;
//                cout<<"ee->v1z: "<<v1.zidx+1<<", ee->v1x: "<<v1.xidx+1<<", ee->v1y: "<<v1.yidx+1<<endl;
//                cout<<"ee->v2z: "<<v2.zidx+1<<", ee->v2x: "<<v2.xidx+1<<", ee->v2y: "<<v2.yidx+1<<endl;
//                cout<<"trig index i: "<<i<<endl;
                //ee->v2z: 28, ee->v2x: 25, ee->v2y: 23
   
                
            }
            //}
		}
        
		if (i % 100000 == 0)
		  OUTPUT_MSG( "reducing boundary 2D: i=" << i <<", trig number=" << myCM.trigNum );
	}
//    cout<<"======================== 12.6 ===================="<<endl;
    if(continue_reduction != 0){
        cout << "Total additional reduction on 2D (should less than 0) = " << total_reduction_effect_2D << endl;
    }

    myCM.setEPersType( low_2D_e2t );    //vector< int > * low_2D_e2t=new vector< int >(myCM.edgeNum,-1);

//	delete boundary_2D;

	OUTPUT_MSG( "boundary_2D all reduced" );
    cout<<"======================== 13: Finish building boundary_2D ===================="<<endl;

    
    
    
//-----------------------------------====== Start building boundary_1D and reduction ======--------------------------------------
    cout<<"======================== 14: Start building boundary_1D ===================="<<endl;
    
    //construct and reduce 1d boundary matrix
    vector< int > * low_1D_v2e= new vector< int >(myCM.vertNum,-1);
    // for each creator vertex, store the edge it is paired to
    
    vector< vector< int > > * boundary_1D=new vector< vector< int > >(myCM.edgeNum, vector< int >()); //first index is col index, each col init to empty
    myCM.buildBoundary1D(boundary_1D);
//    vector<int> tmp_bdry_1D;    int threshold_reduction = 140;
    
    int total_reduction_effect_1D = 0;
//    multiset< VertEdgePair > veQueue;	// robustness pairs
    
    for (i=0;i<myCM.edgeNum;i++){
        
        if ( (* low_2D_e2t)[i] >= 0 ){
            (*boundary_1D)[i].clear();
//            tmptmp++;
            continue;
        }else{
//            tmptmptmp++;
            myassert((* low_2D_e2t)[i] == -1);
            myassert((*boundary_1D)[i].size()==2);
        };
        
        //reduce column i
        low = * (* boundary_1D)[i].rbegin();
        
        while ( ( ! (* boundary_1D)[i].empty() ) && ( (* low_1D_v2e)[low]!=-1 ) ){
            (* boundary_1D)[i]=list_sym_diff((* boundary_1D)[i],(* boundary_1D)[(* low_1D_v2e)[low]]);
//==******************************************** Xudong: Useless ************************************************
//            tmp_bdry_1D = list_sym_diff((* boundary_1D)[i],(* boundary_1D)[(* low_1D_v2e)[low]]);
//            //                    cout<<"======================== 12.33 ===================="<<endl;
//            if(tmp_bdry_1D.size() < threshold_reduction){
//                (* boundary_1D)[i] =  tmp_bdry_1D;                     //boundary_upper[i] = tmp_bdry_v;
//            }
//==******************************************** Xudong: Useless *************************************************
            
            if(! (* boundary_1D)[i].empty()){
                low = * (* boundary_1D)[i].rbegin();
            }
        }
        
        
        
//--=================== separate the continuous reduction =====================
        
//        if (! (* boundary_1D)[i].empty()){
//            myassert(low>=0);
//            myassert((* low_1D_v2e)[low]==-1);
//            (* low_1D_v2e)[low]=i;
//        }
//        
//    }
//        for (i=0;i<myCM.edgeNum;i++){
//--=================== separate the continuous reduction =====================
            
            
            
            
        /************************************ continuous reduction *******************************************/
        
        int before_reduction_size = (* boundary_1D)[i].size();
        
        vector<int> tmp_bdry_v; //MatrixListType tmp_bdry_v;
        int new_pivot = low;
        vector<int>::iterator new_pivot_iter = lower_bound((* boundary_1D)[i].begin(),(* boundary_1D)[i].end(), new_pivot);
        //        MatrixListType::iterator new_pivot_iter = lower_bound(boundary_upper[i].begin(),boundary_upper[i].end(), new_pivot);
        //   pos=lower_bound(myset.begin(),myset.end(),65);  //return  the element >= 65..
        assert(new_pivot_iter == (* boundary_1D)[i].end()-1);
        //        cout << "reducing column " << i << " size = " << boundary_upper[i].size() << endl;
        //        cout<<"======================== 12.3 ===================="<<endl;
        while(new_pivot_iter != (* boundary_1D)[i].begin()){
            new_pivot_iter --;
            //            cout << new_pivot_iter - boundary_upper[i].begin() << " ";
            new_pivot = * new_pivot_iter;
            //            cout<<"======================== 12.31 ===================="<<endl;
            if((* low_1D_v2e)[new_pivot] < BIG_INT){
                // the new_pivot is paired with some column before i, use it for reduction
                assert((* low_1D_v2e)[new_pivot] < i);
                
                if( continue_reduction == 1){
                    //                    (* boundary_2D)[i] = list_sym_diff((* boundary_2D)[i],(* boundary_2D)[(* low_2D_e2t)[low]]);  //->original code.
                    //                    cout<<"======================== 12.32 ===================="<<endl;
                    tmp_bdry_v = list_sym_diff((* boundary_1D)[i], (* boundary_1D)[(* low_1D_v2e)[new_pivot]]);
                    //                    cout<<"======================== 12.33 ===================="<<endl;
                    if(tmp_bdry_v.size() < (* boundary_1D)[i].size()){
                        (* boundary_1D)[i] =  tmp_bdry_v;                     //boundary_upper[i] = tmp_bdry_v;
                        // update the reduction list as well
                        //                        reduction_list[i]=list_sym_diff(reduction_list[i], reduction_list[(* low_2D_e2t)[new_pivot]]);
                    }
                }else if (continue_reduction == 2){
                    (* boundary_1D)[i] = list_sym_diff((* boundary_1D)[i], (* boundary_1D)[(* low_1D_v2e)[new_pivot]]);
                    //                    reduction_list[i]=list_sym_diff(reduction_list[i], reduction_list[(* low_2D_e2t)[new_pivot]]);
                }else{
                    assert(false);
                }
                //                cout<<"======================== 12.34 ===================="<<endl;
            }
            
            new_pivot_iter = lower_bound((* boundary_1D)[i].begin(),(* boundary_1D)[i].end(), new_pivot);
        }
        //        cout<<"======================== 12.4 ===================="<<endl;
        
        int after_reduction_size = (* boundary_1D)[i].size();
        //        //        cout << "reduction effect, size = " << before_reduction_size << " --> " << after_reduction_size << endl;
        //        //        cout << after_reduction_size - before_reduction_size << " ";
        total_reduction_effect_1D += after_reduction_size - before_reduction_size;
        
        /*************************************** continuous reduction ****************************************/
        
        
        
        if (! (* boundary_1D)[i].empty()){
            myassert(low>=0);
            myassert((* low_1D_v2e)[low]==-1);
            (* low_1D_v2e)[low]=i;
            num_e_destroyer++;
            num_v_creator++;
            
            myassert((* boundary_1D)[i].size()==2);
            
            int high =  * (* boundary_1D)[i].begin();
            //reduce high
            while (  (* low_1D_v2e)[high]!=-1 ){
                int edge_high=(*low_1D_v2e)[high];
                (* boundary_1D)[i]=list_sym_diff((* boundary_1D)[i],(* boundary_1D)[edge_high]);
                myassert((* boundary_1D)[i].size()==2);
                high = * (* boundary_1D)[i].begin();
            }
            
            
            //record pair
            vBirth= low;	//creator vertex
            Edge eDestroyer = (* eList)[i];
            vDeath=eDestroyer.v2_order;
            tmp_death=phi3d->data[(* vList)[vDeath].zidx][(* vList)[vDeath].xidx][(* vList)[vDeath].yidx];
            tmp_birth=phi3d->data[(* vList)[vBirth].zidx][(* vList)[vBirth].xidx][(* vList)[vBirth].yidx];
//            tmp_rob=min(fabs(tmp_birth-levelset_val),fabs(tmp_death-levelset_val));
//            if( (tmp_birth<levelset_val) && (tmp_death>levelset_val) && (tmp_rob>rob_thd) ){
            if( (tmp_birth<-780) && (tmp_death>death_1d)  ){
                veQueue.insert(VertEdgePair(low,i,tmp_rob, tmp_birth, tmp_death));
                veQueue_1.insert(VertEdgePair(low,i,tmp_rob, tmp_birth, tmp_death));
                //the component could be killed by either merge or remove
            };
            
        }else{
            myassert(false);
        }
        
        if (i % 100000 == 0)
            OUTPUT_MSG( "reducing boundary 1D: i=" << i <<", edge number=" << myCM.edgeNum );
    }
    myassert(num_v_creator==myCM.vertNum-1);
    myassert(num_e_destroyer==num_v_creator);
    
    if(continue_reduction != 0){
        cout << "Total additional reduction on 1D (should less than 0) = " << total_reduction_effect_1D << endl;
    }
    // 	myassert(num_e_creator+num_e_destroyer==myCM.edgeNum);
    // 	myassert(num_t_destroyer==myCM.trigNum);
    OUTPUT_MSG( "boundary 1D all reduced" );
    cout<<"======================== 15: finish building boundary_1D ===================="<<endl;

    
    
    

    
    //below: xudong show all elements in veQueue and etQueue.
//    int findchord = 0;
//    multiset< VertEdgePair >::iterator ii;
//    int counterve=0;
    //    cout<<"------------ Below is veQueue -----------------"<<endl;
//    for(ii=veQueue.begin(); ii!=veQueue.end(); ii++,counterve++){
//        //        show_ve(*ii);
//        cout<<"================ "<<counterve<<" ================="<<endl;
//        cout<<"vbidx: "<<(*ii).vbidx<<", edidx: "<<(*ii).edidx<<", robustness: "<<(*ii).robustness<<", birth: "<<(*ii).birth<<", death: "<<(*ii).death<<endl;
//        Vertex vv = (*vList)[(*ii).vbidx];
//        cout<<"vv->z: "<<vv.zidx+1<<", vv->x: "<<vv.xidx+1<<", vv->y: "<<vv.yidx+1<<endl;
//        Edge ee = (* eList)[(*ii).edidx];
//        int v1_index=ee.v1_order;
//        int v2_index=ee.v2_order;
//        Vertex v1 = (*vList)[v1_index];
//        Vertex v2 = (*vList)[v2_index];
//        cout<<"ee->v1z: "<<v1.zidx+1<<", ee->v1x: "<<v1.xidx+1<<", ee->v1y: "<<v1.yidx+1<<endl;
//        cout<<"ee->v2z: "<<v2.zidx+1<<", ee->v2x: "<<v2.xidx+1<<", ee->v2y: "<<v2.yidx+1<<endl;
//        
//        if((v1.yidx==34)&&(v2.yidx==33)){
//            findchord = 1;
//        }
//        
//    }
//    cout<<"Total elements in veQueue is: "<<counterve<<endl;
    
    
    
    
    
//    cout<<"------------ Below is etQueue ----------------->>>"<<endl;
//    for(itit=etQueue.begin(); itit!=etQueue.end(); itit++,counterve_2++){
//        show_et(*itit);
//    }
//    cout<<"Total elements in etQueue is: "<<counterve_2<<", "<<etQueue.size()<<endl;
//    cout<<"Total elements in etQueue is: "<<etQueue.size()<<endl;
  
//    multiset< TrigCubePair >::iterator itt;
//    cout<<"------------ Below is tcQueue ----------------->>>"<<endl;
//    for(itt=tcQueue.begin(); itt!=tcQueue.end(); itt++){
//        show_tc(*itt);
//    }
//    cout<<"Total elements in tcQueue is: "<<itt<<endl;
    
//    cout<<"--------------------------------------------- Above is etQueue ----------------------------------------->>>"<<endl;
//    cout<<"The edge order of veQueue element: ";show1dvector(veQueueList);
    //above: xudong show all elements in veQueue and etQueue.



//==================================== fixing topology as followings ========================================
	double dtime, btime;
	int eDeath, eBirth, tDeath;
    
    vector< vector< path > > SortedPath;
    int comp_counter = 0;
//    chords_temp1.push_back(vector< path >());   //xudong: create a storage space for the following stroing-> exist chords_temp1[0], but chords_temp1[0].size()=0.
//    bool chords_temp1_write = false;
//    vector< vector< path > > critchords2d;
//    vector<path> critchords1d;
//    bool critchords1d_write = false;
//    int Nchords_pM;
//    vector< vector < int > > neighborOrder(5, vector< int >(5,-1)); //grammar, can not use "* neighborOrder" when allocate size.
//    vector< vector < int > > * neighborOrder = new vector< vector < int > >(5, vector< int >(5,-1));    // but "* neighborOrder" is allowed when allocate size with "new".
//    int pathiter = 0;
//    int pathiter_2 = 0;
//    int Nvertex_eachComp = 0;
//    int Nvertex_eachComp_2 = 0;
    

#define _d0
    cout<<"======================== 17: starting merge-comp() ===================="<<endl;
    int comp_skipped = 1;
    int counter_xu = 0;

	for(multiset< VertEdgePair >::iterator myveiter1=veQueue.begin(); myveiter1!=veQueue.end(); myveiter1++, comp_counter++){
        
//		if( comp_skipped < ncomp_ub ){
//			comp_skipped ++;
//			continue;
//		};

		vBirth=myveiter1->vbidx;
		eDeath=myveiter1->edidx;
        btime=myveiter1->birth;
        dtime=myveiter1->death;
//#ifdef showinfo
//        cout<<"------------------------ "<<comp_counter<<" --------------------"<<endl;
//#endif
////        pathiter = chords_temp1.size()-1;
//        Nvertex_eachComp = myCM.merge_comp(vBirth, eDeath, btime, dtime, perturbM, critM, n_pix, changeToIntensity,chords_temp_1, chords_temp1_write, chord_ub,chord_lb,phi2d_pM->coord,phi2d_mV->coord,critchords1d,critchords1d_write, chordRange, tol_pMmarker);

//        vector<path> path_temp = myCM.merge_comp(vBirth,min(changeToIntensity_mV, changeToIntensity_pM), eDeath, perturbM, changeToIntensity, counter_xu);
        
        ////-> each chord is from endV1->tmpv1->tmpv2->endV2.
        vector<path> path_temp = myCM.merge_comp(1, min(changeToIntensity_mV, changeToIntensity_pM), eDeath, perturbM, changeToIntensity, counter_xu);
        cout<<comp_counter<<". size of each chord: "<<path_temp.size()<<endl;
        SortedPath.push_back(path_temp);
    }
    cout<<SortedPath.size()<<" Paths(chords) were found"<<endl;
    
    //output vector< vector< path > > SortedPath to critM->data.
    //    for(i=Id_chord; i<Id_chord+1; i++){
//    for(i=0; i<SortedPath.size(); i++){
//        for(j=0; j<SortedPath[i].size(); j++){
//            changeInteCore(0, critM->data, SortedPath[i][j].zcoord, SortedPath[i][j].xcoord, SortedPath[i][j].ycoord, changeToIntensity);
//        }
//    }


    
    
    
    
    
    
    
    
    
    int z_coord, x_coord, y_coord;
    int z_coord_next, x_coord_next, y_coord_next;
    int z_coord_temp, x_coord_temp, y_coord_temp;
    int row, col;

//    vector< vector< int > > SortedPath_id(SortedPath.size(), vector< int >());
    vector< vector< double > > voxel_aveinten(SortedPath.size(), vector< double >());
    for(i=0; i<SortedPath.size(); i++){
        for(j=0; j<SortedPath[i].size(); j++){
            voxel_aveinten[i].push_back(0);
        }
    }
    
    
    double tmpInten;
    int tmpInten_counter;
    double intensity_thre = -400;
    vector< path > temp_chord;
    vector< vector< path > > SortedPath_seg;
    vector< vector< path > > SortedPath_seg_fil;
    vector< vector< path > > SortedPath_seg_range;
    int Id_chord = 8;
    
    if(!SortedPath.empty()){
        cout<<"======================== 18: Checking whether each chord is consecutive ==============="<<endl;
        for(i=0; i<SortedPath.size(); i++){
            if(SortedPath[i].size()<=1) continue;
            for(j=0; j<SortedPath[i].size()-1; j++){
                z_coord = SortedPath[i][j].zcoord;
                x_coord = SortedPath[i][j].xcoord;
                y_coord = SortedPath[i][j].ycoord;
                z_coord_next = SortedPath[i][j+1].zcoord;
                x_coord_next = SortedPath[i][j+1].xcoord;
                y_coord_next = SortedPath[i][j+1].ycoord;
                
                if((abs(z_coord-z_coord_next)==1 && y_coord==y_coord_next && x_coord==x_coord_next)||
                   (z_coord==z_coord_next && abs(y_coord-y_coord_next)==1 && x_coord==x_coord_next)||
                   (z_coord==z_coord_next && y_coord==y_coord_next && abs(x_coord-x_coord_next)==1)){}
                else{cout<<"The "<<i<<"th chord isn't continue !!!!"<<endl;}
                
            }
            
        }
        

        cout<<"======= 18.1 Identifying the chord by drawing ball on it on Id_chord or all chords ====="<<endl;
        // marking each chord.
////        for(i=0; i<SortedPath.size(); i++){
//        for(i=Id_chord; i<Id_chord+1; i++){
//            cout<<"SortedPath[i].size(): "<<SortedPath[i].size()<<endl;
//            for(j=0; j<SortedPath[i].size(); j++){
//                z_coord = SortedPath[i][j].zcoord;
//                x_coord = SortedPath[i][j].xcoord;
//                y_coord = SortedPath[i][j].ycoord;
//                if(0==j){
//                    drawBallOnTips(2,critM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                    //                drawBallOnTips(2,dijM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                }
//                if(40==j){
//                    drawBallOnTips(2,critM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                    //                drawBallOnTips(2,dijM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                }
//                if(80==j){
//                    drawBallOnTips(4,critM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                    //                drawBallOnTips(4,dijM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                }
//                if(120==j){
//                    drawBallOnTips(2,critM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                    //                drawBallOnTips(4,dijM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                }
//                if(460==j){
//                    drawBallOnTips(4,critM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                    //                drawBallOnTips(4,dijM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                }
//                if(500==j){
//                    drawBallOnTips(4,critM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                    //                drawBallOnTips(4,dijM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                }
//                if(SortedPath[i].size()-10 ==j){
//                    drawBallOnTips(4,critM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                    //                drawBallOnTips(4,dijM->data, z_coord, x_coord, y_coord, changeToIntensity);
//                }
//            }
//        }

        
        
        // Create 2d vector voxel_aveinten to store the inten of each voxel on each chord.
        


        cout<<"======================== 19: Calcu pheri intensity of each vocel and wirte into voxel_aveinten ========="<<endl;
        for(i=0; i<SortedPath.size(); i++){
            if(SortedPath[i].size()<=1) continue;
            for(j=0; j<SortedPath[i].size(); j++){
                tmpInten=0;
                tmpInten_counter = 0;
                z_coord = SortedPath[i][j].zcoord;
                x_coord = SortedPath[i][j].xcoord;
                y_coord = SortedPath[i][j].ycoord;
//                cout<<"======================== 23.2 ===================="<<endl;
                if(0 <= z_coord-1){
                    tmpInten += img->data[z_coord-1][x_coord][y_coord];
                    tmpInten_counter++;}
//                 cout<<"======================== 23.3 ===================="<<endl;
                if(z_coord+1 <= l-1){
                    tmpInten += img->data[z_coord+1][x_coord][y_coord];
                    tmpInten_counter++;}
//                cout<<"======================== 23.4 ===================="<<endl;
                if(0 <= x_coord-1){
                    tmpInten += img->data[z_coord][x_coord-1][y_coord];
                    tmpInten_counter++;}
//                cout<<"======================== 23.5 ===================="<<endl;
                if(x_coord+1 <= m-1){
//                    cout<<"== 1 =="<<endl;
//                    cout<<"l: "<<m<<", z_coord+1: "<<z_coord+1<<endl;
//                    cout<<"m: "<<m<<", x_coord+1: "<<x_coord+1<<endl;
//                    cout<<"n: "<<m<<", y_coord+1: "<<y_coord+1<<endl;
//                    cout<<"tmpInten: "<<tmpInten<<endl;
//                    cout<<"intensity: "<<phi3d->data[z_coord][x_coord+1][y_coord]<<endl;
                    tmpInten += img->data[z_coord][x_coord+1][y_coord];
//                    cout<<"== 2 =="<<endl;
                    tmpInten_counter++;
//                    cout<<"== 3 =="<<endl;
                }
//                cout<<"======================== 23.6 ===================="<<endl;
                if(0 <= y_coord-1){
//                    cout<<"== 4 =="<<endl;
                    tmpInten += img->data[z_coord][x_coord][y_coord-1];
//                    cout<<"== 5 =="<<endl;
                    tmpInten_counter++;
//                    cout<<"== 6 =="<<endl;
                }
//                cout<<"======================== 23.7 ===================="<<endl;
                if(y_coord+1 <= n-1){
                    tmpInten += img->data[z_coord][x_coord][y_coord+1];
                    tmpInten_counter++;}
                voxel_aveinten[i][j] = tmpInten/tmpInten_counter;
//                cout<<"======================== 23.8 ===================="<<endl;
  
            }
            
        }
        
//        cout<<"19.1 show voxel_aveinten ===================="<<endl;
        //    //        for(i=0; i<SortedPath.size(); i++){
        //            for(i=Id_chord; i<Id_chord+1; i++){
        //                cout<<"------ "<<i<<" ----"<<endl;
        //                for(j=0; j<SortedPath[i].size(); j++){
        //                    cout<<j<<": "<<voxel_aveinten[i][j]<<endl;
        //                }
        //            }
        
        cout<<"20: Segment each chord ========== Id_chord or all chords =============="<<endl;
        for(i=0; i<voxel_aveinten.size(); i++){
//        for(i=Id_chord; i<Id_chord+1; i++){
            if(voxel_aveinten[i].size()<=1) continue;
            for(j=0; j<voxel_aveinten[i].size(); j++){
                if(voxel_aveinten[i][j] > intensity_thre){
                    temp_chord.push_back(SortedPath[i][j]);
//                    cout<<"temp_chord.size(): "<<temp_chord.size()<<endl;
                }
                else if(!temp_chord.empty()){
//                    cout<<"Push chord happen!!!!"<<endl;
//                    if(temp_chord.size() >= ???){ // set the range of chord length.
                        SortedPath_seg.push_back(temp_chord);
//                    }
                    temp_chord.clear();
//                    if(temp_chord.empty()){cout<<"temp_chord is empty"<<endl;}
                }
                if((j == voxel_aveinten[i].size()-1)&&(!temp_chord.empty())){
//                    cout<<"Push chord happen!!!!"<<endl;
//                    if(temp_chord.size() >= ???){ // set the range of chord length.
                    SortedPath_seg.push_back(temp_chord);
                    temp_chord.clear();
//                    }
                }
                
            }
        }
        cout<<"The side of SortedPath_seg: "<<SortedPath_seg.size()<<endl;
        
    }
    
    

    
    cout<<"21: Filtering by upper and lower bounds of chord length ===================="<<endl;
//    cout<<"Before, The side of SortedPath_seg_fil: "<<SortedPath_seg_fil.size()<<endl;
    // filter the SortedPath_seg with length upper and lower bounds.
    for(i=0; i<SortedPath_seg.size(); i++){
        int chord_length = SortedPath_seg[i].size();
//        cout<<"chord_length: "<<chord_length<<endl;
        if((chord_lb <= chord_length) && (chord_length <= chord_ub)){
//            cout<<"push_back-> chord_lb: "<<chord_lb<<", chord_ub: "<<chord_ub<<", chord_length: "<<chord_length<<endl;
            SortedPath_seg_fil.push_back(SortedPath_seg[i]);
            cout<<"size of sub chord: "<<SortedPath_seg[i].size()<<endl;
        }
    }
    cout<<"After filtering, The side of SortedPath_seg_fil: "<<SortedPath_seg_fil.size()<<endl;
//    show2dchords(SortedPath_seg_fil);
    
    
    
    
    
    
    
    
   
    
    ////---------- Jun 13th, 2018 need to do: 1. inside_mVpM_box() ---------------------------------
#define _SortedPath_seg_fil_range
    cout<<"22: Filtering by inside_mVpM_box() ===================="<<endl;
    int z_bias=0, x_bias=0, y_bias=10;
    SortedPath_seg_range = inside_mVpM_box(z_bias, x_bias, y_bias,SortedPath_seg_fil, chordRange);
    cout<<"After checking the range, The size of SortedPath_seg_range: "<<SortedPath_seg_range.size()<<endl;
//    show2dchords(SortedPath_seg_range);
    
    
    ////------------------------ Draw the SortedPath_seg on dijM -------------------------------------
    my3dMatrix * temp_3d = new my3dMatrix(m,n,l);
    for(i=0; i<SortedPath_seg_range.size(); i++){
//    for(i=0; i<1; i++){
//        vector< vector < vector < double > > > temp_3d;
        for(j=0; j<SortedPath_seg_range[i].size(); j++){
            changeInteCore(0, dijM->data, SortedPath_seg_range[i][j].zcoord, SortedPath_seg_range[i][j].xcoord, SortedPath_seg_range[i][j].ycoord, changeToIntensity);
            changeInteCore(0, temp_3d->data, SortedPath_seg_range[i][j].zcoord, SortedPath_seg_range[i][j].xcoord, SortedPath_seg_range[i][j].ycoord, changeToIntensity);
        }
        
        verifychord->data4d.push_back(temp_3d->data);
        temp_3d->clearTozero();
        
    }
    delete temp_3d;
//    verifychord->data4d.push_back(dijM->data);
//    cout<<"========== Check verifychord->data4d[0] is the same as dijM->data =========="<<endl;
//    datacomp(verifychord->data4d[0],dijM->data);
    cout<<"verifychord->data4d.size(): "<<verifychord->data4d.size()<<endl;

    ////---------- Jun 13th, 2018: 2. do not need to use vertexRangeChecker() ----------------------
    //  since use merge_comp() to find the chords with btime=changeToInten=-900, all the chords
    //  found should be terminated at position with intensity=-900 (touch the muscle or marker).
    //  Hence, do not need to check the distance of each chord voxel and the voxel of pM and mV.
    ////--------------------------------------------------------------------------------------------
    

    
    
     ////---------- Jun 13th, 2018 need to do: 3. make the chord straighter -------------------------
    
   
    

///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    int bias_intensity = 0;
    
     ////-> draw the chords from veQueue_1 to cirtM.
    for(multiset< VertEdgePair >::iterator myveiter2=veQueue_1.begin(); myveiter2!=veQueue_1.end(); myveiter2++,bias_intensity++){
        vBirth=myveiter2->vbidx;
        eDeath=myveiter2->edidx;
        btime=myveiter2->birth;
        dtime=myveiter2->death;
        int vDeath_1 = (* eList)[eDeath].v1_order;
        int vDeath_2 = (* eList)[eDeath].v2_order;

        if(bias_intensity==Id_chord){

//#ifdef showinfo
//            cout<<"------------------------ "<<bias_intensity<<" --------------------"<<endl;
//#endif
//            drawBallOnTips(2,critM->data, (* vList)[vBirth].zidx, (* vList)[vBirth].xidx, (* vList)[vBirth].yidx, changeToIntensity);
//            drawBallOnTips(4,critM->data, (* vList)[vDeath_1].zidx, (* vList)[vDeath_1].xidx, (* vList)[vDeath_1].yidx, changeToIntensity);
//            myCM.merge_comp(vBirth,min(changeToIntensity_mV, changeToIntensity_pM), eDeath, critM, changeToIntensity, counter_xu);    // previous merge_comp().
//            myCM.merge_comp(1, min(changeToIntensity_mV, changeToIntensity_pM), eDeath, critM, changeToIntensity, counter_xu);
        }
        if(bias_intensity==Id_chord){
//
//#ifdef showinfo
//            cout<<"------------------------ "<<bias_intensity<<" --------------------"<<endl;
//#endif
//            drawBallOnTips(2,dijM->data, (* vList)[vBirth].zidx, (* vList)[vBirth].xidx, (* vList)[vBirth].yidx, changeToIntensity);
//            drawBallOnTips(4,dijM->data, (* vList)[vDeath_1].zidx, (* vList)[vDeath_1].xidx, (* vList)[vDeath_1].yidx, changeToIntensity);
//            myCM.merge_comp(1,min(changeToIntensity_mV, changeToIntensity_pM), eDeath, dijM, changeToIntensity, counter_xu);
        }
    }
    
cout<<"======================== 23: write each chord into TXT for SVM ===================="<<endl<<endl;


//    ofstream outf3("SortedChordsForSVM.txt");
//    outf3<<m<<" "<<n<<" "<<l<<endl;
//    outf3<<endl;
//    
//    for(i=0; i<SortedPath_seg_range.size(); i++){    //----> need to change back.
//        //        for(i=0; i<1; i++){
//        for(j=0; j<SortedPath_seg_range[i].size(); j++){
//            outf3<<SortedPath_seg_range[i][j].xcoord<<" "<<SortedPath_seg_range[i][j].ycoord<<" "<<SortedPath_seg_range[i][j].zcoord<<endl;
//            
//        }
//        outf3<<endl;
//    }
//    outf3<<flush;
//    outf3.close();
    
    
//    ofstream outf4("halfSortedChords2.txt");
//    for(i=0; i<halfSortedChords2.size(); i++){    //----> need to change back.
//        //        for(i=0; i<1; i++){
//        for(j=0; j<halfSortedChords2[i].size(); j++){
//            outf4<<halfSortedChords2[i][j].xcoord<<" "<<halfSortedChords2[i][j].ycoord<<" "<<halfSortedChords2[i][j].zcoord<<endl;
//        }
//        outf4<<endl;
//        //            cout<<"halfSortedChords1[i].size(): "<<halfSortedChords1[i].size()<<endl;
//    }
//    
//    outf4<<flush;
//    outf4.close();
    
    delete boundary_1D;
    delete low_1D_v2e;
    delete boundary_2D;
    delete low_2D_e2t;
#endif
	delete vList;
	delete eList;
	delete trigList;
    delete cubeList;
    

	

//    delete boundary_3D;
//    delete low_3D_t2c;

}
        


        


//-----------------************-----------------------------------------------------------------------------------------------------------//
// xudong: the input 3D matrix must be double type, if it's int type, then convert the matrix to be double type.~~------------------------//
//-----------------************-----------------------------------------------------------------------------------------------------------//
const int numInputArgs  = 15;    //const int numInputArgs  = 4;
const int numOutputArgs = 7;
void mexFunction (int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
    cout<<"======================== The followings are from mexFunction: ===================="<<endl;
//->	fstream filestr;
//->	filestr.open (LOG_FILE, fstream::in | fstream::out | fstream::trunc);
//->	filestr << "Debuging Pers_2D" << endl;
//->	filestr.close();
	//mexWarnMsgTxt(LOG_FILE);
    if (nrhs != numInputArgs)
        cout<<"Incorrect number of input arguments"<<endl;
    if (nlhs != numOutputArgs)
        cout<<"Incorrect number of output arguments"<<endl;

//    double rob_thd = 0;
    // rob_thd=40,changeToIntensity_mV = -500;changeToIntensity_pM = -700;Len_chords = 50;chords_temp1_write, 40);chords_temp1_write,150); without changing intensity on original img can find on chords connecting to pM and mV.
//    double levelset_val = 0 ;   // levelset_val: the image value of the levelset (0 in image segmentation), xudong: refer to phi(x)<0 ????
//    double perturb_thd = 0;   // pers_thd: threshold of persistence (only bigger persistence would be recorded
//    int remove_only = 0;
    
    int Nlayer = 3; // if =0, cannot find chords from merge-comp(). xudong: set it = 3.

    const mwSize * img_dim = mxGetDimensions(prhs[0]);  // mxGetDimensions return parameter type as const mwSize *
    int m = img_dim[0]; //# of rows from input 3D img.
    int n = img_dim[1]; //# of cols from input 3D img.
    int l = img_dim[2]; //# of slices from input 3D img.
    cout<<"Dim of original img: m = "<<m<<", n = "<<n<<", l = "<<l<<endl;    //m = 166, n = 0, l = 85
    
    //Xudong:preparing the output 3d matrix.
    mwSize ndim = 3;
    mwSize row_dim = m;
    mwSize col_dim = n;
    mwSize slice_dim = l;
    const mwSize output_3d_dim[3] = {row_dim,col_dim,slice_dim};

    phi3d = new my3dMatrix(m,n,l);  //my3dMatrix * phi3d; //xudong add this line at Line 456.
    phi3d->input1D(mxGetPr(prhs[0]));
    
  	int ncomp_ub = (int)getMatlabScalar( prhs[1] );
  	int chord_ub = (int)getMatlabScalar( prhs[2] );
    int n_pix = (int)getMatlabScalar( prhs[3] ); cout<<"n_pix: "<<n_pix<<endl;// the number of pix around the chords.
    
    double changeToIntensity_mV = (int)getMatlabScalar( prhs[13] );cout<<"changeToIntensity_mV: "<<changeToIntensity_mV<<endl; //= -800 is very good for finding one chord in vePairs.
    double changeToIntensity_pM = (int)getMatlabScalar( prhs[14] );cout<<"changeToIntensity_pM: "<<changeToIntensity_pM<<endl; //= -800 is very good for finding one chord in vePairs.
    
/*
    //--------------------- cover the img with 1500 (white) shell -----------------------------
//    //    changeIntensity(3,phi3d->data, 0, 50, 29, 1500);
//    //    changeIntensity(3,phi3d->data, 0, 50, 28, 1500);
//    //    changeIntensity(3,phi3d->data, 0, 50, 27, 1500);
//    //    changeIntensity(3,phi3d->data, 0, 50, 26, 1500);
//    //    changeIntensity(n_pix,perturbM->data, zcoord, xcoord, ycoord, -1500);
//    int pix = 3;
//    for(int index_l=0;index_l<l;index_l++){
//        for(int index_m=0;index_m<m;index_m++){
//            for(int index_n=0;index_n<n;index_n++){
//                if(index_m!=0 && index_m!=m-1){
//                    index_n = 0; changeInteCore(pix,phi3d->data, index_l, index_m, index_n, 1500);
//                    index_n = n-1; changeInteCore(pix,phi3d->data, index_l, index_m, index_n, 1500);
//                    
//                }
//                else{
//                    changeInteCore(pix,phi3d->data, index_l, index_m, index_n, 1500);
//                }
//            }
//        }
//        
//    }
//  
//    for(int index_m=0;index_m<m;index_m++){
//        for(int index_n=0;index_n<n;index_n++){
//            changeInteCore(0,phi3d->data, 0, index_m, index_n, 1500);
//            changeInteCore(0,phi3d->data, l-1, index_m, index_n, 1500);
//            
//        }
//    }
    //------------------------------------------------------------------------------------------------
    */
    cout<<"pass here 1"<<endl;
    
 
    //---------------------------- Blacking the mitral Valve ----------------------------------
    const mwSize * mitralValve_dim = mxGetDimensions(prhs[4]);  // mxGetDimensions return parameter type as const mwSize *
    int mV_m = mitralValve_dim[0]; //# of rows from input 3D img.
    int mV_n = mitralValve_dim[1]; //# of cols from input 3D img.
    int n_pix_mV = (int)getMatlabScalar( prhs[5] ); // the number of pix around the chords.
    //cout<<"mV_m = "<<mV_m<<", mV_n = "<<mV_n<<", n_pix_mV: "<<n_pix_mV<<endl;    //m = 166, n = 0, l = 85

    mVpM_2dMatrix * phi2d_mV = new mVpM_2dMatrix(mV_m,mV_n);  //xudong add this line.
    phi2d_mV->input1D(mxGetPr(prhs[4]));    //xudong need to write the argument mxGetPr(prhs[0]).

    
    for(int k=0;k<=Nlayer-1;k++){
        for (int i=0;i<(phi2d_mV->coord).size();i++){
            if (3!=(phi2d_mV->coord)[i].size()){ cout<<"Missing data in mitValve marking list"<<endl;continue;}
            changeInteCore(n_pix_mV-k,phi3d->data, (phi2d_mV->coord)[i][2], (phi2d_mV->coord)[i][0], (phi2d_mV->coord)[i][1], changeToIntensity_mV+100*(Nlayer-1-k));
        }
    }
    //changeIntensity(n_pix_mV,phi3d->data, phi2d_mV->coord, -1000, 0, 2);
    //------------------------------------------------------------------------------------------------
//    cout<<"pass here 2"<<endl;


    //------------- Find the consecutive mV coordinates before blacking the tips ------------------------
//    cout<<"before -> size of phi2d_mV: "<<phi2d_mV->coord.size()<<endl; //before -> size of phi2d_mV: 986.
    vector< vector < int > > mV_coord = reconstruct_mV(phi3d->data, phi2d_mV->coord, changeToIntensity_mV+100*(Nlayer-1), changeToIntensity_mV);
//    cout<<"after -> size of phi2d_mV: "<<phi2d_mV->coord.size()<<endl;  //after -> size of phi2d_mV: 15423.
//    cout<<"after -> size of mV_coord: "<<mV_coord.size()<<endl;  //->
    //---------------------------------------------------------------------------------------------------
    

    
    
//    cout<<"pass here 3"<<endl;
    my3dMatrix * tempimg3d = new my3dMatrix(m,n,l);
    //    tempimg3d->copy3d(phi3d);   // copy phi3d->data to tempimg3d->data.
    my3dMatrix * perturbM3d = new my3dMatrix(m,n,l);
    my3dMatrix * critM3d = new my3dMatrix(m,n,l);
    my3dMatrix * dijM3d = new my3dMatrix(m,n,l);
    my4dMatrix * verifychord = new my4dMatrix;
    my3dMatrix * dijk3dimg = new my3dMatrix(m,n,l);  //build a empty 3dimg for dij.
//    cout<<"pass here 6"<<endl;
    
    
    
    //---------------------------- Blacking the tips -------------------------------------------
    const mwSize * papillaryMuscleTips_dim = mxGetDimensions(prhs[6]);  // mxGetDimensions return parameter type as const mwSize *
    int pM_m = papillaryMuscleTips_dim[0]; //# of rows from input 3D img.
    int pM_n = papillaryMuscleTips_dim[1]; //# of cols from input 3D img.
    int n_pix_pM = (int)getMatlabScalar( prhs[7] ); // the number of pix around the chords.
    //cout<<"pM_m = "<<pM_m<<", pM_n = "<<pM_n<<", n_pix_pM: "<<n_pix_pM<<endl;    //m = 166, n = 0, l = 85
    
    mVpM_2dMatrix * phi2d_pM = new mVpM_2dMatrix(pM_m,pM_n);  //xudong add this line.
    phi2d_pM->input1D(mxGetPr(prhs[6]));    //xudong need to write the argument mxGetPr(prhs[0]).
    //cout<<"The position of pM: "<<endl;
//    show2dvector(phi2d_pM->coord);cout<<endl;

    if(n_pix_pM-Nlayer+1 >= 1){
        for(int k=0;k<=Nlayer-1;k++){
            for (int i=(phi2d_pM->coord).size()-2;i<(phi2d_pM->coord).size();i++){
                if (3!=(phi2d_pM->coord)[i].size()){ cout<<"Missing data in mitValve marking list"<<endl;continue;}

                drawBallOnTips(n_pix_pM-k,phi3d->data, (phi2d_pM->coord)[i][2], (phi2d_pM->coord)[i][0], (phi2d_pM->coord)[i][1], changeToIntensity_pM+100*(Nlayer-1-k));

            }
        }
    }
    
    // draw two tips in tempimg3d.
    for(int k=0;k<=Nlayer-1;k++){
        for (int i=(phi2d_pM->coord).size()-2;i<(phi2d_pM->coord).size();i++){
            if (3!=(phi2d_pM->coord)[i].size()){ cout<<"Missing data in mitValve marking list"<<endl;continue;}
            
            drawBallOnTips(n_pix_pM-k,tempimg3d->data, (phi2d_pM->coord)[i][2], (phi2d_pM->coord)[i][0], (phi2d_pM->coord)[i][1], changeToIntensity_pM+100*(Nlayer-1-k));
            
        }
    }
    
    //draw the 7 markers in critM3d.
    for(int k=0;k<=Nlayer-1;k++){
        for (int i=0;i<(phi2d_pM->coord).size()-2;i++){
            if (3!=(phi2d_pM->coord)[i].size()){ cout<<"Missing data in mitValve marking list"<<endl;continue;}
            
            drawBallOnTips(n_pix_pM-k,critM3d->data, (phi2d_pM->coord)[i][2], (phi2d_pM->coord)[i][0], (phi2d_pM->coord)[i][1], changeToIntensity_pM+100*(Nlayer-1-k));
            
        }
    }
    //------------------------------------------------------------------------------------------------

    //------------- Find the consecutive pM coordinates before blacking the tips ------------------------
//    cout<<"before -> size of phi2d_pM: "<<phi2d_pM->coord.size()<<endl; //-> 1
    vector< vector < int > > pM_coord = reconstruct_pM(phi2d_pM->coord, n_pix_pM,l,m,n);
//    cout<<"after -> size of phi2d_pM: "<<phi2d_pM->coord.size()<<endl;  //->
//    cout<<"after -> size of pM_coord: "<<pM_coord.size()<<endl;  //->
    //---------------------------------------------------------------------------------------------------
    
    //    cout<<"The range of pM: "<<endl;
    //    for (int i=0;i<(phi2d_pM->range).size();i++){
    //        cout<<(phi2d_pM->range)[i]<<", ";
    //    }
    //    cout<<endl;
    //    cout<<"The range of mV: "<<endl;
    //    for (int i=0;i<(phi2d_mV->range).size();i++){
    //        cout<<(phi2d_mV->range)[i]<<", ";
    //    }
    //    cout<<endl;

//    cout<<"pass here 5"<<endl;
    double changeToIntensity = (double)getMatlabScalar( prhs[8] );
    //    cout<<"changeToIntensity on chords: "<<changeToIntensity<<endl;
    int chord_lb = (int)getMatlabScalar( prhs[9] );
    int tol = (int)getMatlabScalar( prhs[10] );
    int continue_reduction = (int)getMatlabScalar( prhs[11] );
    double death_1d = (double)getMatlabScalar( prhs[12] );
    
    
    
    
    //------------------- draw mV on img for finding the ends of each chord. ------------------
    my3dMatrix * img = new my3dMatrix(m,n,l);
    img->input1D(mxGetPr(prhs[0]));
    for (int i=0;i<mV_coord.size();i++){
        if (3!=mV_coord[i].size()){ cout<<"Missing data in mitValve marking list"<<endl;continue;}
        changeInteCore(0,img->data, mV_coord[i][2], mV_coord[i][0], mV_coord[i][1], changeToIntensity_mV*2);
    }
    //----------------------------------------------------------------------------------------
    
    
    
    
    
#define _calcPers_Function
    // remove_olny = 1 means only remove(flood) and donot merge. xudong->13:50, original: pair< int, int > betti_numbers=
//    int tol_mVmarker = 30;
//    int tol_pMmarker = 10;
	calcPers(m,
             n,
             l, //xudong: l is # of slices.
             img,
             n_pix,
             dijk3dimg,
             perturbM3d,
             critM3d,
             dijM3d,
             verifychord,

             ncomp_ub,
             
             chord_ub,
             chord_lb,
             changeToIntensity, /*changeToIntensity = -300,*/
             phi2d_pM,
             pM_coord,
             
             phi2d_mV,
             mV_coord,
//             
//             tol_mVmarker,
//             tol_pMmarker,
             continue_reduction,
             death_1d,
             changeToIntensity_mV,
             changeToIntensity_pM);

 
    //-------- Copying recovered chords from perturbM3d, mV and tips to tempimg3d -------------
//    tempimg3d->copy3d_doNotCopyPointWithIntensity0(perturbM3d); // recover the missing chords from perturbM3d (copy the recovered chords from perturbM3d).
    
    //---------------------------------- draw mV on tempimg3d --------------------------------
//    for (int i=0;i<mV_coord.size();i++){
//            if (3!=mV_coord[i].size()){ cout<<"Missing data in mitValve marking list"<<endl;continue;}
//            changeInteCore(0,tempimg3d->data, mV_coord[i][2], mV_coord[i][0], mV_coord[i][1], changeToIntensity_mV);
//    }

    
    //---------------------------------- draw pM on tempimg3d ---------------------------------
//    for (int i=pM_coord.size()-2;i<pM_coord.size();i++){
//        if (3!=pM_coord[i].size()){ cout<<"Missing data in mitValve marking list"<<endl;continue;}
//        drawBallOnTips(0,tempimg3d->data, pM_coord[i][2], pM_coord[i][0], pM_coord[i][1], changeToIntensity_pM);
//    }
//    for (int i=0;i<pM_coord.size()-2;i++){
//        if (3!=pM_coord[i].size()){ cout<<"Missing data in mitValve marking list"<<endl;continue;}
//        drawBallOnTips(0,critM3d->data, pM_coord[i][2], pM_coord[i][0], pM_coord[i][1], changeToIntensity_pM);
//    }
//    cout<<"pM_coord.size(): "<<pM_coord.size()<<endl;
//    cout<<"pass here 7"<<endl;

    //----------------------------------- Return results to Matlab. ------------------------------------
    plhs[0] = mxCreateNumericArray(ndim, output_3d_dim, mxDOUBLE_CLASS, mxREAL);  //xudong 3D case.
    perturbM3d->output1D(mxGetPr(plhs[0]));   //xudong 3D case, give perturbM->data to plhs[0].

    plhs[1] = mxCreateDoubleScalar(n);    // the # of comp for 3D.----------------------------> need to modify.
    plhs[2] = mxCreateDoubleScalar(l);   // the # of holes for 3D.----------------------------> need to modify.
    
    plhs[3] = mxCreateNumericArray(ndim, output_3d_dim, mxDOUBLE_CLASS, mxREAL);
    phi3d->output1D(mxGetPr(plhs[3]));   //Checking the mV and pM on original 3dvolume img.
//    cout<<"pass here 8"<<endl;
    plhs[4] = mxCreateNumericArray(ndim, output_3d_dim, mxDOUBLE_CLASS, mxREAL);  //xudong 3D case.
    critM3d->output1D(mxGetPr(plhs[4]));   //xudong 3D case, give critM3d->data to plhs[0].
    
    plhs[5] = mxCreateNumericArray(ndim, output_3d_dim, mxDOUBLE_CLASS, mxREAL);
    dijM3d->output1D(mxGetPr(plhs[5]));
    
    plhs[6] = mxCreateNumericArray(ndim, output_3d_dim, mxDOUBLE_CLASS, mxREAL);
    tempimg3d->output1D(mxGetPr(plhs[6]));
//    cout<<"Good till now"<<endl;
//    cout<<"pass here 9"<<endl;
    
//    mwSize ndim_4d = 4;
//    mwSize stack_dim = verifychord->data4d.size();
//    const mwSize output_4d_dim[4] = {row_dim,col_dim,slice_dim,stack_dim};
//    plhs[7] = mxCreateNumericArray(ndim_4d, output_4d_dim, mxDOUBLE_CLASS, mxREAL);
//    verifychord->output1D(mxGetPr(plhs[7]));
//    cout<<"pass here 10"<<endl;


    delete phi3d;
    delete img;
    delete perturbM3d;
    delete critM3d;
    delete dijM3d;
    delete tempimg3d;
    delete verifychord;
    delete dijk3dimg;
    
    
    
    delete phi2d_mV;
    delete phi2d_pM;

	return ;
}


