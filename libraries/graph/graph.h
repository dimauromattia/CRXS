#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "TGraph.h"
#include "TGraphErrors.h"



namespace CRACS {
    
    
    
    //! Kinetic Energy to Rigidity
    /*!
     *  \param double Ekin. Kinetic Energy.
     *  \param int Z. Charge number.
     *  \param int A. Mass number.
     *  \return double. R.
     */
    static double Ekin_to_R(double Ekin, int Z, int A){
        double m = A*938; // FIXME for large A
        double E = Ekin + m;
        double p = sqrt(E*E-m*m);
        return p/fabs(1.*Z);
    };
    
    
    //! Rigidity to Kinetic Energy
    /*!
     *  \param double R. Rigidity.
     *  \param int Z. Charge number.
     *  \param int A. Mass number.
     *  \return double. Ekin.
     */
    static double R_to_Ekin(double R, int Z, int A){
        double m = A*938; // FIXME for large A
        double p = R*fabs(1.*Z);
        return sqrt(p*p+m*m)-m;
        
    };
    
    //! Differential Kinetic Energy by Rigidity
    /*!
     *  \param double Ekin. Kinetic Energy.
     *  \param int Z. Charge number.
     *  \param int A. Mass number.
     *  \return double. dEkin/dR.
     */
    static double dEkin_by_dR(double Ekin, int Z, int A){
        double m = A*938; // FIXME for large A
        double E = Ekin + m;
        double p = sqrt(E*E-m*m);
        return p*fabs(1.*Z)/E;
    };
    
    //! Differential Rigidity by Kinetic Energy
    /*!
     *  \param double R.
     *  \param int Z. Charge number.
     *  \param int A. Mass number.
     *  \return double. dR/dEkin.
     */
    static double dR_by_dEkin(double R, int Z, int A){
        double m = A*938; // FIXME for large A
        double p = R*fabs(1.*Z);
        double E = sqrt(p*p+m*m);
        return E/(p*fabs(1.*Z));
    };

    
    /*! \brief Graph class. Store x,y and errors.
     *
     */
    static double fGraphDummy=-1;
    
    class Graph {
    
        
    public:
        //! Constructor.
        /*!
         *  \param std::string filename. (Full path) Filename of the file to read. Calls ReadFile().
         */
        Graph( int reserve = 100 ){ Clear(reserve); };
        //Graph( TGraph* _graph    );
        
        void Clear( int reserve = 100 );
        
        void AddPoint (        double  x, double  y, double  xe=0,           double  ye=0           );
        void SetPoint ( int i, double  x, double  y, double  xe=0,           double  ye=0           );
        void GetPoint ( int i, double &x, double &y, double  &xe,            double  &ye            ) const;
        void GetPoint ( int i, double &x, double &y                                                 ) const;
        
        double Evalf  ( double x                                                                    ) const;
        
        
        void SetName  ( std::string name  ) { fName  = name;  };
        void SetTitle ( std::string title ) { fTitle = title; };
        
        Graph Copy() const;
        
        void  Scale         ( double sx, double sy );
        Graph GetScaledCopy ( double sx, double sy ) const;
        
        Graph GetTiltedCopy ( double percent ) const;
        
        int Size() const { return fN; };
        
        Graph operator+ ( const Graph& right);
        Graph operator- ( const Graph& right);
        Graph operator/ ( const Graph& right);
        
        Graph operator& ( const Graph& right);
        
        
        TGraph*         GetAsTGraph();
        TGraphErrors*   GetAsTGraphErrors();
        
        void WriteToFile(std::string filename, std::string header="");
        
        void ConvertSpectrumFromEkinPerNToRigidity( int Z, int A, double power );
        void ConvertSpectrumFromRigidityToEkinPerN( int Z, int A, double power );
        
    private:
        
        int                 fN;
        std::string         fName;
        std::string         fTitle;
        
        std::vector<double> fX;
        std::vector<double> fY;
        
        std::vector<double> fXe;
        std::vector<double> fYe;
        
    };
    
}

#endif // GRAPH_H
