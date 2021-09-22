#include "iostream"
#include "math.h"

#include <graph.h>
#include <fileTools.h>



//CRACS::Graph::Graph( TGraph* _graph )
//{
//    Clear();
//    for (int n=0; n<_graph->GetN(); n++) {
//        fN++;
//        double x, y;
//        _graph->GetPoint(n, x, y);
//        fX.push_back(x);
//        fY.push_back(y);
//        if (TGraphErrors * _graphError = dynamic_cast<TGraphErrors *>(_graph) ) {
//            double xe = _graphError->GetErrorX(n);
//            double ye = _graphError->GetErrorY(n);
//            fXe.push_back(xe);
//            fYe.push_back(ye);
//        }else{
//            fXe.push_back(0.);
//            fYe.push_back(0.);
//        }
//    }
//}

void CRACS::Graph::Clear( int reserve )
{
    fN = 0;
    
    fX. clear();
    fY. clear();
    fXe.clear();
    fYe.clear();
    
    fX. reserve( reserve );
    fY. reserve( reserve );
    fXe.reserve( reserve );
    fYe.reserve( reserve );
}

void CRACS::Graph::AddPoint( double x, double y, double xe, double ye )
{
    fN++;
    fX. push_back( x  );
    fY. push_back( y  );
    fXe.push_back( xe );
    fYe.push_back( ye );
}

void CRACS::Graph::GetPoint( int i, double &x, double &y, double &xe, double &ye ) const
{
    x  = fX. at(i);
    y  = fY. at(i);
    xe = fXe.at(i);
    ye = fYe.at(i);
}
void CRACS::Graph::GetPoint( int i, double &x, double &y ) const
{
    x  = fX. at(i);
    y  = fY. at(i);
}

void CRACS::Graph::SetPoint( int i, double x,  double y,  double xe,  double ye  )
{
    fX. at(i) = x;
    fY. at(i) = y;
    fXe.at(i) = xe;
    fYe.at(i) = ye;
}

CRACS::Graph CRACS::Graph::Copy() const
{
    Graph copy;
    
    double x, y, xe, ye;
    for ( int i = 0; i<fN; i++ ) {
        GetPoint        ( i, x, y, xe, ye );
        copy.AddPoint   (    x, y, xe, ye );
    }
    return copy;
}

void CRACS::Graph::Scale( double sx, double sy )
{
    double x, y, xe, ye;
    for ( int i = 0; i<fN; i++ ) {
        GetPoint   ( i, x,    y,    xe,    ye    );
        SetPoint   ( i, x*sx, y*sy, xe*sx, ye*sy );
    }
}

CRACS::Graph CRACS::Graph::GetScaledCopy( double sx, double sy ) const
{
    Graph copy = Copy();
    double x, y, xe, ye;
    for ( int i = 0; i<copy.fN; i++ ) {
        copy.GetPoint   ( i, x,    y,    xe,    ye    );
        copy.SetPoint   ( i, x*sx, y*sy, xe*sx, ye*sy );
    }
    return copy;
}

CRACS::Graph CRACS::Graph::GetTiltedCopy( double percent ) const
{
    Graph copy = Copy();
    double x, y, xe, ye;
    for ( int i = 0; i<copy.fN; i++ ) {
        double change = (2./copy.fN*i-1)*percent/100;
        copy.GetPoint   ( i, x,    y,    xe,    ye    );
        copy.SetPoint   ( i, x, y + ye*change, xe, ye );
    }
    return copy;
}

void CRACS::Graph::ConvertSpectrumFromEkinPerNToRigidity( int Z, int A, double power ){
    for (int n=0; n<fN; n++) {
        double EkinPerN,F,xe,ye;
        GetPoint(n, EkinPerN, F, xe, ye);
        double R = Ekin_to_R(EkinPerN*A, Z, A);
        double Fnew = F*pow(R/EkinPerN, power)*dEkin_by_dR(EkinPerN*A,Z,A)/A;
        SetPoint( n, R,  Fnew, xe*R/EkinPerN, ye*Fnew/F);
    }
}

void CRACS::Graph::ConvertSpectrumFromRigidityToEkinPerN( int Z, int A, double power ){
    for (int n=0; n<fN; n++) {
        double R,F, xe, ye;
        GetPoint(n, R, F, xe, ye);
        double EkinPerN = R_to_Ekin(R, Z, A)/A;
        double Fnew = F*pow(EkinPerN/R, power)*dR_by_dEkin(R, Z, A)*A;
        SetPoint( n, EkinPerN, Fnew, xe*EkinPerN/R, ye*Fnew/F);
    }
}

double CRACS::Graph::Evalf(double x) const{
    
    int upper = -1;
    for (int i=0; i<fN; i++) {
        if (fX.at(i)>x) {
            upper = i;
            break;
        }
    }
    
    if(upper<1){
        std::cout << "Warning: Graph interpolation out of range. Return estimate" << std::endl;
        double fl = fY.at(fY.size()-1);
        return fl;
    }
    
    double x0 = log(fX.at(upper-1));
    double x1 = log(fX.at(upper  ));
    
    
    double y0 = log(fY.at(upper-1));
    double y1 = log(fY.at(upper  ));
    
    //std::cout << x0 << "  " << x1 << "  " << y0 << "  " << y1 << "  " << std::endl;
    
    double ret = y0 + (y1-y0)/(x1-x0)*(log(x)-x0);
    //std::cout <<  "     " << ret << std::endl;
    
    return exp(ret);
    
}


CRACS::Graph CRACS::Graph::operator+(const CRACS::Graph& right){
    
    Graph sum;
    
    double xModelLast=0, xModelFirst=0;
    double xMb=0, xMa=0, yMb=0, yMa=0;
    
    right.GetPoint(right.Size()-1, xModelLast, yMb);
    right.GetPoint(0, xModelFirst, yMb);
    right.GetPoint(0, xModelFirst, yMb);
    
    int ModelPoint=0;
    
    for (int pD= 0; pD<fN; pD++) {
        
        double xD=0, yD=0;
        GetPoint(pD, xD, yD);
        if (xD<xModelFirst || xD>xModelLast) {
            continue;
        }
        
        while (ModelPoint<right.Size()-1) {
            if (xMa>xD) {
                break;
            }
            right.GetPoint(ModelPoint, xMb, yMb);
            ModelPoint++;
            right.GetPoint(ModelPoint, xMa, yMa);
            if (xMa>xD) {
                break;
            }
        }
        
        //double yM = yMb + (yMa-yMb)/(xMa-xMb)*(xD-xMb);  //<<<< linear interpolation
        //double yM = yMb + (yMa-yMb)/(log(xMa)-log(xMb))*(log(xD)-log(xMb));  //<<<< log-linear interpolation
        double yM = exp(   log(yMb) + (log(yMa)-log(yMb))/(log(xMa)-log(xMb))*(log(xD)-log(xMb))   );  //<<<< log-log interpolation
        
        
        double yN = yD +yM;
        sum.AddPoint(xD, yN);
    }
    
    return sum;
    
}

CRACS::Graph CRACS::Graph::operator-(const CRACS::Graph& right){
    
    Graph sum;
    
    double xModelLast=0, xModelFirst=0;
    double xMb=0, xMa=0, yMb=0, yMa=0;
    
    right.GetPoint(right.Size()-1, xModelLast, yMb);
    right.GetPoint(0, xModelFirst, yMb);
    right.GetPoint(0, xModelFirst, yMb);
    
    int ModelPoint=0;
    
    for (int pD= 0; pD<fN; pD++) {
        
        double xD=0, yD=0, xDe=0, yDe=0;
        GetPoint(pD, xD, yD, xDe, yDe);
        if (xD<xModelFirst || xD>xModelLast) {
            continue;
        }
        
        while (ModelPoint<right.Size()-1) {
            if (xMa>xD) {
                break;
            }
            right.GetPoint(ModelPoint, xMb, yMb);
            ModelPoint++;
            right.GetPoint(ModelPoint, xMa, yMa);
            if (xMa>xD) {
                break;
            }
        }
        
        //double yM = yMb + (yMa-yMb)/(xMa-xMb)*(xD-xMb);  //<<<< linear interpolation
        //double yM = yMb + (yMa-yMb)/(log(xMa)-log(xMb))*(log(xD)-log(xMb));  //<<<< log-linear interpolation
        double yM = exp(   log(yMb) + (log(yMa)-log(yMb))/(log(xMa)-log(xMb))*(log(xD)-log(xMb))   );  //<<<< log-log interpolation
        
        double yN = yD - yM;
        sum.AddPoint(xD, yN, xDe, yDe);
    }
    
    return sum;
    
}

CRACS::Graph CRACS::Graph::operator/(const CRACS::Graph& right){
    
    Graph ratio;
    
    double xModelLast=0, xModelFirst=0;
    double xMb=0, xMa=0, yMb=0, yMa=0;
    
    right.GetPoint(right.Size()-1, xModelLast, yMb);
    right.GetPoint(0, xModelFirst, yMb);
    right.GetPoint(0, xModelFirst, yMb);
    
    int ModelPoint=0;
    
    for (int pD= 0; pD<fN; pD++) {
        
        double xD=0, yD=0;
        GetPoint(pD, xD, yD);
        if (xD<xModelFirst || xD>xModelLast) {
            continue;
        }
        
        while (ModelPoint<right.Size()-1) {
            if (xMa>xD) {
                break;
            }
            right.GetPoint(ModelPoint, xMb, yMb);
            ModelPoint++;
            right.GetPoint(ModelPoint, xMa, yMa);
            if (xMa>xD) {
                break;
            }
        }
        
        //double yM = yMb + (yMa-yMb)/(xMa-xMb)*(xD-xMb);  //<<<< linear interpolation
        //double yM = yMb + (yMa-yMb)/(log(xMa)-log(xMb))*(log(xD)-log(xMb));  //<<<< log-linear interpolation
        double yM = exp(   log(yMb) + (log(yMa)-log(yMb))/(log(xMa)-log(xMb))*(log(xD)-log(xMb))   );  //<<<< log-log interpolation
        
        double yN = yD / yM;
        ratio.AddPoint(xD, yN);
    }
    
    return ratio;
    
}

CRACS::Graph CRACS::Graph::operator&(const CRACS::Graph& right){
    
    int n = fN + right.Size();
    
    Graph combined(n);
   
    double x, y, xe, ye;
    for ( int i = 0; i<fN; i++ ) {
        GetPoint        ( i, x, y, xe, ye );
        combined.AddPoint   (    x, y, xe, ye );
    }
    for ( int i = 0; i<right.Size(); i++ ) {
        right.GetPoint      ( i, x, y, xe, ye );
        combined.AddPoint   (    x, y, xe, ye );
    }
    
    return combined;
    
}


TGraph* CRACS::Graph::GetAsTGraph(){
    TGraph* g = new TGraph();
    if (fName.length()>0) {
        g->SetName(fName.c_str());
    }
    if (fTitle.length()>0) {
        g->SetTitle(fTitle.c_str());
    }
    for (int i = 0; i<fN; i++) {
        g->SetPoint(i, fX.at(i), fY.at(i));
    }
    return g;
}

TGraphErrors* CRACS::Graph::GetAsTGraphErrors(){
    TGraphErrors* g = new TGraphErrors();
    if (fName.length()>0) {
        g->SetName(fName.c_str());
    }
    if (fTitle.length()>0) {
        g->SetTitle(fTitle.c_str());
    }
    for (int i = 0; i<fN; i++) {
        g->SetPoint         ( i, fX.at(i),  fY.at(i)  );
        g->SetPointError    ( i, fXe.at(i), fYe.at(i) );
    }
    return g;
}


void CRACS::Graph::WriteToFile(std::string filename, std::string header){
    
    std::string toWrite = "";
    if (header!="") {
        toWrite += (header + "\n");
    }
    for (int i=0; i<fN; i++) {
        toWrite += QString("%1 ").arg(fX. at(i)).leftJustified(20).toStdString();
        toWrite += QString("%1 ").arg(fY. at(i)).leftJustified(20).toStdString();
        toWrite += QString("%1 ").arg(fXe.at(i)).leftJustified(20).toStdString();
        toWrite += QString("%1 ").arg(fYe.at(i)).leftJustified(20).toStdString()+"\n";
    }
    CRACS::FileTool::WriteStringToFile(toWrite, filename);
    
}














