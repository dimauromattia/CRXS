#include <math.h>

#include <iostream>
#include <cstdlib>
#include <definitions.h>
#include <fileTools.h>
#include <configHandler.h>
#include <labCrossSectionTabulator.h>

#include <QString>

int CRACS::CS_lab_tab::f_dLog;
int CRACS::CS_lab_tab::f_i_max;
int CRACS::CS_lab_tab::f_j_max;

CRACS::CS_lab_tab::CS_lab_tab(){
    f_dLog        = 30;
    f_i_max =  7*f_dLog;
    f_j_max =  8*f_dLog;
    ReadAll();
}

CRACS::CS_lab_tab* CRACS::CS_lab_tab::fInstance = NULL;

CRACS::CS_lab_tab* CRACS::CS_lab_tab::GetInstance(){
    if(!fInstance)
        fInstance = new CS_lab_tab();
    return fInstance;
}

void CRACS::CS_lab_tab::ReadAll(std::string path){
    
    funOut(CS_lab_tab::ReadAll)
    
    std::string dir = ConfigHandler::GetInstance()->SoftwarePath()+"/data/CS_lab/";
    if (path!="default"){
        dir = path+"/";
    }
    
    fCSfilesnames.clear();
    fCStype.clear();
    fCSfiles.clear();
    fCSparametrization.clear();
    
    
    fCSfilesnames.      push_back("dT_pp_pbar_LAB_diMauro");
    fCStype.            push_back(CS::PP_PBAR);
    fCSparametrization. push_back(CS::DI_MAURO);
    
    fCSfilesnames.      push_back("dT_Hep_pbar_LAB_diMauro");
    fCStype.            push_back(CS::HeP_PBAR);
    fCSparametrization. push_back(CS::DI_MAURO);
    
    fCSfilesnames.      push_back("dT_pHe_pbar_LAB_diMauro");
    fCStype.            push_back(CS::PHe_PBAR);
    fCSparametrization. push_back(CS::DI_MAURO);
    
    fCSfilesnames.      push_back("dT_HeHe_pbar_LAB_diMauro");
    fCStype.            push_back(CS::HeHe_PBAR);
    fCSparametrization. push_back(CS::DI_MAURO);
    
    
    
    fCSfilesnames.      push_back("dT_pp_pbar_LAB_diMauro12");
    fCStype.            push_back(CS::PP_PBAR);
    fCSparametrization. push_back(CS::DI_MAURO12);
    
    fCSfilesnames.      push_back("dT_Hep_pbar_LAB_diMauro12");
    fCStype.            push_back(CS::HeP_PBAR);
    fCSparametrization. push_back(CS::DI_MAURO12);
    
    fCSfilesnames.      push_back("dT_pHe_pbar_LAB_diMauro12");
    fCStype.            push_back(CS::PHe_PBAR);
    fCSparametrization. push_back(CS::DI_MAURO12);
    
    fCSfilesnames.      push_back("dT_HeHe_pbar_LAB_diMauro12");
    fCStype.            push_back(CS::HeHe_PBAR);
    fCSparametrization. push_back(CS::DI_MAURO12);
    
    
    
    fCSfilesnames.      push_back("dT_pp_pbar_LAB_TanNg");
    fCStype.            push_back(CS::PP_PBAR);
    fCSparametrization. push_back(CS::TAN_NG);
    
    fCSfilesnames.      push_back("dT_pHe_pbar_LAB_TanNg");
    fCStype.            push_back(CS::PHe_PBAR);
    fCSparametrization. push_back(CS::TAN_NG);
    
    fCSfilesnames.      push_back("dT_Hep_pbar_LAB_TanNg");
    fCStype.            push_back(CS::HeP_PBAR);
    fCSparametrization. push_back(CS::TAN_NG);
    
    
    
    fCSfilesnames.      push_back("dT_pp_pbar_LAB_DTUNUC_USIN");
    fCStype.            push_back(CS::PP_PBAR);
    fCSparametrization. push_back(CS::DTUNUC);
    
    fCSfilesnames.      push_back("dT_pHe_pbar_LAB_DTUNUC_USIN");
    fCStype.            push_back(CS::PHe_PBAR);
    fCSparametrization. push_back(CS::DTUNUC);
    
    fCSfilesnames.      push_back("dT_Hep_pbar_LAB_DTUNUC_USIN");
    fCStype.            push_back(CS::HeP_PBAR);
    fCSparametrization. push_back(CS::DTUNUC);
    
    fCSfilesnames.      push_back("dT_HeHe_pbar_LAB_DTUNUC_USIN");
    fCStype.            push_back(CS::HeHe_PBAR);
    fCSparametrization. push_back(CS::DTUNUC);
    
    
    
    
    fCSfilesnames.      push_back("dT_pp_pbar_LAB_Duperray");
    fCStype.            push_back(CS::PP_PBAR);
    fCSparametrization. push_back(CS::DUPERRAY);
    
    fCSfilesnames.      push_back("dT_pHe_pbar_LAB_Duperray");
    fCStype.            push_back(CS::PHe_PBAR);
    fCSparametrization. push_back(CS::DUPERRAY);
    
    fCSfilesnames.      push_back("dT_Hep_pbar_LAB_Duperray");
    fCStype.            push_back(CS::HeP_PBAR);
    fCSparametrization. push_back(CS::DUPERRAY);
    
    fCSfilesnames.      push_back("dT_HeHe_pbar_LAB_Duperray");
    fCStype.            push_back(CS::HeHe_PBAR);
    fCSparametrization. push_back(CS::DUPERRAY);
    
    
    
    
    fCSfilesnames.      push_back("dT_pp_pbar_LAB_Korsmeier");
    fCStype.            push_back(CS::PP_PBAR);
    fCSparametrization. push_back(CS::KORSMEIER);
    
    fCSfilesnames.      push_back("dT_pHe_pbar_LAB_Korsmeier");
    fCStype.            push_back(CS::PHe_PBAR);
    fCSparametrization. push_back(CS::KORSMEIER);
    
    fCSfilesnames.      push_back("dT_Hep_pbar_LAB_Korsmeier");
    fCStype.            push_back(CS::HeP_PBAR);
    fCSparametrization. push_back(CS::KORSMEIER);

    
    
    
    fCSfilesnames.      push_back("dT_pp_pbar_LAB_Winkler");
    fCStype.            push_back(CS::PP_PBAR);
    fCSparametrization. push_back(CS::WINKLER);
    
    fCSfilesnames.      push_back("dT_pHe_pbar_LAB_Winkler");
    fCStype.            push_back(CS::PHe_PBAR);
    fCSparametrization. push_back(CS::WINKLER);
    
    fCSfilesnames.      push_back("dT_Hep_pbar_LAB_Winkler");
    fCStype.            push_back(CS::HeP_PBAR);
    fCSparametrization. push_back(CS::WINKLER);
   
    fCSfilesnames.      push_back("dT_HeHe_pbar_LAB_Winkler");
    fCStype.            push_back(CS::HeHe_PBAR);
    fCSparametrization. push_back(CS::WINKLER);
    
    
    fCSfilesnames.      push_back("dT_pp_pbar_LAB_WinklerWithHypWitNbar");
    fCStype.            push_back(CS::PP_PBAR);
    fCSparametrization. push_back(CS::WINKLER_withHYPwithNBAR);
    
    fCSfilesnames.      push_back("dT_pHe_pbar_LAB_WinklerWithHypWitNbar");
    fCStype.            push_back(CS::PHe_PBAR);
    fCSparametrization. push_back(CS::WINKLER_withHYPwithNBAR);
    
    fCSfilesnames.      push_back("dT_Hep_pbar_LAB_WinklerWithHypWitNbar");
    fCStype.            push_back(CS::HeP_PBAR);
    fCSparametrization. push_back(CS::WINKLER_withHYPwithNBAR);
    
    fCSfilesnames.      push_back("dT_HeHe_pbar_LAB_WinklerWithHypWitNbar");
    fCStype.            push_back(CS::HeHe_PBAR);
    fCSparametrization. push_back(CS::WINKLER_withHYPwithNBAR);
    
    
    
    fCSfilesnames.      push_back("dT_pp_pbar_LAB_KR_PPFRAG");
    fCStype.            push_back(CS::PP_PBAR);
    fCSparametrization. push_back(CS::KACHELRIESS);
    
    fCSfilesnames.      push_back("dT_pHe_pbar_LAB_KR_PPFRAG");
    fCStype.            push_back(CS::PHe_PBAR);
    fCSparametrization. push_back(CS::KACHELRIESS);
    
    fCSfilesnames.      push_back("dT_Hep_pbar_LAB_KR_PPFRAG");
    fCStype.            push_back(CS::HeP_PBAR);
    fCSparametrization. push_back(CS::KACHELRIESS);
    
    fCSfilesnames.      push_back("dT_HeHe_pbar_LAB_KR_PPFRAG");
    fCStype.            push_back(CS::HeHe_PBAR);
    fCSparametrization. push_back(CS::KACHELRIESS);
    
    
    
    // DBAR
    fCSfilesnames.      push_back("dT_pp_Dbar_LAB");
    fCStype.            push_back(CS::PP_DBAR);
    fCSparametrization. push_back(CS::DI_MAURO12);
    
    fCSfilesnames.      push_back("dT_ppbar_Dbar_LAB");
    fCStype.            push_back(CS::PPbar_DBAR);
    fCSparametrization. push_back(CS::DEFAULT);
    
    // HEBAR
    
    fCSfilesnames.      push_back("dT_pp_Hebar_LAB");
    fCStype.            push_back(CS::PP_HEBAR);
    fCSparametrization. push_back(CS::DI_MAURO12);
    
    fCSfilesnames.      push_back("dT_ppbar_Hebar_LAB");
    fCStype.            push_back(CS::PPbar_HEBAR);
    fCSparametrization. push_back(CS::DEFAULT);
    
    
    for (int i = 0; i<fCSfilesnames.size(); i++) {
        FileTool* f =  new FileTool(dir+fCSfilesnames.at(i)+".txt");
        f->ExtractNumberTable(f_j_max+2, " ");
        fCSfiles.push_back(f);
        varOut(fCSparametrization.size())
        varOut(fCSfilesnames.size())
    }
    

}

double CRACS::CS_lab_tab::GetCrossSection(double Tn_proj, double Tn_prod, int type, int parametrization){
    FileTool* file = NULL;
    funOut(CS_lab_tab::GetCrossSection)
    varOut(fCSfilesnames.size())
    varOut(fCSfiles.size())
    varOut(fCStype.size())
    for (int i = 0; i<fCSfilesnames.size(); i++) {
        varOut(type)
        varOut(fCStype.at(i))
        varOut(parametrization)
        varOut(fCSparametrization.at(i))
        if (type==fCStype.at(i) && parametrization==fCSparametrization.at(i)) {
            file = fCSfiles.at(i);
            varOut(i)
            varOut(file->GetFileName())
            break;
        }
    }
    if (!file) {
        std::cout << warning_parametrization_not_available << std::endl;
        return 0;
    }
    double id = f_dLog*(   log10( Tn_proj) )+1;
    double jd = f_dLog*( 1+log10( Tn_prod) )+1;
    int i = id;
    int j = jd;
    varOut(file->NumberTable(i  , 0  ))
    varOut(file->NumberTable(i+1, 0  ))
    varOut(file->NumberTable(0  , j  ))
    varOut(file->NumberTable(0  , j+1))
    if (i<1 || i>f_i_max || j<1 || j>f_j_max) {
        return 0;
    }
    //Interpolation
    double v_ll =   file->NumberTable(i  , j  );
    double v_lu =   file->NumberTable(i  , j+1) ;
    double v_ul =   file->NumberTable(i+1, j  ) ;
    double v_uu =   file->NumberTable(i+1, j+1);
    varOut(v_ll)
    varOut(v_lu)
    varOut(v_ul)
    varOut(v_uu)
    if (v_ll>0) v_ll=log(v_ll); else v_ll=-10000;
    if (v_lu>0) v_lu=log(v_lu); else v_lu=-10000;
    if (v_ul>0) v_ul=log(v_ul); else v_ul=-10000;
    if (v_uu>0) v_uu=log(v_uu); else v_uu=-10000;
    double v_l = v_ll + (jd-j)*(v_lu-v_ll);
    double v_u = v_ul + (jd-j)*(v_uu-v_ul);
    double ret =  exp(  v_l +  (id-i)*(v_u-v_l)  );
    varOut(ret);
    return ret;
}


void CRACS::CS_lab_tab::WriteCS(  double (*lab_CS)(double, double), std::string filename, std::string dir, int dLog  ){
    funOut(CS_lab_tab::WriteCS)
    std::cout << std::endl;
    std::string table = "";
    QString line("");
    line += QString("%1 ").arg( "T_proj\\T_prod" ).leftJustified(20);
    for (int j = 0; j<=dLog*8; j++) {
        double Tn_prod = pow(10, 1.*j/dLog-1.);
        line += QString("%1 ").arg( Tn_prod ).leftJustified(20);
    }
    table += line.toStdString();
    int p_last = -1;
    for (int i = 0; i<=dLog*7; i++) {
        int p = i*i/(1.*dLog*7*dLog*7)*100;
        if (p>p_last){
            p_last = p;
            std::cout << "\r" << p << "% done!" <<std::flush;
        }
        varOut(i)
        
        table += "\n";
        double Tn_proj = pow(10, 1.*i/dLog);
        QString line("");
        line += QString("%1 ").arg( Tn_proj ).leftJustified(20);
        for (int j = 0; j<=dLog*8; j++) {
            double Tn_prod = pow(10, 1.*j/dLog-1.);
            double lll = lab_CS(Tn_proj, Tn_prod);
//            if (lll!=lll){
//                out(lll)
//                out(Tn_prod)
//                out(Tn_proj)
//            }
            line += QString("%1 ").arg( lab_CS(Tn_proj, Tn_prod) ).leftJustified(20);
        }
        table += line.toStdString();
    }
    //std::cout << table << std::endl;
    FileTool::WriteStringToFile( table, dir+"/"+filename );
    std::cout << std::endl;
};


void CRACS::CS_lab_tab::WriteCS_lowProd(  double (*lab_CS)(double, double), std::string filename, std::string dir, int dLog  ){
    funOut(CS_lab_tab::WriteCS)
    std::cout << std::endl;
    std::string table = "";
    QString line("");
    line += QString("%1 ").arg( "T_proj\\T_prod" ).leftJustified(20);
    for (int j = 0; j<=dLog*9; j++) {
        double Tn_prod = pow(10, 1.*j/dLog-2.);
        line += QString("%1 ").arg( Tn_prod ).leftJustified(20);
    }
    table += line.toStdString();
    int p_last = -1;
    for (int i = 0; i<=dLog*8; i++) {
        int p = i*i/(1.*dLog*7*dLog*7)*100;
        if (p>p_last){
            p_last = p;
            std::cout << "\r" << p << "% done!" <<std::flush;
        }
        varOut(i)
        
        table += "\n";
        double Tn_proj = pow(10, 1.*i/dLog-1.);
        QString line("");
        line += QString("%1 ").arg( Tn_proj ).leftJustified(20);
        for (int j = 0; j<=dLog*9; j++) {
            double Tn_prod = pow(10, 1.*j/dLog-2.);
            double lll = lab_CS(Tn_proj, Tn_prod);
            if (lll!=lll){
                out(lll)
                out(Tn_prod)
                out(Tn_proj)
            }
            line += QString("%1 ").arg( lab_CS(Tn_proj, Tn_prod) ).leftJustified(20);
        }
        table += line.toStdString();
    }
    //std::cout << table << std::endl;
    FileTool::WriteStringToFile( table, dir+"/"+filename );
    std::cout << std::endl;
    
};



double CRACS::CS_lab_tab::GetInterpolation( double Tn_proj, double Tn_prod ){

    double id = f_dLog*(   log10( Tn_proj) )+1;
    double jd = f_dLog*( 1+log10( Tn_prod) )+1;
    int i = id;
    int j = jd;
    varOut(f_file->NumberTable(i  , 0  ))
    varOut(f_file->NumberTable(i+1, 0  ))
    varOut(f_file->NumberTable(0  , j  ))
    varOut(f_file->NumberTable(0  , j+1))
    if (i<1 || i>f_i_max || j<1 || j>f_j_max) {
        return 0;
    }
    //Interpolation
    double v_ll =   f_file->NumberTable(i  , j  );
    double v_lu =   f_file->NumberTable(i  , j+1) ;
    double v_ul =   f_file->NumberTable(i+1, j  ) ;
    double v_uu =   f_file->NumberTable(i+1, j+1);
    varOut(v_ll)
    varOut(v_lu)
    varOut(v_ul)
    varOut(v_uu)
    if (v_ll>0) v_ll=log(v_ll); else v_ll=-10000;
    if (v_lu>0) v_lu=log(v_lu); else v_lu=-10000;
    if (v_ul>0) v_ul=log(v_ul); else v_ul=-10000;
    if (v_uu>0) v_uu=log(v_uu); else v_uu=-10000;
    double v_l = v_ll + (jd-j)*(v_lu-v_ll);
    double v_u = v_ul + (jd-j)*(v_uu-v_ul);
    double ret =  exp(  v_l +  (id-i)*(v_u-v_l)  );
    varOut(ret);
    return ret;
    
};

void CRACS::CS_lab_tab::ReadCS( std::string filename, std::string dir, int col, int skip_header ){
    out(dir+"/"+filename)
    f_file =    new FileTool(dir+"/"+filename);
    f_file -> ExtractNumberTable(col, " ", true, skip_header);
};

double CRACS::CS_lab_tab::GetInterpolation_Transform( double Tn_proj, double Tn_prod, int nTarg, int nProj, int nProd, double TminProj, double TmaxProj, double TminProd, double TmaxProd, bool changeLoop ){
    funOut(CS_lab_tab::GetInterpolation_Transform)
    double id = (log(Tn_prod/TminProd))/(log(TmaxProd/TminProd))*(nProd);
    double jd = (log(Tn_proj/TminProj))/(log(TmaxProj/TminProj))*(nProj);
    
    int i = id;
    int j = jd;
//    out(Tn_prod)
//    out(i)
    if (i<0 || i>=nProd || j<0 || j>=nProj) {
        return 0;
    }
    varOut(i)
    varOut(j)
    //Interpolation
//    if (changeLoop) {
//        int h = i;
//        i = j;
//        j = h;
//        h = nProj;
//        nProd = nProj;
//        nProj = nProd;
//        double hd = id;
//        id = jd;
//        jd = hd;
//    }
    double v_ll =   f_file->NumberTable(  i   *nProj + j  , nTarg  );
    double v_lu =   f_file->NumberTable(  i   *nProj + j+1, nTarg  );
    double v_ul =   f_file->NumberTable( (i+1)*nProj + j  , nTarg  );
    double v_uu =   f_file->NumberTable( (i+1)*nProj + j+1, nTarg  );
    if (changeLoop) {
        v_ll =   f_file->NumberTable(  i   +nProd *  j   , nTarg  );
        v_lu =   f_file->NumberTable(  i   +nProd * (j+1), nTarg  );
        v_ul =   f_file->NumberTable( (i+1)+nProd *  j   , nTarg  );
        v_uu =   f_file->NumberTable( (i+1)+nProd * (j+1), nTarg  );
        varOut(f_file->NumberTable( i   +nProd *  j, 0  )/Tn_proj)
        varOut(f_file->NumberTable( i   +nProd *  j, 1  )/Tn_prod)
        //out(f_file->NumberTable( i   *nProj + j, 2  ))
    }
    varOut(v_ll)
    varOut(v_lu)
    varOut(v_ul)
    varOut(v_uu)
    if (v_ll>0) v_ll=log(v_ll); else v_ll=-10000;
    if (v_lu>0) v_lu=log(v_lu); else v_lu=-10000;
    if (v_ul>0) v_ul=log(v_ul); else v_ul=-10000;
    if (v_uu>0) v_uu=log(v_uu); else v_uu=-10000;
    double v_l = v_ll + (jd-j)*(v_lu-v_ll);
    double v_u = v_ul + (jd-j)*(v_uu-v_ul);
    double ret =  exp(  v_l +  (id-i)*(v_u-v_l)  );
    return ret;
};



void CRACS::CS_lab_tab::WriteCS_Transform( int col, std::string filename, std::string dir, double factor ){
    funOut(CS_lab_tab::)
    std::cout << std::endl;
    std::string table = "";
    QString line("");
    line += QString("%1 ").arg( "T_proj\\T_prod" ).leftJustified(20);
    for (int j = 0; j<=f_j_max; j++) {
        double Tn_prod = pow(10, 1.*j/f_dLog-1.);
        line += QString("%1 ").arg( Tn_prod ).leftJustified(20);
    }
    table += line.toStdString();
    int p_last = -1;
    for (int i = 0; i<=f_i_max; i++) {
        int p = i*i/(1.*f_i_max*f_i_max)*100;
        if (p>p_last){
            p_last = p;
            std::cout << "\r" << p << "% done!" <<std::flush;
        }
        table += "\n";
        double Tn_proj = pow(10, 1.*i/f_dLog);
        QString line("");
        line += QString("%1 ").arg( Tn_proj ).leftJustified(20);
        for (int j = 0; j<=f_j_max; j++) {
            double Tn_prod = pow(10, 1.*j/f_dLog-1.);
            line += QString("%1 ").arg( factor*GetInterpolation_Transform(Tn_proj, Tn_prod, col)*1e-31 ).leftJustified(20);
        }
        table += line.toStdString();
    }
    //std::cout << table << std::endl;
    FileTool::WriteStringToFile( table, dir+"/"+filename );
    std::cout << std::endl;
};



void CRACS::CS_lab_tab::WriteCS_Transform( int col, std::string filename, std::string dir, int nProj, int nProd, double TminProj, double TmaxProj, double TminProd, double TmaxProd, double factor, bool changeLoop   ){
    funOut(CS_lab_tab::)
    std::string table = "";
    QString line("");
    line += QString("%1 ").arg( "T_proj\\T_prod" ).leftJustified(20);
    for (int j = 0; j<=f_j_max; j++) {
        double Tn_prod = pow(10, 1.*j/f_dLog-1.);
        line += QString("%1 ").arg( Tn_prod ).leftJustified(20);
    }
    table += line.toStdString();
    int p_last = -1;
    for (int i = 0; i<=f_i_max; i++) {
        int p = i*i/(1.*f_i_max*f_i_max)*100;
        if (p>p_last){
            p_last = p;
            std::cout << "\r" << p << "% done!" <<std::flush;
        }
        table += "\n";
        double Tn_proj = pow(10, 1.*i/f_dLog);
        QString line("");
        line += QString("%1 ").arg( Tn_proj ).leftJustified(20);
        for (int j = 0; j<=f_j_max; j++) {
            double Tn_prod = pow(10, 1.*j/f_dLog-1.);
            line += QString("%1 ").arg( GetInterpolation_Transform(Tn_proj, Tn_prod, col, nProj, nProd, TminProj, TmaxProj, TminProd, TmaxProd, changeLoop)*factor ).leftJustified(20);
        }
        table += line.toStdString();
    }
    //std::cout << table << std::endl;
    FileTool::WriteStringToFile( table, dir+"/"+filename );
};

void CRACS::CS_lab_tab::TransformAll(){
    
    std::string dir = ConfigHandler::GetInstance()->SoftwarePath()+"/data/CS_lab/KR/";
    
    ReadCS("KR_all.txt", dir, 6, 0);
    WriteCS_Transform (  2, "dT_pp_pbar_LAB_KR_PPFRAG.txt",     ".", 211, 241, 1e0, 1e7, 1e-1, 1e7, 1e-31*1.0/2.3, true );
    WriteCS_Transform (  3, "dT_Hep_pbar_LAB_KR_PPFRAG.txt",    ".", 211, 241, 1e0, 1e7, 1e-1, 1e7, 1e-31*1.0/2.3, true );
    WriteCS_Transform (  4, "dT_pHe_pbar_LAB_KR_PPFRAG.txt",    ".", 211, 241, 1e0, 1e7, 1e-1, 1e7, 1e-31*1.0/2.3, true );
    WriteCS_Transform (  5, "dT_HeHe_pbar_LAB_KR_PPFRAG.txt",   ".", 211, 241, 1e0, 1e7, 1e-1, 1e7, 1e-31*1.0/2.3, true );
    
    dir = ConfigHandler::GetInstance()->SoftwarePath()+"/data/CS_lab/USIN/";
    
    ReadCS("dSdEpbar_p+HTanNg_p+HeDTUNUC.dat", dir, 2, 60);
    WriteCS_Transform (  0, "dT_pp_pbar_LAB_DTUNUC_USIN.txt",    ".", 1./2  );
    WriteCS_Transform (  1, "dT_pHe_pbar_LAB_DTUNUC_USIN.txt",   ".", 1./2  );
    
    ReadCS("dSdEpbar_He+HHe_DTUNUC2001.dat", dir, 2, 59);
    WriteCS_Transform (  0, "dT_Hep_pbar_LAB_DTUNUC_USIN.txt",   ".", 1./2  );
    WriteCS_Transform (  1, "dT_HeHe_pbar_LAB_DTUNUC_USIN.txt",  ".", 1./2  );
    
    ReadCS("dSdEpbar_p+HHe_Duperray05.dat", dir, 2, 59);
    WriteCS_Transform (  0, "dT_pp_pbar_LAB_Duperray_USIN.txt",   ".", 1./2  );
    WriteCS_Transform (  1, "dT_pHe_pbar_LAB_Duperray_USIN.txt",  ".", 1./2  );
    
    ReadCS("dSdEpbar_He+HHe_Duperray05.dat", dir, 2, 59);
    WriteCS_Transform (  0, "dT_Hep_pbar_LAB_Duperray_USIN.txt",  ".", 1./2  );
    WriteCS_Transform (  1, "dT_HeHe_pbar_LAB_Duperray_USIN.txt", ".", 1./2  );

}

void CRACS::CS_lab_tab::WriteAll(){
    
    out("dT_pp_p_LAB_Anderson")
    WriteCS_lowProd(dT_pp_p_LAB_Anderson,   "dT_pp_p_LAB_Anderson_low.txt"  );
    WriteCS        (dT_pp_p_LAB_Anderson,   "dT_pp_p_LAB_Anderson.txt"  );
    
    out("dT_Hep_pbar_LAB_Winkler")
    WriteCS(dT_Hep_pbar_LAB_Winkler,     "dT_Hep_pbar_LAB_Winkler.txt"    );
    out("dT_Hep_pbar_LAB_WinklerWithHypWitNbar")
    WriteCS(dT_Hep_pbar_LAB_WinklerWitHypWitNbar,     "dT_Hep_pbar_LAB_WinklerWithHypWitNbar.txt"    );
    
    out("dT_pHe_pbar_LAB_Winkler")
    WriteCS(dT_pHe_pbar_LAB_Winkler,     "dT_pHe_pbar_LAB_Winkler.txt"    );
    out("dT_pHe_pbar_LAB_WinklerWithHypWitNbar")
    WriteCS(dT_pHe_pbar_LAB_WinklerWitHypWitNbar,     "dT_pHe_pbar_LAB_WinklerWithHypWitNbar.txt"    );
    
    out("dT_HeHe_pbar_LAB_Winkler")
    WriteCS(dT_HeHe_pbar_LAB_Winkler,     "dT_HeHe_pbar_LAB_Winkler.txt"    );
    out("dT_HeHe_pbar_LAB_WinklerWithHypWitNbar")
    WriteCS(dT_HeHe_pbar_LAB_WinklerWitHypWitNbar,     "dT_HeHe_pbar_LAB_WinklerWithHypWitNbar.txt"    );
    
    out("dT_pp_pbar_LAB_Winkler")
    WriteCS(dT_pp_pbar_LAB_Winkler,     "dT_pp_pbar_LAB_Winkler.txt"    );
    out("dT_pp_pbar_LAB_WinklerWithHypWitNbar")
    WriteCS(dT_pp_pbar_LAB_WinklerWitHypWitNbar,     "dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt"    );

    
    
    out("dT_pp_Dbar_LAB")
    WriteCS(dT_pp_Dbar_LAB,            "dT_pp_Dbar_LAB.txt"            );
    out("dT_ppbar_Dbar_LAB")
    WriteCS(dT_ppbar_Dbar_LAB,         "dT_ppbar_Dbar_LAB.txt"         );
    
    out("dT_pp_Hebar_LAB")
    WriteCS(dT_pp_Hebar_LAB,           "dT_pp_Hebar_LAB.txt"           );
    out("dT_ppbar_Hebar_LAB")
    WriteCS(dT_ppbar_Hebar_LAB,        "dT_ppbar_Hebar_LAB.txt"        );
    
    
    
    
    out("dT_pp_p_LAB_Anderson")
    WriteCS_lowProd(dT_pp_p_LAB_Anderson,   "dT_pp_p_LAB_Anderson_low.txt"  );
    WriteCS        (dT_pp_p_LAB_Anderson,   "dT_pp_p_LAB_Anderson.txt"  );
    
    
   
    
    out("dT_HeHe_pbar_LAB_Duperray")
    WriteCS(dT_HeHe_pbar_LAB_Duperray,   "dT_HeHe_pbar_LAB_Duperray.txt"  );
    out("dT_pHe_pbar_LAB_Duperray")
    WriteCS(dT_pHe_pbar_LAB_Duperray,   "dT_pHe_pbar_LAB_Duperray.txt"  );
    std::system("cp dT_pHe_pbar_LAB_Duperray.txt dT_Hep_pbar_LAB_Duperray.txt");
    
    
    out("dT_HeHe_pbar_LAB_diMauro")
    WriteCS(dT_HeHe_pbar_LAB_diMauro,    "dT_HeHe_pbar_LAB_diMauro.txt"   );
    
    out("dT_HeHe_pbar_LAB_diMauro12")
    WriteCS(dT_HeHe_pbar_LAB_diMauro12,  "dT_HeHe_pbar_LAB_diMauro12.txt" );
    
    
    
    
    out("dT_pp_pbar_LAB_Winkler")
    WriteCS(dT_pp_pbar_LAB_Winkler,     "dT_pp_pbar_LAB_Winkler.txt"    );
    out("dT_pp_pbar_LAB_WinklerWithHypWitNbar")
    WriteCS(dT_pp_pbar_LAB_WinklerWitHypWitNbar,     "dT_pp_pbar_LAB_WinklerWithHypWitNbar.txt"    );
    
    out("dT_Hep_pbar_LAB_diMauro12")
    WriteCS(dT_Hep_pbar_LAB_diMauro12,  "dT_Hep_pbar_LAB_diMauro12.txt" );
    out("dT_pHe_pbar_LAB_diMauro12")
    WriteCS(dT_pHe_pbar_LAB_diMauro12,  "dT_pHe_pbar_LAB_diMauro12.txt" );
    out("dT_pp_pbar_LAB_diMauro12")
    WriteCS(dT_pp_pbar_LAB_diMauro12,   "dT_pp_pbar_LAB_diMauro12.txt"  );
    

    
    out("dT_pp_pbar_LAB_Korsmeier")
    WriteCS(dT_pp_pbar_LAB_Korsmeier,     "dT_pp_pbar_LAB_Korsmeier.txt"    );
    out("dT_pHe_pbar_LAB_Korsmeier")
    WriteCS(dT_pHe_pbar_LAB_Korsmeier,    "dT_pHe_pbar_LAB_Korsmeier.txt"   );
    std::system("cp dT_pHe_pbar_LAB_Korsmeier.txt dT_Hep_pbar_LAB_Korsmeier.txt");
    
    
    out("dT_pp_pbar_LAB_diMauro")
    WriteCS(dT_pp_pbar_LAB_diMauro,     "dT_pp_pbar_LAB_diMauro.txt"    );
    out("dT_pHe_pbar_LAB_diMauro")
    WriteCS(dT_pHe_pbar_LAB_diMauro,    "dT_pHe_pbar_LAB_diMauro.txt"   );
    std::system("cp dT_pHe_pbar_LAB_diMauro.txt dT_Hep_pbar_LAB_diMauro.txt");
    
    
    
    
    out("dT_pp_pbar_LAB_Duperray")
    WriteCS(dT_pp_pbar_LAB_Duperray,    "dT_pp_pbar_LAB_Duperray.txt"   );
    
    
    out("dT_pp_pbar_LAB_TanNg")
    WriteCS(dT_pp_pbar_LAB_TanNg,       "dT_pp_pbar_LAB_TanNg.txt"      );
    out("dT_pHe_pbar_LAB_TanNg")
    WriteCS(dT_pHe_pbar_LAB_TanNg,      "dT_pHe_pbar_LAB_TanNg.txt"     );
    std::system("cp dT_pHe_pbar_LAB_TanNg.txt dT_Hep_pbar_LAB_TanNg.txt");
    
    
    
}

void CRACS::CS_lab_tab::TransformGalprop(std::string filename, int parametrization, int dLog){

    funOut(CS_lab_tab::TransformAllGalpropTransformAllGalprop)
    
    std::string table = "#";
    QString line("");
    line += QString("%1 ").arg( "T/n (projectile"       ).leftJustified(20);
    line += QString("%1 ").arg( "T/n (pbar)"            ).leftJustified(20);
    line += QString("%1 ").arg( "cross section pp"      ).leftJustified(20);
    line += QString("%1 ").arg( "p He"                  ).leftJustified(20);
    line += QString("%1 ").arg( "He p"                  ).leftJustified(20);
    line += QString("%1 ").arg( "He He"                 ).leftJustified(20);
    table += line.toStdString();
    
    for (int i = 0; i<=dLog*7; i++) {
        double Tn_proj = pow(10, 1.*i/dLog);
        for (int j = 0; j<=dLog*8; j++) {
            double Tn_prod = pow(10, 1.*j/dLog-1.);
            
            QString line("\n");
            line += QString("%1 ").arg( Tn_proj ).leftJustified(20);
            line += QString("%1 ").arg( Tn_prod ).leftJustified(20);
            line += QString("%1 ").arg(  GetCrossSection(Tn_proj, Tn_prod, CS::PP_PBAR,   parametrization)  ).leftJustified(20);
            line += QString("%1 ").arg(  GetCrossSection(Tn_proj, Tn_prod, CS::PHe_PBAR,  parametrization)  ).leftJustified(20);
            line += QString("%1 ").arg(  GetCrossSection(Tn_proj, Tn_prod, CS::HeP_PBAR,  parametrization)  ).leftJustified(20);
            line += QString("%1 ").arg(  GetCrossSection(Tn_proj, Tn_prod, CS::HeHe_PBAR, parametrization)  ).leftJustified(20);
            
            table += line.toStdString();
        }
        
    }
    //std::cout << table << std::endl;
    FileTool::WriteStringToFile( table, filename );

}

void CRACS::CS_lab_tab::TransformAllGalprop(){
    
    out(dT_pp_pbar_LAB_WinklerWitHypWitNbar(1000, 100))
    out(dT_pHe_pbar_LAB_WinklerWitHypWitNbar(1000, 100))
    out(dT_Hep_pbar_LAB_WinklerWitHypWitNbar(1000, 100))
    out(dT_HeHe_pbar_LAB_WinklerWitHypWitNbar(1000, 100))
    
    
    TransformGalprop("pbar_CS_diMauro.txt",   CS::DI_MAURO,     30);
    TransformGalprop("pbar_CS_diMauro12.txt", CS::DI_MAURO12,   30);
    TransformGalprop("pbar_CS_Duperray.txt",  CS::DUPERRAY,     30);
    TransformGalprop("pbar_CS_KR.txt",        CS::KACHELRIESS,  30);
    TransformGalprop("pbar_CS_Winkler.txt",   CS::WINKLER,  30);
    TransformGalprop("pbar_CS_Winkler_withHYPwithNBAR.txt",   CS::WINKLER_withHYPwithNBAR,  30);
    
}


void CRACS::CS_lab_tab::test(){
    
//    std::cout << GetCrossSection        ( 50,                10,                 CS::PP_PBAR, CS::DI_MAURO12 ) << "       ";
//    std::cout << dT_pp_pbar_LAB_diMauro ( 50,                10                                              ) << "       ";
//    std::cout << CS::dE_pp_pbar_LAB     ( 50  +fMass_proton, 10+fMass_proton,        1000000, CS::DI_MAURO12 ) << std::endl;
//
//    std::cout << GetCrossSection        (100,                10,                 CS::PP_PBAR, CS::DI_MAURO12 ) << "       ";
//    std::cout << dT_pp_pbar_LAB_diMauro (100,                10                                              ) << "       ";
//    std::cout << CS::dE_pp_pbar_LAB     (100  +fMass_proton, 10+fMass_proton,        1000000, CS::DI_MAURO12 )*1e-31 << std::endl;
  
    out(dT_pp_pbar_LAB_Duperray (20, 5))
    out(dT_pHe_pbar_LAB_Duperray(20, 5))
    out(dT_pHe_pbar_LAB_Duperray(20, 5)/dT_pp_pbar_LAB_Duperray (20, 5))
    std::cout << std::endl;
    out(dT_pp_pbar_LAB_Duperray (100, 10))
    out(dT_pHe_pbar_LAB_Duperray(100, 10))
    out(dT_pHe_pbar_LAB_Duperray(100, 10)/dT_pp_pbar_LAB_Duperray (100, 10))
    std::cout << std::endl;
    out(dT_pp_pbar_LAB_Duperray (100, 20))
    out(dT_pHe_pbar_LAB_Duperray(100, 20))
    out(dT_pHe_pbar_LAB_Duperray(100, 20)/dT_pp_pbar_LAB_Duperray (100, 20))
    std::cout << std::endl;
    
    
};





