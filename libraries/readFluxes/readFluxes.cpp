#include "readFluxes.h"
#include "fileTools.h"
#include "configHandler.h"

#include "iostream"
#include "math.h"
#include "stdlib.h"


// All fluxes are read and saved in MV/MeV, m, s, sr


double CRACS::getBinExpectationValueForPowerLaw2( double binStart, double binEnd, double power  ){
    return (pow(binStart, 2+power) - pow(binEnd, 2+power))*(power+1)/(pow(binStart, 1+power) - pow(binEnd, 1+power))/(power+2);
};

void CRACS::ReadFluxes::Read(){
    
    
    std::string mySoft = CRACS::ConfigHandler::GetInstance()->SoftwarePath();
    
    std::string ams             = mySoft+"/data/AMS02/";
    std::string cream           = mySoft+"/data/CreamII/";
    std::string heao            = mySoft+"/data/HEAO/";
    std::string pamela          = mySoft+"/data/Pamela/";
    std::string voyager         = mySoft+"/data/Voyager/";
    std::string atic            = mySoft+"/data/ATIC02/";
    
    //
    // AMS
    //
    CRACS::Graph protonAMS = ReadAMS(  ams+"proton_ams02.txt"    );
    protonAMS.SetTitle( "Proton AMS; R [GV]; #Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]" );
    fFluxes.    push_back(protonAMS);
    fSpecies.   push_back("proton");
    fExperiment.push_back("ams");
    fZ.         push_back(1);
    fA.         push_back(1);
    fType.      push_back(CRACS::RIGIDITY);
    
    CRACS::Graph heliumAMS = ReadAMS(  ams+"helium_ams02.txt"    );
    heliumAMS.SetTitle( "Helium AMS; R [GV]; #Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]" );
    fFluxes.    push_back(heliumAMS);
    fSpecies.   push_back("helium");
    fExperiment.push_back("ams");
    fZ.         push_back(2);
    fA.         push_back(4);
    fType.      push_back(CRACS::RIGIDITY);
    
    CRACS::Graph lithiumAMS = ReadAMS_R(  ams+"lithium_ams02.txt"    );
    lithiumAMS.SetTitle( "Lithium AMS; R [GV]; #Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]" );
    fFluxes.    push_back(lithiumAMS);
    fSpecies.   push_back("lithium");
    fExperiment.push_back("ams");
    fZ.         push_back(3);
    fA.         push_back(7);
    fType.      push_back(CRACS::RIGIDITY);
    
    CRACS::Graph carbonAMS = ReadAMS_R(  ams+"carbon_ams02.txt"    );
    carbonAMS.SetTitle( "Carbon AMS; R [GV]; #Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]" );
    fFluxes.    push_back(carbonAMS);
    fSpecies.   push_back("carbon");
    fExperiment.push_back("ams");
    fZ.         push_back(6);
    fA.         push_back(12);
    fType.      push_back(CRACS::RIGIDITY);
    
    CRACS::Graph oxygenAMS = ReadAMS_R(  ams+"oxygen_ams02.txt"    );
    oxygenAMS.SetTitle( "Oxygen AMS; R [GV]; #Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]" );
    fFluxes.    push_back(oxygenAMS);
    fSpecies.   push_back("oxygen");
    fExperiment.push_back("ams");
    fZ.         push_back(8);
    fA.         push_back(16);
    fType.      push_back(CRACS::RIGIDITY);
    
    CRACS::Graph nitrogenAMS = ReadAMS_R(  ams+"nitrogen_ams02.txt"    );
    nitrogenAMS.SetTitle( "Oxygen AMS; R [GV];#Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]" );
    fFluxes.    push_back(nitrogenAMS);
    fSpecies.   push_back("nitrogen");
    fExperiment.push_back("ams");
    fZ.         push_back(7);
    fA.         push_back(14);
    fType.      push_back(CRACS::RIGIDITY);
    
    
    CRACS::Graph boronAMS = ReadAMS_Tn(  ams+"boron_ams02.txt"    );
    boronAMS.SetTitle( "Boron AMS;  #Phi E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]" );
    fFluxes.    push_back(boronAMS);
    fSpecies.   push_back("boron");
    fExperiment.push_back("ams");
    fZ.         push_back(5);
    fA.         push_back(11);
    fType.      push_back(CRACS::EKINPERN);
    
    
    CRACS::Graph berylliumAMS = ReadAMS_Tn(  ams+"beryllium_ams02.txt"    );
    berylliumAMS.SetTitle( "Beryllium AMS;  #Phi E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]" );
    fFluxes.    push_back(berylliumAMS);
    fSpecies.   push_back("beryllium");
    fExperiment.push_back("ams");
    fZ.         push_back(4);
    fA.         push_back(9);
    fType.      push_back(CRACS::EKINPERN);
    
    
    CRACS::Graph bovercAMS = ReadAMS_boverc(  ams+"boverc_ams02.txt"    );
    bovercAMS.SetTitle( "B/C AMS; R [GV]; B/C" );
    fFluxes.    push_back(bovercAMS);
    fSpecies.   push_back("boverc");
    fExperiment.push_back("ams");
    fZ.         push_back(5);
    fA.         push_back(11);
    fType.      push_back(CRACS::RATIO);
    
    
    CRACS::Graph pbaroverpAMS = ReadAMS_pbaroverp(  ams+"pbar_ams02.txt"    );
    pbaroverpAMS.SetTitle( "pbar/p AMS (official); R [GV]; pbar/p" );
    fFluxes.    push_back(pbaroverpAMS);
    fSpecies.   push_back("pbaroverp");
    fExperiment.push_back("ams");
    fZ.         push_back(-1);
    fA.         push_back(1);
    fType.      push_back(CRACS::RATIO);
    
    CRACS::Graph pbarAMS = ReadAMS_pbar(  ams+"pbar_ams02.txt"    );
    pbarAMS.SetTitle( "pbar AMS ; R [GV]; #Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]" );
    fFluxes.    push_back(pbarAMS);
    fSpecies.   push_back("pbar");
    fExperiment.push_back("ams");
    fZ.         push_back(-1);
    fA.         push_back(1);
    fType.      push_back(CRACS::RIGIDITY);
   
    
    CRACS::Graph electron_AMS = ReadAMS_lelptons(  ams+"electron_ams02.txt"    );
    pbarAMS.SetTitle( "electron AMS; E [GeV]; #Phi E^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]" );
    fFluxes.    push_back(electron_AMS);
    fSpecies.   push_back("electron");
    fExperiment.push_back("ams");
    fZ.         push_back(-1);
    fA.         push_back(0);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph positron_AMS = ReadAMS_lelptons(  ams+"positron_ams02.txt"    );
    pbarAMS.SetTitle( "positron AMS; E [GeV]; #Phi E^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]" );
    fFluxes.    push_back(positron_AMS);
    fSpecies.   push_back("positron");
    fExperiment.push_back("ams");
    fZ.         push_back(1);
    fA.         push_back(0);
    fType.      push_back(CRACS::EKINPERN);
    
    
    //
    // PAMELA
    //
    CRACS::Graph protonPamela = ReadPamela_pHe(  pamela+"proton_pamela.txt"    );
    protonPamela.SetTitle( "Proton Pamela; R [GV]; #Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]" );
    fFluxes.    push_back(protonPamela);
    fSpecies.   push_back("proton");
    fExperiment.push_back("pamela");
    fZ.         push_back(1);
    fA.         push_back(1);
    fType.      push_back(CRACS::RIGIDITY);
    
    CRACS::Graph pbarPamela = ReadPamela_pHe(  pamela+"pbar_pamela.txt"    );
    pbarPamela.SetTitle( "Pbar Pamela; R [GV]; #Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]" );
    fFluxes.    push_back(pbarPamela);
    fSpecies.   push_back("pbar");
    fExperiment.push_back("pamela");
    fZ.         push_back(1);
    fA.         push_back(-1);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph heliumPamela = ReadPamela_pHe(  pamela+"helium_pamela.txt"    );
    heliumPamela.SetTitle( "helium Pamela; R [GV]; #Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]" );
    fFluxes.    push_back(heliumPamela);
    fSpecies.   push_back("helium");
    fExperiment.push_back("pamela");
    fZ.         push_back(2);
    fA.         push_back(4);
    fType.      push_back(CRACS::RIGIDITY);
    
    
    CRACS::Graph carbonPamela = ReadPamela(  pamela+"carbon_pamela.txt"    );
    carbonPamela.SetTitle( "Carbon Pamela; R [GV]; #Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]" );
    fFluxes.    push_back(carbonPamela);
    fSpecies.   push_back("carbon");
    fExperiment.push_back("pamela");
    fZ.         push_back(6);
    fA.         push_back(12);
    fType.      push_back(CRACS::RIGIDITY);
    //
    //  CREAM
    //
    CRACS::Graph protonCream = ReadCream_PandHe(  cream+"proton_cream.txt"    );
    protonCream.SetTitle( "Proton flux, data form Cream; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(protonCream);
    fSpecies.   push_back("proton");
    fExperiment.push_back("cream");
    fZ.         push_back(1);
    fA.         push_back(1);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph heliumCream = ReadCream_PandHe(  cream+"helium_cream.txt"    );
    heliumCream.SetTitle( "Helium flux, data form Cream; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(heliumCream);
    fSpecies.   push_back("helium");
    fExperiment.push_back("cream");
    fZ.         push_back(2);
    fA.         push_back(4);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph carbonCREAM = ReadCream(  cream+"carbon_cream.txt"    );
    carbonCREAM.SetTitle( "Carbon flux, data form Cream; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(carbonCREAM);
    fSpecies.   push_back("carbon");
    fExperiment.push_back("cream");
    fZ.         push_back(6);
    fA.         push_back(12);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph nitrogenCREAM = ReadCream(  cream+"nitrogen_cream.txt"    );
    nitrogenCREAM.SetTitle( "Nitrogen flux, data form Cream; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(nitrogenCREAM);
    fSpecies.   push_back("nitrogen");
    fExperiment.push_back("cream");
    fZ.         push_back(7);
    fA.         push_back(14);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph oxygenCREAM = ReadCream(  cream+"oxygen_cream.txt"    );
    oxygenCREAM.SetTitle( "Oxygen flux, data form Cream; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(oxygenCREAM);
    fSpecies.   push_back("oxygen");
    fExperiment.push_back("cream");
    fZ.         push_back(8);
    fA.         push_back(16);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph neonCREAM = ReadCream(  cream+"neon_cream.txt"    );
    neonCREAM  .SetTitle( "Neon flux, data form Cream; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(neonCREAM);
    fSpecies.   push_back("neon");
    fExperiment.push_back("cream");
    fZ.         push_back(10);
    fA.         push_back(20);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph magnesiumCREAM = ReadCream(  cream+"magnesium_cream.txt"    );
    magnesiumCREAM.SetTitle( "Magnesium flux, data form Cream; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(magnesiumCREAM);
    fSpecies.   push_back("magnesium");
    fExperiment.push_back("cream");
    fZ.         push_back(12);
    fA.         push_back(24);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph siliconCREAM = ReadCream(  cream+"silicon_cream.txt"    );
    siliconCREAM.SetTitle( "Silicon flux, data form Cream; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(siliconCREAM);
    fSpecies.   push_back("silicon");
    fExperiment.push_back("cream");
    fZ.         push_back(14);
    fA.         push_back(28);
    fType.      push_back(CRACS::EKINPERN);
    
    
    
    //
    // HEAO
    //
    CRACS::Graph nitrogenHEAO = ReadHEAO( heao+"nitrogen_heao.txt"    );
    nitrogenHEAO.SetTitle( "Nitrogen flux, data form HEAO; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(nitrogenHEAO);
    fSpecies.   push_back("nitrogen");
    fExperiment.push_back("heao");
    fZ.         push_back(7);
    fA.         push_back(14);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph oxygenHEAO = ReadHEAO( heao+"oxygen_heao.txt"    );
    oxygenHEAO.SetTitle( "Oxygen flux, data form HEAO; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(oxygenHEAO);
    fSpecies.   push_back("oxygen");
    fExperiment.push_back("heao");
    fZ.         push_back(8);
    fA.         push_back(16);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph neonHEAO = ReadHEAO( heao+"neon_heao.txt"    );
    neonHEAO.SetTitle( "Neon flux, data form HEAO; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(neonHEAO);
    fSpecies.   push_back("neon");
    fExperiment.push_back("heao");
    fZ.         push_back(10);
    fA.         push_back(20);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph magnesiumHEAO = ReadHEAO( heao+"magnesium_heao.txt"    );
    magnesiumHEAO.SetTitle( "Magnesium flux, data form HEAO; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(magnesiumHEAO);
    fSpecies.   push_back("magnesium");
    fExperiment.push_back("heao");
    fZ.         push_back(12);
    fA.         push_back(24);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph siliconHEAO = ReadHEAO( heao+"nitrogen_heao.txt"    );
    siliconHEAO.SetTitle( "Silicon flux, data form HEAO; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(siliconHEAO);
    fSpecies.   push_back("silicon");
    fExperiment.push_back("heao");
    fZ.         push_back(14);
    fA.         push_back(28);
    fType.      push_back(CRACS::EKINPERN);
    
    //
    //Voyager
    //
    CRACS::Graph protonVoyager = ReadVoyager(  voyager+"proton_voyager_EkinPerNuc.txt"    );
    protonVoyager.SetTitle( "Proton Voyager; Ekin [GeV]; #Phi Ekin/n^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]" );
    fFluxes.    push_back(protonVoyager);
    fSpecies.   push_back("proton");
    fExperiment.push_back("voyager");
    fZ.         push_back(1);
    fA.         push_back(1);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph heliumVoyager = ReadVoyager(  voyager+"helium_voyager_EkinPerNuc.txt");
    heliumVoyager.SetTitle( "Helium Voyager; Ekin [GeV]; #Phi Ekin/n^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]" );
    fFluxes.    push_back(heliumVoyager);
    fSpecies.   push_back("helium");
    fExperiment.push_back("voyager");
    fZ.         push_back(2);
    fA.         push_back(4);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph carbonVoyager = ReadVoyager(  voyager+"carbon_voyager_EkinPerNuc.txt");
    carbonVoyager.SetTitle( "Helium Voyager; Ekin [GeV]; #Phi Ekin/n^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]" );
    fFluxes.    push_back(carbonVoyager);
    fSpecies.   push_back("carbon");
    fExperiment.push_back("voyager");
    fZ.         push_back(6);
    fA.         push_back(12);
    fType.      push_back(CRACS::EKINPERN);
    
    //
    // ATIC 02
    //
    CRACS::Graph protonATIC = ReadAtic( atic+"proton_atic.txt"    );
    protonATIC.SetTitle( "Proton flux, data form Atic; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(protonATIC);
    fSpecies.   push_back("proton");
    fExperiment.push_back("atic");
    fZ.         push_back(1);
    fA.         push_back(1);
    fType.      push_back(CRACS::EKINPERN);
    
    CRACS::Graph heliumATIC = ReadAtic( atic+"helium_atic.txt"    , 4);
    heliumATIC.SetTitle( "Helium flux, data form Atic; E_{kin}/n [GeV]; #Phi (E_{kin}/n)^{2.7} [GeV^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
    fFluxes.    push_back(heliumATIC);
    fSpecies.   push_back("helium");
    fExperiment.push_back("atic");
    fZ.         push_back(2);
    fA.         push_back(4);
    fType.      push_back(CRACS::EKINPERN);
    
    
}

CRACS::Graph CRACS::ReadFluxes::GetSpectrum(std::string species, std::string experiment, CRACS::spectrumType type){
    CRACS::Graph g;
    for (int i = 0; i<fSpecies.size(); i++) {
        if (fSpecies.at(i)==species) {
            if (fExperiment.at(i)==experiment) {
                if (fType.at(i)==type) {
                    g = fFluxes.at(i);
                    g.SetName(fSpecies.at(i)+"_"+fExperiment.at(i));
                    g.Scale(1./1000, pow(1./1000, 1.7));
                    if (type==RIGIDITY) {
                        g.SetTitle(fSpecies.at(i)+" "+fExperiment.at(i)+";R [GV]; #Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]");
                    }else if (type==EKINPERN){
                        g.SetTitle(fSpecies.at(i)+" "+fExperiment.at(i)+";E_{kin}/n [GeV/n]; #Phi (E_{kin}/n)^{2.7} [(GeV/n)^{2.7} m^{-2} sr^{-1} s^{-1} GeV^{-1}]");
                    }
                    
                    return g;
                }else{
                    g = fFluxes.at(i);
                    g.SetName(fSpecies.at(i)+"_"+fExperiment.at(i));
                    if (fType.at(i)==CRACS::RATIO) {
                        std::cout << "Warning: there is no way to convert ratios" << std::endl;
                        continue;
                    }
                    if (fType.at(i)==CRACS::RIGIDITY) {
                        g.ConvertSpectrumFromRigidityToEkinPerN(fZ.at(i), fA.at(i), 2.7);
                        g.SetTitle(fSpecies.at(i)+" "+fExperiment.at(i)+";E_{kin}/n [GeV/n]; #Phi (E_{kin}/n)^{2.7} [(GeV/n)^{2.7} m^{-2} sr^{-1} s^{-1} (GeV/n)^{-1}]");
                    }else{
                        g.ConvertSpectrumFromEkinPerNToRigidity(fZ.at(i), fA.at(i), 2.7);
                        g.SetTitle(fSpecies.at(i)+" "+fExperiment.at(i)+"R [GV]; #Phi R^{2.7} [GV^{2.7} m^{-2} sr^{-1} s^{-1} GV^{-1}]");
                    }
                    g.Scale(1./1000, pow(1./1000, 1.7));
                    
                    return g;
                }
            }
        }
    }
    std::cout << "Warning: Spectrum not found. Return empty graph for:" << std::endl;
    std::cout << species    << std::endl;
    std::cout << experiment << std::endl;
    
    return g;
}



CRACS::Graph CRACS::ReadFluxes::ReadAMS(std::string filename){
    
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(10, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        
        double R = getBinExpectationValueForPowerLaw2(f.NumberTable(row, 0), f.NumberTable(row, 1), -2.7);
        double F = f.NumberTable(row, 2)*f.NumberTable(row, 9);
        double helper = 0;
        for (int col = 3; col<9; col++) {
            helper += pow(f.NumberTable(row, col), 2);
        }
        double sigma = sqrt(helper)*f.NumberTable(row, 9);
        graph.AddPoint(R*1000, F*pow(R*1000, 2.7)/1000, 0, sigma*pow(R*1000, 2.7)/1000);
    }
    
    return  graph;
    
}

CRACS::Graph CRACS::ReadFluxes::ReadAMS_pbar(std::string filename){
    
   CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(14, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        
        double R = getBinExpectationValueForPowerLaw2(f.NumberTable(row, 0), f.NumberTable(row, 1), -2.7);
        double F = f.NumberTable(row, 10)*f.NumberTable(row, 13);
        double helper = 0;
        for (int col = 11; col<=12; col++) {
            helper += pow(f.NumberTable(row, col), 2);
        }
        double sigma = sqrt(helper)*f.NumberTable(row, 13);
        graph.AddPoint(R*1000, F*pow(R*1000, 2.7)/1000, 0, sigma*pow(R*1000, 2.7)/1000);
    }
    
    return  graph;
    
}

CRACS::Graph CRACS::ReadFluxes::ReadAMS_pbaroverp(std::string filename){
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(14, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        double R = getBinExpectationValueForPowerLaw2(f.NumberTable(row, 0), f.NumberTable(row, 1), -2.7);
        double F = f.NumberTable(row, 3)*f.NumberTable(row, 9);
        double helper = 0;
        for (int col = 4; col<=7; col++) {
            helper += pow(f.NumberTable(row, col), 2);
        }
        double sigma = sqrt(helper)*f.NumberTable(row, 9);
        graph.AddPoint(R*1000., F, 0, sigma);
    }
    return  graph;
}


CRACS::Graph CRACS::ReadFluxes::ReadAMS_R(std::string filename){
    
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(6, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        double R = f.NumberTable(row, 0);
        double F = f.NumberTable(row, 1);
        double sigma = f.NumberTable(row, 4);
        graph.AddPoint(R*1000, F*pow(1000, 2.7)/1000, 0, sigma*pow(1000, 2.7)/1000);
    }
    return  graph;
    
}

CRACS::Graph CRACS::ReadFluxes::ReadAMS_Tn(std::string filename){
    
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(6, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        double E        = f.NumberTable(row, 0);
        double flux     = f.NumberTable(row, 1);
        double sigma    = f.NumberTable(row, 4);
        graph.AddPoint( E*1000, pow(1000, 2.7)*flux/1000,  0, pow(1000, 2.7)*sigma/1000 );
    }
    return  graph;
    
}


CRACS::Graph CRACS::ReadFluxes::ReadAMS_lelptons(std::string filename){
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(6, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        double Ekin = f.NumberTable(row, 0);
        double F    = f.NumberTable(row, 3);
        double sigma= f.NumberTable(row, 4);
        graph.AddPoint(Ekin*1000, F*pow(Ekin*1000, 2.7)/1000, 0, sigma*pow(Ekin*1000, 2.7)/1000);
    }
    
    return  graph;
    
}


                       
CRACS::Graph CRACS::ReadFluxes::ReadAMS_boverc(std::string filename){
    
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(9, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        double R = getBinExpectationValueForPowerLaw2(f.NumberTable(row, 0), f.NumberTable(row, 1), -2.7);
        double F = f.NumberTable(row, 2);
        double helper = 0;
        for (int col = 3; col<9; col++) {
            helper += pow(f.NumberTable(row, col), 2);
        }
        double sigma = sqrt(helper);
        graph.AddPoint(R*1000, F, 0, sigma);
    }
    return  graph;
    
}
                       


CRACS::Graph CRACS::ReadFluxes::ReadPamela(std::string filename){
    
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(8, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        double R = f.NumberTable(row, 2);
        double F = f.NumberTable(row, 4)*f.NumberTable(row, 7);
        double sigma = sqrt( pow(f.NumberTable(row, 5), 2) + pow(f.NumberTable(row, 6), 2) )*f.NumberTable(row, 7);
        graph.AddPoint(R*1000, F*pow(R*1000, 2.7)/1000, 0, sigma*pow(R*1000, 2.7)/1000);
    }
    return  graph;
    
}

CRACS::Graph CRACS::ReadFluxes::ReadPamela_pHe(std::string filename){
    
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(4, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        double R = f.NumberTable(row, 0);
        double F = f.NumberTable(row, 1);
        double sigma = sqrt( pow(f.NumberTable(row, 2), 2) + pow(f.NumberTable(row, 3), 2) );
        graph.AddPoint(R*1000, F*pow(R*1000, 2.7)/1000, 0, sigma*pow(R*1000, 2.7)/1000);
    }
    return  graph;
    
}

CRACS::Graph CRACS::ReadFluxes::ReadAtic(std::string filename, int A){
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(3, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        double Ekin = f.NumberTable(row, 0);
        double F    = f.NumberTable(row, 1);
        double sigma = f.NumberTable(row, 2);
        graph.AddPoint(Ekin/A*1000, F*pow(Ekin/A*1000, 2.7)/1000*A, 0, sigma*pow(Ekin/A*1000, 2.7)/1000*A);
    }
    return  graph;
}


CRACS::Graph CRACS::ReadFluxes::ReadCream_PandHe(std::string filename){
    
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(6, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        double EperN    = getBinExpectationValueForPowerLaw2(f.NumberTable(row, 0), f.NumberTable(row, 1), -2.7);
        double F        = f.NumberTable(row, 2)*f.NumberTable(row, 5);
        double sigma = sqrt( pow(f.NumberTable(row, 3), 2) + pow(f.NumberTable(row, 4), 2) )*f.NumberTable(row, 5);
        graph.AddPoint(  EperN*1000, F*pow(EperN*1000, 2.7)/1000, 0, sigma*pow(EperN*1000, 2.7)/1000  );
    }
    return  graph;
    
}

CRACS::Graph CRACS::ReadFluxes::ReadCream(std::string filename){
    
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(8, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        double EperN    = f.NumberTable(row, 2);
        double F        = f.NumberTable(row, 4)*f.NumberTable(row, 7);
        double sigma = sqrt( pow(f.NumberTable(row, 5), 2) + pow(f.NumberTable(row, 6), 2) )*f.NumberTable(row, 7);
        graph.AddPoint( EperN*1000, F*pow(EperN*1000, 2.7)/1000, 0, sigma*pow(EperN*1000, 2.7)/1000 );
    }
    
    return  graph;
    
}




CRACS::Graph CRACS::ReadFluxes::ReadHEAO(std::string filename){
    
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(4, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        double EperN    = f.NumberTable(row, 1);
        double F        = f.NumberTable(row, 2);
        double sigma    = f.NumberTable(row, 3);
        graph.AddPoint(EperN*1000, F*pow(EperN*1000, 2.7)/1000, 0, sigma*pow(EperN*1000, 2.7)/1000);
    }
    
    return  graph;
    
}



CRACS::Graph CRACS::ReadFluxes::ReadVoyager(std::string filename){
    
    CRACS::Graph graph;
    CRACS::FileTool f(filename);
    f.ExtractNumberTable(6, " ", true);
    for (int row = 0; row<f.NumberTableGetNrows(); row++) {
        double E        = f.NumberTable(row, 0);
        double flux     = f.NumberTable(row, 1);
        double sigma    = f.NumberTable(row, 4);
        graph.AddPoint( E*1000, pow(E*1000, 2.7)*flux/1000,  0, pow(E*1000, 2.7)*sigma/1000 );
    }
    return  graph;
    
}

