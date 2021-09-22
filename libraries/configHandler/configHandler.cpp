#include "configHandler.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

#include <QString>
#include <QStringList>

CRACS::ConfigHandler* CRACS::ConfigHandler::fInstance = NULL;

CRACS::ConfigHandler* CRACS::ConfigHandler::GetInstance(){
    if(!fInstance)
        fInstance = new ConfigHandler();
    return fInstance;
}

void CRACS::ConfigHandler::SetInput(int argc, char *argv[], std::string description){
    
    std::string s = "";
    for(int i=1;i<argc;i++)
        s = s + argv[i]+" ";
    
    fOptionString = s;
    fName = argv[0];
    fDescription = description;
    fOk = true;
    fPositionOptions = 0;
    
    ParseOption();

    //Read some standard options

    std::string mySoft = "";
    if (const char* env_p = std::getenv("CRACS")) {
        mySoft = env_p;
    }else{
        std::cout << "Error: environment variable CRACS not set." << std::endl;
    }

    fSoftwarePath  = mySoft;
    AddOption        (  "softwarePath",      fSoftwarePath,           "Full path of of the software package location. Default: $CRACS"                         );
};


void CRACS::ConfigHandler::AddOptionPosition(std::string &value, std::string description, bool mandatory){
    fPositionOptionsDescription.push_back(description);
    fPositionOptionsMandotory.push_back(mandatory);
    if (!fOk){
        return;
    }
    
    if (fPositionOptions<fPositionOptionsValue.size()) {
        value = fPositionOptionsValue.at(fPositionOptions);
        fPositionOptions++;
    }else if(mandatory){
        value = "not defined";
        fOk = false;
        fError += "To few positional arguments. \n";
    }
};

void CRACS::ConfigHandler::AddOption(std::string option, std::string& value, std::string description){
    fOptionsDescription.push_back(description);
    fOptionsExpected.push_back(option);
    if (!fOk)
        return;
    for (int i = 0; i<fOptions.size(); i++) {
        if (fOptions.at(i)==option) {
            value = fOptionValues.at(i);
        }
    }
};

void CRACS::ConfigHandler::AddOptionInt(std::string option, int& value, std::string description){
    fOptionsDescription.push_back(description);
    fOptionsExpected.push_back(option);
    if (!fOk)
        return;
    for (int i = 0; i<fOptions.size(); i++) {
        if (fOptions.at(i)==option) {
            bool toInt = true;
            QString sToInt(fOptionValues.at(i).c_str());
            int v = sToInt.toInt(&toInt);
            if (!toInt) {
                fOk = false;
                fError += "Cannot convert argument of  '--" + option + " to int!'. \n";
            }
            value = v;
        }
    }
};

void CRACS::ConfigHandler::AddOptionDouble(std::string option, double& value, std::string description){
    fOptionsDescription.push_back(description);
    fOptionsExpected.push_back(option);
    if (!fOk)
        return;
    for (int i = 0; i<fOptions.size(); i++) {
        if (fOptions.at(i)==option) {
            bool toDoub = true;
            QString sToDoub(fOptionValues.at(i).c_str());
            double v = sToDoub.toDouble(&toDoub);
            if (!toDoub) {
                fOk = false;
                fError += "Cannot convert argument of  '--" + option + " to double!'. \n";
            }
            value = v;
        }
    }
};

void CRACS::ConfigHandler::AddOptionTrue(std::string option, bool& stored, std::string description){
    fOptionsDescription.push_back(description);
    fOptionsExpected.push_back(option);
    if (!fOk)
        return;
    for (int i = 0; i<fOptions.size(); i++) {
        if (fOptions.at(i)==option) {
            stored = true;
            return;
        }
    }
};

void CRACS::ConfigHandler::ParseOption(){
    
    QString option(fOptionString.c_str());
    if (option.indexOf("-h")>-1) {
        fOk = false;
        return;
    }
    
    QStringList listPos =  option.split(" ", QString::SkipEmptyParts);
    for (int i = 0; i<listPos.size(); i++) {
        if (listPos.at(i).left(1) == "-")
            break;
        fPositionOptionsValue.push_back( (listPos.at(i)).toStdString() );
    }
    
    
    QStringList list =  option.split("--", QString::SkipEmptyParts);
    int i = 0;
    if (option.left(2)=="--") {
        i=0;
    }else{
        i=1;
    }
    for (; i<list.size(); i++) {
        QString option = list.at(i);
        QStringList optionList = option.split(" ", QString::SkipEmptyParts);
        
        fOptions.push_back(optionList.at(0).toStdString());
        if (optionList.size()==2) {
            fOptionValues.push_back(optionList.at(1).toStdString());
            
        }else{
            fOptionValues.push_back("");
        }
        if (optionList.size()>2) {
            fOk = false;
            fError += ("Option '--" + option.toStdString() + "' as incorrect pattern.\n");
        }
        
    }
    
    
};

void CRACS::ConfigHandler::CheckOption(){
    if (fOk && fPositionOptions<fPositionOptionsValue.size()) {
        fOk = false;
        fError += "To many positional arguments. \n";
    }
    
    if (!fOk) {
        PrintHelp();
        exit(EXIT_FAILURE);
    }
}


void CRACS::ConfigHandler::PrintHelp(){
    std::cout << fError << std::endl;
    std::string use = fName + "  ";
    for (int i = 0; i<fPositionOptionsDescription.size(); i++) {
        if (fPositionOptionsMandotory.at(i)) {
            use += QString("VALUE_%1 ").arg(i+1).toStdString();
        }else{
            use += QString("[VALUE_%1] ").arg(i+1).toStdString();
        }
    }
    
    for (int i = 0; i<fOptionsExpected.size(); i++) {
        use += QString("[--%1 VALUE] ").arg(fOptionsExpected.at(i).c_str()).toStdString();
    }
    
    std::cout << "Use: " << use << std::endl;
    
    
    std::cout << fDescription << std::endl;
    std::cout << std::endl;
    
    for (int i = 0; i<fPositionOptionsDescription.size(); i++) {
        std::cout << "  " << "VALUE_" << i+1 << " \t\t" << fPositionOptionsDescription.at(i) <<  std::endl;
    }
    
    for (int i = 0; i<fOptionsExpected.size(); i++) {
        std::cout << "  " << "--" << fOptionsExpected.at(i) << " \t\t" << fOptionsDescription.at(i) <<  std::endl;
    }
    
    
    
    
}
