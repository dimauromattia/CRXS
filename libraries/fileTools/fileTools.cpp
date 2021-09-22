#include <fileTools.h>
#include <iostream>

CRACS::FileTool::FileTool(std::string filename){
    QString f(filename.c_str());
    QStringList pathList = f.split("/");
    fFileName = pathList.at(pathList.size()-1).toStdString();
    
    fFilePath = "";
    if (f.left(1)=="/")
        fFilePath +="/";
    for (int i = 0; i<pathList.size()-1; i++) {
        fFilePath += pathList.at(i).toStdString() + "/";
    }
    fReadingOK = ReadFile();
        
    BackToOriginal();
};

CRACS::FileTool::~FileTool(){

}

bool CRACS::FileTool::ReadFile(){
    
    QFile file(QString((fFilePath+fFileName).c_str()));
    if (!file.open(QIODevice::ReadOnly)) {
        std::cout << "File "<< fFileName <<" could not be opened!" << std::endl;
        return false;
    }
    QTextStream in(&file);
    
    while (!in.atEnd()) {
        fOriginal.push_back( in.readLine() );
    }
    file.close();
    return true;
};

bool CRACS::FileTool::Exists(std::string filename){
    
    QFile file(  QString( (filename).c_str() )  );
    if (!file.open(QIODevice::ReadOnly)) {
        std::cout << "File " << filename << " could not be opened!" << std::endl;
        return false;
    }
    file.close();
    return true;
};

void CRACS::FileTool::BackToOriginal(){
    fModified.clear();
    for (int i = 0; i<fOriginal.size(); i++) {
        fModified.push_back(fOriginal.at(i));
    }

    fNumberTable.clear();
    fNumberTableAvailable = false;
    fNumberTableColumns = 0;
    fNumberTableRows = 0;
};

void CRACS::FileTool::FindAndReplace(std::string find, std::string replace, bool caseSensitiv){
    
    for (int i = 0; i<fModified.size(); i++) {
        if (caseSensitiv) {
            fModified.at(i).replace(find.c_str(), replace.c_str());
        }else{
            fModified.at(i).replace(find.c_str(), replace.c_str(), Qt::CaseInsensitive);
        }
    }

};

bool CRACS::FileTool::ExtractNumberTable(int NumberOfColumn, std::string separator, bool skipEmptyParts, int lineStart, int lineStop){
    
    for (int line = 1; line<=fModified.size(); line++) {
        
        if (line<lineStart) {
            continue;
        }
        if (line>lineStop && lineStop!=-1) {
            break;
        }
        
        QString lineString = fModified.at(line-1);
        while (lineString.left(1)==" ") {
            lineString = lineString.mid(1);
        }
        
        
        
        if (lineString =="")
            continue;
        if (lineString.left(1)=="*")
            continue;
        if (lineString.left(1)=="#")
            continue;
        if (lineString.left(2)=="//")
            continue;
        
        
        QStringList list;
        if (skipEmptyParts) {
            list = lineString.split(separator.c_str(), QString::SkipEmptyParts);
        }else{
            list = lineString.split(separator.c_str());
        }
        
        if (list.size()!=NumberOfColumn){
            std::cout << "Error extracting NumberTable: \nDo not find the wanted columns of entries in line " << line << "." << std::endl;
            std::cout << "Find:  " << list.size() << "instead of " << NumberOfColumn << "!" << std::endl;
            
            std::cout << "Line content:  " << lineString.toStdString() << std::endl;
            fNumberTableAvailable = false;
            return false;
        }
        
        std::vector<double> entry;
        bool conversionWorked = true;
        bool noError = true;
        for (int n = 0; n<NumberOfColumn; n++) {
            //std::cout << list.at(n).toStdString() << std::endl;
            entry.push_back(list.at(n).toDouble(&conversionWorked));
            noError = noError || !conversionWorked;
        }
        if (!noError){
            std::cout << "Error extracting NumberTable: \nConversion to double did not work in line " << line << "." << std::endl;
            fNumberTableAvailable = false;
            return false;
        }
        fNumberTable.push_back(entry);
    }
    fNumberTableRows = fNumberTable.size();
    fNumberTableColumns = NumberOfColumn;
    fNumberTableAvailable = true;
    return true;
};

bool CRACS::FileTool::WriteFile(std::string filename, std::string path){
    std::string f = "";
    if (path=="SameAsOriginal") {
        f += fFilePath;
    }else{
        f += path;
    }
    f+=filename;
    //std::cout << f << std::endl;
    if (f==(fFilePath+fFileName)) {
        std::cout << "Error writing file: If you want to overwrite file use OverwriteFile." << std::endl;
        return false;
    }
    
    
    QFile file(f.c_str());
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        std::cout << "Error writing file: Can't open file " << f << "."<< std::endl;
        return false;
    }
    QTextStream out(&file);
    for (int i=0; i<fModified.size(); i++) {
        out << fModified.at(i) << "\n";
    }
    
    file.close();
    return true;
};


bool CRACS::FileTool::WriteNumberTable(std::string filename, std::string path){
    std::string f = "";
    if (path=="SameAsOriginal") {
        f += fFilePath;
    }else{
        f += path;
    }
    f+=filename;
    //std::cout << f << std::endl;
    if (f==(fFilePath+fFileName)) {
        std::cout << "Error writing file: If you want to overwrite file use OverwriteFile." << std::endl;
        return false;
    }
    
    
    QFile file(f.c_str());
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        std::cout << "Error writing file: Can't open file " << f << "."<< std::endl;
        return false;
    }
    QTextStream out(&file);
    for (int i=0; i<fNumberTableRows; i++) {
        for (int j=0; j<fNumberTableColumns; j++) {
            out << QString("%1").arg(NumberTable(i, j)).rightJustified(30);
        }
        out << "\n";
    }
    
    file.close();
    return true;
};


bool CRACS::FileTool::OverwriteFile(){
    
    QFile file( (fFilePath+fFileName).c_str() );
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        std::cout << "Error writing file: Can't open file " << fFilePath+fFileName << "."<< std::endl;
        return false;
    }
    QTextStream out(&file);
    for (int i=0; i<fModified.size(); i++) {
        out << fModified.at(i) << "\n";
    }
    
    file.close();
    return true;
};


double CRACS::FileTool::NumberTable(int row, int column){
    if (!fNumberTableAvailable) {
        std::cout << "Error returning Value: NumberTable not available -> return -1" << std::endl;
        return -1;
    }
    if (column<0 || column>= fNumberTableColumns) {
        std::cout << "Error returning Value: Cannot acces your column -> return -1" << std::endl;
        return -1;
    }
    if (row<0 || row>= fNumberTableRows) {
        std::cout << "Error returning Value: Cannot acces your row -> return -1" << std::endl;
        return -1;
    }
    return fNumberTable.at(row).at(column);
};

void CRACS::FileTool::SetNumberTable(double val, int row, int column){
    if (!fNumberTableAvailable) {
        std::cout << "Error setting Value: NumberTable not available" << std::endl;
        return;
    }
    if (column<0 || column>= fNumberTableColumns) {
        std::cout << "Error setting Value: Cannot acces your column" << std::endl;
        return;
    }
    if (row<0 || row>= fNumberTableRows) {
        std::cout << "Error setting Value: Cannot acces your row" << std::endl;
        return;
    }
    fNumberTable.at(row).at(column) = val;
};

//TGraphErrors* FileReader::GetErrorGraph(std::string name, std::string title, int column_x, int column_y, int column_xe, int column_ye, double xe, double ye){
//    
//    if (column_x<-1 || column_y<-1 || column_xe<-1 || column_ye<-1) {
//        std::cout << "Error returning Graph: Cannot acces your column -> return nullptr" << std::endl;
//        return static_cast<TGraphErrors*>(nullptr);
//    }
//    if (column_x>=fColumns || column_y>=fColumns || column_xe>=fColumns || column_ye>=fColumns) {
//        std::cout << "Error returning Graph: Cannot acces your column -> return nullptr" << std::endl;
//        return static_cast<TGraphErrors*>(nullptr);
//    }
//    
//    TGraphErrors* graph = new TGraphErrors();
//    graph->SetNameTitle(name.c_str(), title.c_str());
//    for (unsigned int n = 0; n<fData.size(); n++) {
//        graph->SetPoint(n, (column_x==-1)? n : fData.at(n).at(column_x), (column_y==-1)? n : fData.at(n).at(column_y));
//        graph->SetPointError(n, (column_xe==-1)? xe : fData.at(n).at(column_xe), (column_ye==-1)? ye : fData.at(n).at(column_ye) );
//    }
//    
//    return graph;
//};

std::vector<double> CRACS::FileTool::NumberTableColumn(int column){
    if (!fNumberTableAvailable) {
        std::cout << "Error returning Column: NumberTable not available -> return empty vector" << std::endl;
        std::vector<double> empty;
        return empty;
    }
    if (column<0 || column>= fNumberTableColumns) {
        std::cout << "Error returning Column: Cannot acces your column -> return empty vector" << std::endl;
        std::vector<double> empty;
        return empty;
    }
    
    std::vector<double> columnVec;
    columnVec.reserve(fNumberTableRows);
    for (unsigned int n = 0; n<fNumberTableRows; n++) {
        columnVec.push_back(fNumberTable.at(n).at(column));
    }
    
    return columnVec;
};

std::vector<double> CRACS::FileTool::NumberTableRow(int row){
    if (!fNumberTableAvailable) {
        std::cout << "Error returning Row: NumberTable not available -> return empty vector" << std::endl;
        std::vector<double> empty;
        return empty;
    }
    if (row<0 || row>= fNumberTableRows) {
        std::cout << "Error returning Row: Cannot acces your row -> return empty vector" << std::endl;
        std::vector<double> empty;
        return empty;
    }
    
    return fNumberTable.at(row);
};





void CRACS::FileTool::PrintModified(){
    for (int i = 0; i<fModified.size(); i++) {
        std::cout << fModified.at(i).toStdString() << std::endl;
    }
};

void CRACS::FileTool::PrintOriginal(){
    for (int i = 0; i<fOriginal.size(); i++) {
        std::cout << fOriginal.at(i).toStdString() << std::endl;
    }
};

void CRACS::FileTool::PrintNumberTable(){
    std::cout << "Rows: " << fNumberTableRows << std::endl;
    std::cout << "Cols: " << fNumberTableColumns << std::endl;
    for (int i = 0; i<fNumberTableRows; i++) {
        for (int j = 0; j<fNumberTableColumns; j++) {
            //std::cout << "(" << i << "," << j << ")" << std::endl;
            std::cout << fNumberTable.at(i).at(j) << "\t";
        }
        std::cout << std::endl;
    }
};



bool CRACS::FileTool::WriteStringToFile(std::string string, std::string filename, std::string path){
    std::string f = "";
    if (path!="") {
        f += path;
    }
    f+=filename;
    QFile file( (f).c_str() );
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        std::cout << "Error writing file: Can't open file " << f << "."<< std::endl;
        return false;
    }
    QTextStream out(&file);
    out << string.c_str();
    
    file.close();
    return true;
};



int CRACS::FileTool::GetNumberOfLines(){
    return fModified.size();
}

std::string CRACS::FileTool::GetLine(int line){
    return fModified.at(line).toStdString();
}


void CRACS::FileTool::ReplaceLine(std::string replace, int line){
    fModified.at(line) = QString(replace.c_str());
}

void CRACS::FileTool::AddLine(std::string add, int line){
    if (line==-1) {
        fModified.push_back(QString(add.c_str()));
        return;
    }
    std::vector<QString> helper = fModified;
    fModified.clear();
    fModified.reserve(helper.size()+1);
    for (int i = 0; i<helper.size(); i++) {
        if (i==line) {
            fModified.push_back( QString(add.c_str()) );
        }
        fModified.push_back(helper.at(i));
    }
}
































