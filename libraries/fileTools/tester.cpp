#include "iostream"
#include "fileTools.h"

int main(int argc, char *argv[])
{
    
    CRACS::FileTool fileTool("/Users/michaelkorsmeier/MasterThesisSoftware/libraries/fileReader/test.txt");
    
    
    fileTool.PrintOriginal();
    
    fileTool.FindAndReplace("Test", "Auto");
    
    fileTool.PrintModified();
    
    fileTool.ExtractNumberTable(3, "\t", false);
    std::cout << "Number Table: " << std::endl;
    fileTool.PrintNumberTable();
    
    fileTool.WriteFile("test2");
    fileTool.OverwriteFile();
    
    
}
