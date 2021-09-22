#include "iostream"
#include "configHandler.h"

int main(int argc, char *argv[])
{
    
    
    ConfigHandler config(argc ,argv, "Programm Description");
    std::string option1 = "d1";
    config.AddOptionPosition(option1, "This is Option1.", true);
    std::string option2 = "d2";
    config.AddOptionPosition(option2, "This is Option2.", true);
    std::string file = "df";
    config.AddOption("file",file, "This is Option file." );
    config.CheckOption();
    
    std::cout << option1 << std::endl;
    std::cout << option2 << std::endl;
    std::cout << file << std::endl;
    return -1;
}