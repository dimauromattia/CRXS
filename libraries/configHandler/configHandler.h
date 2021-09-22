#ifndef CONFIGHANDLER_H
#define CONFIGHANDLER_H

#include <string>
#include <vector>
#include <functional>
#include <iostream>


namespace CRACS {
    
    
    
    /*! \brief ConfigHander Class. This singlton class parses options handed over to a program.
     *
     *  The ConfigHandler recognizes opions as follows:
     *
     *  exampleProgramm VALUE_POSTION_OPTION1 VALUE_POSTION_OPTION2 --option1 VALUE1 --option2 VALUE2
     *
     *  Example
     *  //FIXME to be updated
     *  \include configHandler.cpp
     */

    class ConfigHandler{
        
    public:

        static ConfigHandler* GetInstance();

        //! Set input.
        /*!
         *  \param int argc. Number of arguments from argv.
         *  \param char *argv[]. Programm name and options.
         *  \param std::string descripion. Short description of the programm. Displayed on help page (-h, --h).
         */
        void SetInput(int argc, char *argv[], std::string description);

        //! Adds positional argument.
        /*!
         *  \param std::string& value. Value of position argument.
         *  \param bool mandatory. If ture the ConfigHandler insists on this Option.
         *  \param std::string descripion. Short description of the option. Displayed on help page.
         */
        void AddOptionPosition(std::string& value, std::string description="", bool mandatory=false);
        
        //! Adds optinal argument: "--option VALUE". VALUE is a std::string.
        /*!
         *  \param std::string option. Name of the Option.
         *  \param std::string& value. Value of the Option as string.
         *  \param std::string descripion. Short description of the option. Displayed on help page.
         */
        void AddOption(std::string option, std::string &value, std::string description="");
        
        //! Adds optinal argument: "--option VALUE".VALUE is an int.
        /*!
         *  \param std::string option. Name of the Option.
         *  \param int& value. Value of the Option as int.
         *  \param std::string descripion. Short description of the option. Displayed on help page.
         */
        void AddOptionInt(std::string option, int &value, std::string description="");
        
        //! Adds optinal argument: "--option VALUE". VALUE is a double.
        /*!
         *  \param std::string option. Name of the Option.
         *  \param double& value. Value of the Option as double.
         *  \param std::string descripion. Short description of the option. Displayed on help page.
         */
        void AddOptionDouble(std::string option, double &value, std::string description="");
        
        //! Adds optinal argument: "--option". True if stored
        /*!
         *  \param std::string option. Name of the Option.
         *  \param bool& stored. Ture if this option is used.
         *  \param std::string descripion. Short description of the option. Displayed on help page.
         */
        void AddOptionTrue(std::string option, bool &stored, std::string description="");
        
        //! Checks all options. Run after adding.
        /*!
         */
        void CheckOption();
           
        std::string SoftwarePath()         { return fSoftwarePath;          };       
        
    private:
        static ConfigHandler* fInstance;
        std::string fSoftwarePath;

        //! Constructor. Private because this is a singelton.
        ConfigHandler(){};
        
        //! Function parsing the options. You need to execute this before you can access arguments.
        /*!
         */
        void ParseOption();
        
        
        std::string fOptionString;                              /// String with all options.
        std::string fName;                                      /// String with name of programm.
        std::string fDescription;                               /// Programm description.
        std::string fError;                                     /// String with error message.
        
        std::vector<std::string>    fPositionOptionsDescription;    /// Positional arguments description.
        std::vector<bool>           fPositionOptionsMandotory;      /// Positional option mandatory.
        std::vector<std::string>    fPositionOptionsValue;          /// Positional Value.
        int                         fPositionOptions;               /// Current Positional argument number.
        
        std::vector<std::string> fOptions;                      /// Option of optinal arguments.
        std::vector<std::string> fOptionValues;                 /// Value of optinal arguments.
        
        std::vector<std::string> fOptionsExpected;              /// Expected optional arguments.
        std::vector<std::string> fOptionsDescription;           /// Optional arguments description.
        
        
        //! Private function printing help text.
        /*!
         */
        void PrintHelp();
        
        bool fOk;
        
        
    };
    

}



//end CONFIGHANDLER_H

#endif
