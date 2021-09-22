#ifndef FILETOOL_H
#define FILETOOL_H

#include <iostream>
#include <string>
#include <vector>

#include "QFile"
#include "QTextStream"
#include <QString>
#include <QStringList>



namespace CRACS {
    
    /*! \brief File Tool Class to read/write file. Extract tables, ...
     *
     */
    class FileTool {
    private:
        
        std::string fFileName;                  /// Name of the file for reading.
        std::string fFilePath;                  /// Path of the file for reading.
        std::vector<QString> fOriginal;         /// Vector of the original file content by line.
        std::vector<QString> fModified;         /// Vector of the Modified file content by line.
        
        bool fNumberTableAvailable;                         /// Bool if table with numbers was read in.
        std::vector< std::vector<double> > fNumberTable;    /// 2D Vector with table data. Order: (line, column).
        int fNumberTableColumns;                            /// Number of columns in number data table.
        int fNumberTableRows;                               /// Number of rows in number data table
        
        bool fReadingOK;
        
        //! Reads the file. Called by constructor
        /*!
         *  \return bool. Ture if reading was succesfull
         */
        bool ReadFile();
        
    public:
        //! Constructor.
        /*!
         *  \param std::string filename. (Full path) Filename of the file to read. Calls ReadFile().
         */
        FileTool(std::string filename);
        FileTool(){};
        ~FileTool();
        
        
        //! Static function to check whether a file exists.
        /*!
         *  \param std::string filename. (Full path) Filename of the file to check.
         */
        static bool Exists(std::string filename);
        
        //! Writes file (fModified)
        /*!
         *  \param std::string filename. Name of the file to write.
         *  \param std::string path. Path for the file to write. Default same as original file.
         *  \return bool. Ture if writing was succesfull.
         */
        bool WriteFile(std::string filename, std::string path="SameAsOriginal");
        //! Overwrites original file (with fModified)
        /*!
         *  \return bool. Ture if writing was succesfull.
         */
        bool OverwriteFile();
        
        bool WriteNumberTable(std::string filename, std::string path="SameAsOriginal");
        
        std::string GetFileName(){return fFileName;};
        
        
        //! FindAndReplace. Searches <find> and replaces with <replace> in fModified
        /*!
         *  \param std::string find. String to find.
         *  \param std::string replace. String to replace.
         *  \param bool caseSensitiv. Replace only case sensitiv. Default=true.
         */
        void FindAndReplace(std::string find, std::string replace, bool caseSensitiv=true);
        
        
        //! BackToOriginal restores the original stat. fModified is replaced by fOriginal. fNumberTable is removed;
        /*!
         */
        void BackToOriginal();
        
        //! Prints original file
        /*!
         */
        void PrintOriginal();
        
        //! Prints modified file
        /*!
         */
        void PrintModified();
        
        
        //! Returns number of lines.
        /*!
         */
        int GetNumberOfLines();
        //! Returns line as string.
        /*!
         *  \param int line.
         */
        std::string GetLine(int line);
        
        //! Replaces line.
        /*!
         *  \param std::string replace.
         *  \param int line. Is replaced.
         */
        void ReplaceLine(std::string replace, int line);
        //! Adds line.
        /*!
         *  \param std::string add.
         *  \param int line. Add add in this line
         */
        void AddLine(std::string add, int line);
        
        
        bool IsFileOk(){return fReadingOK;};
        
        
        //! ExtractNumberTable extracts a table from this file
        /*!
         *  You may specify a start and stop line
         *
         *  \param int NumberOfColumn. Number of columns of your table
         *  \param std::string separator. Seperator as string (e.g. "\t", ";", ...)
         *  \param bool skipEmptyParts. Use true if there seperator apears more than once between tow values (eg. value1;;value2). Default: false
         *  \param int lineStart. Start line of table, use -1 for not specified. Default: -1.
         *  \param int lineStop. Stop line of table, use -1 for not specified. Default: -1.
         *  \return bool. True if extraction successful.
         */
        bool ExtractNumberTable( int NumberOfColumn, std::string separator, bool skipEmptyParts=true, int lineStart = -1, int lineStop = -1 );
        
        // If error columns are -1 the constants xe and ye are used as error
        // If x or y column are -1 use 0,1,2,3,...
        // TGraphErrors* NumberTable_Graph(std::string name, std::string title, int column_x, int column_y, int column_xe=-1, int column_ye=-1, double xe = 0, double ye = 0);
        
        //! NumberTableColumn
        /*!
         *  Get Colum from number table.
         *
         *  \param int column. Column number (counting starts at 0).
         *  \return std::vector<double> . Column as double vector.
         */
        std::vector<double> NumberTableColumn(int column);
        
        //! NumberTableRow
        /*!
         *  Get Row from number table.
         *
         *  \param int row. Row number (counting starts at 0).
         *  \return std::vector<double> . Row as double vector.
         */
        std::vector<double> NumberTableRow(int row);
        
        //! NumberTable
        /*!
         *  Get Row from number table.
         *
         *  \param int row. Row number (counting starts at 0).
         *  \param int column. Column number (counting starts at 0).
         *  \return double . Value of position (row, column).
         */
        double NumberTable(int row, int column);
        
        void   SetNumberTable(double val, int row, int column);
        
        //! NumberTableGetNrows
        /*!
         *
         *  \return int . Number of Rows.
         */
        int NumberTableGetNrows(){return fNumberTableRows;};
        //! NumberTableGetNcolumns
        /*!
         *
         *  \return int . Number of Columns.
         */
        int NumberTableGetNcolumns(){return fNumberTableColumns;};
        
        
        //! Prints number table
        /*!
         */
        void PrintNumberTable();
        
        //! WriteStringToFile
        /*!
         *  Writes a std::string to file.
         *
         *  \param std::string string. String to be written.
         *  \param std::string filename. Filename.
         *  \param std::string path. Path. Default: "".
         *  \return bool . True if writing was successful.
         */
        static bool WriteStringToFile(std::string string, std::string filename, std::string path="");
    };

}

#endif // FILETOOL_H
