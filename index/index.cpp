#include "index.h"
#include "cgi.h"

#include <iostream>
#include <algorithm>

//using namespace IndivExome;

int main()
{
  
  
  //
  // Get CGI variables from POST
  //
  
  char **cgivars = getcgivars() ;
  
  
  //
  // Core parameters
  //

  // Project path
  std::string project_path = "";

  // Password
  std::string pwd = "(if required)";
  
  // Use defaults Default locus-set ("refseq", "symbol")

  std::string loc_set = PLINKSeq::DEFAULT_LOC_GROUP();
  std::string gene_symbol = PLINKSeq::DEFAULT_GENE_SYMBOL();


  // Filters
  
  std::set<std::string> var_inc;
  std::set<std::string> var_ex;
  std::set<std::string> var_req;
  
  std::set<std::string> varset_inc;
  std::set<std::string> varset_ex;
  std::set<std::string> varset_req;
  
  std::set<std::string> loc_req;
  std::set<std::string> reg_req;
  

  // Start HTML

  std::cout << "Content-type: text/html\n\n"
	    << "<html><head><title>Individual Exome Viewer (IndEx)</title>"
	    << "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />"
	    << "<title>IndEx</title>";
  
  // CSS
  
  std::cout << "<style media=\"screen\" type=\"text/css\">"
	    << "table { border-collapse: collapse; border: 1px solid #666666; font: normal 11px verdana, arial, helvetica, sans-serif;"
	    << "color: #363636; background: #f6f6f6; text-align:left;}"
	    << "caption {"
	    << "text-align: center;"
	    << "font: bold 16px arial, helvetica, sans-serif;"
	    << "background: transparent;"
	    << "padding:6px 4px 8px 0px;"
	    << "color: #CC00FF;"
	    << "   text-transform: uppercase;"
	    << "}"
	    << "thead, tfoot {"
	    << "  text-align:left;"
	    << "height:30px;"
	    << "}"
	    << "thead th, tfoot th {"
	    << "padding:5px;"
	    << "}"
	    << "table a {"
	    << "color: #333388;"
	    << "  text-decoration:underline;"
	    << "}"
	    << "table a:hover {"
	    << "  text-decoration:underline;"
	    << "}"
	    << "tr.odd {"
	    << "background: #f171f1;"
	    << "}"
	    << "tbody th, tbody td {"
	    << "padding:5px;"
	    << "}"
	    << "</style>";
  
            std::cout << "</head><body>";
	    
	   
  /// -------------- END OF HEAD ----------------

  
  exit(0);
}

