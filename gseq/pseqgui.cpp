
#include "cgi.h"
#include "helper.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

void report_error( const std::string & errmsg );
void finished();
std::string streampseq( const std::string & syscall , std::string (*f)(const std::string &, void * p ) = NULL , void * p = NULL  );


//output functions
std::string f_tabulate( const std::string & res , void * p = NULL );
std::string f_pre( const std::string & res , void * p = NULL );

typedef std::string (*ofunc_t)( const std::string & , void * );

int main()
{
  
  const std::string pseq = "/psych/genetics/pseq/dist/pseq ";
  
  //
  // Get CGI variables from POST
  //
  
  char **cgivars ;  int i ;  
  cgivars = getcgivars() ;


  //
  // Initiate HTML page
  //

  std::cout << "Content-type: text/html\n\n"
	    << "<html><head><title>GSEQ</title>"
	    << "</head><body><font face=\"arial\">";
  
  
  //
  // Set up different PSEQ functions to have different browser outut modes
  //

  std::map<std::string,ofunc_t> omode;
  
  omode["version"] = f_pre;
  omode["v-freq"] = f_tabulate;

  
  //
  // Core values were are expecting from the previous page
  //
  
  std::string project = ".";
  std::string pcomm   = "";
  std::string cmdline = "";
  std::string action  = "list-commands";
  std::string value   = "";
  std::string value2  = "";
  std::string currdir = ""; 

  for (i=0; cgivars[i]; i+= 2)
    {      
      std::string str = cgivars[i];      
      if      ( str == "project" ) project = cgivars[i+1];
      else if ( str == "cmdline" ) cmdline = cgivars[i+1];
      else if ( str == "action"  ) action  = cgivars[i+1];
      else if ( str == "value"   ) value   = cgivars[i+1];
      else if ( str == "value2"  ) value2  = cgivars[i+1];
      else if ( str == "pcomm"   ) pcomm   = cgivars[i+1];
    }
  

  //
  // Free all CGI variables
  //
  
  for (i=0; cgivars[i]; i++) free(cgivars[i]) ;
  free(cgivars) ;
  

  //
  // File browser
  //
  
  if ( action == "file-browser" ) 
    {
      //      std::
    }

  //
  // Attach a new project 
  //

  if ( action == "attach-project" )
    {
      project = value;
    }


  //
  // Check project exists
  //
  
  bool using_project = project != ".";
  
  if ( ! fileExists( project ) ) 
    report_error( "project file not found" );
  
  
  
    
  //
  // Show .history (or get job number)
  //
  
  int job = 0;

  if ( action == "show-history" || action == "run-pseq" )
    {
      
      std::map<std::string,std::string> history;
      std::map<std::string,std::string> comment;
      std::map<std::string,std::string> ofile;
      std::map<std::string,std::string> status;
      
      if ( using_project && fileExists( project + ".history" ) )
	{
	  std::ifstream IN1( (project + ".history").c_str() , std::ios::in );
	  while ( ! IN1.eof() )
	    {
	      std::string l;
	      getline(IN1,l);
	      if ( l == "" ) continue;	      
	      
	      std::vector<std::string> tok = char_split( l , '\t' );
	      if ( tok.size() == 5 && tok[0] == "_COMM" ) 
		{
		  history[ tok[1] ] = tok[2];
		  comment[ tok[1] ] = tok[3];
		  ofile[ tok[1] ] = tok[4];
		  int j = atoi( tok[1].c_str() );
		  if ( j > job ) job = j;
		}
	      if ( tok.size() >= 3 && tok[0] == "_STATUS" ) status[ tok[1] ] = tok[2];
	    }
	  IN1.close();
	}
      else if ( action == "show-history" )
	report_error( "could not find project history" );

      // advance to next job

      ++job;

      //
      // Display history
      //
      
      if ( action == "show-history" )
	{
	  std::cout << "<table border=1 width=100%><tr><th>Job</th><th>Status</th><th>Comment</th><th>Output file</th><th>Command</th></tr>";
	  std::map<std::string,std::string>::iterator h = history.begin();
	  while ( h != history.end() )
	    {
	      std::cout << "<tr><td>" << h->first << "</td>";
	      if ( status.find( h->first ) == status.end() ) std::cout << "<td><em><font color=\"gray\">in flight</font></em></td>";
	      else std::cout << "<td>" << status[ h->first ] << "</td>";
	      std::cout << "<td>" << comment[ h->first ] << "</td>";
	      std::cout << "<td>" << ofile[ h->first ] << "</td>";	      
	      std::cout << "<td>" << h->second << "</td></tr>";
	      ++h;
	    }
	  std::cout << "</table>";
	  
	  finished();
	}
    }

   

  //
  // What is current action?
  //
  
  // known commands:

  //  list-commands 
  //  run-pseq 
  //  select-mask
  //  add-cmd  
  
  //
  // Assume we've been given something to add to the cmd line 
  //

  if ( action == "add-cmd" )
    {
      pcomm = value;
    }

  if ( action == "add-mask" || action == "add-opt" )
    {      
      cmdline += " " + value;
      if ( value2 != "" ) cmdline += "=" + value2 + " ";
    }
  
  if ( action == "add-arg" )
    {
      cmdline += " " + value;
      if ( value2 != "" ) cmdline += " " + value2 + " ";      
    }

  
  //
  // Show current command
  //
  
  // extract mask options, arguments and options from the current command -line
  // (merge additional spaces, allow quotes )

  std::vector<std::string> spl = quoted_char_split( cmdline , ' ' , true );  
  
  bool mask = false;
  bool arg = false;
  bool opt = false;
  
  std::string mask_str = "";
  std::string arg_str = "";
  std::string opt_str = "";

  for (int i=0; i<spl.size(); i++)
    {
      if ( spl[i] == "." ) continue;

      if ( spl[i] == "--options" ) { opt = true; continue; }
      else if ( spl[i] == "--mask" ) { mask = true; continue; }
      else if ( spl[i].substr(0,2) == "--" ) { arg = true; mask = opt = false; } 
      
      if ( mask ) mask_str += " " + spl[i];
      else if ( opt ) opt_str += " " + spl[i];
      else arg_str += " " + spl[i];

    }
  

  //
  // Command to clear mask element
  //

  if ( action == "clear-mask" ) 
    {
      mask_str = "";
    }


  //
  // Reformulate command line: mask, args, options
  //
  
  cmdline = "";
  if ( mask_str != "" ) cmdline += " --mask" + mask_str;
  if ( arg_str != "" ) cmdline += " " + arg_str;
  if ( opt_str != "" ) cmdline += " " + opt_str;
  



  //
  // Run a completed PSEQ command ? 
  //
  
  if ( action == "run-pseq" )
    {

      std::string syscall = pseq + " " + project + " " + pcomm + " " + cmdline;


      //
      // Are we tracking job completeion in a history file?
      //
            
      if ( using_project ) 
	syscall = pseq + " " + project + " " + pcomm + " --history " + project + ".history " + int2str(job) + " " + cmdline;

      
      //
      // Write output to a file
      //

      if ( value != "" ) 
	{
	  syscall += " > " + value;
	  std::cout << "<em>redirecting output to [ " << value << " ]<br> ";
	  streampseq( syscall ); // do not capture output
	}


      //
      // ...or otherwise, display in browser
      //
      
      else 
	{
	  // (for now, do not allow any auxiliary object to be passed to the output
	  // function, but that can be added later: i.e. form is:
	  //  streampseq( syscall , f_tabulate , &aux_data );
	  
	  std::map<std::string,ofunc_t>::iterator f = omode.find( pcomm );	  
	  std::cout << streampseq( syscall , f == omode.end() ? f_pre : f->second );     	  
	  
	}

      // add to project history

      if ( using_project ) 
	{
	  std::ofstream O1( ( project + ".history").c_str() , std::ios_base::out | std::ios::app );
	  O1 << "_COMM\t" 
	     << job << "\t"     // job counter
	     << pcomm << " " << cmdline << "\t" // cmdline
	     << value2 << "\t"  // comment
	     << value << "\n";  // output file
	  
	  O1.close();
	} 

      finished();
    }


  //
  // Select a mask item and add to the cmdline 
  //

  if ( action == "select-mask" ) 
    {
      std::string syscall = pseq + " " + project + " " + "masks";
      std::string r = streampseq( syscall );
      
      std::string curr_group = "";
      std::vector<std::string> l = char_split( r , '\n' );
      for (int i=0; i<l.size(); i++)
	{
	  // expecting 4 tab-delimt cols
	  std::vector<std::string> row = char_split( l[i] , '\t' );
	  if ( row.size() != 4 ) continue;
	  if ( curr_group != row[0] ) 
	    {
	      if ( curr_group != "" ) std::cout << "</table>";
	      std::cout << "<h3>" << row[0] << "</h3>";
	      curr_group = row[0];
	      std::cout << "<table border=0 width=\"70%\">";
	    }
	  
	  
	  std::cout << "<tr><td valign=\"top\" width=\"15%\"><tt>" << row[1] << "</tt></td>";

	  std::cout << "<td align=\"top\" valign=\"center\" width=\"25%\">" 
		    << "<form action=\"./gseq.cgi\">"
		    << "<input type=\"hidden\" name=\"action\" value=\"add-mask\">"
		    << "<input type=\"hidden\" name=\"project\" value=\"" << project << "\">"
		    << "<input type=\"hidden\" name=\"cmdline\" value=\"" << cmdline << "\">"
		    << "<input type=\"hidden\" name=\"pcomm\" value=\"" << pcomm << "\">"
		    << "<input type=\"hidden\" name=\"value\" value=\"--mask " << row[1] << "\">";
	  if ( row[2] != "flag" ) 
	    std::cout << "<input name=\"value2\" size=\"25\">";
	  std::cout << " <input type=\"submit\" value=\"Set\">"
		    << "</form>";	  

	  std::cout << "<td valign=\"top\" width=\"60%\">" << row[3] << "</td>";

	  std::cout << "</tr>";
	  
	}
      std::cout << "</table>";
      finished();
    }

  
  //
  // Primary table
  //


  std::cout << "<table width=100%><tr><td align=left><b><em>GSEQ command table</em></b></td><td align=right>";
  
  if ( project != "." ) 
    {
      std::cout << " Project: <em>" << project << "</em>";
    }
  else
    {
      std::cout << "<em>no project attached</em>";
    }
  std::cout << "</td></tr></table><hr>";
  
  

  std::cout << "<table width=100% border=2 cellspacing=0 cellpadding=1 bordercolor=#aaaaaa>";
  std::cout << "<tr><td width=\"40%\" bgcolor=#EEEEFF  valign=\"top\">";
  
  std::cout << "<p><a href=\"gseq.cgi?project=" << project 
	    << "&cmdline=" << cmdline 
	    << "&action=list-commands&value=root\">back</a> &nbsp;&nbsp;&nbsp; (to initial command-table)</p>";
  
  //
  // Either we've selected a command already, or still choosing:
  //

  
  
  if ( action == "list-commands" )
    {      
      
      // default value under list-commands:
      if ( value == "" ) value = "root";

      std::string syscall = "../client/pseq-dev " + project + " commands --name " + value ;
      
      std::vector<std::string> resl = char_split( streampseq( syscall ) , '\n' );
      
      if ( resl.size() > 0 ) 
	{

	  std::cout << "<table border=0><tr>";


	  for (int i=0; i<resl.size(); i++)
	    {
	      std::vector<std::string> k = char_split( resl[i] , '\t' );
	      if ( k.size() == 3 ) 
		{
		  
		  if ( k[0] == "GROUP" ) 
		    {
		      std::cout << "<td width=30%><a href=\"gseq.cgi?project=" << project << "&action=list-commands"
				<< "&value=" << k[1] 
				<< "&cmdline=" << cmdline << "\">"
				<< k[1] << "</a></td><td> " << k[2] << "</td></tr>";
		    }
		  else if ( k[0] == "COMM" ) 
		    {
		      std::cout << "<tr width=30%><td><a href=\"gseq.cgi?project=" << project << "&action=add-cmd&value=" 
				<< k[1] << "&cmdline="<<cmdline<<"\">" << k[1] << "</a></td><td> " << k[2] << "</td></tr>";
		    }
		  
		}
	    }
	  
	  std::cout << "</table>";
	} 
    }


  //
  // ...or we've already selected a command; give verbose description and run button
  //
  
  if ( pcomm != "" ) 
    {
      
      std::cout << "<p><b>Command: <font color=\"darkblue\">" << pcomm << "</font></b></p>"; 
      
      std::cout << "<form action=\"./gseq.cgi\">"
		<< "<input type=\"submit\" value=\"Run PSEQ\">"
		<< "<input type=\"hidden\" name=\"project\" value=\"" << project << "\">"
		<< "<input type=\"hidden\" name=\"cmdline\" value=\"" << cmdline << "\">"
		<< "<input type=\"hidden\" name=\"pcomm\" value=\"" << pcomm << "\">"
		<< "<input name=\"value2\" size=20> Comment<br>"
		<< "<input name=\"value\" size=20> Output to file<br>"
		<< "<input type=\"hidden\" name=\"action\" value=\"run-pseq\">"
		<< "</form>";
      
    }


  //
  // Right panel will now list arg/opt choices
  //
  
  std::cout << "</td><td bgcolor=#EEEEFF valign=\"top\" width=\"70%\" >";
  
  // Get accepted arguments and display
  
  std::string syscall = pseq + " " + project + " commands --name " + pcomm ;    

  std::vector<std::string> resl = char_split( streampseq( syscall ) , '\n' );
  
  //
  // Get accepted options and display
  //
  
  std::cout << "<h4>Mask "
	    << "(<a href=\"gseq.cgi?project=" << project 
	    << "&action=select-mask" 
	    << "&cmdline=" << cmdline 
	    << "&pcomm=" << pcomm 
	    << "\">add</a>)</h4>";
  
  std::cout << mask_str << "<br>";

  
  //
  // List chosen/available arguments, if we've selected a command
  //
  
  if ( pcomm != "" )
    {

      std::cout << "<h4>Arguments</h4>" << arg_str << "<br>";
      
      for (int i=0; i<resl.size(); i++)
	{
	  std::vector<std::string> d = char_split( resl[i] , '\t' );
	  if ( d.size() != 4 ) continue;
	  if ( d[0] != "ARG" ) continue;
	  
	  std::cout << "<form action=\"./gseq.cgi\">" << d[2] << " = "		    
		    << "<input type=\"hidden\" name=\"project\" value=\"" << project << "\">"
		    << "<input type=\"hidden\" name=\"cmdline\" value=\"" << cmdline << "\">"
		    << "<input type=\"hidden\" name=\"pcomm\" value=\"" << pcomm << "\">"		    
		    << "<input type=\"hidden\" name=\"action\" value=\"add-arg\">";
	  if ( d[3] == "flag" ) 
	    std::cout << "<input type=\"hidden\" name=\"value\" value=\"\">";
	  else
	    std::cout << "<input name=\"value\" value=\"\">";
	  std::cout << "<input type=\"submit\" value=\"add\">"
		    << "</form><br>";
	}
      
      //
      // List chosen / available options 
      //
      
      std::cout << "<h4>Options</h4> " << opt_str << "<br>";
      
      for (int i=0; i<resl.size(); i++)
	{
	  std::vector<std::string> d = char_split( resl[i] , '\t' );
	  if ( d.size() != 4 ) continue;
	  if ( d[0] != "OPT" ) continue;	  
	  std::cout << " opt = " << resl[i] << "<br>";
	}
      
    }
  
  // Close the righthand-side panel and whole table
  std::cout << "</td></tr></table><hr>";
  

  //
  // Direct command line input
  //
  
  // Reset button

  std::cout << "<table border=0 width=100%><tr>";

  std::cout << "<td width=\"20%\"><form action=\"./gseq.cgi\">"
	    << " <input type=\"submit\" value=\"Reset\">"
	    << "<input type=\"hidden\" name=\"project\" value=\"" << project << "\">"
	    << "<input type=\"hidden\" name=\"cmdline\" value=\"\">"
	    << "<input type=\"hidden\" name=\"pcomm\" value=\"\">"
	    << " <input type=\"hidden\" name=\"action\" value=\"list-commands\">"
	    << " <input type=\"hidden\" name=\"value\" value=\"root\">"
	    << "</form></td>";

  // Attach a new project

  std::cout << "<td width=\"20%\"><form action=\"./gseq.cgi\">"
	    << "<input type=\"submit\" value=\"Attach new project\"><br>"
	    << "<input type=\"hidden\" name=\"action\" value=\"attach-project\">"
	    << "<input name=\"value\" size=30 value=\"\">"
	    << "</form></td>";
  
  // New command (but keep arguments)

  std::cout << "<td width=\"20%\"><form action=\"./gseq.cgi\">"
	    << " <input type=\"submit\" value=\"New command\">"
	    << "<input type=\"hidden\" name=\"project\" value=\"" << project << "\">"
	    << "<input type=\"hidden\" name=\"cmdline\" value=\"" << cmdline << "\">"
	    << " <input type=\"hidden\" name=\"action\" value=\"list-commands\">"
	    << " <input type=\"hidden\" name=\"value\" value=\"root\">"
	    << "</form></td>";


  // Add a mask element

  std::cout << "<td width=\"20%\"><form action=\"./gseq.cgi\">"
	    << " <input type=\"submit\" value=\"Clear Mask\">"
	    << "<input type=\"hidden\" name=\"project\" value=\"" << project << "\">"
	    << "<input type=\"hidden\" name=\"cmdline\" value=\"" << cmdline << "\">"
	    << "<input type=\"hidden\" name=\"pcomm\" value=\"" << pcomm << "\">"
	    << " <input type=\"hidden\" name=\"action\" value=\"clear-mask\">"
	    << "</form></td>";
  
  // Show history

  std::cout << "<td  width=\"20%\"><form action=\"./gseq.cgi\">"
	    << " <input type=\"submit\" value=\"History\">"
	    << "<input type=\"hidden\" name=\"project\" value=\"" << project << "\">"
	    << "<input type=\"hidden\" name=\"cmdline\" value=\"" << cmdline << "\">"
	    << "<input type=\"hidden\" name=\"pcomm\" value=\"" << pcomm << "\">"
	    << " <input type=\"hidden\" name=\"action\" value=\"show-history\">"
	    << " <input type=\"hidden\" name=\"value\" value=\"select-mask\">"
	    << "</form></td>";
  
  std::cout << "</tr></table><hr>";

  //
  // Generate current form/representation of command line
  //
  
  
  finished();
}


void report_error( const std::string & errmsg )
{
  std::cout << "<p>" << errmsg << "</p>";
  std::cout << "</body></html>";
  exit(0);
}

void finished()
{
  std::cout << "</body></html>";  
  exit(0);  
}


//
// Function to run PSEQ via a system call, and collect the output in some 
// user-defined way, e.g. to file, to HTML page, tabularize or annotate
// certain output types, etc
//

std::string streampseq( const std::string & syscall , std::string (*f)(const std::string & , void * ) , void * aux )
{
  FILE * pseqstream;
  char buff[ 512 ];
  if ( ! ( pseqstream = popen( syscall.c_str() , "r" ) ) )
    report_error( "unable to popen( pseq )" );  
  std::string res = "";
  while ( fgets( buff , sizeof(buff) , pseqstream ) != NULL )
    {
      res += buff;
      // TODO: add check to stop process if output exceeds some threshol
    }

  pclose( pseqstream );   

  return f ? f(res,aux) : res;
}

std::string f_pre( const std::string & res , void * aux )
{
  return "<pre>" + res + "</pre>";
}

std::string f_tabulate( const std::string & res , void * hdr )
{

  std::string hstr1 = hdr && *(bool*)hdr ? "<th>" : "<td>";
  std::string hstr2 = hdr && *(bool*)hdr ? "</th>" : "</td>";

  // assume equal # of cols 
  int col = char_split( res.substr(0,res.find("\n") ) , '\t' ).size();
  if ( col == 0 ) return "";

  // header 
  std::vector<std::string> lines = char_split( res , '\n' );

  std::cout << "<table border=1><tr>";  
  for (int r=0;r<lines.size(); r++)
    {
      std::vector<std::string> row = char_split( lines[r] , '\t' );         
      if ( row.size() == col ) 
	{
	  std::cout << "<tr>";
	  for (int c=0; c<col; c++)
	    std::cout << ( r==0 ? hstr1 : "<td>" ) << row[c] << ( r==0 ? hstr2 : "</td>") ;      
	  std::cout << "</tr>";
	}
    }        
  std::cout << "</table>";
  return "";
}
