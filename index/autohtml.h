#ifndef __PSEQ_AUTOHTML_H__
#define __PSEQ_AUTOHTML_H__

#include <map>
#include <vector>
#include <string>

namespace HTMLHelper { 

  void header( const std::string & title ) 
  {     

    std::cout << "Content-type: text/html\n\n";
    
    std::cout << "<html lang=\"en\">";
    
    std::cout << "<head><title>" << title << "</title>"
	      << "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=iso-8859-1\" />";
  }
  
  void add_css1()
  {
    
    
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
    
  }



  void init_combocheckbox()
  {

    std::cout << "<style type=\"text/css\">"
	      << " .checklist { border: 1px solid #ccc; list-style: none; height: 10em; overflow: auto; width: 16em; } "
	      << " .checklist, .checklist li { margin-left: 0; padding: 0; } "
	      << " .checklist label { display: block; padding-left: 25px; text-indent: -25px; } "
	      << " .checklist label:hover, .checklist label.hover { background: #777; color: #fff; } "
	      << " * html .checklist label { height: 1%; } "
      
	      << " .cl1 { font-size: 0.9em; width: 90%; height: 15em; } "
	      << " .cl1 .alt { background: #f5f5f5; } "
	      << " .cl1 input { vertical-align: middle; } "
	      << " .cl1 label:hover, .cl1 label.hover { background: #ddd; color: #000; } "
	      << " </style>";
    
    // javascript

    std::cout << " <script type=\"text/javascript\">";
      
    // addLoadEvent: Add event handler to body when window loads

    std::cout << "function addLoadEvent(func) {"
	      << " var oldonload = window.onload; "
    
	      << " if (typeof window.onload != \"function\") { "
	      << " window.onload = func; "
	      << " } else { "
	      << " window.onload = function () { "
	      << " oldonload(); "
	      << " func(); } } } ";

    // Functions to run when window loads 

    std::cout << " addLoadEvent(function () { initChecklist(); }); ";
    
    // initChecklist: Add :hover functionality on labels for IE

    std::cout << " function initChecklist() {"
	      << " if (document.all && document.getElementById) { "
	      << " var lists = document.getElementsByTagName(\"ul\"); " 
	      << " for (i = 0; i < lists.length; i++) { "
	      << "    var theList = lists[i]; "      		    
	      << " if (theList.className.indexOf(\"checklist\") > -1) { "
	      << " var labels = theList.getElementsByTagName(\"label\"); "
      
	      << " for (var j = 0; j < labels.length; j++) { "
	      << " var theLabel = labels[j]; "
	      << " theLabel.onmouseover = function() { this.className += \" hover\"; }; "
	      << " theLabel.onmouseout = function() { this.className = this.className.replace(\" hover\", \"\"); }; "
	      << " } } } } } "
	      << "</script>";

  }

  void end_header() 
  {
    std::cout << "</head> <body> ";
  }

  void end_page( const std::string & msg = "" )
  {
    std::cout << "<b><em>" << msg << "</em></b></body></html>";
    exit(0);
  }


   void add_combocheckbox( const std::string & name , const std::map<std::string,std::string> & recs , const std::set<std::string> & sel ) 
   { 
  
     std::cout << "<ul class=\"checklist cl1\">"; 
  
     std::map<std::string,std::string>::const_iterator i = recs.begin(); 
  
     bool even = false; 

     while ( i != recs.end() ) 
       { 
	 
	 if ( sel.find( i->first ) != sel.end() )	
	   {
	     if ( even )  
	       { 
		 std::cout << "<li class=\"alt\"><label for=\""  
			   << i->first  
			   << "\"><input id=\""  
			   << i->first  
			   << "\" name=\""  
			   << name 
			   << "\" type=\"checkbox\"  CHECKED > "  
			   << i->second << "</label></li>"; 
		 even = false; 
	       } 
	     else 
	       { 
		 std::cout << "<li><label for=\""  
			   << i->first  
			   << "\"><input id=\""  
			   << i->first  
			   << "\" name=\""  
			   << name 
			   << "\" type=\"checkbox\" CHECKED > "  
			   << i->second << "</label></li>"; 
		 even = true; 
	       } 
	   }
 	++i; 

       } 
     
     // repeat for unselected items

     i = recs.begin(); 
     
     while ( i != recs.end() ) 
       { 	 
	 
         if ( sel.find( i->first ) == sel.end() )	
	   {
	     if ( even )  
	       { 
		 std::cout << "<li class=\"alt\"><label for=\""  
			   << i->first  
			   << "\"><input value=\""  
			   << i->first  
			   << "\" name=\""  
			   << name 
			   << "\" type=\"checkbox\" /> "  
			   << i->second << "</label></li>"; 
		 even = false; 
	       } 
	     else 
	       { 
		 std::cout << "<li><label for=\""  
			   << i->first  
			   << "\"><input value=\""  
			   << i->first  
			   << "\" name=\""  
			   << name 
			   << "\" type=\"checkbox\" /> "  
			   << i->second << "</label></li>"; 
		 even = true; 
	       } 
	   }
 	++i; 

       } 

  
     std::cout << "</ul>"; 
   } 



  void add_combocheckbox( const std::string & name , const std::set<std::string> & recs , const std::set<std::string> & sel )
  {
    
    std::cout << "<ul class=\"checklist cl1\">";
    
    // first add selected items, then unselected
    
    std::set<std::string>::const_iterator i = recs.begin();
    
    bool even = false;

    while ( i != recs.end() )
      {
	if ( sel.find( *i ) != sel.end() )	
	  {
	    if ( even ) 
	      {
		std::cout << "<li class=\"alt\"><label for=\"" 
			  << *i
			  << "\"><input value=\"" 
			  << *i
			  << "\" name=\"" 
			  << name
			  << "\" type=\"checkbox\" CHECKED > " 
			  << *i << "</label></li>";
		even = false;
	      }
	    else
	      {
		std::cout << "<li><label for=\"" 
			  << *i
			  << "\"><input value=\"" 
			  << *i
			  << "\" name=\"" 
			  << name
			  << "\" type=\"checkbox\" CHECKED > " 
			  << *i << "</label></li>";
		even = true; 
	      }
	  }

	++i;

      }
    
    // now unselected items
    i = recs.begin();
    
    while ( i != recs.end() )
      {
	if ( sel.find( *i ) == sel.end() )
	  {
	    if ( even ) 
	      {
		std::cout << "<li class=\"alt\"><label for=\"" 
			  << *i
			  << "\"><input value=\"" 
			  << *i
			  << "\" name=\"" 
			  << name
			  << "\" type=\"checkbox\" /> " 
			  << *i << "</label></li>";
		even = false;
	      }
	    else
	      {
		std::cout << "<li><label for=\"" 
			  << *i
			  << "\"><input value=\"" 
			  << *i
			  << "\" name=\"" 
			  << name
			  << "\" type=\"checkbox\" /> " 
			  << *i << "</label></li>";
		even = true;
	      }
	  }
	++i;

      }


    std::cout << "</ul>";
  }


  void form( const std::string & name , const std::string & cgi )
  {
    std::cout << "<form name=\"" << name << "\" action=\"" << cgi << "\" method=\"GET\"> ";
  }

  void end_form()
  {
    std::cout << "</form>";
  }

  
  void init_filterlist()
  {    
    // requires that 'filterlist.js' is in run-path, i.e. same place as 'index'
    std::cout << " <script type=\"text/javascript\" src=\"filterlist.js\"></script>";        
  }


  void add_hidden_value( const std::string & key , const std::string & val )
  {
    std::cout << "<input type=\"hidden\" NAME=\"" << key << "\" value=\"" << val << "\">";
  }


  void add_filterlist( const std::string & form , const std::string & name , const std::set<std::string> & recs )
  {

    std::cout << "<select name=\"" << name << "\" size=5>";

    std::set<std::string>::iterator i = recs.begin();
    
    while ( i != recs.end() )
      {
	std::cout << "<option>" << *i ; 
	++i;
      }
    std::cout << "</select><br>";
    
    std::cout << "<SCRIPT TYPE=\"text/javascript\">"
	      << "var myfilter_" << name << " = new filterlist(document." << form << "." << name << ");"
	      << "</script>";
    
    //
    // Note: if you have a very large select list,
    // you should deactivate the real-time filtering
    // on the INPUT field below - remove the onKeyUp attribute.
    //

    std::cout << "<A HREF=\"javascript:myfilter_fl1.reset()\" TITLE=\"Clear the filter\">Clear</A>"
	      << "<A HREF=\"javascript:myfilter_fl1.set('^A')\" TITLE=\"Show items starting with A\">A</A>"
	      << "<A HREF=\"javascript:myfilter_fl1.set('^B')\" TITLE=\"Show items starting with B\">B</A>"
	      << "<A HREF=\"javascript:myfilter_fl1.set('^C')\" TITLE=\"Show items starting with C\">C</A>";


    std::cout << "<INPUT NAME=regexp onKeyUp=\"myfilter_" << name << ".set(this.value)\">";

    //    std::cout << "<INPUT TYPE=button onClick="myfilter.set(this.form.regexp.value)\" value=\"Filter\"> "
    ///* <INPUT TYPE=button onClick="myfilter.reset();this.form.regexp.value=''" value="Clear"> */
    
    std::cout << "<br>";

/* <INPUT TYPE=checkbox NAME="toLowerCase" */
/* onClick="myfilter.set_ignore_case(!this.checked)"> */
/* Case-sensitive */

  }

}

#endif;
