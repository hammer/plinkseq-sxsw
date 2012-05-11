
#include "lines.h"
#include <sstream>

std::vector<std::string> Pseq::Helper::Lines::names;

std::string Pseq::Helper::Lines::preamble() const
{

  std::stringstream  ss;

  ss << "<script type=\"text/javascript\">"
     << "var canvas_" << name << ";"
     << "var ctx_" << name << ";"
     << "var textHeight_" << name << " = 12;";
  

  // define all "inlink" variables
  std::map<std::string,LinkBox>::const_iterator ii = linkmap.begin();
  while ( ii != linkmap.end() )
    {
      ss << "var inLink_" << name << "_" << ii->first << ";";
      ++ii;
    }  

  // on-init draw() function (called from body)
  ss << "function canvasdraw_" << name << "(){"
     << " canvas_" << name <<" = document.getElementById(\"" << name << "\");"
     << "if( canvas_" << name << ".getContext ){"
     << "ctx_" << name << "=canvas_" << name << ".getContext(\"2d\");";
  
  // Font
  ss << "ctx_" << name << ".font='12px sans-serif';";
  
  // add draw events
  for (int e=0; e<elements.size(); e++)
    ss << elements[e] ;
  
  //add mouse listeners
  ss << "canvas_" << name << ".addEventListener(\"mousemove\", on_mousemove_" << name << ", false);"
     << "canvas_" << name << ".addEventListener(\"click\", on_click_" << name << ", false);";
  
  ss << "} }";


  // define the mouse-handler 
  
  ss << "function on_mousemove_" << name << "(ev) {"
     << "var x, y;"
    
     << "if (ev.layerX || ev.layerX == 0 ) { " // for firefox
     << "x = ev.layerX;"
     << "y = ev.layerY;"
     << "}"
     << "x-=canvas_" << name << ".offsetLeft;"
     << "y-=canvas_" << name << ".offsetTop;";

  // check all links
  ss << "if (false) { var d=false; }"; // never get called, but allows all other events as else if () 
  
  std::map<std::string,LinkBox>::const_iterator jj = linkmap.begin();
  while ( jj != linkmap.end() )
    {
      
      if ( jj->second.x2 < 0 ) 
	{
	  ss << "else if (  x >= " << jj->second.x1 
	     << " && x <= " << jj->second.x1 << " + ctx_" << name << ".measureText(\"" << jj->second.text << "\").width "  
	     << " && y >= " << jj->second.y1 << " - textHeight_" << name << " " 
	     << " && y <= " << jj->second.y1 
	     << " ) {"
	     << "document.body.style.cursor = \"pointer\";";
	}
      else
	{
	  ss << "else if (  x >= " << jj->second.x1 
	     << " && x <= " << jj->second.x2 
	     << " && y >= " << jj->second.y1 
	     << " && y <= " << jj->second.y2 
	     << " ) {"
	     << "document.body.style.cursor = \"pointer\";";
	}

      std::map<std::string,LinkBox>::const_iterator kk = linkmap.begin();
      while ( kk != linkmap.end() ) 
	{ 
	  if ( kk != jj ) { ss << "inLink_" << name << "_" << kk->first << "=false;"; } 
	  ++kk; 
	}       
      ss << "inLink_" << name << "_" << jj->first << "=true; }";
      ++jj;
    }

  // no links hit
  ss << "else { document.body.style.cursor = \"\"; ";
  jj = linkmap.begin();
  while ( jj != linkmap.end() ) 
    { 
      ss << "inLink_" << name << "_" << jj->first << "=false;"; 
      ++jj; 
    } 
  ss << " } "
     << " } ";

  // finally, define on_click actions
  ss << "function on_click_" << name << "(e) { ";
  ss << "if ( false ) { var d = false ; } ";
  
  std::map<std::string,LinkBox>::const_iterator kk = linkmap.begin();
  while ( kk != linkmap.end() ) 
    { 
      ss << "else if (inLink_" << name << "_" << kk->first << ") { window.location = \"" << kk->second.url << "\"; }";
      ++kk;
    }
  ss << " } ";

  ss << "</script>";

  return ss.str();
}


std::string Pseq::Helper::Lines::place() const
{
  std::stringstream ss;
  ss << " <canvas id=\"" << name << "\" width=\"" << xs << "\" height=\"" << ys << "\">Canvas not supported.</canvas>";
  return ss.str();
}
 
     
void Pseq::Helper::Lines::box( const int x1, const int y1, const int x2, const int y2, 
			       const std::string & bg  ,
			       const std::string & fg ,
			       const int width , 
			       const std::string & n ,
			       const std::string & link )
{
  std::stringstream ss;
  
  ss << "ctx_" << name << ".fillStyle = " << colmap[ bg ] << ";"
     << "ctx_" << name << ".strokeStyle = " << colmap[ fg ] << ";"
     << "ctx_" << name << ".lineWidth = \"" << width << "\";"
     << "ctx_" << name << ".fillRect(" << x1 << "," << y1 << "," << x2-x1 << "," << y2-y1 << ");"
     << "ctx_" << name << ".strokeRect(" << x1 << "," << y1 << "," << x2-x1 << "," << y2-y1 << ");";

  elements.push_back( ss.str() );  
  if ( n != "" ) linkmap[ n ] = LinkBox( x1,y1,x2,y2,link );
}
      
void Pseq::Helper::Lines::line( const int x1 , const int y1 , const int x2 , const int y2 , 
				const std::string & fg , const int width )
				
				
{

  std::stringstream ss;
  
  ss << "ctx_" << name << ".strokeStyle = " << colmap[ fg ] << ";"
     << "ctx_" << name << ".lineWidth = \"" << width << "\";"     
     << "ctx_" << name << ".beginPath();"
     << "ctx_" << name << ".moveTo(" << x1 << "," << y1 << ");"
     << "ctx_" << name << ".lineTo(" << x2 << "," << y2 << ");"
     << "ctx_" << name << ".closePath();"
     << "ctx_" << name << ".stroke();";

  elements.push_back( ss.str() );  
  
    
}

void Pseq::Helper::Lines::text( const int x1, const int y1 , 
				const std::string & text , 
				const std::string & fg ,
				const std::string & n ,
				const std::string & link )
{
  
  std::stringstream ss; 
  
  ss << "ctx_" << name << ".fillStyle = " << colmap[ fg ] << ";"     
     << "ctx_" << name << ".fillText(\""<<text<<"\"," << x1 << "," << y1 << ");";

  elements.push_back( ss.str() );

  if ( n != "" ) linkmap[ n ] = LinkBox( x1 , y1 , text, link );  
}


