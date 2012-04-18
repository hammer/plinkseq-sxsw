#ifndef __PSEQ_HTML_LINES_H__
#define __PSEQ_HTML_LINES_H__

#include <string>
#include <sstream>
#include <map>
#include <vector>

namespace Pseq {
  namespace Helper { 
    
    struct LinkBox { 

      LinkBox() { }
 
      LinkBox(const int x1, const int y1 , const int x2 , const int y2 , const std::string & url )
      : x1(x1), y1(y1), x2(x2), y2(y2), url(url) { } 

      LinkBox(const int x1, const int y1 , const std::string & text , const std::string & url )
      : x1(x1), y1(y1), x2(-1), y2(-1), text(text), url(url) { } 
      
      int x1, x2, y1, y2;
      std::string url;
      std::string text;
    };
    
    class Lines {

    public:
      
      Lines() { } 
      
      void setobj( const std::string & n , const int x , const int y ) 
      {
	name = n; xs = x; ys = y;
	names.push_back( name );
	
	colmap[ "black" ] = "\"rgb(0,0,0)\"";  
	colmap[ "white" ] = "\"rgb(255,255,255)\"";  
	colmap[ "gray" ]  = "\"rgb(100,100,100)\"";  
	
	colmap[ "red" ]   = "\"rgb(200,000,000)\"";  
	colmap[ "green" ] = "\"rgb(000,200,000)\"";  
	colmap[ "blue" ] = "\"rgb(000,000,200)\"";  

	colmap[ "lightred" ]   = "\"rgb(250,150,150)\"";  
	colmap[ "lightgreen" ] = "\"rgb(150,250,150)\"";  
	colmap[ "lightblue" ] = "\"rgb(150,150,250)\"";  

	colmap[ "yellow" ] = "\"rgb(000,200,200)\"";  
	colmap[ "pink" ] = "\"rgb(200,100,100)\"";  
      }
      
      Lines( const std::string & n , const int x , const int y ) 
	{
	  setobj(n,x,y);
	}

            
      // keep track of all different canvases on a page
      
      static std::vector<std::string> names;
      static std::string loader() 
      { 
	std::stringstream ss;
	ss << "<script type=\"text/javascript\">"
	   << "function canvasdraw() { ";
	for (int i=0;i<names.size();i++)
	  ss << "canvasdraw_" << names[i] << "();";
	ss << "} </script>";	
	return ss.str();
      }
      
      std::string preamble() const;

      std::string place() const;
      
      void box( const int, const int, const int , const int , 
		const std::string & bg = "white" ,
		const std::string & fg = "black" ,
		const int width = 1 ,
		const std::string & name = "" ,
		const std::string & link = "" );
      
      void line( const int , const int , const int , const int , 
		 const std::string & fg = "black" , const int = 1);		 
      
      void text( const int , const int , 
		 const std::string & text , 
		 const std::string & fg = "black" ,
		 const std::string & name = "" ,
		 const std::string & link = "" );

    private:
      
      int xs, ys;  // dimension of this canvas
      std::string name; // and its name
      std::map<std::string,LinkBox> linkmap;
      std::map<std::string,std::string> colmap;
      std::vector<std::string> elements;
    };
  };
};

#endif
