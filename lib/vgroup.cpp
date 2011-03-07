#include "vgroup.h"
#include "gstore.h"

using namespace std;
using namespace Helper;

extern GStore * GP;

void VariantGroup::force_add(Variant & v)
{
  vars.push_back(v);
}


void VariantGroup::add(Variant &v ) 
{ 
  
  if ( is_complete ) return;
  
  // If the first variant added, will obviously be okay
  
  if ( vars.size() == 0 ) 
    {		      
      gname = v.meta.get1_string( PLINKSeq::META_GROUP() );
      
      if ( mask.named_grouping() && gname == "" ) 
	is_complete = true;
      vars.push_back(v);
      return;
    }
  
  // If a named grouping has been specified, does this
  // new variant conform?
  
  if ( mask.named_grouping() )
    {
      if ( gname != v.meta.get1_string( PLINKSeq::META_GROUP() ) )
	{
	  is_complete = true;
	  return;
	}
      vars.push_back(v);
      return;
    }
 
  // We should not have reached here, as some grouping function 
  // should have been specified in the mask in order to use 
  // a VariantGroup. 
  
  is_complete = true;
  
  return;
  
}



int VariantGroup::n_individuals() const
{
  return size() > 0 ? vars[0].size() : 0;
}


void VariantGroup::clear(Variant & v)
{ 

  vars.clear();
  vars.push_back(v);
  
  gname = v.meta.get1_string( PLINKSeq::META_GROUP() );

  if ( mask.named_grouping() && gname == "" )
    is_complete = true;
  else
    is_complete = false;
  
  return;
}


void VariantGroup::clear()
{ 
  vars.clear(); 
  gname = "";
  is_complete = false;
}

string VariantGroup::dump( bool vmeta , 
			   bool vexpand , 
			   bool geno , 
			   bool gmeta , 
			   bool transpose , 
			   bool rarelist , 
			   bool show_phenotype )
{


  stringstream ss;
  
  // Overall group details
  
  ss << "NAME=[" << name() << "]\t" 
     << "N(V)=" << size() << "\t"
     << "N(I)=" << n_individuals() << "\n";

  if ( size() == 0 ) 
    {
      return ss.str();
    }

  //
  // Rare list
  //

  if ( rarelist ) 
    {
      
      std::map<int,std::string> istr;
      std::map<int,std::string> vstr;
      
      std::map<int,int> icnt;
      std::map<int,int> vcnt;

      for (int i=0; i < size(); i++)
	{

	  // Get the minor allele
	  int c = 0 , n = 0;
	  bool altmin = var(i).n_minor_allele( c , n );
	  
	  for (int j=0; j < n_individuals(); j++ )
	    {	      
	      const Genotype & g = (*this)(i,j);
	      
	      if ( g.minor_allele( altmin ) ) 
		{
		  
		  if ( istr[j] != "" ) istr[j] += ",";
		  istr[j] += var(i).coordinate() + ";" + var(i).consensus.label( (*this)(i,j) );
		  
		  if ( vstr[i] != "" ) vstr[i] += ",";
		  vstr[i] += ind(j)->id() + ";" + var(i).consensus.label( (*this)(i,j) );
		  
		  ++icnt[j];
		  ++vcnt[i];
		  
		  if ( gmeta ) 
		    {
		      std::stringstream s2;
		      s2 << (*this)(i,j).meta; 
		      istr[j] += ";" + s2.str(); 
		      vstr[i] += ";" + s2.str(); 
		    }
		}
	    }
	}


      //
      // Create rarelist display
      //
      
      std::map<int,std::string>::iterator ii = istr.begin();
      while ( ii != istr.end() )
	{
	  ss << "I\t" 
	     << name() << "\t"	 
	     << ind( ii->first )->id() << "\t";

	  if ( show_phenotype )
	    {
	      if ( GP->phmap.type() == PHE_DICHOT )
		{
		  if ( ind( ii->first )->affected() == CASE ) ss << "CASE\t";
		  else if ( ind( ii->first)->affected() == CONTROL ) ss << "CONTROL\t";
		  else ss << ".\t";
		}
	    }

	  
	  ss << icnt[ ii->first ] << "\t"
	     << "\t"
	     << ii->second 
	     << "\n";	    
	  ++ii;
	}

      std::map<int,std::string>::iterator iv = vstr.begin();
      while ( iv != vstr.end() )
	{
	  ss << "V\t" 
	     << name() << "\t"
	     << var( iv->first ) << "\t"
	     << var( iv->first ).reference() << "\t"
	     << vcnt[ iv->first ] << "\t"
	     << iv->second 
	     << "\n";	    
	  ++iv;
	}
      
      return ss.str();
    }
  
  

  //
  // Per-variant summary
  //
  
  for (int i=0; i< size(); i++)
    {
      
      if ( ! ( transpose && !vmeta ) )
	{
	  ss << "V" << i+1 << "\t";

	  ss << var(i) << "\t";
	  
	  // what about sample meta-information?
	  
	  if ( vmeta ) 
	    {
	      if ( vexpand ) 
		ss << "\n" << var(i).meta.display();
	      else
		ss << var(i).meta;   
	    }
	  ss << "\n";
	}
    }

  if ( ! geno ) return ss.str();
  
  ss << "\n";




  //
  // Show genotypes 
  //

  if ( transpose ) 
    {
      // show var x ind
      
      // indiv. header line

      ss << "V";
      for (int j=0; j < n_individuals(); j++ )
	ss << "\t" << ind(j)->id();	
      ss << "\n";
      
      // Variant rows
      for (int i=0; i< size(); i++)
	{	    

	  ss << var(i) << "\t" << var(i).reference();

	  for (int j=0; j < n_individuals(); j++ )
	    {
	      ss << "\t" << var(i).consensus.label( (*this)(i,j) ) ;
	      if ( gmeta ) 
		ss << " " << (*this)(i,j).meta;
	    }
	  ss << "\n";
	}
    }
  else
    {
      // ind x var

      // var. header line

      ss << "I";
      for (int i=0; i< size(); i++)
	ss << "\tV" << i+1;

      ss << "\n";

      // Indiv rows
      for (int j=0; j < n_individuals(); j++ )
	{ 
	  ss << ind(j)->id();
	  for (int i=0; i< size(); i++)
	    {	    
	      ss << "\t" << var(i).consensus.label( (*this)(i,j) ) ;
	      if ( gmeta ) 
		ss << " " << (*this)(i,j).meta;
	    }
	  ss << "\n";
	}
    }
      
  ss << "\n";
  
  return ss.str();
}



std::string VariantGroup::coordinate() const 
{
  if ( size() == 0 ) return "NA";
  int chr = var(0).chromosome();
  int bpmin = var(0).position();
  int bpmax = var(0).stop();
  for (int i=1; i<size(); i++)
    {
      if ( var(i).chromosome() != chr ) return "NA";
      if ( var(i).position() < bpmin ) 
	bpmin = var(i).position();
      if ( var(i).stop() > bpmax ) 
	bpmax = var(i).stop();      
    }
  return Helper::chrCode( chr ) + ":" + Helper::int2str( bpmin ) + ".." + Helper::int2str( bpmax ); 
}

int VariantGroup::midposition() const
{
  if ( size() == 0 ) return -1;
  int chr = var(0).chromosome();
  int bpmin = var(0).position();
  int bpmax = var(0).stop();
  for (int i=1; i<size(); i++)
    {
      if ( var(i).chromosome() != chr ) return -1;
      if ( var(i).position() < bpmin ) 
	bpmin = var(i).position();
      if ( var(i).stop() > bpmax ) 
	bpmax = var(i).stop();      
    }
  
  return (int)(bpmin + (bpmax-bpmin)/2 );
  
}

int VariantGroup::span() const
{
  if ( size() == 0 ) return -1;
  int chr = var(0).chromosome();
  int bpmin = var(0).position();
  int bpmax = var(0).stop();
  for (int i=1; i<size(); i++)
    {
      if ( var(i).chromosome() != chr ) return -1;
      if ( var(i).position() < bpmin ) 
	bpmin = var(i).position();
      if ( var(i).stop() > bpmax ) 
	bpmax = var(i).stop();      
    }
  return bpmax-bpmin+1;
}
