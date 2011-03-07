#ifndef __INDIVIDUAL_H__
#define __INDIVIDUAL_H__

#include <string>
#include <vector>
#include <ostream>

#include "meta.h"
#include "defs.h"

class Measure;

class Individual 
{ 

  // Corresponds to DB id
  uint64_t indiv_id;
  
  // Corresponds to investigator IDs
  std::string id1, f_id, i_id;
  
  sType sex_code;
  std::string pat_id, mat_id;

  Individual * patp, * matp;
  
  // Meta-information (public) to store phenotypes

  // Slots for commonly used attributes 

  bool    is_missing;
  
  affType aff;
  
  double  phenotype;    
  
  int     grp;
  
  void construct()
    {
      indiv_id = 0;
      id1 = ".";
      f_id = i_id = ".";  
      pat_id = mat_id = ".";
      patp = matp = NULL;
      is_missing = false;
      sex_code = UNKNOWN_SEX;

      aff = UNKNOWN_PHE;
      phenotype = 0;
      grp = 0;
      
    }

 public:  
  
  Individual() 
    {
      construct();
    }

  Individual(std::string i, std::string f, std::string ind) 
    {
      construct();
      id1 = i;
      f_id = f;
      i_id = ind;
    }

  Individual(std::string i) 
    {
      construct();
      id1 = i;
    }

  Individual(uint64_t n, std::string i) 
    {
      construct();
      indiv_id = n;
      id1 = i;
   }

  uint64_t idx() const { return indiv_id; }

  std::string id() const { return id1; }
  void id(const std::string s) { id1 = s; } 
  
  std::string fid() const { return f_id; }
  std::string iid() const { return i_id; }

  void idx(uint64_t n) { indiv_id = n; }
  
  void fid(std::string s) { f_id = s; } 
  void iid(std::string t) { i_id = t; } 
  
  void pat( const std::string & s) { pat_id = s; }
  void mat( const std::string & s) { mat_id = s; }
  void sex( const sType s ) { sex_code = s; }
  void sex( const std::string & s ) 
  {
    if ( s == "1" ) sex_code = MALE;
    else if ( s == "2" ) sex_code = FEMALE;
    else if ( s == "0" ) sex_code = UNKNOWN_SEX;
    else if ( s .size() > 0 )
      {
	std::string ss = s.substr(0,1);
	if ( ss == "m" || ss == "M" ) sex_code = MALE;
	else if ( ss == "f" || ss == "F" ) sex_code = FEMALE;
      }
    else      
      sex_code = UNKNOWN_SEX;
  }

  std::string father() const { return pat_id ; }
  std::string mother() const { return pat_id ; }

  Individual * pat() const { return patp; }
  Individual * mat() const { return matp; }

  void affected( const affType d) { aff = d; }
  void affected( const bool b) { if (b) aff = CASE; else aff = CONTROL; }
  affType affected() const { return missing() ? UNKNOWN_PHE : aff; }
  
  void group( const std::string & l ) 
    { 
      // handle missing codes
      if ( l == "" || l == "." ) 
	{	  
	  grp = -1;
	  return;
	}

      std::map<std::string,int>::iterator j = factor.find(l);
      if ( j == factor.end() )
	{
	  grp = factor.size() + 1;
	  factor[l] = grp;
	  factor2[grp] = l;
	}
      else
	grp = j->second;
    }
  
  int group() const { return grp; }
  
  std::string group_label() const 
    { return Individual::label( group() ); }
  

  static std::string label(int i) 
    {
      std::map<int,std::string>::iterator j = factor2.find(i);
      return j == factor2.end() ? "." : j->second;
    }

  void qt( const double q ) { phenotype = q; }
  double qt() const { return phenotype; }

  void missing(const bool b) { is_missing = b; }
  bool missing() const { return is_missing; }
  bool included() const { return ! is_missing; }

  int sex() const 
    { 
      if ( sex_code == MALE ) return 1;
      else if ( sex_code == FEMALE ) return 2;
      else return 0;
    }

  void pat(Individual * p) { patp = p; }
  void mat(Individual * m) { matp = m; }

  bool operator< (const Individual & b) const
    { return id1 < b.id1 ; }

  bool operator== (const Individual & b) const
    { return id1 == b.id1 ; }

  MetaInformation<IndivMeta> meta;    

  friend std::ostream & operator<<( std::ostream & out, const Individual & p)
  { 
    out << p.id() << "\t" 
	<< p.fid() << "\t" 
	<< p.iid() << "\t" 
	<< p.missing() << "\t"
	<< p.sex() << "\t" 
	<< p.pat() << "\t" 
	<< p.mat() << "\t"
	<< p.meta;
    
    return out;
  }


 private:
  
  static std::map<std::string,int> factor;

  static std::map<int,std::string> factor2;

};


#endif
