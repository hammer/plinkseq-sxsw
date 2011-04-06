#ifndef __PHMAP_H__
#define __PHMAP_H__

#include <string>
#include <map>
#include <set>

#include "defs.h"
#include "individual.h"
#include "matrix.h"

class IndDBase;
class Individual;
class IndividualMap;

// The PhenotypeMap is a helper class responsible for knowing and 
// setting the current phenotype, and for all transactions with the
// INDDB.  This class takes
// care of the fact that VARDB and INDDB might not contain perfectly matching
// lists of individuals

// The main access to individuals (and their phenotype information)
// will be through the IndividualMap, however.


class PhenotypeMap {
  
 public:
  
  PhenotypeMap(IndDBase * inddb = NULL );
  
  ~PhenotypeMap();
 
  
  // Total number of individuals from database, for whom we have
  // phenotype information
  
  int size() const
  {
    return pmap.size();
  }

  bool exists( const std::string & id ) const
  { 
    return pmap.find( id ) != pmap.end(); 
  }

  
  // Return a pointer to an individual

  Individual * ind( const std::string & id ) const
  {
    std::map<std::string,Individual*>::const_iterator i = pmap.find(id);
    return i == pmap.end() ? NULL : i->second;
  }
  

  // 
  // Add individuals to map
  //

  Individual * new_individual( const std::string & );
  void begin();
  void commit();

  //
  // Prune phenotype map to be aligned to an IndividualMap
  //

  void align( const std::set<std::string> & );

  // 
  // Clean up
  //

  void reset();

  
  // 
  // Phenotype related functions
  //

  int set_phenotype( const std::string & phenotype );

  int make_phenotype( const std::string & make_phenotype );

  int set_strata( const std::string & s );

  std::map<std::string,int> summarise_phenotype( const std::string & phenotype );

  std::map<std::string,int> summarise_phenotype();
  
  pType type() const { return phenotype_type; }

  //int attach_covariates( const std::string & );
  
  Data::Matrix<double> covariates( const std::vector<std::string> & c , const IndividualMap & indmap );


  //
  // Queuries
  //

  bool phenotype_set() const { return phenotype_type != PHE_NONE; }

  std::string phenotype() const { return phenotype_name; }
  
  bool strata_set() const { return use_strata; } 
  
  std::string strata() const { return strata_name; }
  
  //
  // Directly load
  //

  void direct_load( const std::string & f , const std::string & l );

  
  //
  // Display functions
  //

  friend std::ostream & operator<<( std::ostream & out , PhenotypeMap & m )
    {
      std::map<std::string, Individual* >::const_iterator i = m.pmap.begin();
      while ( i != m.pmap.end() )
	{
	  out << i->first << "\t"
	      << *(i->second) << "\n"; 
	  ++i;
	}      
      return out;
    }
  

 private:


  //
  // A phenotype map must be attached to an INDDB when constructed
  //
  
  IndDBase * inddb;
  
  
  //
  // Pointers to individuals given a unique ID
  //
  
  std::map<std::string,Individual*> pmap;
  
  
  //
  // Phenotype information
  //

  std::string phenotype_name;

  pType phenotype_type;

  bool use_strata;
  
  std::string strata_name;

};



#endif
