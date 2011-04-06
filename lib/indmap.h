
#ifndef __INDMAP_H__
#define __INDMAP_H__

#include <string>
#include <vector>
#include <ostream>

#include "meta.h"
#include "defs.h"

class VarDBase;
class Mask;
class PhenotypeMap;
class Individual;

class IndividualMap {
  
 public:

  IndividualMap() 
    {
      reset();
    }
  
  int populate( VarDBase & , PhenotypeMap & , const Mask & );

  int populate( const std::vector<std::string> & );
  
  void reset();
  
  inline bool flat() const 
  {
    return unflattering ? false : is_flat;
  }
  
  /// Over-ride usual 'flat' designation
  void force_unflat( const bool b ) 
  {
    unflattering = b;
  }
  
  
  inline bool multi_sample() const 
    {
      return unflattering ? true : is_multi_sample;
    }


  inline bool single_sample() const
    {
      return unflattering ? false : ! is_multi_sample;
    }
  
  inline int get_slot(int f, int i) const
  {	    
    // consensus 1:1 mapping?
    if ( f == 0 ) return i;
    // retrive from file
    int2 p(f,i);
    std::map<int2,int_string_pair>::const_iterator k = ialign.find( p );
    return k == ialign.end() ? -1 : k->second.i;
  }
  
  inline int size() const
  {
    return n_uniq;
  }
  
  int n_samples() const 
  {
    return wsample.size();
  }

  std::set<int> samples() const 
  {
    std::set<int> w;
    std::map<int,std::map<int,int> >::const_iterator i = wsample.begin();
    while ( i != wsample.end() )
      {
	w.insert( i->first );
	++i;
      }
    return w;
  }

  inline int size(const int f) const
    {    
      if ( f == 0 ) return size();
      std::map<int,std::map<int,int> >::const_iterator i = wsample.find( f );
      if ( i == wsample.end() ) return 0;
      return i->second.size();
    }

  inline bool unique(const int i) const
  { 
    return uniq[i] != int2(-1,-1); 
  }
  
  inline int2 unique_mapping( const int i ) const 
  { 
    return uniq[i]; 
  }
  
  inline std::set<int2> multiple_mapping( const int i ) const 
  { 
    return mult[i]; 
  }
  
  inline int sample_remapping(int f, int j) const
  {
    std::map<int,std::map<int,int> >::const_iterator i = wsample.find( f );
    if ( i == wsample.end() ) return -1;
    std::map<int,int>::const_iterator k = i->second.find(j);
    return k == i->second.end() ? -1 : k->second;
  }
  
  /// Primary access to individuals
  
  Individual * operator() (const int i) const 
    {
      return indiv[i];
    }

  inline Individual * ind(const int i) const
  {
    return indiv[i];
  }
  
  inline Individual * ind( const std::string & s) const
  {
    std::map< std::string,int >::const_iterator i = id2pos.find( s );
    if ( i == id2pos.end() ) return NULL;
    return ind( i->second );  
  }
  
  inline int ind_n( const std::string & s ) const 
  {
    std::map< std::string,int >::const_iterator i = id2pos.find( s );
    return i == id2pos.end() ? -1 : i->second; 
  }
  
  inline std::vector<std::string> ind_id() const
  {
    return idvec;
  }
  
  inline std::string ind_id(const int i) const
  {
    return idvec[i];
  }
  
  
  /// Indicate which fileset a given individual belongs to, if unique (else 0);
  
  inline int sample(const int i) const
    {
      int2 k = uniq[i];
      return k.p1 != -1 ? k.p1 : 0 ;
    }
  
  


  //
  // Display functions
  //
  

  friend std::ostream & operator<<( std::ostream & out , const IndividualMap & a ) 
    {
      
      out << "Primary mapping\n";      
      std::map<int2,int_string_pair>::const_iterator i = a.ialign.begin();
      while ( i != a.ialign.end() )
	{
	  out << i->first.p1 << " , " 
	      << i->first.p2 << " --> "
	      << i->second.i << " , " 
	      << i->second.s << "\n";
	  ++i;
	}
      
      out << "Reverse mapping\n";     
      out << a.uniq.size() << " " << a.mult.size() << "\n";
      for (int j=0; j<a.uniq.size(); j++)
	{
	out << j << "\t"
		  << a.uniq[j].p1 << " , "
		  << a.uniq[j].p2 << "\t"
		  << a.mult[j].size() << " { ";
	std::set<int2>::iterator k = a.mult[j].begin();
	while ( k != a.mult[j].end() )
	  {
	    out << k->p1 << ", " 
		<< k->p2 << " ";
	    ++k;
	  }
	out << "}\n";
      }

    return out;

  }

  void construct_wsint_vector();
  std::vector<int> * svar2consensus(const int f);

 private:
  
  //
  // Private functions
  // 

  void insert( const int file_id, 
	       const int supposed_j , 
	       const int slot_j, 
	       const int slot_k, 
	       const std::string & id );
    
  std::map<int,std::string> map_slot_to_id() const;

  // map of final svar slots -> consensus slot in convenient
  // std::vector<int> form
  std::map<int,std::vector<int> > wsint;


  //
  // Private data-members
  // 


  // Because N will typically be small, we track sample information 
  // in a few different, partially redundant ways here, to make 
  // common modes of access more convenient


  // For a given variant database containing multiple files, this is
  // constructed by taking all files and aligning the individuals.
  // Depending on the merge-mode, we either assume that people
  // across files with the same ID are dupes, or else we append with
  // a unique file ID. We are then left with a mapping of
  
  // {file-i,slot-j} --> { string ID , slot-k }
  
  // where slot-k indicates the position in the
  // file
  
  std::map<int2,int_string_pair> ialign;
  
  // Within-sample mapping (i.e. post filters)
  std::map<int,std::map<int,int> > wsample;
  // and a stable set of vectors vector
  
  // map from the 0..N-1 to the individual file, or file(s) (int2 = f/i pairs)
  // these are populated by vardb.populate_individual_alignment()
  
  std::vector<int2> uniq;       
  std::vector<std::set<int2> > mult;

  
  // Keep track of IDs of included individuals
  
  std::set< std::string > ids;
  
  std::map< std::string, int > id2pos;

  // Pointers to actual individuals (in the MAP if present) (0..n-1)
  
  std::vector<Individual*> indiv;
  
  // And, for convenience, the actual string IDs too
  
  std::vector< std::string> idvec;


  // And also do the same for within-sample 

  std::map<int, std::map<int, Individual*> > sample_indiv;
  std::map<int, std::map<int, std::string> > sample_idvec;

  // Track number of unique individuals, a common lookup; 
  //  This will correspond to the size of indiv* and 
  //  
  
  int n_uniq;

  /// Does this mapping a simple, or 'flat' one?
  
  bool is_flat;
  bool is_multi_sample;
  bool unflattering; // over-ride (i.e. force to run as a multi-sample variant)
  
};


#endif
