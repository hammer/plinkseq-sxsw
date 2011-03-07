#ifndef __VGROUP_H__
#define __VGROUP_H__

#include "variant.h"
#include "individual.h"
#include "mask.h"
#include <map>

/*!
  @class VariantGroup
  @brief A group of variants
  $Header$
*/


class VariantGroup {    

 public:

  VariantGroup(Mask & mask) : mask(mask) 
    { 
      fileset = 0;
      is_complete = false;
    } 
  
  
  VariantGroup(const VariantGroup & g)
    : mask( g.mask )
    {
	vars = g.vars;
	gname = g.gname;
	fileset = g.fileset;
	is_complete = g.is_complete;
    }
  
  VariantGroup & operator= (const VariantGroup & g)
    {
      if ( this != &g)
	{
	  
	  mask = g.mask;

	  // Manually copy other elements	  
	  vars = g.vars;
	  gname = g.gname;
	  fileset = g.fileset;
	  is_complete = g.is_complete;

	}
      return *this;
    }

  /*!
    
    @return Total number of variants in the group 
  */

  int n_variants() const { return vars.size(); }
  

  /*!
    This function is identical to n_variants()
    @return Total number of variants in the group 
  */
  
  int size() const { return vars.size(); }


  /*! 
    Return chromosomal location (span) of group
    @return Textual representation of position, NA if >1 chr
  */

  std::string coordinate() const ;

  int midposition() const;

  int span() const;


  /*!
    @return Number of individuals in this group
  */
  
  int n_individuals() const;
    
  
  
  
  std::string name() const { return gname; } 

  void name(const std::string n) { gname = n; }


    //
    // Access to underlying variant data
    //


    Variant & var(const int i) { return vars[i]; }

    Variant & operator()(const int i) { return vars[i]; }

    const Variant & var(const int i) const { return vars[i]; }    
    
    const Variant & operator()(int i) const { return vars[i]; } 
    
    const Genotype & operator()(const int i, const int j) const 
    { return vars[i](j); }
    

    

    //
    // Access to genotype data
    //

    const Genotype & geno(int i,int j) const 
      { 
	return vars[i](j); 
      }    


    // Convenience features to get to individuals, assuming that
    // all variants in the group share the same individual-structure
    // (they will, if created internally)

    Individual * ind(const int i)  const
    {
      return vars.size() == 0 ? NULL : vars[0].ind(i); 
    }

    
    /*!
      Attempt to add a variant to this group; if this 
      breaks the grouping rules specified in 'mask', 
      then set finished flag to true      
      @param v Reference to a variant object
    */

    void add(Variant &v );

    /*!
      Force an addition to the group, not taking note of 
      any grouping rules, and not setting any flags. Useful
      if one has an already-constructed set of Variants that
      one wants to treat as a VariantGroup. Also, no checking for 
      uniqueness.
      @param v Reference to a variant object
    */
    
    void force_add(Variant&v);


    /*!
      Test whether the last attempt to add a variant to this group
      failed, i.e. which indicates that the group is now complete, 
      as variants are sort by group
      @return True if the group is now complete 
    */
    
    bool complete() const { return is_complete; } 


    /*!
      Completely clear/reset the group
    */

    void clear();
    
    
    /*!
      Clear/reset the group, and start a new one with a single variant
      @param v Reference to a variant object
    */

    void clear(Variant & v);

    
    
    friend std::ostream & operator<<( std::ostream & out , const VariantGroup & v)
      {

	if ( v.size() == 0 ) 
	  {
	    out << "-NULL-\n";
	    return out;
	  }
	
	// Overall group details
	
	out << "NAME=[" << v.name() << "]\t" 	    
	    << "N(V)=" << v.size() << "\t"	    
	    << "N(I)=" << v.n_individuals() << "\n";
	
	// Per-variant summary
	
	for (int i=0; i<v.size(); i++)
	  {
	    out << "  " 
		<< v.var(i) << "\t" 
		<< v.var(i).meta << "\n";	    
	  }
		
	return out;
      }

    std::string dump(bool vmeta = true , 
		     bool vexpand = false , 
		     bool geno = true ,
		     bool gmeta = false , 
		     bool transpose = false , 
		     bool rarelist = false , 
		     bool show_phenotype = false );


 private:
    
    // Main variant store
    std::vector<Variant> vars;
    
    // Misc. properties of group
    std::string gname;
    int fileset;
    bool is_complete;

    //
    // References to key other structures
    //

    Mask & mask;

};



#endif
