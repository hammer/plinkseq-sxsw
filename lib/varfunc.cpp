#include "plinkseq/gstore.h"
#include "plinkseq/varfunc.h"


std::vector<bool> VarFunc::missing_genotype_mask( const Variant & v )
{
  std::vector<bool> b( v.size() , false );
  for (int i=0; i<v.size(); i++) if ( v(i).null() ) b[i] = true;
  return b;
}


std::vector<bool> VarFunc::missing_genotype_mask( const VariantGroup & vars )
{
  // case-wise missingness
  const int n = vars.n_individuals();
  std::vector<bool> b( n , false );
  for (int v=0; v<vars.size(); v++)
    {
      const Variant & var = vars(v);
      for (int i=0; i<n; i++) 
	if ( var(i).null() ) b[i] = true;
    }
  return b;
}


Data::Vector<double> VarFunc::alternate_allele_count( const Variant & v )
{
  Data::Vector<double> d( v.size() );
  for (int i=0; i<v.size(); i++)
    {
      if ( v(i).null() ) 
	{
	  d.set_elem_mask( i );
	}
      else
	d(i) = v(i).allele_count( );
    }
  return d;
}
