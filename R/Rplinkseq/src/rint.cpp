#include "rint.h"
#include "plinkseq/annot.h"
#include <string>


GStore * gp;
extern GStore * GP;
bool R_project_attached;
  
std::string Rversion()
{
  std::string vmaj = R_MAJOR;
  std::string vmin = R_MINOR;  
  return vmaj + "-" + vmin;
}


void R_init_Rplinkseq(DllInfo *info)
{
  gp = new GStore;    
  gp->R_mode( true );
}


SEXP Rattach(SEXP d)
{
  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 
  std::string s = CHAR(STRING_ELT(d, 0));  
  gp->vardb_attach(s);
  return(R_NilValue);
}

void R_error( const std::string & s )
{
  Rf_error( s.c_str() ) ; 
} 

void R_warning( const std::string & s )
{
  Rf_warning( s.c_str() ) ; 
} 

void R_flush_plog()
{
  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return; } 
  Rprintf( "want to print more\n" );
  //  Rprintf( plog.R_flush().c_str() );
}

SEXP Rvardb_new(SEXP d)
{
  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 
  std::string s = CHAR(STRING_ELT(d, 0));  
  gp->vardb_new(s);
  return(R_NilValue);
}
 

SEXP Rdettach()
{
  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 
  gp->vardb_dettach();
  return(R_NilValue);
}


template<class T> 
SEXP Rmeta( const MetaInformation<T> & m )
{
  
  std::vector<std::string> metaKeys = m.keys();   
  int nm = metaKeys.size();
  
  SEXP vmlist_names;
  PROTECT(vmlist_names = allocVector( STRSXP, nm ));    
  
  SEXP vmlist;
  PROTECT(vmlist = allocVector( VECSXP, nm )); 
  
  for(int i = 0; i < nm; i++)   
    {
      
      SET_STRING_ELT(vmlist_names,i,mkChar( metaKeys[i].c_str() )); 
      
      mType mt = MetaInformation<T>::type( metaKeys[i] );
      
      SEXP vmval;
      
      switch (mt) {
      case META_INT :
	{		
	  int s = m.size( metaKeys[i] );
	  PROTECT(vmval = allocVector( INTSXP, s ));		
	  std::vector<int> d = m.get_int( metaKeys[i] );
	  for (int j =0; j<s; j++)
	    INTEGER(vmval)[j] = d[j];
	  break;
	}	    
      case META_FLOAT :
	{		
	  int s = m.size( metaKeys[i] );
	  PROTECT(vmval = allocVector( REALSXP, s ));		
	  std::vector<double> d = m.get_double( metaKeys[i] );
	  for (int j =0; j<s; j++)
	    REAL(vmval)[j] = d[j];
	  break;
	}	    
      case META_BOOL :
	{		
	  int s = m.size( metaKeys[i] );
	  PROTECT(vmval = allocVector( INTSXP, s ));		
	  std::vector<int> d = m.get_int( metaKeys[i] );
	  for (int j =0; j<s; j++)
	    INTEGER(vmval)[j] = (int)d[j];
	  break;
	}	    
      default :
	{
	  int s = m.size( metaKeys[i] );
	  PROTECT(vmval = allocVector( STRSXP , s ));		
	  std::vector<std::string> d = m.get_string( metaKeys[i] );
	  for (int j =0; j<s; j++)
	    SET_STRING_ELT( vmval,j,mkChar( d[j].c_str() ) );
	}
	
      }
      
      SET_VECTOR_ELT( vmlist, i, vmval ); 
      
      // Clean-up value
      UNPROTECT(1);
    }
  
  // Attach the vector of names
  
  setAttrib(vmlist, R_NamesSymbol, vmlist_names); 
  
  UNPROTECT(2);
  
  return( vmlist );
}



SEXP Rvariant_group(VariantGroup & v, Rdisplay_options & opt) 
{
  
  // Take a list of variants and make an R object that 
  // corresponds to that variant group
  

  //
  // Attach variants to group
  //
  

  SEXP vlist;
  PROTECT(vlist = allocVector( VECSXP, v.size() ));   
  for (int i=0; i<v.size(); i++)
    SET_VECTOR_ELT(vlist, i, Rvariant(v.var(i) , opt) ); 
  
  
  //
  // Header information
  //
  
  const int sz = 4;
  std::vector<std::string> names(sz);
  names[0] = "NAME"; 
  names[1] = "NV"; 
  names[2] = "NI"; 
  names[3] = "VAR"; 

  SEXP list_names;  
  PROTECT(list_names = allocVector(STRSXP,sz));    
  for(int i = 0; i < sz; i++)
    SET_STRING_ELT(list_names,i,mkChar(names[i].c_str())); 
  

 // 
 // Store variantgroup meta-information, and variant list in 
 // a final list
 //
 
 
 SEXP gname;
 PROTECT(gname = allocVector(STRSXP, 1));
 SET_STRING_ELT(gname, 0, mkChar( v.name().c_str() ) ); 

 SEXP n_tot;
 PROTECT(n_tot = allocVector(INTSXP, 1));
 INTEGER(n_tot)[0] = v.size();

 SEXP n_ind;
 PROTECT(n_ind = allocVector(INTSXP, 1));
 INTEGER(n_ind)[0] = v.n_individuals();

  
 // Attach values to a list

 SEXP list;
 PROTECT(list = allocVector( VECSXP, sz )); 
 SET_VECTOR_ELT(list, 0, gname  ); 
 SET_VECTOR_ELT(list, 1, n_tot  ); 
 SET_VECTOR_ELT(list, 2, n_ind  ); 
 SET_VECTOR_ELT(list, 3, vlist   ); 
 
 setAttrib(list, R_NamesSymbol, list_names); 
 
 UNPROTECT(6);
 
 return list;

}


SEXP Rvariant(Variant & v , Rdisplay_options & opt )
{

  int sz = 6;
  
  // Main list names

  SEXP list_names;
  
  std::vector<std::string> names(6);
    
  names[0] = "CHR";
  names[1] = "BP1"; 
  names[2] = "BP2";
  names[3] = "ID";
  names[4] = "NS";
  names[5] = "META";
  
  if ( opt.show_consensus ) 
    {
      names.push_back( "CON" );
      ++sz;
    }
  
  
  if ( opt.show_multi_sample ) 
    {
      std::set<int> fset = v.unique_files();
      std::set<int>::iterator i = fset.begin();
      while ( i != fset.end() )
	{
	  names.push_back( "S" + Helper::int2str( *i ) ); 
	  ++sz;
	  ++i; 
	}
    }
  
  
  
  // 
  // Store variantgroup meta-information, and variant list in 
  // a final list
  //
  

  // Make primary list object
  
  SEXP vlist;
  PROTECT(vlist = allocVector( VECSXP, sz ) );
  
  PROTECT(list_names = allocVector(STRSXP,sz));    
  for(int i = 0; i < sz; i++)   
    SET_STRING_ELT(list_names,i,mkChar(names[i].c_str())); 
  
  
  SEXP chr;
  PROTECT(chr = allocVector(INTSXP, 1));
  INTEGER(chr)[0] = v.chromosome();
  
  SEXP pos;
  PROTECT(pos = allocVector(INTSXP, 1));
  INTEGER(pos)[0] = v.position();
  
  SEXP stop;
  PROTECT(stop = allocVector(INTSXP, 1));
  INTEGER(stop)[0] = v.stop();
  
  SEXP vname;
  PROTECT(vname = allocVector(STRSXP, 1));
  SET_STRING_ELT(vname, 0, mkChar( v.name().c_str() ) ); 
  
  SEXP n_file;
  PROTECT(n_file = allocVector(INTSXP, 1));
  INTEGER(n_file)[0] = v.n_samples();
  
  SEXP metas;
  PROTECT( metas = Rmeta( v.meta ) );
  
  
  //
  // Attach values
  //
  
  SET_VECTOR_ELT(vlist, 0, chr    );
  SET_VECTOR_ELT(vlist, 1, pos    ); 
  SET_VECTOR_ELT(vlist, 2, stop   ); 
  SET_VECTOR_ELT(vlist, 3, vname  ); 
  SET_VECTOR_ELT(vlist, 4, n_file ); 
  SET_VECTOR_ELT(vlist, 5, metas  ); 
  
  
  //
  // Attach SampleVariants
  //
  
  int s = 6;
  
  // Consensus
  if ( opt.show_consensus ) 
    {
      SET_VECTOR_ELT(vlist, s++, Rsample_variant( -1 , v, opt ) );
    }
  
  
  // Individual sample variants
  
  if ( opt.show_multi_sample && v.multi_sample() )
    {
      // creating a problem?

      const int ns = v.n_samples();

      for (int sample = 0; sample < ns ; sample ++ )
	{	  

	  SET_VECTOR_ELT( vlist, 
 			  s++, 
 			  Rsample_variant( sample , v , opt ) );
	  	  
	}

    }
  
  setAttrib( vlist, R_NamesSymbol, list_names ); 
  
  UNPROTECT(8);

  return vlist; 
  
}
  


SEXP Rsample_variant( const int si , Variant & parent , Rdisplay_options & opt )
{
  
  SampleVariant & v = parent.sample( si );
  
  bool is_consensus = v.fileset() == 0 ; 
  
  const int nind = parent.size( si );
  
  // indicates we will need to cherry pick a sample from the consensus.
  const bool subset = parent.flat() && ! is_consensus;
  std::vector<int> imask;  
  if ( subset ) imask = parent.indiv_mask( v.fileset() );
  
  
  // Construct an R list object to represent this variant, it's
  // meta-information; its genotypes and their meta-information
  

  // Fileset ID
  
  SEXP file_id;
  PROTECT(file_id = allocVector(INTSXP, 1));
  INTEGER(file_id)[0] = v.fileset();
  
  // Reference allele

  SEXP ref;
  PROTECT(ref = allocVector(STRSXP, 1));
  SET_STRING_ELT(ref, 0, mkChar( v.reference().c_str() ) ); 
  
  // Alternate allele

  SEXP alt;
  PROTECT(alt = allocVector(STRSXP, 1));
  SET_STRING_ELT(alt, 0, mkChar( v.alternate().c_str() ) ); 
     
  // Quality score
    
  SEXP quality;
  PROTECT(quality = allocVector(REALSXP, 1));
  if ( v.quality() < 0 ) 
    REAL(quality)[0] = NA_REAL;
  else
    REAL(quality)[0] = v.quality();

  // Filter-strings 
  
  SEXP filter;
  std::vector< std::string > keys = v.meta_filter.keys();
  PROTECT(filter = allocVector(STRSXP, keys.size() ));  
  for (int k=0; k<keys.size(); k++)
    SET_STRING_ELT(filter, k, mkChar( v.meta_filter.get1_string( keys[k] ).c_str() ) ); 
  
  
  // Variant meta-information
  // (a list itself)

  std::vector< std::string > metaKeys = v.meta.keys();
  
  int nm = metaKeys.size();
  
  SEXP vmlist_names;
  PROTECT(vmlist_names = allocVector( STRSXP, nm ));    
  
  // Creating a value list with nm vector elements:

  SEXP vmlist;
  PROTECT(vmlist = allocVector( VECSXP, nm )); 
    
  for(int i = 0; i < nm; i++)   
    {
      
      // Add variant meta-field name label
      
      SET_STRING_ELT(vmlist_names,i,mkChar( metaKeys[i].c_str() )); 
      
      // Add appropriately-typed value
      
      mType mt = MetaInformation<VarMeta>::type( metaKeys[i] );
      
      SEXP vmval;
      
      switch (mt) {
      case META_INT :
	{		
	  int s = v.meta.size( metaKeys[i] );
	  PROTECT(vmval = allocVector( INTSXP, s ));		
	  std::vector<int> d = v.meta.get_int( metaKeys[i] );
	  for (int j =0; j<s; j++)
	    INTEGER(vmval)[j] = d[j];
	  break;
	}	    
      case META_FLOAT :
	{		
	  int s = v.meta.size( metaKeys[i] );
	  PROTECT(vmval = allocVector( REALSXP, s ));		
	  std::vector<double> d = v.meta.get_double( metaKeys[i] );
	  for (int j =0; j<s; j++)
	    REAL(vmval)[j] = d[j];
	  break;
	}	    
      case META_BOOL :
	{		
	  int s = v.meta.size( metaKeys[i] );
	  PROTECT(vmval = allocVector( INTSXP, s ));		
	  std::vector<int> d = v.meta.get_int( metaKeys[i] );
	  for (int j =0; j<s; j++)
	    INTEGER(vmval)[j] = (int)d[j];
	  break;
	}	    
      default :
	{
	  int s = v.meta.size( metaKeys[i] );
	  PROTECT(vmval = allocVector( STRSXP , s ));		
	  std::vector<std::string> d = v.meta.get_string( metaKeys[i] );
	  for (int j =0; j<s; j++)
	    SET_STRING_ELT( vmval,j,mkChar( d[j].c_str() ) );
	}
	
      }
      
      SET_VECTOR_ELT( vmlist, i, vmval ); 
      
      // Clean-up value
      UNPROTECT(1);
    }
  
  // Attach the vector of names
  
  setAttrib(vmlist, R_NamesSymbol, vmlist_names); 


  
  //
  // Genotypes, and genotype meta-information
  //
  
  SEXP geno;
  
  if ( opt.show_genotypes )
    {

      // Note -- we first have to tally up all the meta-fields for this 
      //         variant -- this is v. inefficient, but do for now
      
      
      // Number of genotype meta-fields (if showing any)
      
      std::set<std::string> mk;
      
      if ( opt.show_genotype_metainformation )
	{
	  
	  for (int i = 0 ; i < nind ; i++)
	    {
	      
	      const int idx = subset ? imask[i] : i ;
	      
	      Genotype * g = parent.genotype( si , idx );
	      
	      if ( g ) 
		{
		  std::vector<std::string> keys = g->meta.keys();
		  for (int k=0; k<keys.size(); k++) mk.insert( keys[k] ); 
		}
	    }
	  
	}

            
      // Geno :  ind1 ind2 ... indN
      // M1   :  ind1 ind2 ... indN
      // M2   :  ind1 ind2 ... indN
      // ...
      
      
      PROTECT(geno = allocVector( VECSXP , 1 + mk.size() ));
	
      //
      // Genotype labels for calls and meta-fields
      //
      
      SEXP glist_names;
      PROTECT(glist_names = allocVector( STRSXP, 1 + mk.size() ));    
      SET_STRING_ELT( glist_names , 0 , mkChar("GT" ) ) ;
      

      //
      // Actual genotype calls
      //
	
      SEXP g_calls;
      PROTECT(g_calls = allocVector( INTSXP, nind ));
      
      for (int i = 0 ; i < nind ; i++)
	{
	  
	  const int idx = subset ? imask[i] : i ;	 
	  
	  Genotype * g = parent.genotype( si , idx );
	  
	  if ( ! g ) 
	    {
	      INTEGER(g_calls)[i] = NA_INTEGER;
	    }
	  else
	    {
	      
	      // For now, just add allele-count code	
	      if ( g->more() || g->null() )
		{
		  INTEGER(g_calls)[i] = NA_INTEGER;
		}
	      else
		{
		  INTEGER(g_calls)[i] = g->allele_count( );
		}
	    }
	}
      
      
      SET_VECTOR_ELT( geno, 0, g_calls  ); 

      UNPROTECT(1); // g_calls
      

      //
      // Genotype meta-information
      //
      
      if ( opt.show_genotype_metainformation )
	{
	  
	  std::set<std::string>::iterator k = mk.begin();
	  
	  int x = 1;
	  
	  while ( k != mk.end() )
	    {

	      //
	      // Handles types:  int/bool -> integer
	      //                 double   -> float
	      //                 string   -> string
	      //  ? flag
	      
	      // Handles len:  1  = vector
	      //               >1 = matrix
	      //               -1 = list
	      
	      meta_index_t midx = MetaInformation<GenMeta>::field( *k );
	      
	      mType mt = midx.mt;	      
	      int len = midx.len;
	      
	      if ( mt == META_INT || mt == META_BOOL ) 
		{
		  SEXP g_meta;
		  
		  if ( len == 1 ) // vector
		    {
		      PROTECT(g_meta = allocVector( INTSXP, nind ));
		      for (int j=0; j<nind; j++) 
			{
			  const int idx = subset ? imask[j] : j ;
			  Genotype * g = parent.genotype( si , idx );
			  if ( g && g->meta.has_field( *k ) )
			    INTEGER(g_meta)[j] = g->meta.get1_int( *k ) ;
			  else
			    INTEGER(g_meta)[j] = NA_INTEGER;
			}
		      SET_VECTOR_ELT(geno, x , g_meta  ); 
		      UNPROTECT(1); // g_meta
		    }
		  else if ( len > 1 ) // matrix
		    {
		      PROTECT(g_meta = allocMatrix( INTSXP, nind , len ));
		      for (int j=0; j<nind; j++) 
			{
			  const int idx = subset ? imask[j] : j ;
			  Genotype * g = parent.genotype( si , idx );
			  if ( g && g->meta.has_field( *k ) )
			    {
			      std::vector<int> vec = g->meta.get_int( *k ) ;
			      for (int z=0; z<len; z++)
				INTEGER(g_meta)[j + nind * z] = vec[z];			      
			    }
			  else
			    {
			      for (int z=0; z<len; z++)
				INTEGER(g_meta)[j + nind * z] = NA_INTEGER;
			    }
			}		      
		      SET_VECTOR_ELT(geno, x , g_meta  ); 
		      UNPROTECT(1); // g_meta
		    }
		  else  // list (variable-numbered arrays)
		    {
		      PROTECT(g_meta = allocVector( VECSXP, nind ));
		      for (int j=0; j<nind; j++) 
			{
			  const int idx = subset ? imask[j] : j ;
			  Genotype * g = parent.genotype( si , idx );
			  if ( g && g->meta.has_field( *k ) )
			    {
			      std::vector<int> vec = g->meta.get_int( *k ) ;
			      SEXP list1;
			      PROTECT( list1 = allocVector( INTSXP , vec.size() ) );
			      for (int z=0; z<vec.size(); z++)
				INTEGER(list1)[z] = vec[z];	
			      SET_VECTOR_ELT( g_meta, j , list1 );
			      UNPROTECT(1); // list1
			    }
			}	      
		      SET_VECTOR_ELT(geno, x , g_meta  ); 
		      UNPROTECT(1); // g_meta
		    }
		}
	      else if ( mt == META_FLOAT ) 
		{
		  SEXP g_meta;

		  if ( len == 1 ) // vector
		    {
		      PROTECT(g_meta = allocVector( REALSXP, nind ));
		      for (int j=0; j<nind; j++) 
			{
			  const int idx = subset ? imask[j] : j ;
			  Genotype * g = parent.genotype( si , idx );
			  if ( g && g->meta.has_field( *k ) )
			    REAL(g_meta)[j] = g->meta.get1_double( *k ) ;
			  else
			    REAL(g_meta)[j] = NA_REAL;
			}
		      SET_VECTOR_ELT(geno, x , g_meta  ); 
		      UNPROTECT(1); // g_meta
		    }
		  else if ( len > 1 ) // matrix  
		    {
		      PROTECT(g_meta = allocMatrix( REALSXP, nind , len ));
		      for (int j=0; j<nind; j++) 
			{
			  const int idx = subset ? imask[j] : j ;
			  Genotype * g = parent.genotype( si , idx );
			  if ( g && g->meta.has_field( *k ) )
			    {
			      std::vector<double> vec = g->meta.get_double( *k ) ;
			      for (int z=0; z<len; z++)
				REAL(g_meta)[j + nind * z] = vec[z];			      
			    }
			  else
			    {
			      for (int z=0; z<len; z++)
				REAL(g_meta)[j + nind * z] = NA_REAL;
			    }
			}		      
		      SET_VECTOR_ELT(geno, x , g_meta  ); 
		      UNPROTECT(1);
		    }
		  else  // list
		    {
		      PROTECT(g_meta = allocVector( VECSXP, nind ));
		      for (int j=0; j<nind; j++) 
			{
			  const int idx = subset ? imask[j] : j ;
			  Genotype * g = parent.genotype( si , idx );
			  if ( g && g->meta.has_field( *k ) )
			    {
			      std::vector<double> vec = g->meta.get_double( *k ) ;
			      SEXP list1;
			      PROTECT( list1 = allocVector( REALSXP , vec.size() ) );
			      for (int z=0; z<vec.size(); z++)
				REAL(list1)[z] = vec[z];	
			      SET_VECTOR_ELT( g_meta, j , list1 );
			      UNPROTECT(1);
			    }
			}		      
		      SET_VECTOR_ELT(geno, x , g_meta  ); 
		      UNPROTECT(1);
		    }
		}
	      else 
		{
		  
		  SEXP g_meta;
		  
		  if ( len == 1 ) // vector
		    {
		      PROTECT(g_meta = allocVector( STRSXP, nind ));
		      for (int j=0; j<nind; j++) 
			{
			  const int idx = subset ? imask[j] : j ;
			  Genotype * g = parent.genotype( si , idx );
			  if ( g && g->meta.has_field( *k ) )
			    SET_STRING_ELT( g_meta, j , mkChar( g->meta.get1_string( *k ).c_str() ) );
			  else
			    SET_STRING_ELT( g_meta, j , NA_STRING );			  
			}
		      SET_VECTOR_ELT(geno, x , g_meta  ); 
		      UNPROTECT(1);
		    }
		  else if ( len > 1 ) // matrix
		    {
		      PROTECT(g_meta = allocMatrix( STRSXP, nind , len ));
		      for (int j=0; j<nind; j++) 
			{
			  const int idx = subset ? imask[j] : j ;
			  Genotype * g = parent.genotype( si , idx );
			  if ( g && g->meta.has_field( *k ) )
			    {
			      std::vector<std::string> vec = g->meta.get_string( *k ) ;
			      for (int z=0; z<len; z++)
				SET_STRING_ELT( g_meta, j+nind*z , mkChar( vec[z].c_str() ) );			      
			    }
			  else
			    {
			      for (int z=0; z<len; z++)
				SET_STRING_ELT( g_meta, j+nind*z , NA_STRING );			      			      
			    }
			}		      
		      SET_VECTOR_ELT(geno, x , g_meta  ); 
		      UNPROTECT(1);
		    }
		  else  // list
		    {
		      PROTECT(g_meta = allocVector( VECSXP, nind ));
		      for (int j=0; j<nind; j++) 
			{
			  const int idx = subset ? imask[j] : j ;
			  Genotype * g = parent.genotype( si , idx );
			  if ( g && g->meta.has_field( *k ) )
			    {
			      std::vector<std::string> vec = g->meta.get_string( *k ) ;
			      SEXP list1;
			      PROTECT( list1 = allocVector( STRSXP , vec.size() ) );
			      for (int z=0; z<vec.size(); z++)
				SET_STRING_ELT( list1 , z , mkChar( vec[z].c_str() ) );	
			      SET_VECTOR_ELT( g_meta, j , list1 );
			      UNPROTECT(1);
			    }
			}		      
		      SET_VECTOR_ELT(geno, x , g_meta  ); 
		      UNPROTECT(1);
		    }
		}

	      
	      SET_STRING_ELT( glist_names , x , mkChar( (*k).c_str() ) ) ;
	      

	      //
	      // Next genotype meta-field
	      //

	      ++x;	
	      ++k;
	    }
	}
      
      setAttrib(geno, R_NamesSymbol, glist_names); 

    }
  

    //
    // Make final variant list:
    //


    SEXP list_names;
    int sz = 6;    
    std::vector< std::string > names(6);
    names[0] = "FSET";
    names[1] = "REF";
    names[2] = "ALT";
    names[3] = "QUAL";
    names[4] = "FILTER";
    names[5] = "META";

    if ( opt.show_genotypes ) 
      {
	++sz;
	names.push_back( "GENO" );
      }
    
    PROTECT(list_names = allocVector(STRSXP,sz));    
    
    for(int i = 0; i < sz; i++)   
      SET_STRING_ELT(list_names,i,mkChar(names[i].c_str())); 
    
  
    // Creating a value list with sz vector elements:
    
    SEXP list;
    
    PROTECT(list = allocVector( VECSXP, sz )); 
    
    // Attach values
    
    SET_VECTOR_ELT(list, 0,  file_id ); 
    SET_VECTOR_ELT(list, 1,  ref     ); 
    SET_VECTOR_ELT(list, 2,  alt     ); 
    SET_VECTOR_ELT(list, 3,  quality ); 
    SET_VECTOR_ELT(list, 4,  filter  ); 
    SET_VECTOR_ELT(list, 5,  vmlist  ); 
    if ( opt.show_genotypes )
      SET_VECTOR_ELT(list, 6,  geno    ); 
    
    setAttrib(list, R_NamesSymbol, list_names); 
    

    //
    // Free up protected resources
    //
    
    UNPROTECT(9);

    if ( opt.show_genotypes ) 
      UNPROTECT(2); // geno, glist_names    
    
    return list;
    
}


struct Rdata {
  Rdisplay_options * opt;
  SEXP fncall;
  SEXP rho;
};


bool R_filter_func( SEXP fn )
{
  Helper::halt("Not implemented R filter functions yet");
}


//
// R iterate functions()
//


void R_accumulate_func( Variant & v , void * p )
{      
  R_CheckUserInterrupt();
  std::vector<Variant> * d = (std::vector<Variant>*)p;
  d->push_back(v);
}

void R_group_accumulate_func( VariantGroup & v , void * p )
{
  R_CheckUserInterrupt();
  std::vector<VariantGroup> * d = (std::vector<VariantGroup>*)p;
  d->push_back(v);
  VariantGroup vcopy = v;
}

void R_iterate_func( Variant & v , void * p )
{
  R_CheckUserInterrupt();

  Rdata * d = (Rdata*)p;
  
  // Call R function, with variant and with result object
  
  // Bind appropriate parameters
  
  SEXP var;
  
  PROTECT( var = Rvariant( v , *(d->opt) ) );
  
  SETCADR( d->fncall, var );
  

  // Evalute function

  eval(d->fncall, d->rho);

  
  UNPROTECT(1);    
  
}


void R_group_iterate_func( VariantGroup & v , void * p )
{
  
  R_CheckUserInterrupt();

  Rdata * d = (Rdata*)p;
  
  // Call R function, with variant and with result object
  
  // Bind appropriate parameters
  
  SEXP var;
  
  PROTECT( var = Rvariant_group( v , *(d->opt)) );
  
  SETCADR( d->fncall, var );
  
  // Evalute function
  
  eval(d->fncall, d->rho);
  
  UNPROTECT(1);    
  
}


SEXP R_getListElement(SEXP list, char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

std::map<std::string,std::vector<std::string> > R_getListMap(SEXP list)
{
  std::map<std::string,std::vector<std::string> > m;
  SEXP names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++)
    { 
      std::string key = CHAR(STRING_ELT(names, i));	
      SEXP elmt = VECTOR_ELT(list, i);
      for (int j = 0 ; j < length(elmt); j++ )
	m[key].push_back( CHAR(STRING_ELT(names, i)) );	
    }
  return m;    
}

std::set<std::string> R_getListNames(SEXP list)
{
  SEXP names = getAttrib(list, R_NamesSymbol);
  std::set<std::string> n;
  for (int i = 0; i < length(list); i++)
    n.insert( CHAR(STRING_ELT(names, i)) );
  return n;
}




Mask R_make_mask(SEXP r)
{
  
  //
  // Use PLINK/SEQ textual specification to construct mask
  //
  
  if ( length(r) < 1 ) return Mask("");    
  std::string s = CHAR(STRING_ELT(r, 0));  	

  Mask m(s);

  gp->register_mask( m );
  return m;
}


struct R_aux_getmeta 
{ 
  std::vector<std::string> * nmeta_int;
  std::vector<std::string> * nmeta_float;
  std::vector<std::string> * nmeta_string;  
  std::vector<std::string> * nmeta_filter;  
  std::set<std::string> * nmeta_special;  

  std::vector<std::vector<int> > * dmeta_int;
  std::vector<std::vector<double> > * dmeta_float;
  std::vector<std::vector<std::string> > * dmeta_string;
  std::vector<std::vector<bool> > * dmeta_filter;

  std::vector<std::vector<bool> > * mmeta_int;
  std::vector<std::vector<bool> > * mmeta_float;
  std::vector<std::vector<bool> > * mmeta_string;
  std::vector<std::vector<bool> > * mmeta_filter;
  int ne;

};


void R_aux_getmeta_func( Variant & var , void * p )
{
  
  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return; } 

  R_aux_getmeta * d = (R_aux_getmeta*)p;
  
  d->ne++;

  for (int i=0;i< d->nmeta_int->size(); i++)
    {

      if ( d->nmeta_special->find( (*d->nmeta_int)[i] ) != d->nmeta_special->end() )
	{
	  // must be POS
	  (*d->dmeta_int)[i].push_back( var.position() );
	  (*d->mmeta_int)[i].push_back( true );
	}
      else
	{
	  MetaInformation<VarMeta> * m = var.meta.has_field( (*d->nmeta_int)[i] ) ? & var.meta 
	    : ( var.consensus.meta.has_field( (*d->nmeta_int)[i]) ? &var.consensus.meta : NULL ) ;
	  if ( m ) 
	    {
	      (*d->dmeta_int)[i].push_back( m->get1_int( (*d->nmeta_int)[i] ) );
	      (*d->mmeta_int)[i].push_back( true );
	    }
	  else
	    (*d->mmeta_int)[i].push_back( false );
	}
    }

  for (int i=0;i< d->nmeta_float->size();i++)
    {

      if ( d->nmeta_special->find( (*d->nmeta_float)[i] ) != d->nmeta_special->end() )
	{
	  // must be QUAL
	  (*d->dmeta_float)[i].push_back( var.consensus.quality() );
	  (*d->mmeta_float)[i].push_back( true );
	}
      else
	{
	  MetaInformation<VarMeta> * m = var.meta.has_field( (*d->nmeta_float)[i] ) ? & var.meta 
	    : ( var.consensus.meta.has_field( (*d->nmeta_float)[i]) ? &var.consensus.meta : NULL ) ;
	  
	  if ( m )
	    {
	      (*d->dmeta_float)[i].push_back( m->get1_double( (*d->nmeta_float)[i] ) );
	      (*d->mmeta_float)[i].push_back( true );
	    }
	  else
	    (*d->mmeta_float)[i].push_back( false );
	}
    }

  for (int i=0;i<d->nmeta_string->size();i++)
    {

      if ( d->nmeta_special->find( (*d->nmeta_string)[i] ) != d->nmeta_special->end() )
	{
	  if ( (*d->nmeta_string)[i] == "REF" )
	    (*d->dmeta_string)[i].push_back( var.reference() );
	  else if ( (*d->nmeta_string)[i] == "ALT" )
	    (*d->dmeta_string)[i].push_back( var.alternate() );
	  else if ( (*d->nmeta_string)[i] == "CHROM" )
	    (*d->dmeta_string)[i].push_back( Helper::chrCode( var.chromosome() ) );
	  else // must be ID
	    (*d->dmeta_string)[i].push_back( var.name() );

	  (*d->mmeta_string)[i].push_back( true );
	}
      else
	{
	  MetaInformation<VarMeta> * m = var.meta.has_field( (*d->nmeta_string)[i] ) ? & var.meta 
	    : ( var.consensus.meta.has_field( (*d->nmeta_string)[i]) ? &var.consensus.meta : NULL ) ;
	  
	  if ( m )
	    {
	      (*d->dmeta_string)[i].push_back( m->get1_string( (*d->nmeta_string)[i] ) );
	      (*d->mmeta_string)[i].push_back( true );
	      
	    }
	  else
	    (*d->mmeta_string)[i].push_back( false );
	}
    }


  // FILTERs

  for (int i=0;i< d->nmeta_filter->size(); i++)
    (*d->dmeta_filter)[i].push_back( var.consensus.meta_filter.has_field( (*d->nmeta_filter)[i]) );
  
}


SEXP Rgetmeta(SEXP meta, SEXP rmask, SEXP rho)
{

  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 

  // Which variant meta-fields to return?

  int n = length(meta);
  if ( n == 0 ) return R_NilValue;  
  
  std::vector<std::string> nmeta_int;
  std::vector<std::string> nmeta_float;
  std::vector<std::string> nmeta_string;
  std::vector<std::string> nmeta_filter;
  std::set<std::string> nmeta_special;
  
  std::vector<std::vector<int> > dmeta_int;
  std::vector<std::vector<double> > dmeta_float;
  std::vector<std::vector<std::string> > dmeta_string;
  std::vector<std::vector<bool> > dmeta_filter;

  std::vector<std::vector<bool> > mmeta_int;
  std::vector<std::vector<bool> > mmeta_float;
  std::vector<std::vector<bool> > mmeta_string;
  std::vector<std::vector<bool> > mmeta_filter;


  for (int i=0;i<n;i++)
    {
      std::string m = CHAR( STRING_ELT(meta,i));
      mType mt = MetaInformation<VarMeta>::type( m );
      
      if ( mt == META_INT ) 
	{
	  nmeta_int.push_back( m );
	  dmeta_int.resize( nmeta_int.size() );
	  mmeta_int.resize( nmeta_int.size() );
	}
      else if ( mt == META_FLOAT )
	{
	  nmeta_float.push_back( m );
	  dmeta_float.resize( nmeta_float.size() );	  
	  mmeta_float.resize( nmeta_float.size() );
	}
      else if ( mt == META_TEXT )
	{
	  nmeta_string.push_back( m );
	  dmeta_string.resize( nmeta_string.size() );
	  mmeta_string.resize( nmeta_string.size() );
	}
      else
	{
	  // Is this a filter?
  
	  mType mt = MetaInformation<VarFilterMeta>::type( m );
 	  if ( mt != META_UNDEFINED )
 	    {
	      nmeta_filter.push_back( m );
	      dmeta_filter.resize( nmeta_filter.size() );	      
 	    }
	  else
	    {
	      // A special code from VCF header?
	      
	      if ( m == "POS" ) 
		{
		  nmeta_int.push_back(m);		 
		  nmeta_special.insert(m);
		  dmeta_int.resize( nmeta_int.size() );
		  mmeta_int.resize( nmeta_int.size() );
		}
	      else if ( m == "QUAL" )
		{
		  nmeta_float.push_back(m);
		  nmeta_special.insert(m);
		  dmeta_float.resize( nmeta_float.size() );
		  mmeta_float.resize( nmeta_float.size() );
		}
	      else if ( m == "CHROM" || m == "REF" || m == "ALT" )
		{
		  nmeta_string.push_back(m);
		  nmeta_special.insert(m);
		  dmeta_string.resize( nmeta_string.size() );
		  mmeta_string.resize( nmeta_string.size() );
		}
	    }
 	}
    }


  // Mask
  
  Mask mask;
  if ( length(rmask) > 0 ) 
    mask = R_make_mask(rmask);
  
  // Set up auxiliary objects
  
  R_aux_getmeta aux;
  aux.nmeta_int = &nmeta_int;
  aux.nmeta_float = &nmeta_float;
  aux.nmeta_string = &nmeta_string;
  aux.nmeta_filter = &nmeta_filter;
  aux.nmeta_special = &nmeta_special;
  
  aux.dmeta_int = &dmeta_int;
  aux.dmeta_float = &dmeta_float;
  aux.dmeta_string = &dmeta_string;
  aux.dmeta_filter = &dmeta_filter;
  
  aux.mmeta_int = &mmeta_int;
  aux.mmeta_float = &mmeta_float;
  aux.mmeta_string = &mmeta_string;

  aux.ne = 0;

  //
  // Iterate
  //

  gp->vardb.iterate( R_aux_getmeta_func , &aux , mask );
  

  //
  // Generate R return values
  //

  // total # of cols being returned:
  n = nmeta_int.size() + nmeta_float.size() + nmeta_string.size() + nmeta_filter.size();

  SEXP ret, retnames;
  
  PROTECT(ret = allocVector(VECSXP, n)); 
  PROTECT(retnames = allocVector(STRSXP, n)); 
  
  int le = 0;

  /* fill list and names */ 
  
  for (int i=0;i<nmeta_filter.size();i++)
    {
      SEXP l1;
      PROTECT( l1 = allocVector( INTSXP , aux.ne ) );
      int p = 0;
      for (int j=0;j<aux.ne;j++)
	INTEGER(l1)[j] = dmeta_filter[i][p++];  // no missing codes for filter
      SET_STRING_ELT( retnames, le , mkChar( nmeta_filter[i].c_str() ) );
      SET_VECTOR_ELT( ret , le++ , l1 );
      UNPROTECT(1); // l1
    }

  for (int i=0;i<nmeta_int.size();i++)
    {
      SEXP l1;
      PROTECT( l1 = allocVector( INTSXP , aux.ne ) );
      int p = 0;
      for (int j=0;j<aux.ne;j++)
	{
	  if ( mmeta_int[i][j] )
	    INTEGER(l1)[j] = dmeta_int[i][p++];
	  else
	    INTEGER(l1)[j] = NA_INTEGER;
	}
      SET_STRING_ELT( retnames, le , mkChar( nmeta_int[i].c_str() ) );
      SET_VECTOR_ELT( ret , le++ , l1 );
      UNPROTECT(1); // l1
    }


  for (int i=0;i<nmeta_float.size();i++)
    {
      SEXP l1;
      PROTECT( l1 = allocVector( REALSXP , aux.ne ) );
      int p = 0;
      for (int j=0;j<aux.ne;j++)
	{
	  if ( mmeta_float[i][j] )
	    REAL(l1)[j] = dmeta_float[i][p++];
	  else
	    REAL(l1)[j] = NA_REAL;
	}
      SET_STRING_ELT( retnames, le , mkChar( nmeta_float[i].c_str() ) );
      SET_VECTOR_ELT( ret , le++ , l1 );
      UNPROTECT(1); // l1
    }


  for (int i=0;i<nmeta_string.size();i++)
    {
      SEXP l1;
      PROTECT( l1 = allocVector( STRSXP , aux.ne ) );
      int p = 0;
      for (int j=0;j<aux.ne;j++)
	{
	  if ( mmeta_string[i][j] )
	    SET_STRING_ELT( l1 , j , mkChar( dmeta_string[i][p++].c_str() ) );
	  else
	    SET_STRING_ELT( l1 , j , NA_STRING );
	}
      SET_STRING_ELT( retnames, le , mkChar( nmeta_string[i].c_str() ) );
      SET_VECTOR_ELT( ret , le++ , l1 );
      UNPROTECT(1); // l1
    }
  
  setAttrib( ret , R_NamesSymbol, retnames );		

  UNPROTECT(2);
  return ret;
}




SEXP Riterate(SEXP fn, SEXP rmask, SEXP ret, SEXP rho)
{    
  
  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 
  
  // 
  // Do we want to accumulate and return the list of variants?
  //
  
  // ret = 0  implies we will iterate over a function, and not explicitly 
  //          return any values
  
  // ret > 0  implies we will rturn up to 'ret' values (i.e. make sure we do 
  //          not explode memory in R
  
  bool return_vars = false;
  
  if ( length(ret)>0 && INTEGER(ret)[0] != 0 ) 
    return_vars = true;
  
  if ( ! return_vars )
    {
      if ( ! isFunction(fn) ) 
	Helper::halt("'fn' must be a function");
    }
  
  if ( ! isEnvironment(rho)) 
    Helper::halt("'rho' should be an environment");
  


  //
  // Set-up R function call
  //
  
  SEXP R_fcall, ans;
  
  PROTECT(R_fcall = lang2(fn, 
			  R_NilValue)); 
  
  // Object to pass to C/C++ function, that points
  // to the R function and an R object that can 
  // collect any results


  Rdisplay_options opt;
  
  Rdata * d = new Rdata;
  d->opt = &opt;
  d->fncall = R_fcall;
  d->rho    = rho;

  std::vector<Variant> vars;
  std::vector<VariantGroup> varGroups;
  
  //
  // Construct mask
  //
  
  Mask mask;
  
  if ( length(rmask) > 0 ) 
    mask = R_make_mask(rmask);
  

  if ( length(ret) > 0 && INTEGER(ret)[0] != 0 ) 
    mask.limit( INTEGER(ret)[0]);
  
  

  //
  // Set up individual-map
  //
  
  gp->indmap.populate( gp->vardb , gp->phmap , mask );



  //
  // Iterate over all variants
  //
  
  if ( return_vars ) 
    {
      if ( mask.any_grouping() )	  
	gp->vardb.iterate( R_group_accumulate_func , &varGroups , mask );
      else
	gp->vardb.iterate( R_accumulate_func , &vars , mask );
      
    }
  else if ( mask.any_grouping() )
    {
      gp->vardb.iterate( R_group_iterate_func , d , mask );
    }
  else
    {
      gp->vardb.iterate( R_iterate_func , d , mask );
    }    
  


  //
  // Clean up...
  //
  
  UNPROTECT(1);        
  
  delete d;  
  
  
  //
  // Return list of variants, or nothing, as specified
  //
  
  if ( return_vars )
    {
      
      // Return list of groups, or list of single variants?
      
      if ( mask.any_grouping() )
	{
	  SEXP rvars;
	  PROTECT( rvars = allocVector( VECSXP, varGroups.size() ));
	  for (int j=0; j<varGroups.size(); j++)
	    {
	    SET_VECTOR_ELT( rvars , j , Rvariant_group( varGroups[j] , opt ) );
	    }
	  UNPROTECT(1);
	  return(rvars);
	}
      else // return list of single variants
	{
	  SEXP rvars;
	  PROTECT( rvars = allocVector( VECSXP, vars.size() ));
	  for (int j=0; j<vars.size(); j++)
	    SET_VECTOR_ELT( rvars , j , Rvariant( vars[j] , opt ) );
	  UNPROTECT(1);
	  return(rvars);
	}
    }
  else
    return(R_NilValue);
  
}



SEXP Rhdr_list_int(int f)
{

  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 

  std::vector<std::map<std::string,std::string> > hdrs = gp->vardb.fetch_headers( f );    
  SEXP rhdrs;
  PROTECT( rhdrs = allocVector( VECSXP, hdrs.size() ));
  
  for (int i = 0 ; i < hdrs.size() ; i++)
    {
      // Expect either 1 or two key->value pairs here
      SEXP rhdr;
      PROTECT( rhdr = allocVector( STRSXP,  hdrs[i].size()) );
      
      std::map<std::string,std::string>::iterator j = hdrs[i].begin();
      int k = 0;
      while ( j != hdrs[i].end() )
	{	    
	  std::string key = j->first;
	  std::string value = j->second;
	  SET_STRING_ELT( rhdr ,       k , mkChar( value.c_str() ) );	    
	  ++k;
	  ++j;
	}
	
      // Pull altogether in a list
      SET_VECTOR_ELT( rhdrs , i , rhdr );	
      UNPROTECT(1);		
    }
  
    // Names and list
    UNPROTECT(1);    
    return(rhdrs);
}


SEXP Rmeta_list_int( int f )
{

  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 

  std::vector<std::map<std::string,std::string> > metas = GP->vardb.fetch_metatypes( f );
  
  SEXP rmetas;
  
  PROTECT( rmetas = allocVector( VECSXP, metas.size() ));

  for (int i = 0 ; i < metas.size() ; i++)
    {
      
      SEXP rmetas2;
      PROTECT( rmetas2 = allocVector( VECSXP, metas[i].size() ));
      
      SEXP lname;
      PROTECT( lname = allocVector( STRSXP, metas[i].size() ) );
      
      int j = 0;
      std::map<std::string,std::string>::iterator k = metas[i].begin();
      while ( k != metas[i].end() )
	{	    
	  std::string key = k->first;
	  std::string value = k->second;

	  SEXP rmeta;
	  PROTECT( rmeta = allocVector( STRSXP, 1) );
	  SET_STRING_ELT( rmeta, 0 , mkChar( value.c_str() ) );
	  
	  SET_STRING_ELT( lname, j , mkChar( key.c_str() ) );
	  
	  // Pull altogether in a list
	  SET_VECTOR_ELT( rmetas2 , j , rmeta );
	  
	  UNPROTECT(1);	
	  
	  ++j;
	  ++k;
	}
      
      setAttrib(rmetas2 , R_NamesSymbol, lname );		
      SET_VECTOR_ELT( rmetas , i , rmetas2 );
      UNPROTECT(2);	
      
    }
  
  // Names and list
  UNPROTECT(1);    
  return(rmetas);  
  
}


SEXP Rind_list(SEXP rmask, SEXP phe)
{

  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 
  
  // Given a mask, and a list of phenotype names, return the individuals (as
  // they will be ordered given a var.fetch() or var.iterate()
  
  Mask mask;
  
  if ( length(rmask) > 0 ) 
    mask = R_make_mask(rmask);
  
  //
  // Set up individual-map
  //
  
  gp->indmap.populate( gp->vardb , gp->phmap , mask );
  
  const int n = gp->indmap.size();
  
  const int p = length(phe);
    
  SEXP indiv;
  PROTECT( indiv = allocVector( VECSXP, p+1 ) );
  
  // Labels
  
  SEXP list_names;
  PROTECT(list_names = allocVector( STRSXP, p+1 ));    
  SET_STRING_ELT(list_names, 0 , mkChar( "ID" ) );
  
    
  //
  // IDs
  //
  
  SEXP idstr;
  PROTECT( idstr = allocVector( STRSXP, n ) );
  for (int i = 0 ; i < n ; i++)
    SET_STRING_ELT( idstr , i , mkChar( gp->indmap.ind(i)->id().c_str() ) );   
  SET_VECTOR_ELT( indiv , 0 , idstr );
  
  
  //
  // Phenotypes
  //
  
  for (int j = 0 ; j < p ; j++)
    {
      
      std::string pname = CHAR( STRING_ELT(phe, j) );
      
      SET_STRING_ELT(list_names, j+1 , STRING_ELT(phe, j) );
      
      
      //
      // Try to load this phenotype
      //
      
      if ( gp->phmap.set_phenotype( CHAR(STRING_ELT(phe,j ) ) ) == 0 ) continue;
	

      //
      // Determine type
      //
      
      mType mt = MetaInformation<IndivMeta>::type( pname );

      if ( mt == META_UNDEFINED ) continue;
      
      if ( mt == META_BOOL )
	{
	  SEXP phe1;
	  PROTECT( phe1 = allocVector( INTSXP, n ) );
	  for (int i=0; i<n; i++)
	    {
	      if ( ! gp->indmap.ind(i)->meta.has_field( pname ) )
		INTEGER(phe1)[i] = NA_INTEGER;
	      else 
		INTEGER(phe1)[i] = gp->indmap.ind(i)->meta.get1_bool( pname );
	    }
	  SET_VECTOR_ELT( indiv, j+1 , phe1 );
	  UNPROTECT(1);	    
	}
      else if ( mt == META_INT )
	{
	  SEXP phe1;
	  PROTECT( phe1 = allocVector( INTSXP, n ) );
	  for (int i=0; i<n; i++)
	    {
	      if ( ! gp->indmap.ind(i)->meta.has_field( pname ) )
		INTEGER(phe1)[i] = NA_INTEGER;
	      else 
		INTEGER(phe1)[i] = gp->indmap.ind(i)->meta.get1_int( pname );
	    }
	  SET_VECTOR_ELT( indiv, j+1 , phe1 );
	  UNPROTECT(1);	    
	}
      else if ( mt == META_FLOAT )
	{
	  SEXP phe1;
	  PROTECT( phe1 = allocVector( REALSXP, n ) );
	  for (int i=0; i<n; i++)
	    {
	      if ( ! gp->indmap.ind(i)->meta.has_field( pname ) )
		REAL(phe1)[i] = NA_REAL;
	      else 
		REAL(phe1)[i] = gp->indmap.ind(i)->meta.get1_double( pname );
	    }
	  SET_VECTOR_ELT( indiv, j+1 , phe1 );
	  UNPROTECT(1);
	  
	}
      else if ( mt == META_TEXT )
	{
	  SEXP phe1;
	  PROTECT( phe1 = allocVector( STRSXP, n ) );
	  for (int i=0; i<n; i++)
	    {
	      if ( ! gp->indmap.ind(i)->meta.has_field( pname ) )
		SET_STRING_ELT(phe1, i , NA_STRING );
	      else 
		SET_STRING_ELT(phe1, i , mkChar( gp->indmap.ind(i)->meta.get1_string( pname ).c_str() ) );
	    }
	  SET_VECTOR_ELT( indiv, j+1 , phe1 );
	  UNPROTECT(1);
	  
	}
      
    }
  
  
  
  // Add labels
  
  setAttrib( indiv, R_NamesSymbol, list_names );
      
  UNPROTECT(3);	
  
  return(indiv);
}



SEXP Rind_pedlist( SEXP rmask )
{
  
  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 
  
  // Given a mask, return the individuals (as they will be ordered
  // given a var.fetch() or var.iterate()) an their associated
  // pedigree info
  
  Mask mask;
  
  if ( length(rmask) > 0 ) 
    mask = R_make_mask(rmask);
  
  //
  // Set up individual-map
  //
  
  gp->indmap.populate( gp->vardb , gp->phmap , mask );
  
  const int n = gp->indmap.size();
  
  // ID, FID, IID, PAT, MAT, SEX
  
  SEXP indiv;
  PROTECT( indiv = allocVector( VECSXP, 6 ) );
  
  // Labels
  
  SEXP list_names;
  PROTECT(list_names = allocVector( STRSXP, 6 ));    
  SET_STRING_ELT(list_names, 0 , mkChar( "ID" ) );
  SET_STRING_ELT(list_names, 1 , mkChar( "FID" ) );
  SET_STRING_ELT(list_names, 2 , mkChar( "IID" ) );
  SET_STRING_ELT(list_names, 3 , mkChar( "PAT" ) );
  SET_STRING_ELT(list_names, 4 , mkChar( "MAT" ) );
  SET_STRING_ELT(list_names, 5 , mkChar( "SEX" ) );
  
    
  //
  // IDs
  //
  
  SEXP idstr;
  PROTECT( idstr = allocVector( STRSXP, n ) );
  for (int i = 0 ; i < n ; i++)
    SET_STRING_ELT( idstr , i , mkChar( gp->indmap.ind(i)->id().c_str() ) );   
  SET_VECTOR_ELT( indiv , 0 , idstr );
  
  SEXP fidstr;
  PROTECT( fidstr = allocVector( STRSXP, n ) );
  for (int i = 0 ; i < n ; i++)
    SET_STRING_ELT( fidstr , i , mkChar( gp->indmap.ind(i)->fid().c_str() ) );   
  SET_VECTOR_ELT( indiv , 1 , fidstr );

  SEXP iidstr;
  PROTECT( iidstr = allocVector( STRSXP, n ) );
  for (int i = 0 ; i < n ; i++)
    SET_STRING_ELT( iidstr , i , mkChar( gp->indmap.ind(i)->iid().c_str() ) );   
  SET_VECTOR_ELT( indiv , 2 , iidstr );

  SEXP patstr;
  PROTECT( patstr = allocVector( STRSXP, n ) );
  for (int i = 0 ; i < n ; i++)
    SET_STRING_ELT( patstr , i , mkChar( gp->indmap.ind(i)->father().c_str() ) );   
  SET_VECTOR_ELT( indiv , 3 , patstr );

  SEXP matstr;
  PROTECT( matstr = allocVector( STRSXP, n ) );
  for (int i = 0 ; i < n ; i++)
    SET_STRING_ELT( matstr , i , mkChar( gp->indmap.ind(i)->mother().c_str() ) );   
  SET_VECTOR_ELT( indiv , 4 , matstr );

  SEXP sexint;
  PROTECT( sexint = allocVector( INTSXP, n ) );
  for (int i = 0 ; i < n ; i++)
    {
      if ( gp->indmap.ind(i)->sex() == 1 ) 
	INTEGER(sexint)[i] = 1;
      else if ( gp->indmap.ind(i)->sex() == 2 ) 
	INTEGER(sexint)[i] = 2;
      else
	INTEGER(sexint)[i] = NA_INTEGER;
    }
  SET_VECTOR_ELT( indiv , 5 , sexint );


  //
  // Add labels
  //

  setAttrib( indiv, R_NamesSymbol, list_names );
      
  UNPROTECT(8);	
  
  return(indiv);
}




SEXP Rfile_list()
{
  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 
    
  // Files contain
  //  1) Index ID
  //  2) Name
  //  3) Headers
  //  4) Meta-types
  
  std::map<int,std::string> files = gp->vardb.fetch_files();
  
  SEXP rfiles;
  PROTECT( rfiles = allocVector( VECSXP, files.size() ));
  
  int j = 0;
  std::map<int,std::string>::iterator i = files.begin();
  
  while ( i != files.end() ) 
    {
      SEXP rfile;
      
      PROTECT( rfile = allocVector( VECSXP, 4) );
      
      // ID
      SEXP tmp = allocVector(INTSXP, 1);
      INTEGER(tmp)[0] = i->first;
      SET_VECTOR_ELT( rfile, 0 , tmp );
      
      // Name	
      SEXP fname;
      PROTECT(fname = allocVector(STRSXP, 1));
      SET_STRING_ELT(fname, 0, mkChar( i->second.c_str() ) ); 
      SET_VECTOR_ELT( rfile , 1 , fname );
      UNPROTECT(1);
      
      // Headers
      SET_VECTOR_ELT( rfile , 2 , Rhdr_list_int( i->first ) );
      
      // Meta-information
      SET_VECTOR_ELT( rfile , 3 , Rmeta_list_int( i->first ) );

      //
      // Set list names
      //
      
      SEXP list_names;
      PROTECT(list_names = allocVector( STRSXP, 4 ));    
      SET_STRING_ELT(list_names, 0 , mkChar( "IDX" ) );
      SET_STRING_ELT(list_names, 1 , mkChar( "FNAME" ) );
      SET_STRING_ELT(list_names, 2 , mkChar( "HDR" ) );
      SET_STRING_ELT(list_names, 3 , mkChar( "META" ) );
      setAttrib(rfile, R_NamesSymbol, list_names );	
      
      // Pull altogether in a list
      SET_VECTOR_ELT( rfiles , j , rfile );
      
      // Clear list and names
      UNPROTECT(2);	
      
      ++j;
      ++i;
    }
  
  
  UNPROTECT(1);
  
  return(rfiles);
}


SEXP Rhdr_list(SEXP f) 
{ 
  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 
  return Rhdr_list_int( INTEGER(f)[0] ); 
}

SEXP Rmeta_list(SEXP f) 
{ 
  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 
  if ( length(f) != 1 ) return R_NilValue;
  std::string fn = CHAR(STRING_ELT( f , 0) );    
  int fi = gp->vardb.file_tag( fn );
  return fi ? Rmeta_list_int( fi ) : R_NilValue ; 
}



////////////////////////////////////////////////////// 
//                                                  //
// Locus database interface                         //
//                                                  //
////////////////////////////////////////////////////// 



SEXP Rvardb_attach(SEXP filename) 
{ 
  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 
  gp->vardb.attach( CHAR(STRING_ELT(filename, 0)) );    
  return(R_NilValue);
}

SEXP Rvardb_dettach() 
{ 
  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 
  gp->vardb.dettach(); 
  return(R_NilValue);
}


SEXP Rset_project(SEXP n)
{
  std::string proj_name = CHAR(STRING_ELT(n, 0));
  R_project_attached = gp->set_project( FileMap::tilde_expansion( proj_name ) );
  if ( R_project_attached ) Rprintf("successfully attached project\n");
  else Rprintf( ("failed: could not attach " + proj_name + "\n").c_str() );
  return( R_project_attached ? R_NilValue : R_NilValue );
}

SEXP Rsummary()
{

  if ( ! R_project_attached ) { plog.warn( "no project attached" ); return( R_NilValue); } 

  // split on line return, no empty lines
  
  std::vector<std::string> s = Helper::char_split( gp->summary(false) , '\n' , false );
  
  SEXP summ;
  PROTECT( summ = allocVector( STRSXP, s.size() ));  
  for (int i=0; i<s.size(); i++)
    {
      // tmp fix
      std::string t = Helper::search_replace( s[i] , "\t" , "  " );
      SET_STRING_ELT( summ , i , mkChar( t.c_str() ) );
    }
  // TODO: make into a list object to return instead of print to STDOUT
  //  R_flush_plog();
  
  UNPROTECT(1);
  return( summ ) ;  
}


SEXP Rdirect_load_vcf( SEXP rfilename , SEXP m , SEXP retn )
{
  
  if ( length(rfilename) != 1 ) return(R_NilValue);    
  if ( length(m) != 1 ) return(R_NilValue);    
  if ( length(retn) != 1 ) return(R_NilValue);    
  
  std::string filename = CHAR(STRING_ELT(rfilename, 0)); 
  filename = FileMap::tilde_expansion( filename );
  std::string mask = CHAR(STRING_ELT(m, 0)); 
  int ret = INTEGER(retn)[0];
  
  // This will make iterate() function read directly from the VCF, not
  // from the VARDB
  
  // Create a temporary GStore, so that everything here is
  // self-contained and does not impact any other attached projects
  
  GStore gtmp;
  GStore * orig = gp;
  gp = GP = &gtmp;
  
  gtmp.single_file_mode( true );
  gtmp.vardb.attach( ":memory:" );
  gtmp.inddb.attach( ":memory:" );
  
  mask += " ex-vcf=" + filename;
  if ( ret != 0 ) mask += " limit=" + Helper::int2str( ret );

  // Variants, first populate vars[] 
  Mask pmask( mask );

  // Reader headers (necessary? we repeat the process here)
  gtmp.register_mask( pmask );

  gtmp.vardb.vcf_iterate_read_header( pmask );

  
  // construct an object that contains
  // a) VCF meta-fields
  // b) IDs from header
  // c) list of variants (with any filters applied)
  
  // Simply get first rows from file directly
  std::vector<std::string> meta_hdrs;
  std::vector<std::string> indiv_ids;
  
  InFile vcf( filename );
  while ( 1 ) 
    {
      if ( vcf.eof() ) break;
      std::string l = vcf.readLine();      
      if ( l == "" ) continue;      
      if ( l.substr(0,2) == "##" ) meta_hdrs.push_back( l );
      else if ( l.substr(0,4) == "#CHR" ) 
	{
	  std::vector<std::string> tok = Helper::char_split( l , '\t' );
	  for (int i=9; i<tok.size(); i++) indiv_ids.push_back( tok[i] );
	}
      else break;  // end of headers
    }
  vcf.close();
  

  // List of 3 lists
  //  -- meta-headers
  //  -- IDs
  //  -- list of variants
    
  // Meta-headers
  SEXP mhdr;
  PROTECT( mhdr = allocVector( STRSXP,meta_hdrs.size() ));
  for (int i=0; i<meta_hdrs.size(); i++)
    SET_STRING_ELT( mhdr, i , mkChar( meta_hdrs[i].c_str() ) );
  
  // IDs
  SEXP idstr;
  PROTECT( idstr = allocVector( STRSXP, indiv_ids.size() ) );
  for (int i = 0 ; i < indiv_ids.size() ; i++)
    SET_STRING_ELT( idstr , i , mkChar( indiv_ids[i].c_str() ) );
    

  // Get variants
  std::vector<Variant> vars;   
  IterationReport rep = gtmp.vardb.iterate( R_accumulate_func , &vars , pmask );    
  plog << rep.report() << "\n";


  // Compile as R list
  Rdisplay_options opt;
  SEXP rvars;
  PROTECT( rvars = allocVector( VECSXP, vars.size() ));
  for (int j=0; j<vars.size(); j++)
    SET_VECTOR_ELT( rvars , j , Rvariant( vars[j] , opt ) );
  
  SEXP vcflist;
  PROTECT( vcflist = allocVector( VECSXP, 3 ) );
  SET_VECTOR_ELT( vcflist , 0 , mhdr );
  SET_VECTOR_ELT( vcflist , 1 , idstr );
  SET_VECTOR_ELT( vcflist , 2 , rvars );
    
  SEXP vcflist_names;
  PROTECT(vcflist_names = allocVector( STRSXP, 3 ));  
  SET_STRING_ELT(vcflist_names, 0 , mkChar( "META" ) );
  SET_STRING_ELT(vcflist_names, 1 , mkChar( "ID" )   );	    
  SET_STRING_ELT(vcflist_names, 2 , mkChar( "VAR" )  );  
  setAttrib( vcflist, R_NamesSymbol, vcflist_names );

  // release
  UNPROTECT( 5 );

  gp = GP = orig;
  
  return vcflist;


}


//
// Make a set of variants, given a set of regions
//

SEXP Rvardb_make_set(SEXP id, SEXP name)
{
//   std::string grp = CHAR(STRING_ELT(id, 0));
//   gp->vardb_make_set( grp , CHAR(STRING_ELT(name, 0)) );    
  return(R_NilValue);  
}


//
// Load GTF file
//


SEXP Rlocdb_attach(SEXP name)
{
    gp->locdb.attach( CHAR(STRING_ELT(name, 0)) );    
    return(R_NilValue);
}

SEXP Rlocdb_load_gtf( SEXP file, SEXP name)
{
    gp->locdb.load_GTF( CHAR(STRING_ELT(file, 0)) , CHAR(STRING_ELT(name, 0)) );   
    gp->locdb.index();
    return(R_NilValue); 
}

SEXP Rlocdb_load_alias( SEXP file , SEXP name )
{
//     string grp = CHAR(STRING_ELT(id, 0));
//     uint64_t grp_id = gp->locdb.lookup_group_id( grp );    
//     gp->locdb.
    return(R_NilValue); 
}

SEXP Rlocdb_collapse_subregions( SEXP id , SEXP name )
{
    gp->locdb.merge( CHAR(STRING_ELT(id, 0)), CHAR(STRING_ELT(name, 0)) );   
    return(R_NilValue); 
}

SEXP Rfetch_set_names( SEXP x, SEXP y)
{
  std::string loc_group = CHAR(STRING_ELT(x, 0));
  std::string set_group = CHAR(STRING_ELT(y, 0));
  std::vector<std::string> r = gp->locdb.fetch_set_names( loc_group, set_group );
  return Rmake_string_vector( r );
}

SEXP Rmake_string_vector( std::vector<std::string> & r )
{
  SEXP v;
  PROTECT( v = allocVector( STRSXP, r.size()) );
  for (int i=0;i<r.size(); i++)
    SET_STRING_ELT(v, i, mkChar( r[i].c_str() ) );     
  UNPROTECT(1);
  return v;
}

SEXP Rfetch_set_members( SEXP x, SEXP y, SEXP z)
{
  std::string loc_group = CHAR(STRING_ELT(x, 0));
  std::string set_group = CHAR(STRING_ELT(y, 0));
  std::string set_name = CHAR(STRING_ELT(z, 0));
  std::vector<std::string> r = gp->locdb.fetch_set_members( loc_group, set_group , set_name );
  return Rmake_string_vector( r );
}


SEXP Rfetch_regions(SEXP g)
{

  std::string grp = CHAR(STRING_ELT(g, 0));
  uint64_t grp_id = gp->locdb.lookup_group_id( grp );    
  std::set<Region> reg = gp->locdb.get_regions( grp_id );

  // Construct similar R objects
  
  SEXP rregs;
  PROTECT( rregs = allocVector( VECSXP, reg.size()) );
    
  int j = 0;
  std::set<Region>::iterator i = reg.begin();

  while ( i != reg.end() )
    {
      
      const Region & r = *i;
      
      // ID
      // NAME
      // CHR
      // BP1, BP2
      // GROUP
      // META
      // List of sub-regions
      
      SEXP rreg;
      PROTECT( rreg = allocVector( VECSXP, 8 ) );
	       
      // Basic region informaiton: 

      SEXP id;
      PROTECT(id = allocVector(INTSXP, 1));
      INTEGER(id)[0] = r.id;
      
      SEXP vname;
      PROTECT(vname = allocVector(STRSXP, 1));
      SET_STRING_ELT(vname, 0, mkChar( r.name.c_str() ) ); 
      
      SEXP chr;
      PROTECT(chr = allocVector(INTSXP, 1));
      INTEGER(chr)[0] = r.start.chromosome();
      
      SEXP bp1;
      PROTECT(bp1 = allocVector(INTSXP, 1));
      INTEGER(bp1)[0] = r.start.position();
	
      SEXP bp2;
      PROTECT(bp2 = allocVector(INTSXP, 1));
      INTEGER(bp2)[0] = r.stop.position();
	
      SEXP group;
      PROTECT(group = allocVector(INTSXP, 1));
      INTEGER(group)[0] = r.group;
	

      //
      // Meta information
      //

      SEXP metas;
      PROTECT( metas = Rmeta( r.meta ) );

      // Subregions?
	
      SEXP subregs;
      PROTECT( subregs = allocVector( VECSXP, r.subregion.size()) );

      for (int s = 0 ; s <  r.subregion.size(); s++ )
 	{
	  
	  SEXP rsreg;
	  PROTECT( rsreg = allocVector( VECSXP, 5 ) );
	  
	  // Basic region informaiton: 
	  
	  SEXP sid;
	  PROTECT(sid = allocVector(INTSXP, 1));
	  INTEGER(sid)[0] = r.subregion[s].id;
	
	  SEXP svname;
	  PROTECT(svname = allocVector(STRSXP, 1));
	  SET_STRING_ELT(svname, 0, mkChar( r.subregion[s].name.c_str() ) ); 
	     
	  SEXP sbp1;
	  PROTECT(sbp1 = allocVector(INTSXP, 1));
	  INTEGER(sbp1)[0] = r.subregion[s].start.position();
	    
	  SEXP sbp2;
	  PROTECT(sbp2 = allocVector(INTSXP, 1));
	  INTEGER(sbp2)[0] = r.subregion[s].stop.position();
	    
	  SEXP smetas;
	  PROTECT( smetas = Rmeta( r.subregion[s].meta ) );


	  SET_VECTOR_ELT( rsreg , 0 , sid );
	  SET_VECTOR_ELT( rsreg , 1 , svname );
	  SET_VECTOR_ELT( rsreg , 2 , sbp1 );
	  SET_VECTOR_ELT( rsreg , 3 , sbp2 );
	  SET_VECTOR_ELT( rsreg , 4 , smetas );
	    
	  SEXP slist_names;
	  PROTECT(slist_names = allocVector( STRSXP, 5 ));    
	  SET_STRING_ELT(slist_names, 0 , mkChar( "IDX" ) );
	  SET_STRING_ELT(slist_names, 1 , mkChar( "NAME" ) );	    
	  SET_STRING_ELT(slist_names, 2 , mkChar( "BP1" ) );
	  SET_STRING_ELT(slist_names, 3 , mkChar( "BP2" ) );
	  SET_STRING_ELT(slist_names, 4 , mkChar( "META" ) );
	  setAttrib(rsreg, R_NamesSymbol, slist_names );
	    
	  SET_VECTOR_ELT( subregs , s , rsreg ); 

	  UNPROTECT( 7 );

 	}	


	// Add labels

        SEXP list_names;
        PROTECT(list_names = allocVector( STRSXP, 8 ));    
        SET_STRING_ELT(list_names, 0 , mkChar( "IDX" ) );
        SET_STRING_ELT(list_names, 1 , mkChar( "NAME" ) );
        SET_STRING_ELT(list_names, 2 , mkChar( "CHR" ) );
        SET_STRING_ELT(list_names, 3 , mkChar( "BP1" ) );
        SET_STRING_ELT(list_names, 4 , mkChar( "BP2" ) );
        SET_STRING_ELT(list_names, 5 , mkChar( "GROUP" ) );
        SET_STRING_ELT(list_names, 6 , mkChar( "META" ) );
        SET_STRING_ELT(list_names, 7 , mkChar( "SUBREG" ) );
        setAttrib(rreg, R_NamesSymbol, list_names );

	// Add elements
	
	SET_VECTOR_ELT( rreg , 0 , id     );
	SET_VECTOR_ELT( rreg , 1 , vname  );
	SET_VECTOR_ELT( rreg , 2 , chr    );
	SET_VECTOR_ELT( rreg , 3 , bp1    );
	SET_VECTOR_ELT( rreg , 4 , bp2    );
	SET_VECTOR_ELT( rreg , 5 , group  );
	SET_VECTOR_ELT( rreg , 6 , metas  );
	SET_VECTOR_ELT( rreg , 7 , subregs );
	
	SET_VECTOR_ELT( rregs , j , rreg );

	UNPROTECT(10);
	
	++j;
	++i;
    }

    UNPROTECT(1);

    return( rregs );
    
}


SEXP Rlocdb_summary()
{
  plog << gp->locdb.summary(false);
  return(R_NilValue);
}


////////////////////////////////////////////////////// 
//                                                  //
// Reference database                               //
//                                                  //
////////////////////////////////////////////////////// 


SEXP Rseqdb_loadFASTA( SEXP filename )
{
  if ( length(filename) != 1 ) return(R_NilValue);    
  std::string s = CHAR(STRING_ELT(filename, 0));  
  std::map<std::string,std::string> m;
  gp->seqdb.loadFASTA( s , m );   
  return(R_NilValue); 
}

SEXP Rseqdb_attach( SEXP filename )
{
  if ( length(filename) != 1 ) return(R_NilValue);    
  std::string s = CHAR(STRING_ELT(filename, 0));  
  gp->seqdb.attach( s );   
  Annotate::init();
  return(R_NilValue); 
}

SEXP Rseqdb_lookup( SEXP pos )
{
  if ( ! isVector( pos ) ) return(R_NilValue);
  if ( length( pos ) != 3 ) return(R_NilValue);
  
  int chr = INTEGER(pos)[0];
  int bp1 = INTEGER(pos)[1];
  int bp2 = INTEGER(pos)[2];
  
  SEXP seq;
  PROTECT( seq = allocVector( STRSXP, 1 ));    
  SET_STRING_ELT(seq, 0 , mkChar( gp->seqdb.lookup(chr,bp1,bp2).c_str() ) );
  UNPROTECT(1);
  
  return seq;
}

SEXP Rseqdb_annotate_load( SEXP loc_id )
{
  if ( length( loc_id ) != 1 ) return(R_NilValue);
  std::string trans_id = CHAR(STRING_ELT(loc_id, 0));    
  Annotate::setDB( LOCDB );
  if ( ! Annotate::set_transcript_group( trans_id ) ) plog.warn( "trouble attaching 'refseq' group from LOCDB" );    
  return(R_NilValue);
}

SEXP Rseqdb_annotate( SEXP pos , SEXP alleles )
{
  
  if ( length( pos ) != 2 ) return(R_NilValue);
  if ( length( alleles ) != 2 ) return(R_NilValue);
  
  int chr = INTEGER(pos)[0];
  int bp = INTEGER(pos)[1];
  
  std::string ref = CHAR(STRING_ELT(alleles, 0));  
  std::string alt = CHAR(STRING_ELT(alleles, 1));  
    
  Variant v( "Var" , chr, bp );

  v.consensus.reference( ref );
  v.consensus.alternate( alt );
  
  Annotate::annotate( v );
    
  return(R_NilValue);
  
}

SEXP Rrefdb_load(SEXP x)
{
  return(R_NilValue); 
}


SEXP Rrefdb_attach(SEXP x)
{
  std::string s = CHAR(STRING_ELT(x, 0));  
  gp->refdb_attach(s);
  return(R_NilValue); 
}

SEXP Rrefdb_summary()
{
  plog << gp->refdb.summary(false);
  return(R_NilValue); 
}

SEXP Rrefdb_lookup(SEXP pos, SEXP y)
{    
  if ( length( pos ) != 2 ) return(R_NilValue);
  if ( length( y ) != 1 ) return(R_NilValue);
  
  int chr = INTEGER(pos)[0];
  int bp = INTEGER(pos)[1];
  
  Variant v( "Var" , chr, bp );

  std::string g = CHAR(STRING_ELT(y, 0));
  int gid = gp->refdb.lookup_group_id( g );
  if ( gid == 0 ) return(R_NilValue);

  RefVariant rv = gp->refdb_lookup( v, gid );

  if ( rv.observed() )
    plog << rv << "\t" << rv.meta << "\n";
  else
    plog << "NA\n";

  return(R_NilValue); 
  
}


SEXP Rrefdb_index_lookup(SEXP x)
{    
  return(R_NilValue); 
}

