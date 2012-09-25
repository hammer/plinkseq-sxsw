
#include "views.h"
#include "plinkseq/em.h"
#include "plinkseq/dose.h"
#include "util.h"
#include "func.h"

#include "plinkseq.h"

#include <iostream>
#include <cmath>

extern GStore g;
extern Pseq::Util::Options args;

void f_view( Variant & v , void * p )
{
  
  Out & pout = Out::stream( "vars" );
  
  OptVView * opt = (OptVView*)p;
  
  if ( opt->simple ) 
  {
    pout << v.coordinate() << "\t" 
	 << v.consensus.reference() << "\t"
	 << v.consensus.alternate() << "\n";
    return;
  }
  
  // do not show now: accumulate for a ind x var display of (rare) genotypes:
  if ( opt->mview ) 
    {
      opt->vars.force_add( v );
      return;
    }
  
  pout << v << "\t"
       << v.name() << "\t"
       << v.consensus << "\t"      
       << "." << "\t" 
       << v.n_samples() << "\t";

  pout << v.print_meta_filter( ";" );

  if ( opt->vmeta )
    {
      if ( opt->vexpand ) pout << "\n" << v.meta.display();
      else pout << "\t" << v.meta ;      

      if ( opt->vexpand ) pout << "\n" << v.consensus.meta.display();
      else pout << "\t" << v.consensus.meta ;      
    }

  pout << "\n";
  
  //
  // Show separaetly specific sample information? 
  //
  
  if ( opt->show_samples )
    {
      
      const int ns = v.n_samples();
      
      for (int s = 0 ; s < ns; s++ )
	{

	  const SampleVariant & sample = v.sample(s);
	  
	  pout << v << "\t"
	       << v.name() << "\t"
	       << sample << "\t"	       
	       << sample.file_name() << "\t"	    
	       << sample.filter();
	  
	  if ( opt->vmeta ) 
	    {
	      if ( opt->vexpand ) pout << "\n" << sample.meta.display();
	      else pout << "\t" << sample.meta ;
	    }
	  
	  pout << "\n";
	  
	}      
    }
  

  //
  // Show genotypes and genotype meta-information? One row per individual
  //
  
  if ( opt->geno )
    {

      bool altmin = true;

      if ( opt->show_only_minor ) 
	altmin = v.n_minor_allele();
		
      const int n = v.size();
      
      for ( int i=0; i<n; i++)
	{

	  // show null and/or major genotypes?

	  if ( (!opt->show_nonmissing_geno) && v(i).null() ) continue;
	  if ( opt->show_only_alt   && ! v(i).nonreference() ) continue;
	  if ( opt->show_only_minor && v(i).minor_allele_count( altmin ) <= 0 ) continue;
	  
	  pout << v.ind( i )->id() << "\t";
	  
	  pout << v.sample_label( i, "," ) << "\t";
	  
	  //
	  // Optionally, phenotype
	  //

	  if ( g.phmap.type() == PHE_DICHOT )
	    {
	      if ( v.ind(i)->affected() == CASE ) pout << "CASE\t";
	      else if ( v.ind(i)->affected() == CONTROL ) pout << "CONTROL\t";
	      else pout << ".\t";
	    }
	  else if ( g.phmap.type() == PHE_QT )
	    pout << v.ind(i)->qt() << "\t";
	  else if ( g.phmap.type() == PHE_FACTOR )
	    pout << v.ind(i)->group_label() << "\t";

	  //
	  // Genotype
	  //

	  pout << v.label( i , "," );

	  //
	  // Gentype meta-information
	  //
	  
	  if ( opt->gmeta )
	    pout << "\t[" << v.gmeta_label(i) << "]";

	  
	  pout << "\n";
	}
    }



}


void g_view( VariantGroup & vars , void * p )
{

  OptGView * opt = (OptGView*)p;
  Out & pout = Out::stream( "groups" );  

  pout << vars.dump( opt->vmeta , 
		     opt->vexpand , 
		     opt->geno , 
		     opt->gmeta , 
		     opt->transpose , 
		     opt->rarelist , 
		     opt->show_phenotype ) << "\n";
  
}


void f_view_tped( Variant & v , void * p )
{

  Out & tped = Out::stream( "tped" );
  
  tped << v.chromosome() << "\t";
  
  if ( v.name() == "." )
    tped << "var_" 
	 << Helper::chrCode( v.chromosome() ) << "_" 
	 << v.position() << "\t";
  else
    tped << v.name() << "\t";
  
  tped  << "0\t"
	<< v.position();
  
  int n = v.size();
  
  for (int i=0;i<n;i++)
    tped << "\t" << v.print_PED( v(i) );

  tped << "\n";

}

void f_view_lik( Variant & v , void * p )
{
  
  //
  // output genotype likelihoods from GENO_LIK field ("GL" or "PL" by default
  // in GATK, specified in lib/defs.*)
  //

  const double EPS = 1e-5;

  GStore * g = (GStore*)p;
  
  int n = v.size();

  Out & pout = Out::stream( "lik" );
  
  pout << v << "\t" 
       << v.reference() << "\t"
       << v.alternate();
  
  for (int i = 0 ; i < n ; i++)
    {
      
      MetaInformation<GenMeta> & genotype = v(i).meta;
      
      std::vector<double> l;
      
      // Extract either PL (or GL, if no PL)

      bool using_gl = false;
      
      if ( genotype.has_field( PLINKSeq::META_GENO_PHRED() ) )
	{

	  std::vector<int> t = genotype.get_int( PLINKSeq::META_GENO_PHRED() );
	
	  if ( t.size() == 3 ) 
	    {
	      l.resize(3,0);
	      for (int i=0;i<3;i++) l[i] = t[i];
	    }
	}
      else if ( genotype.has_field( PLINKSeq::META_GENO_LIK() ) )
	{
	  l = genotype.get_double( PLINKSeq::META_GENO_LIK() );
	  using_gl = true;
	}
      
      

      // Currently, only handle bi-allelic cases


      if ( l.size() == 3 ) 
	{

	  // subtract out maximum log-likelihood, so that exp(max)==1
	  // for a decent scaling (we can set really small values to 0
	  // then and avoid over/under flow)
	  
	  double g0 = l[0];
	  double g1 = l[1];
	  double g2 = l[2];
	  
	  if ( ! using_gl ) 
	    {
	      g0 = -g0/10.0;
	      g1 = -g1/10.0;
	      g2 = -g2/10.0;
	    }
	  
	  double mx = g0 > g1 ? g0 : g1 ;
	  mx = g2 > mx ? g2 : mx;
	  
	  g0 -= mx;
	  g1 -= mx;
	  g2 -= mx;
	  
	  g0 = exp(g0);
	  g1 = exp(g1);
	  g2 = exp(g2);
	  
	  if ( ! Helper::realnum( g0 ) ) g0 = 0;
	  if ( ! Helper::realnum( g1 ) ) g1 = 0;
	  if ( ! Helper::realnum( g2 ) ) g2 = 0;
	  
	  double s = 1 / ( g0 + g1 + g2 );
	  
	  g0 *= s;
	  g1 *= s;
	  g2 *= s;
	  
	  // For a slightly smaller file, set all trivial likelihoods
	  // to exactly 0 (?okay)
	  
	  if ( g0 < EPS ) g0 = 0;
	  if ( g1 < EPS ) g1 = 0;
	  if ( g2 < EPS ) g2 = 0;
	  
	  if ( ! ( Helper::realnum( g0 ) && Helper::realnum( g1 ) && Helper::realnum( g2 ) ) )
	    {
	      pout << "\t0.333\t0.333\t0.333";
	    }
	  else
	    pout << "\t" << g0 
		 << "\t" << g1 
		 << "\t" << g2;
	}   
      else
	{
	  
	  // do we have an observed genotype?  if so, make a dummy likelihood
	  const Genotype & g = v(i);
	  

	  // flat ratio ~ missing data point
	  if ( g.null() ) 	
	    {	      
	      pout << "\t0.333\t0.333\t0.333";
	    }
	  else
	    {
	      int ac = g.allele_count( );
	      if ( ac == 0 ) pout << "\t1\t0\t0";
	      else if ( ac == 1 ) pout << "\t0\t1\t0";
	      else if ( ac == 2 ) pout << "\t0\t0\t1"; 
	      else pout << "\t0.333\t0.333\t0.333";
	    }
	}	    

    }
  
  pout << "\n"; 
   
}



void f_view_matrix( Variant & v , void * p )
{


  Out & pout = Out::stream( "matrix" );
  
  // output simple 0/1/2/NA calls 
  
  pout << v << "\t" 
       << v.reference() << "\t"
       << v.alternate();
  
  int n = v.size(); 

  for (int i = 0 ; i < n ; i++)
    {
	if ( v(i).null() ) 
	    pout << "\tNA";
	else
	    pout << "\t" << v(i).score();
    }
  pout << "\n";
  
}


void f_view_meta_matrix( Variant & v , void * p )
{

  Out & pout = Out::stream( "matrix" );

  std::map<int,std::string> * samples = (std::map<int,std::string>*)p;
  

  // Core variant information
  
  pout << v.chromosome() << "\t" 
       << v.position() << "\t"
       << v.name();
  
  
  int nelem_static = v.consensus.meta.n_visible_keys_static();
  
  int nelem_nonstatic = v.consensus.meta.n_visible_keys_nonstatic();

  // Need to separate out 'static' (Variant) and 'non-static'
  // (SampleVariant/Consensus) information
  
  // Static, population/variant-level meta-fields (if any)
  if ( nelem_static ) pout << "\t" << v.meta.display_row_static("NA");
  
  // Consensus information

  if ( nelem_nonstatic ) pout << "\t" << v.consensus.meta.display_row_nonstatic("NA");  

  // Filters, if only a single sample
  if ( samples->size() == 1 ) 
    {
      std::string f = v.consensus.meta_filter.display_row("NA");
      if ( f != "" ) pout << "\t" << f;
    }

  // Per-sample information

  if ( samples->size() > 1 )
    {
        
      // Samples (even if only 1 sample, loop through once and dump
      // meta_filter info from SampleVariant
            
      std::map<int,std::string>::iterator i = samples->begin(); 
      
      while ( i != samples->end() )
	{
	  
	  // note; NULL is returned if unique_svar_slot() return -1 (not seen)
	  
	  SampleVariant * s = ( ! v.multi_sample() ) ? &(v.consensus) :  v.psample( v.unique_svar_slot( i->first ) );
	  
	  if ( ! s ) 
	    {
	      
	      // If sample for this variant not present, insert missing values
	      
	      int nelem = MetaInformation<VarMeta>::n_visible_keys_nonstatic();
	      for (int i=0; i<nelem; i++) pout << "\tNA";
	      
	      nelem = MetaInformation<VarFilterMeta>::n_visible_keys();
	      for (int i=0; i<nelem; i++) pout << "\tNA";
	      
	    }
	  else
	    {
	      
	      std::string m = s->meta.display_row_nonstatic( "NA" );
	      std::string f = s->meta_filter.display_row("NA"); 

	      // Variant INFO tags
	      if ( m != "" ) pout << "\t" << m;
	      
	      // Variant FILTER tags
	      if ( f != "" ) pout << "\t" << f;
	      
	    }
	  ++i;
	}
    }

  // next row
  pout << "\n";
}



void f_view_var_meta_matrix( Variant & v , void * p )
{

  Out & pout = Out::stream( "matrix" );

  std::string name = *(std::string*)p; 
  pout << v << "\t";
  for (int i=0; i<v.size(); i++)
    {
      pout << ( v(i).meta.has_field(name) ? 
		v(i).meta.as_string(name) : 
		"NA" ) << "\t";
    }
  pout << "\n";
}



void f_view_gene_matrix( VariantGroup & vars , void * p )
{

  Out & pout = Out::stream( "matrix" );
  
  OptGMatrix * opt = (OptGMatrix*)p;
  
  // for each gene
  
  const int nv = vars.size();
  
  std::vector<bool> altmin( vars.size() , false );
  for (int v = 0 ; v < nv ; v++)
    altmin[v] = vars(v).n_minor_allele();
  
  std::string s = vars.name() + "\t" + Helper::int2str( nv );  
  bool variation = false;
  const int ni = g.indmap.size();

  for (int i = 0 ; i < ni ; i++)
    {

      bool valid = false;
      int cnt = 0;
      
      for (int v = 0 ; v < nv ; v++)
	{
	  // need to see at least 1 non-null variant
	  if ( ! vars(v,i).null() ) 
	    valid = true;
	  
	  // we are searching for minor alleles:
	  if ( vars.geno(v,i).minor_allele( altmin[v] ) )
	    {
	      cnt += vars.geno(v,i).minor_allele_count( altmin[v] );
	      variation = true;
	    }
	}
      
      if ( valid ) 
	s += "\t" + Helper::int2str( opt->collapse_01 ? ( cnt ? 1 : 0 ) : cnt  );
      else
	s += "\tNA";
    }
  
  if( ! ( opt->hide_zero_variance && ! variation ) )
    {
      pout << s;    
      pout << "\n";
    }

}



void f_view_gene_meta_matrix( VariantGroup & vars , void * p )
{
  
  Out & pout = Out::stream( "matrix" );
  
  // Output a region x individual matrix in which a genotype-level meta-field
  // is summarised for each individual across the region
  
  OptGMetaMatrix * opt = (OptGMetaMatrix*)p;
  
  // for each gene
  
  const int nv = vars.size();
  
  bool valid = true;
  
  // for now, just assume mean.. 

  if ( ! opt->show_mean ) Helper::halt("only mean option enabled for now");
  
  std::string s = vars.name() + "\t" + Helper::int2str( nv );  
  
  bool variation = false;
  const int ni = vars.n_individuals();
  
  for (int i = 0 ; i < ni ; i++)
    {

      bool valid = false;

      double numer = 0 , denom = 0;
      int cnt = 0;
      
      for (int v = 0 ; v < nv ; v++)
	{	  
	  if ( vars.geno(v,i).meta.has_field( opt->name ) )
	    {
	      valid = true;
	      ++denom;
	      numer += vars.geno(v,i).meta.get1_int( opt->name ) ;
	    }
	}
      

      /// Calculate 

      double value = 0;

      if ( opt->show_mean && valid ) 
	value = numer / denom; 
      
      /// Display 
      
      if ( valid ) 
 	s += "\t" + Helper::dbl2str( value ) ;
      // 	s += "\t" + ( opt->show_mean || opt->show_flag ) ? Helper::int2str( value ) : value ;
      else
	s += "\tNA";
      
    }
    
  pout << s << "\n";

}


std::map<XQC_varid,std::set<int> >::iterator XQCstats::flush( std::map<XQC_varid,std::set<int> >::iterator n )
{
  
  Out & pout = Out::stream( "vfreq" );

  XQC_varid pos = n->first;
  std::vector<std::string> & d = data[pos];
  std::vector<int> & k = datak[pos];
  
  if ( d.size() > 0 )
    {
      
      pout.data_group( d[0] );
      for (int i=1; i<d.size(); i++)
	pout.data( d[i] , k[i]);      
      
      pout.data( neighbours[ pos ].size() ); 
	
      // end of data for this variant
      pout.print_data_group();    
    }
    
  data.erase( data.find( pos ) );
  datak.erase( datak.find( pos ) );
  
  neighbours.erase( n++ );
  return n;
} 


void XQCstats::flush() 
{
  
  Out & pout = Out::stream( "vfreq" );
  
  std::map<XQC_varid, std::vector<std::string> >::iterator i = data.begin();
  while ( i != data.end() )
    {	
      std::vector<std::string> & d = i->second;
      std::vector<int> & k = datak[ i->first ];
      
      pout.data_group( d[0] );
      for (int j=1; j<d.size(); j++)
	pout.data( d[j] , k[j]);
      
      pout.data( neighbours[ i->first ].size() , 0 );
	
      pout.print_data_group();
      
      ++i;
    }
  data.clear();
  datak.clear();
  neighbours.clear();
}



void f_extra_qc_metrics( Variant & var , void * p )
{
    
  XQCstats * x = (XQCstats*)p;
  
  XQC_varid pos = XQC_varid( var );
  
  x->add( pos , var , 0 ); // first entry will be data_group()
  x->add( pos , var.chromosome() , 0 );
  x->add( pos , var.position() , 0 );
  x->add( pos , var.reference() , 0 );
  x->add( pos , var.alternate() , 1 );

  x->add( pos , var.print_meta_filter(",") , 0 );
  x->add( pos , var.consensus.quality() == -1 ? "." : Helper::dbl2str( var.consensus.quality() ) , 0 );

  bool ti = var.transition();
  bool tv = var.transversion();
  if ( ! ( ti || tv ) ) 
    x->add( pos , "." , 0 );
  else
    x->add( pos , ti , 0 );
  
  //
  // Frequency information for alternate allele(s)
  //

  int c     = 0; // minor allele
  int c_tot = 0; // total counts	  
  double maf;
  bool altmin = var.n_minor_allele( &c , &c_tot , &maf );

  int naa, nab, nbb;
  double hwe = Helper::hwe( var , &naa , &nab , &nbb ) ;
  int t = naa + nab + nbb;

  // number of obs. genotypes
  x->add( pos , (double)t / (double)var.size() , 0 );
  

  // INFO score for imputed / soft-called data?
  if ( Genotype::using_soft_calls )
    {
      double aaf = 0;
      double rsq = DoseReader::Rsq( var , &aaf );
      if ( rsq < 0 ) 
	{
	  x->add( pos , "NA" , 1 );
	  x->add( pos , "NA" , 1 );
	}
      else
	{
	  x->add( pos , rsq , 1 );
	  x->add( pos , aaf , 1 );
	}

    }


  //
  // Frequencies based on hard-calls
  //
  
  // MAC
  if ( c_tot == 0 ) 
    x->add( pos , "NA" , 1 ); 
  else
    x->add( pos , c , 1 ); 

  // MAF
  if ( c_tot == 0 ) 
    x->add( pos , "NA" , 1 ); 
  else
    x->add( pos , maf , 1 ); 

  
  // flag T is REF is minor allele
  x->add( pos , ! altmin , 1 );
  

  // Hardy-Weinberg equilibrium, -log10(p), and heterozygote prop.
  x->add( pos ,  hwe <= 0 ? 0 : hwe , 1 );

  // Heterozygosity
  if ( t>0 ) x->add( pos ,  (double)nab/(double)(naa+nab+nbb) , 1 );
  else x->add( pos , "NA" , 1 );
  
  
  //
  // Note -- the caller for this function will have ensured that the
  // EM-caller was run prior to coming here; thus we can just look up the
  // statistics here
  //
  
  if ( x->em_stats )
    {
      double ent_all, ent_alt;  
      
      var.em.entropy( ent_all , ent_alt );
      
      x->add( pos , var.em.mean_max_posterior() , 1 );
      x->add( pos , ent_all , 1 );
      
      if ( Helper::realnum( ent_alt ) )
	x->add( pos , ent_alt , 1 );
      else
	x->add( pos , "NA" , 1 );
      
      x->add( pos , var.em.frequency() , 1 );
      
    }


  // End of basic QC info for this variant; but we cannot finish up until 
  // we get info on the number of nearby variants; thus just track and 
  // we will print later ( via flush() )

  
  //
  // Dump contents if new are on a new chromosome
  //
  
  
  if ( x->curr_chr != var.chromosome() )
    {
      if ( x->curr_chr != -1 ) 
	x->flush();
      x->curr_chr = var.chromosome();
    }
  else
    {
      //
      // Look at all other SNPs in the library -- if this SNP is within 100bases 
      // then add
      //
      
      std::map<XQC_varid,std::set<int> >::iterator n = x->neighbours.begin();
      while ( n != x->neighbours.end() ) 
	{
	  
	  if ( abs( var.position() - n->first.bp1 ) <= 100 ) 
	    {
	      n->second.insert( var.position() );
	      ++n;
	    }
	  else
	    {
	      n = x->flush( n );
	    }
	}
    }

  
  // this will be populated later (N SNPs)
  x->neighbours[ pos ];
  

}


//
// Dump contents of group from LOCDB
//

bool Pseq::LocDB::loc_view( const std::string & group , const std::vector<std::string> & alias , const bool meta , const bool subregions )
{
  
  Out & pout = Out::stream( "loci" );

  LocDBase * db = g.resolve_locgroup( group ) ;
  if ( ! db ) return false;

  int gid = db->lookup_group_id( group );
  if ( gid == 0 ) return false;
  std::set<Region> loc = db->get_regions( gid ); 
  
  std::set<Region>::iterator i = loc.begin();
  while ( i != loc.end() ) 
    {
      pout << i->name << "\t"
	   << i->coordinate() << "\t"
	   << i->altname << "\t"
	   << db->alias( i->name , false );

      if ( meta ) 
	pout << "\t" << i->meta;
      pout << "\n";

      if ( subregions ) 
	{
	  for (int s=0; s<i->subregion.size(); s++)
	    {
	      const Subregion & sub =  i->subregion[s];
	      pout << "_S" << " " << i->name << "\t" 
		   << sub;
	      if ( meta ) 
		pout << "\t" << sub.meta ;
	      pout << "\n"; 
	    }
	}

      ++i;
    }
  return true;
}



bool Pseq::SeqDB::loc_translate( const std::string & group )
{

  Out & pout = Out::stream( "loci" );
  Out & pbase = Out::stream( "bstats" );

  bool verbose = args.has( "verbose" );
  
  LocDBase * db = g.resolve_locgroup( group ) ;
  
  if ( ! db ) return false;
  
  int gid = db->lookup_group_id( group );
  
  if ( gid == 0 ) return false;
  
  Annotate::setDB( db , &g.seqdb );
  
  std::set<Region> loc = db->get_regions( gid ); 
  
  //
  // Header row
  //

  pout << "ID\tPOS\tNEXON\tNAA\tAASEQ\n";
  
  if ( verbose )
    {
      pbase << "GENE" << "\t"
	    << "POS" << "\t"
	    << "EXON" << "\t"
	    << "STRAND" << "\t"
	    << "A\tC\tG\tT\tN";
      
      std::map<std::string,std::string>::iterator pp = Annotate::aa.begin();
      while ( pp != Annotate::aa.end() )
	{
	  pbase << "\t" << pp->second ;
	  ++pp;
	}
      pbase << "\n";
    }
  
  
  //
  // Process each locus
  //

  std::set<Region>::iterator i = loc.begin();
  while ( i != loc.end() ) 
    {
      
      std::string aa = Annotate::translate_reference( *i , verbose );
      
      pout << i->name << "\t"
	   << i->coordinate() << "\t"
	   << i->subregion.size() << "\t"
	   << aa.size() << "\t";
      
      if ( aa.size() == 0 ) pout << ".";
      
      for (int p=0; p<aa.size(); p++)
	{
	  if ( p > 0 && p % 10 == 0 ) pout << " ";
	  pout << aa[p];
	}
      
      pout << "\n";

      ++i;
    }
  return true;
}


void f_simple_counts( Variant & var , void * p )
{
  
  OptSimpleCounts * data = (OptSimpleCounts*)p;

  Out & pout = Out::stream( data->qt_pheno ? "means" : ( data->genotypes ? "gcounts" : "counts" ) );
  
  std::string gene = ".";
  std::vector<std::string> genevec;
  std::string annot1 = ".";
  std::string annotfull = ".";
  std::string protein = ".";
  
  //
  // Do we have a gene attached, or do we have to look it up? 
  //
  
  if ( data->apply_annot ) 
    {
      
      if ( var.meta.has_field( PLINKSeq::META_GENE() ) )
	{
	  gene = var.meta.get1_string( PLINKSeq::META_GENE() );  
	}
      
      //
      // Do we have a functional annotation attached?
      //
      
      if ( ! var.meta.has_field( PLINKSeq::META_ANNOT() ) )
	{
	  
	  Annotate::annotate( var );

	  annot1 = var.meta.get1_string( PLINKSeq::ANNOT() );	  
	  annotfull = var.meta.as_string( PLINKSeq::ANNOT_TYPE() , "," );
	  gene = var.meta.as_string( PLINKSeq::ANNOT_GENE() , "," );
	  genevec = var.meta.get_string( PLINKSeq::ANNOT_GENE() );
	  protein = var.meta.as_string( PLINKSeq::ANNOT_PROTEIN() , "," );
	  
	}
      else
	annot1 = var.meta.as_string( PLINKSeq::META_ANNOT() , "," );
    }

  
  // Output

  pout << Helper::chrCode( var.chromosome() ) << ":" << var.position() ;

  

  // Alleles, or genotype counts
  
  if ( data->genotypes )
    {
      // print ref/alt allele(s)
      pout << "\t" << var.reference() << "/" << var.alternate();

      std::map<std::string,int> gcnta;
      std::map<std::string,int> gcntu;
      
      if ( data->dichot_pheno )
	{
	  gcnta = var.genotype_counts( CASE );
	  gcntu = var.genotype_counts( CONTROL );

	  pout << "\t";
	  std::map<std::string,int>::iterator ii = gcnta.begin();
	  while ( ii != gcnta.end() )
	    {
	      if ( ii != gcnta.begin() ) pout << " ";
	      pout << ii->first << "=" << ii->second << ":" << gcntu[ ii->first ] ;
	      ++ii;
	    }	  
	}
      else // no phenotype
	{

	  gcntu = var.genotype_counts( UNKNOWN_PHE );
      
	  pout << "\t";
	  std::map<std::string,int>::iterator ii = gcntu.begin();
	  while ( ii != gcntu.end() )
	    {
	      if ( ii != gcntu.begin() ) pout << " ";
	      pout << ii->first << "=" << ii->second ;
	      ++ii;
	    }
	}

    }
  else
    {
      int case_count = 0 , case_tot = 0; 
      int control_count = 0 , control_tot = 0;
      
      bool ma, control_ma;
      
      // QT means
      double refmean = 0 , altmean = 0;
      double refsd = 0   , altsd   = 0;
      int    refn  = 0   , altn    = 0;

      if ( data->dichot_pheno ) 
	{
	  ma = var.n_minor_allele( &case_count , &case_tot , NULL , CASE );
	  control_ma = var.n_minor_allele( &control_count , &control_tot , NULL , CONTROL );
	  // flip alleles?
	  if ( control_ma != ma ) control_count = control_tot - control_count; 
	}
      else if ( data->qt_pheno )
	{

	  // 
	  //	pout << "\tREFMEAN\tREFSD\tREFN"
	  //  << "\tALTMEAN\tALTSD\tALTN";
	  	  
	  const int n = var.size();
	  
	  for (int i = 0; i < n; i++)
	    {	  	  
	      if ( ! ( var.ind(i)->missing() || var(i).null() ) )
		{
		  // always get 'alternate' allele
		  int ac = var(i).minor_allele_count( true );
		  
		  if ( ac ) 
		    {
		      ++altn;
		      altmean += var.ind(i)->qt();
		    }
		  else
		    {
		      ++refn;
		      refmean += var.ind(i)->qt();
		    }
		}
	    }
	  
	  refmean /= (double)refn;
	  altmean /= (double)altn;

	  for (int i = 0; i < n; i++)
	    {	  	  
	      if ( ! ( var.ind(i)->missing() || var(i).null() ) )
		{
		  // always get alternate allele 
		  int ac = var(i).minor_allele_count( true );
		  
		  if ( ac ) 
		    {
		      altsd += ( var.ind(i)->qt() - altmean ) * ( var.ind(i)->qt() - altmean );
		    }
		  else
		    {
		      refsd += ( var.ind(i)->qt() - refmean ) * ( var.ind(i)->qt() - refmean );		      
		    }
		}
	    }
	  
	 
	  altsd /= (double)(altn-1);
	  refsd /= (double)(refn-1);
	  
	}
      else
	ma = var.n_minor_allele( &case_count , &case_tot );

      // print REF / ALT / Minor
      
      pout << "\t" << var.reference() << "/" << var.alternate() ;

      // print minor/major allele(s)
      
      pout << "\t" << ( ma ? var.alternate()  : var.reference() ) ;            
      
      if ( data->dichot_pheno )
	pout << "\t" << case_count
	     << "\t" << control_count ;
      else if ( data->qt_pheno )
	{

	  if ( refn >= 1 )
	    pout << "\t" << refmean;
	  else
	    pout << "\tNA";

	  if ( refn >= 2 )
	    pout << "\t" << refsd;
	  else
	    pout << "\tNA";
	  
	  pout << "\t" << refn;

	  // -----------------------------

	  if ( altn >= 1 )
	    pout << "\t" << altmean;
	  else
	    pout << "\tNA";

	  if ( altn >= 2 )
	    pout << "\t" << altsd;
	  else
	    pout << "\tNA";
	  
	  pout << "\t" << altn;
	  
	}
      else
	pout << "\t" << case_count ;
      
      if ( data->dichot_pheno )
	pout << "\t" << case_tot 
	     << "\t" << control_tot;
      else if ( ! data->qt_pheno )
	pout << "\t" << case_tot;
    }


  if ( data->apply_annot ) 
    {
      pout << "\t" << annot1
	   << "\t" << ( gene == "" ? "." : gene ) ;

      if ( data->apply_full_annot )
	{
	  pout << "\t" << annotfull	       
	       << "\t" << protein;
	  std::set<std::string> aliases;
	  for (int s=0;s<genevec.size();s++) aliases.insert( g.locdb.alias( genevec[s] ) );
	  
	  std::set<std::string>::iterator i = aliases.begin();
	  while ( i != aliases.end() )
	    {
	      if ( i != aliases.begin() ) pout << ","; else pout << "\t";
	      pout << *i ;
	      ++i;
	    }
	}
    }
  
  
  if ( data->meta.size() )
    {
      std::set<std::string>::iterator i = data->meta.begin();
      while ( i != data->meta.end() )
	{
	  if ( var.meta.has_field( *i ) )
	    pout << "\t" << var.meta.as_string( *i );
	  else if ( var.consensus.meta.has_field( *i) )
	    pout << "\t" << var.consensus.meta.as_string( *i );
	  else
	    pout << "\t.";
	  ++i;
	}
    }

  if ( data->show_filter ) 
    pout << "\t" << var.print_meta_filter(";");
  
  pout << "\n";
  
}



bool Pseq::VarDB::header_VCF( const bool show_meta , 
			      const bool show_header , 
			      Mask & mask )
{

  Out & pout = Out::stream( "indiv" );

  // Assumes we are in single-VCF mode, or else do nothing
  if ( ! g.single_file_mode() ) return false;
  
  std::string filename = mask.external_vcf_filename();
  Helper::checkFileExists( filename );
  
  IndividualMap imap;
  VarDBase tmpdb( imap );
  tmpdb.attach( ":memory:" );
  
  File vcffile( filename , VCF );
  VCFReader v( &vcffile , "" , &tmpdb ,  NULL );

  tmpdb.begin();

  bool show_dichot_pheno = g.phmap.type() == PHE_DICHOT;
  bool show_qt = g.phmap.type() == PHE_QT;

  while ( 1 ) 
    { 
      
      VCFReader::line_t l = v.parseLine( );
      
      if ( l == VCFReader::VCF_EOF ) break;      
      if ( l == VCFReader::VCF_INVALID ) continue;
      
      if ( l == VCFReader::VCF_META ) {  } 
      
      if ( l == VCFReader::VCF_HEADER && show_header ) 
	{
	  int n = imap.populate( tmpdb , g.phmap , mask );
	  for (int i=0; i<n; i++)
	    {	      
	      pout << imap(i)->id();
	      if ( show_dichot_pheno ) pout << "\t" << imap(i)->affected();
	      else if ( show_qt ) pout << "\t" << imap(i)->qt();   
	      pout << "\n";
	    }
	}
      
      if ( l == VCFReader::VCF_VARIANT ) break; // assume no more meta-info

    }
  
  if ( show_meta ) 
    {
      pout << MetaInformation<VarMeta>::list_fields("META_VARIANT")
	   << MetaInformation<VarFilterMeta>::list_fields("META_FILTER")
	   << MetaInformation<GenMeta>::list_fields("META_GENOTYPE");
      
    }

  return true;

}


void g_loc_view( VariantGroup & vars , void * p )
{
  Out & pout = Out::stream( "loci" );
  pout << vars.name() << "\t" 
       << vars.coordinate() << "\n";
}
