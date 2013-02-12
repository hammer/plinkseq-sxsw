#include "pyplinkseqint.h"
#include "plinkseq.h"
#include "boost/lexical_cast.hpp"
#include <string>
#include <map>
#include <vector>

GStore *g;

void Py_init_Pyplinkseq() {
  g = new GStore;
  return;
}


std::string Py_gstore_version() {
  std::map<std::string, std::string> g_version = g->version();
  std::map<std::string, std::string>::iterator iter;
  std::string strToReturn;

  for (iter = g_version.begin(); iter != g_version.end(); ++iter) {
    strToReturn.append(iter->first);
    strToReturn.append("=");
    strToReturn.append(iter->second);
    strToReturn.append("\n");
  }

  return strToReturn.c_str();
}

void Py_set_project(std::string project) {
  std::string project_str(project);
  g->set_project(project_str);
}

std::string Py_summary() {
  return g->summary(false);
}

/*Want to access P-value, I statistic, and name of test.
 *
 * void gene_accumulate_assoc( VariantGroup & vars , void *p ){
	Py_geneassocGroup geneassocgroup;

	std::vector<Py_test> tests;
	for( int i = 0; i < )


}*/
void group_accumulate_func( VariantGroup &v, void *p ){
  //Variant Group
	Py_variantGroup vargroup;
	vargroup.NV = v.n_variants();
	vargroup.SIZE = v.size();
	vargroup.COORD = v.coordinate();
	vargroup.SPAN = v.span();
	vargroup.MIDPOS = v.midposition();
	vargroup.NAME = v.name();
	vargroup.NIND = v.n_individuals();
	std::vector<Py_variant> vars;
	for(int i = 0; i < v.n_variants(); i++){
		// Add new variant to the vector of variant groups.
		Variant vn = v.var(i);
		Py_variant var;
		var.CHR = vn.chromosome();
		var.BP1 = vn.position();
		var.BP2 = vn.stop();
		var.ID = vn.name();
		var.NS = vn.n_samples();

		// Sample variant
		SampleVariant &sv = vn.sample(-1);
		Py_sample_variant svar;
		svar.FSET = sv.fileset();
		svar.REF = sv.reference();
		svar.ALT = sv.alternate();
		svar.QUAL = sv.quality();

		// Genotype
		Py_genotype py_geno;
		std::vector<int> gt;
		int nind = vn.size(-1);
		for (int i = 0; i < nind; i++) {
			Genotype* geno = vn.genotype(-1, i);
		    if (geno->more() || geno->null()) {
		    	gt.push_back(-1); // Code NA as -1
		    } else {
		    	gt.push_back(geno->allele_count());
		    }
		  }
		 py_geno.GT = gt;
		 svar.GENO = py_geno;
		 var.CON = svar;

		vars.push_back(var);

	}
	vargroup.VARS = vars;
	std::vector<Py_variantGroup>* d = (std::vector<Py_variantGroup>*)p;
	d->push_back(vargroup);
}




void accumulate_func(Variant &v, void *p) {
  // Variant
  Py_variant var;
  var.CHR = v.chromosome();
  var.BP1 = v.position();
  var.BP2 = v.stop();
  var.ID = v.name();
  var.NS = v.n_samples();

  // Sample variant
  SampleVariant &sv = v.sample(-1);
  Py_sample_variant svar;
  svar.FSET = sv.fileset();
  svar.REF = sv.reference();
  svar.ALT = sv.alternate();
  svar.QUAL = sv.quality();

  // Genotype
  Py_genotype py_geno;
  std::vector<int> gt;
  int nind = v.size(-1);
  for (int i = 0; i < nind; i++) {
    Genotype* geno = v.genotype(-1, i);
    if (geno->more() || geno->null()) {
      gt.push_back(-1); // Code NA as -1
    } else {
      gt.push_back(geno->allele_count());
    }
  }
  py_geno.GT = gt;
  svar.GENO = py_geno;
  var.CON = svar;

  // Add new variant to the vector of variants
  std::vector<Py_variant>* d = (std::vector<Py_variant>*)p;
  d->push_back(var);
}


std::vector<Py_variantGroup> Py_iterateGroup(std::string mask){
	Mask m(mask, "", true, false);
	g->register_mask(m);
	std::vector<Py_variantGroup> varGroups;
	g->vardb.iterate(group_accumulate_func, &varGroups, m);
	return varGroups;

}

std::vector<Py_locGroup> Py_locview(std::string group){
	Py_locGroup pyloc;
	LocDBase * db = g->resolve_locgroup( group );
	std::vector<Py_locGroup> pylocvec;
    int gid = db->lookup_group_id( group );
	std::set<Region> loc = db->get_regions( gid );
	std::set<Region>::iterator i = loc.begin();
  while ( i != loc.end() ){
	  pyloc.NAME = i->name;
	  pyloc.COORD = i->coordinate();
	  pyloc.ALIAS = i->altname;

	  pylocvec.push_back(pyloc);
	  }
  return pylocvec;
}

std::vector<Py_variant> Py_iterate(std::string mask, int limit) {
  Mask m(mask, "", true, false);
  g->register_mask(m);
  m.limit(limit);
  std::vector<Py_variant> vars;
  g->vardb.iterate(accumulate_func, &vars, m);
  return vars;
}

Py_individual_map Py_ind_list(std::string mask , std::string pheno ) {

  Mask m(mask, "", true, false);
  g->register_mask(m);
  g->indmap.populate(g->vardb, g->phmap, m);
  const int n = g->indmap.size();

  
  // Labels 
  std::vector<std::string> ids;
  for (int i = 0; i < n; i++) {
    ids.push_back(g->indmap.ind(i)->id());
  }
  Py_individual_map ind_map;

  ind_map.ID = ids;

  // Phenotypes MR
  std::vector<std::string>  phe1;
  std::string pname = pheno;

  // Try to load
  //  if( g->phmap.set_phenotype(pname) == 0) break ; 
      
  // Determine Type
  mType mt = MetaInformation<IndivMeta>::type( pname );

  //  if( mt == META_UNDEFINED ) break ; 
      
  if( mt == META_BOOL ){
    //	std::vector<std::string> phe1;
    for( int i=0; i < n; i++)
      {
	if( ! g->indmap.ind(i)->meta.has_field( pname ) )
	  phe1.push_back("-9");
	else
	  phe1.push_back( boost::lexical_cast<std::string>(g->indmap.ind(i)->meta.get1_bool( pname )));
      }


  }
  else if ( mt == META_INT )
    {
      // std::vector<std::string> phe1;
      for( int i=0; i < n; i++)
	{
	  if( ! g->indmap.ind(i)->meta.has_field( pname ) )
	    phe1.push_back("-9");
	  else
	    phe1.push_back(boost::lexical_cast<std::string>(g->indmap.ind(i)->meta.get1_int( pname )));
	}

      
    }
  else if (mt == META_FLOAT )
    {
      //  std::vector<std::string> phe1;
      for( int i=0; i < n; i++)
	{
	  if( ! g->indmap.ind(i)->meta.has_field( pname ) )
	    phe1.push_back("-9");
	  else
	    phe1.push_back(boost::lexical_cast<std::string>(g->indmap.ind(i)->meta.get1_double( pname )));
	}	  
    }
  else if ( mt == META_TEXT )
    {
      //	  std::vector<std::string> phe1;
      for( int i=0; i < n; i++)
	{
	  if( ! g->indmap.ind(i)->meta.has_field( pname ) )
	    phe1.push_back("-9");
	  else
	    phe1.push_back(boost::lexical_cast<std::string>(g->indmap.ind(i)->meta.get1_string( pname ).c_str()));
	}
    }
  Py_Phenotype phendata;
  phendata.PHENOTYPE = phe1;
  phendata.LABELS = pname; 
  ind_map.PHENO = phendata; 
  return ind_map;
}


void Py_refdbattach(std::string filename)
{
	g->refdb.attach(filename);
}

void Py_seqdbattach(std::string filename)
{
  g->seqdb.attach(filename);
  Annotate::init();

}

void Py_locdbattach(std::string filename)
{
  g->locdb.attach(filename);
  
}


void Py_annotate_load(std::string loc_id)
{
  Annotate::setDB( LOCDB );
  Annotate::set_transcript_group( loc_id );
   
}

void Py_locdb_load_gtf(std::string filename , std::string name )
{
  g->locdb.load_GTF(filename , name );
  g->locdb.index();
  
}

void Py_locdb_collapse_subregions(std::string id , std::string name )
{
  g->locdb.merge( id , name );
  g->locdb.index();
  
}



void Py_protdbattach( std::string name ){

  g->protdb.attach( name );


}



std::string Py_protdb_summary()
{
  return g->protdb.summary();
}
std::string Py_locdb_summary() 
{
  return g->locdb.summary(false);
}

std::string Py_seqdb_summary() 
{
  return g->seqdb.summary(false);
}

std::string Py_refdb_summary() 
{
  return g->refdb.summary(false);
}



std::set<Py_Feature> Py_protdb_fetch( std::string db , std::string transcript ){
  Py_ProtFeatureSet pys;

  g->protdb.attach( db );
  std::set<Feature> f  = g->protdb.fetch( transcript );
  std::set<Feature>::iterator ii = f.begin();
  while(ii != f.end() ){
    Py_Feature fs;
    fs.SOURCE_ID = ii->source_id ; 
    fs.FEATURE_ID = ii->feature_id;
    fs.FEATURE_NAME = ii->feature_name;
    fs.PROTEIN_ID = ii->protein_id;
    fs.PSTART = ii->pstart;
    fs.PSTOP = ii->pstop;
    fs.MSTR = ii->mstr;
    fs.CHR = ii->chr;
    fs.GSTART = ii->gstart;
    fs.GSTOP = ii->gstop;
    pys.add( transcript , fs );
    ++ii;
  };
  return pys.get( transcript );
}


/*
Py_ProtFeatureSet Py_protdb_lookup( int chr, int bp, std::string ref, std::string alt, std::string info){
  Variant v("Var", chr, bp);
  v.consensus.reference( ref );
  v.consensus.alternate( alt );
  std::set<Py_ProtFeatureSet> s;
  //  std::set<Py_Feature> f;
  ProtFeatureSet pfs =  g->protdb.lookup( v );
  std::map<std::string, std::set<Feature> >::iterator ii = pfs.feat.begin();
  while(ii != pfs.feat.end() )
    {
      std::set<Feature>::iterator jj = ii->second.begin();
      while(jj != ii->second.end() ){
      Py_Feature f;
      f.SOURCE_ID = jj->source_id ; 
      f.FEATURE_ID = jj->feature_id;
      f.FEATURE_NAME = jj->feature_name;
      f.PROTEIN_ID = jj->protein_id;
      f.PSTART = jj->pstart;
      f.PSTOP = jj->pstop;
      f.MSTR = jj->mstr;
      f.CHR = jj->chr;
      f.GSTART = jj->gstart;
      f.GSTOP = jj->gstop;

      s.add( f );
      ++jj;
      }
      ++ii;
    }
  return s;
}


*/
std::string Py_seqdb_annotate( int chr , int bp , std::string ref , std::string alt , std::string info , std::string transcript )
{

  Variant v("Var" , chr, bp);
  v.consensus.reference( ref );
  v.consensus.alternate( alt );
  Annotate::annotate( v );
  if( info == "annot" ){
    return v.meta.get1_string( PLINKSeq::ANNOT() );
  }
  else if( info == "func" ){
    return v.meta.as_string( PLINKSeq::ANNOT_TYPE() , ",");
    
  }
  else if( info == "genomic" ){
    return v.meta.as_string( PLINKSeq::ANNOT_CHANGE() , ",");
    
  }
  else if( info == "codon" ){
    return v.meta.as_string( PLINKSeq::ANNOT_CODON() , ",");
    
  }

  else if( info == "transcript" ){
    return v.meta.as_string( PLINKSeq::ANNOT_GENE() , ",");
    
  }
  else if( info == "alias"){
 /* allow alias to be scanned */
     return  g->locdb.alias( transcript , false ) ;
  }
  else if( info == "summary" ){
    return v.meta.as_string( PLINKSeq::ANNOT_SUMMARY() , ",");
    
  }

  else if( info == "protein" ){
    return v.meta.as_string( PLINKSeq::ANNOT_PROTEIN() , ",");

  }
  else if( info == "details"){
    return v.meta.as_string( PLINKSeq::ANNOT_DETAILS() , ",");

  }
  else{
    return v.meta.as_string( PLINKSeq::META_ANNOT(), ",");
    
  }
  
}
