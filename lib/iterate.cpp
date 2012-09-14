
#include "plinkseq/vardb.h"
#include "plinkseq/gstore.h"

extern GStore * GP;

void add_requires_excludes( std::string & query , const Mask & mask );


//
// Individual and group variant iteration wrappers
//


IterationReport VarDBase::iterate( void (*f)(Variant&, void *) ,
				   void * data ,	
				   Mask & mask )
{
  if ( mask.any_grouping() ) 
    {
      Helper::halt( "grouping mask specified with a non-group command" );      
      return IterationReport( false );
    }
  
  if ( f ) 
    {
      IterationReport report = generic_iterate( f, NULL , data , mask );
      plog >> report.report() ;
      return report;
    }

  Helper::halt( "internal error: no function pointer for generic_iterate() " );

  return IterationReport( false );
  
}

IterationReport VarDBase::iterate( void (*f)(VariantGroup &, void *) ,
				   void * data ,	
				   Mask & mask )
{
  if ( !mask.any_grouping() ) 
    {
      Helper::halt( "no grouping mask specified with a group command" );
      return IterationReport( false );
    }
  // ensure only a single include-group
  mask.ensure_single_include_group();
  
  if ( f ) 
    {
      IterationReport report = generic_iterate( NULL , f, data , mask );
      plog >> report.report();
      return report;
    }

  Helper::halt( "internal error: no function pointer for generic_iterate() " );
  return IterationReport( false );
}


void f_group_variants(Variant & var, void * vars)
{
  ((VariantGroup*)vars)->add( var );
}


//
// Generic iteration function
//


IterationReport VarDBase::generic_iterate( void (*f)(Variant&, void *) ,
					   void (*g)(VariantGroup&, void *) ,
					   void * data ,
					   Mask & mask ) 
{

  

  //
  // Does everything look set up correctly?
  //
  
  if ( mask.invalid() ) return IterationReport( false );
  

  //
  // Upfront, determine if any variant/genotype masks present
  //

  mask.determine_variant_mask();
  
  mask.determine_genotype_mask();


  //
  // Get access to Mask in all functions
  //

  if      (  mask.load_genotype_data()   && (!mask.load_genotype_meta())  &&  mask.load_variant_meta()    ) fetch_mode = NO_GMETA;
  else if ( (!mask.load_genotype_data()) && (!mask.load_genotype_meta())  &&  mask.load_variant_meta()    ) fetch_mode = ONLY_VMETA;
  else if ( mask.load_genotype_data()    && (!mask.load_genotype_meta())  &&  (!mask.load_variant_meta()) ) fetch_mode = ONLY_GENO;
  else if ( (!mask.load_genotype_data()) && (!mask.load_genotype_meta())  &&  (!mask.load_variant_meta()) ) fetch_mode = ONLY_CORE;
  else fetch_mode = ALL;


  //
  // Is this looking at the VARDB fully, or just a single VCF/BCF?
  //

  if ( mask.external_vcf_iteration() )
    {
      if ( mask.any_grouping() ) 
	{
	  plog.warn("cannot specify any grouping in single-VCF mode");
	  return 0;
	}
      
      if ( GP->single_file_bcf() )
	return bcf_iterate( f , data , mask );
      else
	return vcf_iterate( f , data , mask );
      
    }
  

  //
  // Do we have a valid VARDB?
  //
  
  if ( ! attached() )  return IterationReport( false );


  //
  // Track how many variants we evaluate/process 
  //
  
  IterationReport rep( true , mask.any_grouping() , mask.variant_limit() );



  //
  // Create crafted query
  //
  
  
  // and some other settings from the Mask: should overlapping
  // sample-variants be merged to make a single variant?
  
  merge_mode = mask.mergemode();
  

  // downcode mode for multi-allelic markers
  
  downcode_mode = mask.downcode();
  

  //
  // If this is a group query, use the individual level 
  // function to accumulate variants into a group
  //
  
  void * gdata = NULL;
  
  VariantGroup vars(mask);
  
  
  
  
  //
  // Swap around functions if grouping is specified
  //

  if ( mask.any_grouping() )
    {

      f = f_group_variants;

      gdata = data;

      data = &vars;

    }
  
  
  
  //
  // Attach the appropriate databases
  //
  
  if  ( mask.loc_any() )
    {
      if ( locdb_attached ) sql.query(" DETACH DATABASE locdb; " );      
      sql.query(" ATTACH \"" + mask.locdb_name() + "\" AS locdb; " );
      locdb_attached = true;
    }

//   if  ( mask.seg_any() )
//     sql.query(" ATTACH \"" + mask.segdbname() + "\" AS segdb; " );
  
  if  ( mask.ref_any() )
    {
      if ( refdb_attached ) sql.query(" DETACH DATABASE refdb; " );
      sql.query(" ATTACH \"" + mask.refdb_name() + "\" AS refdb; " );
      refdb_attached = true;
    }
  


  //
  // Create a temp database in memory
  //

  if ( mask.build_temporary_db() ) 
    build_temporary_db( mask );


  
  //
  // Register new, variable-length meta-information fields
  //
  
  if ( mask.var() || mask.var_append() )
    {
      MetaInformation<VarMeta>::field( PLINKSeq::META_VSET(), META_TEXT , -1 , "Variant set name");
      MetaInformation<VarMeta>::field( PLINKSeq::META_VGRP(), META_INT  , -1 , "Variant set group");
    }
  
  if ( mask.loc() || mask.loc_append() )
    {
      MetaInformation<VarMeta>::field( PLINKSeq::META_LSET(), META_TEXT , -1 , "Locus name");
      MetaInformation<VarMeta>::field( PLINKSeq::META_LGRP(), META_INT  , -1 , "Locus group");
    }
  
  if ( mask.loc_set() || mask.loc_set_append() )
    {
      MetaInformation<VarMeta>::field( PLINKSeq::META_LSSET(), META_TEXT , -1 , "Locus set name");
      MetaInformation<VarMeta>::field( PLINKSeq::META_LSGRP(), META_INT  , -1 , "Locus set group");
    }
  
  
  //
  // Build the specific compound query, with masks and appends as appropriate
  //
  
  // INCLUDES
  //       var, loc, seg, locset, segset, reg

  // REQUIRES
  //       var, loc, set, reg, file
  
  // EXCLUDES
  //       var, loc, reg, file

  // APPENDS 
  //       var, loc, ref
  
  // GROUPS
  //       var, loc, locset, file
  
  
  
  //
  // SELECT ...
  //
  
  bool multiple_queries = false;
  
  std::string query = "";
  
  
  if ( mask.var() ) 
    {
      multiple_queries = true; 


//       SELECT v.var_id, v.file_id, v.name, v.chr, v.bp1 , v.bp2 , v.offset , 1 , 0 , x.set_id FROM variants AS v 
//       INNER JOIN ( SELECT * FROM set_data WHERE set_id IN ( 1 , 3 ) ) AS x ON x.var_id == v.var_id ;

//       if ( mask.var_exceptions() )
//         {
//           query += " INNER JOIN tmp.var AS x "
//             " ON v.var_id == x.var_id ";
//         }
//       else if ( mask.var() )
//         {

// 		       + mask.var_include_string() + " ) ) AS x "
//             " ON x.var_id == v.var_id ";
//         }

      //// ---
      
 //   query += "SELECT v.var_id, v.file_id, v.name, v.chr, v.bp1 , v.bp2 , v.offset , 1 , x.name, x.group_id "; 
      query += "SELECT v.var_id, v.file_id, v.name, v.chr, v.bp1 , v.bp2 , v.offset , 1 , 0 , x.set_id "; 
      
      query += " FROM variants AS v ";
      
      if ( mask.var_exceptions() )
	{	
	  query += " INNER JOIN tmp.var AS x "
	    " ON v.var_id == x.var_id ";
	}
      else if ( mask.var() )
	{
	  query += " INNER JOIN "
	    "   ( SELECT * FROM set_data WHERE set_id IN ( " 
	    + mask.var_include_string() + " ) ) AS x "
	    " ON x.var_id == v.var_id ";	    
	}
      
      // Grouping? There will only be a single include here, so no need to check group
      // (based on a varset.group)

       if ( mask.named_grouping() )
 	{
	  query += " AND v.var_id IN ( SELECT var_id FROM set_data WHERE set_id IN "
	    " ( SELECT set_id FROM sets WHERE name == :grp_value ) ) ";
 	}

    

      // Requires/excludes performed here
      
      add_requires_excludes( query , mask );
      
    }
  
  
  if ( mask.loc() ) 
    {
      
      if ( multiple_queries ) 
	query += " UNION ALL ";
      
      multiple_queries = true; 
      
      query += "SELECT v.var_id, v.file_id, v.name, v.chr, v.bp1 , v.bp2, v.offset , 2 , x.name, x.group_id "; 
      
      query += " FROM variants AS v ";
      
      if ( mask.named_grouping() ) 
	{
	  
	  // Gouping? There will only be a single include here, so no need to check group
	  
	  if ( mask.loc_exceptions() )
	    {	
	      query += " INNER JOIN ( SELECT * FROM tmp.loc WHERE group_id == " 
		+ Helper::int2str( mask.group_set() ) + " AND name == :grp_value ) AS x "
		" ON v.chr == x.chr AND v.bp1 BETWEEN x.bp1 AND x.bp2 ";	
	    }
	  else if ( mask.loc() )
	    {
	      
	      query += " INNER JOIN ( SELECT * FROM locdb.loci WHERE group_id == " 
		+ Helper::int2str( mask.group_set() ) + " AND name == :grp_value ) AS x "
		"  ON x.group_id IN ( " + mask.loc_include_string() + " ) "
		" AND v.chr == x.chr AND v.bp1 BETWEEN x.bp1 AND x.bp2 "; 
	    }

	}
      else
	{
	  if ( mask.loc_exceptions() )
	    {	
	      query += " INNER JOIN tmp.loc AS x "
		" ON v.chr == x.chr AND v.bp1 BETWEEN x.bp1 AND x.bp2 ";	
	    }
	  else if ( mask.loc() )
	    {
	      query += " INNER JOIN locdb.loci AS x "
		"  ON x.group_id IN ( " + mask.loc_include_string() + " ) "
		" AND v.chr == x.chr AND v.bp1 BETWEEN x.bp1 AND x.bp2 ";		
	    }
	}
      

      // Requires/excludes performed here
      
      add_requires_excludes( query , mask );

    }
  
  
    //
    // Mask from locus-database?
    //
    
    if ( mask.loc_set() ) 
      {

	if ( multiple_queries ) 
	  query += " UNION ALL ";

	multiple_queries = true; 

	query += "SELECT v.var_id, v.file_id, v.name, v.chr, v.bp1 , v.bp2, v.offset , 3 , x.name, x.group_id "; 
	
	query += " FROM variants AS v ";
	

	if ( mask.named_grouping() ) 
	  {
	    query += " INNER JOIN ( SELECT * FROM tmp.locset WHERE group_id == "
	      + Helper::int2str( mask.group_set() ) + " AND name == :grp_value ) AS x "
	      " ON v.chr == x.chr AND v.bp1 BETWEEN x.bp1 AND x.bp2 ";		    
	  }
	else 
	  {
	    query += " INNER JOIN tmp.locset AS x "
	      " ON v.chr == x.chr AND v.bp1 BETWEEN x.bp1 AND x.bp2 ";	
	  }
	
	
	// Requires/excludes performed here
	
	add_requires_excludes( query , mask );

      }
    
    
    //
    // Specific regions
    //
        
    if ( mask.reg() ) 
      {

	if ( multiple_queries ) 
	  query += " UNION ALL ";
	
	multiple_queries = true; 
	
	query += "SELECT v.var_id, v.file_id, v.name, v.chr, v.bp1 , v.bp2, v.offset , 0 , NULL , NULL  "; 
	
	query += " FROM variants AS v ";
	
	if ( mask.reg() )
	  {      
	    query += 
	      " INNER JOIN tmp.regin AS x "
	      "     ON v.chr == x.chr "
	      "    AND v.bp1 BETWEEN x.bp1 AND x.bp2 ";
	  }
		

	// Requires/excludes performed here
	
	add_requires_excludes( query , mask );

      }


    //
    // Named variant IDs
    //
        
    if ( mask.id() ) 
      {

	if ( multiple_queries ) 
	  query += " UNION ALL ";
	
	multiple_queries = true; 
	
	query += "SELECT v.var_id, v.file_id, v.name, v.chr, v.bp1 , v.bp2, v.offset , 0 , NULL , NULL  "; 
	
	query += " FROM variants AS v ";
	
	query += 
	  " INNER JOIN tmp.vidin AS x "
	  "     ON v.name == x.name ";
	
	// Requires/excludes performed here
	
	add_requires_excludes( query , mask );

      }


    
    //
    // If no includes, implies all variants 
    //

    if ( ! ( mask.var() || mask.reg() || mask.loc() || mask.id() || mask.loc_set() || mask.var_set() ) ) 
      {
	
	query += "SELECT v.var_id, v.file_id, v.name, v.chr, v.bp1 , v.bp2, v.offset, 0 , NULL , NULL "; 
	query += " FROM variants AS v ";
	
	// Requires/excludes performed here
	
	add_requires_excludes( query , mask );
	
      }
    

    

    //
    // ORDER BY ...
    //
    
    query += " ORDER BY v.chr , v.bp1, v.bp2 , v.var_id ";
    
    
    //
    // END OF QUERY 
    //

    query += " ; ";


    //
    // Prepare this query and execute; note-- a particular variant may
    // be represented over multiple rows (i.e. if joined with distinct
    // meta-information from other tables). Thus keep looping until we
    // are sure this is the variant.
    //
    
    //    std::cout << "Q = [" << query << "]\n";

 

    //
    // If we are grouping, need to specify the outer loop
    //
    
    sqlite3_stmt * s0 = NULL;

    if ( mask.named_grouping() ) 
      {
	s0 = sql.prepare( "SELECT grp FROM tmp.grp" );
      }


    //
    // Compile the core query
    //

    sqlite3_stmt * s = sql.prepare( query );
    

    
    //////////////////////////////////////////////////////////////////////////////////////
    
    //
    // Step through results of core query (s); if a group-based query
    // has been specified then wrap the core query within that (s0)
    //


    while ( 1 ) 
      {
	
	std::string grp_name = "";

	if ( mask.named_grouping() ) 
	  {
	    
	    vars.clear();
	    
	    if ( sql.step(s0) )
	      {		
		grp_name = sql.get_text( s0 , 0 ); 		
		sql.bind_text( s , ":grp_value" , grp_name );
		vars.name( grp_name );
	      }
	    else
	      break; // all done

	  }
	
	
	//
	// Process core query
	//

	bool expecting_new_variant = true;

	bool any_data = false;
		
	Variant var;

	
	// Track ID and position 

	SampleVariant * sample = NULL;
	

	//
	// Read first row
	//


	if ( sql.step(s) ) 
	  {	
	      
	      sample = &(construct( var , s , &indmap ));
	      
	      sample->decode_BLOB_alleles();
	      
	      if ( mask.named_grouping () ) 
		  var.meta.add( PLINKSeq::META_GROUP() , grp_name );
	      
	      addMetaFields( var , s, mask ); 
	      expecting_new_variant = false;
	      any_data = true;
	      
	  }
	else 
	  {

	    // Handle empty groups if desired

	    if ( mask.named_grouping() && mask.process_empty_groups() )
	      {
		g( vars , gdata );
	      }

	  }

	
	
	//
	// Iterate through all other rows
	//
	
	if ( any_data ) 
	  while ( sql.step(s) ) 
	    {	      
	      

	      Variant nextrow;

	      SampleVariant & next_sample = construct( nextrow, s , &indmap );
	      
	      if ( mask.named_grouping () ) 
		nextrow.meta.add( PLINKSeq::META_GROUP() , grp_name );
	      
	      
	      // 
	      // Does this row contain just more meta-information for exactly the 
	      // same original variant
	      //
	      
	      if ( sample->index() == next_sample.index() ) 
		{
		  addMetaFields( var , s, mask );
		  continue;
		}
	      
	      
	      // 
	      // Does this row contain a distinct SampleVariant, but
	      // that should be attached to the same genomic Variant?
	      // (i.e. overlapping reference position)
	      //
	      
	      // will always be ordered by start position, so test
	      // whether the next variant starts after the current
	      // set. then update the span of the current variant
	      
	      // NOTE: we will need to come up with a strategy to handle very
	      // large events (e.g. Mb CNVs) that will otherwise bring many
	      // distinct variants into one mega-variant, which will typically
	      // be unwieldy and not what is desired
	      
	      bool same_variant; 
	      

	      // ensure that we know the REF and ALT codes at this stage

	      next_sample.decode_BLOB_alleles();
	      

	      if ( merge_mode == MERGE_MODE_NONE ) 
	      {
		  
		// checks		
		// std::cout << "nr " << nextrow.position() << " " << nextrow.reference() << " " << nextrow.alternate() << "\n";
		// std::cout << "s  " << var.position() << " " << sample->reference() << " " << sample->alternate() << "\n";
		// std::cout << "svS  " << sample->reference() << " " << sample->alternate() << "\n";
		// std::cout << "svNS  " << next_sample.reference() << " " << next_sample.alternate() << "\n";
		
		// requires an exact match of REF and ALT
		
		same_variant = nextrow.chromosome() == var.chromosome() 
		  && nextrow.position() == var.position() 
		  && sample->reference() == nextrow.reference() 
		  && sample->alternate() == nextrow.alternate() ;
		
	      }
	      else if ( merge_mode == MERGE_MODE_EXACT )
	      {

		  // here the REF and ALT alleles must all have the same size

		  same_variant = nextrow.chromosome() == var.chromosome() 
		      && nextrow.position() == var.position() 
		      && nextrow.stop() == var.stop();
		  
		  if ( same_variant && var.alternate() != nextrow.alternate() )		      
		    {
		      // also check that the ALT alleles have a similar size
		      std::set<int> sz;		      
		      std::vector<std::string> a1 = Helper::char_split( sample->alternate() , ',' );
		      std::vector<std::string> a2 = Helper::char_split( nextrow.alternate() , ',' );
		      for (int i=0;i<a1.size();i++) sz.insert( a1[i].size() );
		      for (int i=0;i<a2.size();i++) sz.insert( a2[i].size() );
		      if ( sz.size() != 1 ) same_variant = false; 							
		  }
		  
	      }
	      else // merge anything that overlaps (or attempt to)
	      {
		  same_variant = nextrow.chromosome() == var.chromosome() && nextrow.position() <= var.stop();
	      }
	      
	     		  
	      
	      if ( same_variant )
		{
		  
		  // if allowing non-exact matches, update spans and keep track of
		  // individual SV base positions
		  
		  if ( merge_mode == MERGE_MODE_ANY_OVERLAP )
		    {
		      
		      if ( nextrow.stop() > var.stop() ) 
			var.stop( nextrow.stop() );
		      
		      // also, ensure that all allele codes start from the same position
		      
		      next_sample.offset = nextrow.position() - var.position();
		      
		    }
		  
		    
		  // track what the latest SampleVariant is in 'sample'
		    
		  sample = &(var.add( next_sample ));
		  
		  continue;
		}
	      
	      
	      
	      //
	      // If we are here, it means that nextrow contains a
	      // genomically different variant (non-overlapping)
	      //
	      
	      //
	      // Did we observe any overlapping variants within the same file? If so, 
	      // we will, for this variant only, force the individual-alignment to be
	      // non-flat/multi-sample, so that we appropriately track the multiple 
	      // instances of SampleVariant within the same file (which we assume will
	      // not happen too often).  This will include when "overlapping" variants
	      // are being merged though, e.g. a SNP in an indel, that have been 
	      // coded separately.  In either case, the variant will bed contained 
	      // within a single file
	      //
	      
	      if ( var.infile_overlap() ) 
		{		  
		  indmap.force_unflat( true );
		}

	      //
	      // In single-variant (non-group) iteration, just apply the function
	      // In group-iteration mode, the single-variant function has been 
	      // adapted to add this variant to the group, potentially, so run
	      // it in any case
	      //

	      if ( eval_and_call( mask, &indmap , var , f , data ) ) 
		{
		  indmap.force_unflat( false );
		  if ( ! rep.accepted_variant() )
		    break;
		}
	      else
		{
		  indmap.force_unflat( false );
		  rep.rejected_variant();
		}		


	      //
	      // In single-variant iteration mode, we are done now
	      //
	      
	      if ( ! mask.any_grouping() ) 
		{

		  if ( rep.processed() % 100 == 0 ) 
		    plog.counter2( Helper::int2str( rep.processed() ) + " accepted variants" );
		  
		  // Update our current 'variant'
		  // and add the meta-fields
		
		  var = nextrow;
		  
		  sample = var.psample( 0 ) ;
		  
		  addMetaFields( var , s, mask );	
		  
		  continue;
		}
	      
	      
	      //
	      //  Group-iteration ... do we now want to call g-function?
	      //
	      
 	      if ( vars.complete() )
 		{	    
		  
		  if ( rep.groups() % 10 == 0 ) 
		    plog.counter2( Helper::int2str( rep.groups() ) + " accepted variant-groups" );
		  
 		  // If the above is true, it means that the last 
 		  // variant could not be added to the group. Therefore, 
 		  // process the group as-is, then begin a new group, 
 		  // starting with the most recent variant, that could
 		  // not be added to the group
		  
 		  g( vars , gdata );

		  rep.processed_group();

 		  vars.clear(var);
 		}
	      
	      
	      var = nextrow;

	      sample = var.psample( 0 );

	      addMetaFields( var , s, mask );
	      
	    } // next row
	
	
	//
	// Finished iterating -- process the last variant (if any)
	//
	
	if ( any_data && ( ! expecting_new_variant ) && ( ! rep.reached_limit() ) )
	  {
	    
	    // Individual variant function call

	    if ( var.infile_overlap() ) 
	      indmap.force_unflat( true );
	    
	    if( eval_and_call( mask, &indmap, var , f , data ) )
	      rep.accepted_variant(); // no more variants at this point anyway
	    else
	      rep.rejected_variant();
	    
	    indmap.force_unflat( false );

	    // Call group-iteration function (might need to 
	    // call twice)
	    
	    if ( mask.any_grouping() )
	      {
		if ( ! vars.complete() )
		  {
		    g( vars , gdata );
		    rep.processed_group();
		  }
		else
		  {

		    g( vars , gdata );
		    rep.processed_group();

		    vars.clear(var);

		    g( vars , gdata );
		    rep.processed_group();

		  }

		if ( rep.groups() % 10 == 0 ) 
		  plog.counter2( Helper::int2str( rep.groups() ) + " accepted variant-groups" );
		
	      }
	    
	  }
	
	
	//
	// Clean up any temp/attached databases and other resources
	//
	
	sql.reset( s );
	

	//
	// If we are not grouping, then we are all done now
	//
	
	if ( ! mask.named_grouping() ) 
	  break;
    
      }
    
    
    if ( mask.build_temporary_db() )
      {
	tmpdb_attached = false;
	sql.query("DETACH DATABASE tmp; ");    
      }
    
    if ( mask.loc() )
      {
	locdb_attached = false;
	sql.query("DETACH DATABASE locdb; ");
      }

    

    //
    // Return report of how many variants considered, etc
    //
    
    return rep;

}




bool VarDBase::eval_and_call( Mask & mask, 
			      IndividualMap * align, 
			      Variant & var , 
			      void (*f)(Variant&, void *) ,
			      void * data )
{
  
  
  //
  // Evaluate any variant filters (on a per-SampleVariant basis)
  //
  
  const int n = var.n_samples();

  
  bool good = false;
  
  std::vector<int> svar_rlist; // any SampleVariants that need to be removed
  
  int ngood = 0; // track # of passing sampels
      
  for ( int s = 0 ; s < n ; s++ ) 
    {
      
      SampleVariant & svar = var.sample( s );
      
      SampleVariant * target = ( ! align->multi_sample() ) ? &(var.consensus) : &svar ;
      
      SampleVariant * vmeta_target = MetaMeta::force_consensus() ? &(var.consensus) : target ;
      
      SampleVariant * genotype_target = align->flat() ? &(var.consensus) : &svar ;
      
      bool sample_okay = true;

      // exclude any empty sets? -- but would this impact summary handles of VCFs?  yes...
      //  thus, not a good idea to implement, as empty VCFs are a handy device
      
      
      // Basic stuff (allele encoding)
      
      svar.decode_BLOB_basic( target ); 

      
      //
      // Is this a populated file, but from which nobody was actually
      // included?  In that case, we should treat this as if it does
      // not exist.  In the case of an initially empty VCF
      // (i.e. sites-only VCF) then N=0 is no longer a critieron for
      // exclusion, however.
      //
      
      if ( align->size( svar.fileset() ) == 0 && ! mask.site_only( svar.fileset() ) ) sample_okay = false;
      
      if ( ! sample_okay ) { svar_rlist.push_back( s - svar_rlist.size() ); continue; } 


      //
      // Attach any extra meta-information that belongs to this sample-variant
      //

      if ( mask.attach_meta() ) 
	{	  
	  if ( mask.attach_all_meta() ) 
	    {
	      attach_indep_metadata( svar.index() , *vmeta_target , var , NULL );
	    }
	  else
	    {
	      std::set<std::string> s = mask.attach_meta_fields();
	      attach_indep_metadata( svar.index() , *vmeta_target , var , &s );
	    }
	}
      
      
      // These filters are based on QUAL and FILTER fields -- these
      // will travel with target, not vmeta_target 
      
      if ( ! mask.eval_filters( *target ) ) 
	{	  
	  if ( mask.fail_on_sample_variant() ) return false;
	  else sample_okay = false;	
	}
      
      if ( ! sample_okay ) { svar_rlist.push_back( s - svar_rlist.size() ); continue; } 
      
      //
      // At this point, we now have to create the full sample-variant with all
      // meta-information and genotype information. Apply any variant-level
      // meta-information filters at this step also
      //
      
      if ( ! svar.decode_BLOB_vmeta( &mask , &var , vmeta_target ) )
	{
	  if ( mask.fail_on_sample_variant() ) return false;
	  else sample_okay = false;
	}
      
      if ( ! sample_okay ) 
	{ 
	  svar_rlist.push_back( s - svar_rlist.size() ); 
	  continue; 
	} 
      

      //
      // Assign all genotypes (and genotype meta-information) to
      // either a normal sample variant, or directly to the consensus
      // sample variant, if a 'flat' alignment.  This function will 
      // also evaluate any include="expressions" that depend on 
      // genotpe data (i.e. include a g()) and thus we can reject a 
      // SampleVariant at this level because of that
      //
      


      if ( ! svar.decode_BLOB_genotype( align , &mask , &var, &svar , target , genotype_target ) )
	{
	  if ( mask.fail_on_sample_variant() ) return false;
	  else sample_okay = false;
	}

      if ( ! sample_okay ) 
	{ 
	  svar_rlist.push_back( s - svar_rlist.size() ); 
	  continue; 
	} 
      
      //
      // Apply vmeta filters (include statements) if they depend on genotype (
      // "g()" functions ), as we wouldn't have been able to earlier
      //

      if ( mask.filter_expression() && mask.filter_expression_requires_genotypes() ) 
	{
	  if ( ! mask.calc_filter_expression( var , *vmeta_target, *genotype_target ) ) 
	    {
	      if ( mask.fail_on_sample_variant() ) return false;
	      else sample_okay = false;
	    }
	}
      

      //
      // Here, we must have at least 1 good sample
      //
      
      if ( sample_okay )
	{
	  good = true;
	  ++ngood;
	}


      //
      // If this particular sample not okay, need to remove from Variant (based on SampleVariant filters), if 
      // dealing with multiple-SampleVariants in a Variant (adjusting for shrinkage of list as first items are 
      // removed)
      //
      
      if ( ! sample_okay ) svar_rlist.push_back( s - svar_rlist.size() );
            
      // next SampleVariant
    } 


  //
  // Did we observe at least one good sample?
  //

  if ( ! good ) return false;

  
  //
  // Evaluate fail.nfile (i.e. that we did not fail on more than nfailed
  // samples, and that we have at least 'ngood' good samples)
  //

  if ( ! mask.test_fail_on_sample_variant( var.n_samples() - ngood , ngood  ) ) return false;


  //
  // Do we need to remove any individual SampleVariants?
  //

  if ( svar_rlist.size() > 0 ) 
    {
      for ( int i = 0 ; i < svar_rlist.size(); i++)
	var.remove( svar_rlist[i] );      
    }

  

  // 
  // If everything looks good with the individual SampleVariants, now 
  // create a consensus set, from which we can calculate a basic allele
  // frequency
  //

  var.make_consensus( align );



  //
  // Variant/file filters
  //
  
  if ( ! mask.eval_obs_file_filter( var ) ) return false;
  
  if ( ! mask.eval_alt_file_filter( var ) ) return false;

  if ( ! mask.eval_file_count( var ) ) return false;


  //
  // Frequency filter, etc
  //

  if ( ! var.frequency_filter(&mask) ) return false;

  if ( ! mask.polymorphism_filter( var ) ) return false;

  if ( mask.hwe_filter() && ! mask.hwe_filter( Helper::hwe( var ) ) ) return false;

  if ( ! var.null_filter(&mask) ) return false;

  if ( ! var.case_control_filter(&mask) ) return false;

 
  // 
  // And also add appends now, as these may be used by the filter function
  //
  
  if ( mask.var_append() )
    {
      mask.vardb_pointer()->append_metainformation( var , mask.appended_var() );      
    }
  
  if ( mask.loc_append() )
    {
      mask.locdb_pointer()->append_metainformation( var , mask.appended_loc() );      
    }
  
//   if ( mask.seg_append() )
//     {
//       mask.segdb_pointer()->append_metainformation( var , mask.appended_seg() );      
//     }
  
  if ( mask.loc_set_append() )
    {
      Helper::halt( "locset.append not functional yet" );
    }
    

  
  //
  // Appends from reference-database
  //
  
  if ( mask.ref_append() )
    {

      // Either a) just appended information (name/value) 
      // and meta-information, or b) if 1 or more 
      // are also 'required', then check that also
      
      const std::set<int> r = mask.appended_ref();
      const std::set<int> inc = mask.included_ref();
      const std::set<int> req = mask.required_ref();
      const std::set<int> exc = mask.excluded_ref();
      
      const bool use_inc = inc.size() > 0 ;
      const bool use_req = req.size() > 0 ;
      const bool use_exc = exc.size() > 0 ;
      
      if ( ! ( use_inc || use_req || use_exc ) )
	{
	  std::set<int>::iterator i = r.begin();
	  while ( i != r.end() )
	    mask.refdb_pointer()->annotate( var , *(i++) );
	}
      else
	{
	  bool inc1 = false;
	  std::set<int>::iterator i = r.begin();
	  std::set<int> fnd;
	  while ( i != r.end() )
	    {
	      bool b = mask.refdb_pointer()->annotate( var , *i );
	      
	      if ( use_req && b && req.find( *i ) != req.end() ) fnd.insert(*i); 
	      if ( use_exc && b && exc.find( *i ) != exc.end() ) return false;
	      if ( use_inc && b && inc.find( *i ) != inc.end() ) inc1 = true;
	      ++i;
	    }
	  	   
	  if ( fnd.size() != req.size() ) return false;
	  if ( use_inc && (!inc1) && (!use_req) ) return false;
	}
    }
  
  
  //
  // Filter on reference/appended meta-information (i.e. all Variant-level information)
  //
  
  if ( mask.var_filter_expression() )
    {
      if ( ! mask.var_calc_filter_expression( var ) ) return false; 
    }

  
  //
  // Apply any user-defined functions
  //

  if ( mask.func() ) 
    {
      if ( ! mask.eval( var ) ) return false;
    }
  

  //
  // If EM phasing requested, perform here
  //

  if ( mask.EM_caller() ) 
    {
      var.em.load(var);
      var.em.estimate();
      var.em.call( mask.EM_threshold() );
    }
  
  
  //
  // Otherwise, if we are here, we are happy with this variant, so
  // perform the call, optionally about to handle multi-allelic markers 
  // via a downcoding procedure
  //

  
  if ( var.multiallelic() && downcode_mode != DOWNCODE_MODE_NONE ) 
    {      

      // need to perform some variety of downcoding:

      if ( downcode_mode == DOWNCODE_MODE_ALL_ALT )
	{
	  // perform a once-off downcoding (of the consensus variant only)	  
	  var.consensus.collapse_alternates( &var );	  
	  f( var , data );
	}
      else if ( downcode_mode == DOWNCODE_MODE_EACH_ALT )
	{
	  int na = var.n_alleles();
	  for (int k=0; k<na; k++)
	    {
	      Variant var2 = var;
	      var2.consensus.collapse_alternates( &var2 , k );	      
	      f( var2 , data );
	    }
	}
    }
  else // just call the function as is
    {      
      f( var , data );
    }
  

    
  return true;

}


void add_requires_excludes( std::string & query , const Mask & mask )
{
  
  bool donewhere = false;
  
  if ( mask.files() )
    {	    
      query += " WHERE v.file_id IN ( " + mask.files_include_string() + " ) ";
      donewhere = true;
    }
  if ( mask.xfiles() )
    {
      query += donewhere ? " AND " : " WHERE ";
      query += " v.file_id NOT IN ( " + mask.files_exclude_string() + " ) ";
      donewhere = true;
    }
  if ( mask.requires() ) 
    {
      query += donewhere ? " AND " : " WHERE ";
      query += " v.var_id IN (SELECT var_id FROM tmp.require) ";
      donewhere = true;
    }
  if ( mask.excludes() )
    {
      query += donewhere ? " AND " : " WHERE ";
      query += " v.var_id NOT IN (SELECT var_id FROM tmp.exclude) ";
    }

}




void VarDBase::build_temporary_db( Mask & mask )
{      

  if ( tmpdb_attached ) 
    sql.query(" DETACH DATABASE tmp ; " );
  sql.query(" ATTACH \":memory:\" AS tmp ; " );
  tmpdb_attached = true;

  // Construct and populate the necessary temporary tables
  
  if ( mask.var_exceptions() )
    {
      
      sql.query("CREATE TABLE tmp.var "
		" ( var_id    INTEGER , "
		"   name      VARCHAR(10) , "
		"   group_id  INTEGER );");
      
      sql.query("CREATE TABLE tmp.var2( name VARCHAR(20) ); ");
      
      sqlite3_stmt * stmt_tmp_insert = 
	sql.prepare("INSERT OR IGNORE INTO tmp.var(var_id,name,group_id) "
		    " values( :var_id,:name,:grp) ; ");
      
      sqlite3_stmt * stmt_tmp_varset_iterate = 
	sql.prepare("SELECT sd.var_id,sm.name "
		    "  FROM set_members AS sm, set_data AS sd "
		    " WHERE sm.group_id == :grp "
		    "   AND sm.name IN ( SELECT name FROM tmp.var2 ) "
		    "   AND sd.set_id == sm.set_id ; " );
      
      sqlite3_stmt * stmt_tmp_insert0 = 
	sql.prepare( " INSERT INTO tmp.var2 ( name ) values ( :name ) ; " );
      
      
      
      //
      // Consider each group 
      //
      
      sql.begin();
      
      const std::set<int> & vset = mask.included_var();
      std::set<int>::iterator i =  vset.begin();
      while ( i != vset.end() )
	{
	      
	  std::set<std::string> names = mask.subset_var(*i);
	  std::set<std::string>::iterator j = names.begin();
	  
	  sql.query("DELETE FROM tmp.var2;");
	  while ( j != names.end() )
	    {
	      sql.bind_text( stmt_tmp_insert0 , ":name" , *j );
	      sql.step( stmt_tmp_insert0 );
	      sql.reset( stmt_tmp_insert0 );
	      ++j;
	    }
	  
	  std::set<std::string> ignore_names = mask.skip_var(*i);
	  
	  sql.bind_int( stmt_tmp_varset_iterate , ":grp" , *i );
	  while ( sql.step( stmt_tmp_varset_iterate ) )
	    {		    
	      std::string n = sql.get_text( stmt_tmp_varset_iterate , 1 );
	      
	      // Explicit exclusion of this locus?
	      if ( ignore_names.find(n) != ignore_names.end() ) { continue; }
	      
	      uint64_t var_id = sql.get_int64( stmt_tmp_varset_iterate , 0 );
	      
	      sql.bind_int64( stmt_tmp_insert, ":var_id" , var_id ); 
	      sql.bind_text( stmt_tmp_insert, ":name" , n ); 
	      sql.bind_int( stmt_tmp_insert, ":grp" , *i ); 
	      sql.step( stmt_tmp_insert );
	      sql.reset( stmt_tmp_insert );
	    }
	  sql.reset( stmt_tmp_varset_iterate );
	  ++i;
	}
      
      sql.commit();
      
      sql.finalise( stmt_tmp_insert );
      sql.finalise( stmt_tmp_insert0 );
      sql.finalise( stmt_tmp_varset_iterate );	    	    
      
    }
  
  
  
  //
  // Any named loci ? 
  //
  
  if ( mask.loc_exceptions() )
    {
      
      sql.query("CREATE TABLE tmp.loc "
		" ( chr       INTEGER, "
		"   bp1       INTEGER, "
		"   bp2       INTEGER, "
		"   name      VARCHAR(10) , "
		"   group_id  INTEGER  ) ; ");
      
      
      sqlite3_stmt * stmt_tmp_insert = 
	sql.prepare("INSERT INTO tmp.loc(chr,bp1,bp2,name,group_id) "
		    " values( :chr,:bp1,:bp2,:name,:grp) ; ");
      
      sqlite3_stmt * stmt_tmp_locus_iterate = 
	sql.prepare("SELECT chr,bp1,bp2,name "
		    " FROM locdb.loci "
		    " WHERE group_id == :grp ; ");
      
      //
      // Consider each group
      //
      
      sql.begin();
      const std::set<int> & lset = mask.included_loc();
      std::set<int>::iterator i =  lset.begin();
      while ( i != lset.end() )
	{
	  
	  std::set<std::string> names = mask.subset_loc(*i);
	  std::set<std::string> skip_names = mask.skip_loc(*i);
	  
	  sql.bind_int( stmt_tmp_locus_iterate , ":grp" , *i );
	  
	  while ( sql.step( stmt_tmp_locus_iterate ) )
	    {		    
	      
	      std::string n = sql.get_text( stmt_tmp_locus_iterate , 3 );
		  
	      //
	      // Implicit or explicit exclusion of this locus?
	      //
	      
	      if ( skip_names.find(n) != skip_names.end() ) { continue; }
	      if ( names.size() > 0 
		   && names.find(n) == names.end() ) { continue; }
	      
	      
	      int chr = sql.get_int( stmt_tmp_locus_iterate , 0 );
	      int bp1 = sql.get_int( stmt_tmp_locus_iterate , 1 );
	      int bp2 = sql.get_int( stmt_tmp_locus_iterate , 2 );
	      
	      sql.bind_int( stmt_tmp_insert, ":chr" , chr ); 
	      sql.bind_int( stmt_tmp_insert, ":bp1" , bp1 ); 
	      sql.bind_int( stmt_tmp_insert, ":bp2" , bp2 ); 
	      sql.bind_int( stmt_tmp_insert, ":grp" , *i );
	      sql.bind_text( stmt_tmp_insert, ":name" , n ); 
	    
	      sql.step( stmt_tmp_insert );
	      sql.reset( stmt_tmp_insert );
	      
	    }
	  sql.reset( stmt_tmp_locus_iterate );
	  ++i;
	}
      
      sql.commit();
      
      sql.finalise( stmt_tmp_insert );
      sql.finalise( stmt_tmp_locus_iterate );

      // sql.query("CREATE INDEX tmp.ind1 ON loc(chr,bp1,bp2); ");
      // sql.query("CREATE INDEX tmp.ind2 ON loc(name); ");
      
    }
  
      
  //
  // Any user-specified regions? 
  //
  
      
  // Either inclusion (implicit exclusion)...
      
  if ( mask.reg() )
    {
      
      //
      // Region-includes
      //
      
      sql.query(" CREATE TABLE tmp.regin "
		" ( chr INTEGER, bp1 INTEGER, bp2 INTEGER ); " );
      
      stmt_tmp_insert = 
	sql.prepare( "INSERT INTO tmp.regin(chr,bp1,bp2) "
		     " values ( :chr, :bp1, :bp2 ) ; " );
      
      const std::set<Region> & rset = mask.included_reg();
      std::set<Region>::iterator r = rset.begin();
      while ( r != rset.end() )
	{
	  sql.bind_int( stmt_tmp_insert , ":chr" , r->start.chromosome() );
	  sql.bind_int( stmt_tmp_insert , ":bp1", r->start.position() );
	  sql.bind_int( stmt_tmp_insert , ":bp2" , r->stop.position() );
	  sql.step( stmt_tmp_insert );
	  sql.reset( stmt_tmp_insert );
	  ++r;
	}
      
      sql.finalise( stmt_tmp_insert );
	  
    }
  

  //
  // Region-requires
  //
  
  if ( mask.rreg() )
    {
      
      sql.query(" CREATE TABLE tmp.regreq "
		" ( chr INTEGER, bp1 INTEGER, bp2 INTEGER ); " );
      
      stmt_tmp_insert = 
	sql.prepare( "INSERT INTO tmp.regreq(chr,bp1,bp2) "
		     " values ( :chr, :bp1, :bp2 ) ; " );
      
      const std::set<Region> & rset = mask.required_reg();
      std::set<Region>::iterator r = rset.begin();
      while ( r != rset.end() )
	{
	  sql.bind_int( stmt_tmp_insert , ":chr" , r->start.chromosome() );
	  sql.bind_int( stmt_tmp_insert , ":bp1", r->start.position() );
	  sql.bind_int( stmt_tmp_insert , ":bp2" , r->stop.position() );
	  sql.step( stmt_tmp_insert );
	  sql.reset( stmt_tmp_insert );
	  ++r;
	}
      sql.finalise( stmt_tmp_insert );
      
    }
  
  
  //
  // Region-excludes
  //
  
  if ( mask.xreg() )
    {
      
      sql.query(" CREATE TABLE tmp.regex "
		" ( chr INTEGER, bp1 INTEGER, bp2 INTEGER ); " );
      
      stmt_tmp_insert = 
	sql.prepare( "INSERT INTO tmp.regex(chr,bp1,bp2) "
		     " values ( :chr, :bp1, :bp2 ) ; " );
      
      const std::set<Region> & rset = mask.excluded_reg();
      std::set<Region>::iterator r = rset.begin();
      while ( r != rset.end() )
	{
	  sql.bind_int( stmt_tmp_insert , ":chr" , r->start.chromosome() );
	  sql.bind_int( stmt_tmp_insert , ":bp1", r->start.position() );
	  sql.bind_int( stmt_tmp_insert , ":bp2" , r->stop.position() );
	  sql.step( stmt_tmp_insert );
	  sql.reset( stmt_tmp_insert );
	  ++r;
	}
      sql.finalise( stmt_tmp_insert );
      
    }
  

  //
  // Variant ID based includes/requires/excludes
  //


  // ID-includes
      
  if ( mask.id() )
    {

      sql.query(" CREATE TABLE tmp.vidin ( name VARCHAR(20) ); " );
      
      stmt_tmp_insert = 
	sql.prepare( "INSERT INTO tmp.vidin(name) values (:name); ");
      
      const std::set<std::string> & nset = mask.included_id();
      std::set<std::string>::iterator n = nset.begin();
      while ( n != nset.end() )
	{
	  sql.bind_text( stmt_tmp_insert , ":name" , *n );
	  sql.step( stmt_tmp_insert );
	  sql.reset( stmt_tmp_insert );
	  ++n;
	}
      
      sql.finalise( stmt_tmp_insert );
	  
    }
  

  
  // ID-requires
  
  if ( mask.rid() )
    {
      
      sql.query(" CREATE TABLE tmp.vidreq ( name VARCHAR(20) ); " );
      
      stmt_tmp_insert = 
	sql.prepare( "INSERT INTO tmp.vidreq(name) values (:name); ");
      
      const std::set<std::string> & nset = mask.required_id();
      std::set<std::string>::iterator n = nset.begin();
      while ( n != nset.end() )
	{
	  sql.bind_text( stmt_tmp_insert , ":name" , *n );
	  sql.step( stmt_tmp_insert );
	  sql.reset( stmt_tmp_insert );
	  ++n;
	}
      
      sql.finalise( stmt_tmp_insert );
	  
    }


  // ID-excludes
  
  if ( mask.xid() )
    {
      
      sql.query(" CREATE TABLE tmp.videx ( name VARCHAR(20) ); " );
      
      stmt_tmp_insert = 
	sql.prepare( "INSERT INTO tmp.videx(name) values (:name); ");
      
      const std::set<std::string> & nset = mask.excluded_id();
      std::set<std::string>::iterator n = nset.begin();
      while ( n != nset.end() )
	{
	  sql.bind_text( stmt_tmp_insert , ":name" , *n );
	  sql.step( stmt_tmp_insert );
	  sql.reset( stmt_tmp_insert );
	  ++n;
	}
      
      sql.finalise( stmt_tmp_insert );
	  
    }


  
  // 
  // Locus sets?
  //
  
  if ( mask.loc_set() )
    {
      
      sql.query("CREATE TABLE tmp.locset "
		" ( chr       INTEGER, "
		"   bp1       INTEGER, "
		"   bp2       INTEGER, "
		"   name      VARCHAR(10) , "
		"   group_id  INTEGER  ) ; ");
      
      
      sqlite3_stmt * stmt_tmp_insert = 
	sql.prepare("INSERT INTO tmp.locset(chr,bp1,bp2,name,group_id) "
		    " values( :chr,:bp1,:bp2,:name,:grp) ; ");
      
      sqlite3_stmt * stmt_tmp_extract_pathway = 
	sql.prepare("SELECT l.chr,l.bp1,l.bp2,sm.name,sm.group_id "
		    " FROM locdb.loci AS l , locdb.set_data AS sd ,locdb.set_members AS sm "
		    " WHERE   sm.group_id IN ( " + mask.loc_set_include_string() + " ) "
		    "     AND sm.set_id == sd.set_id "
		    "     AND l.loc_id == sd.loc_id ; " );
      
      sql.begin();
      
      while ( sql.step( stmt_tmp_extract_pathway ) )
	{		    
	  
	  int chr = sql.get_int( stmt_tmp_extract_pathway , 0 );
	  int bp1 = sql.get_int( stmt_tmp_extract_pathway , 1 );
	  int bp2 = sql.get_int( stmt_tmp_extract_pathway , 2 );
	  std::string n = sql.get_text( stmt_tmp_extract_pathway , 3 );
	  int g = sql.get_int( stmt_tmp_extract_pathway , 4 );
	  
	  if ( mask.insert_locset( g , n ) ) 
	    {		    
	      sql.bind_int( stmt_tmp_insert, ":chr" , chr ); 
	      sql.bind_int( stmt_tmp_insert, ":bp1" , bp1 ); 
	      sql.bind_int( stmt_tmp_insert, ":bp2" , bp2 ); 
	      sql.bind_int( stmt_tmp_insert, ":grp" , g );
	      sql.bind_text( stmt_tmp_insert, ":name" , n ); 
	      
	      sql.step( stmt_tmp_insert );
	      sql.reset( stmt_tmp_insert );
	    } 
	}
      
      sql.reset( stmt_tmp_extract_pathway );
      
      sql.commit();
      
      sql.finalise( stmt_tmp_insert );
      sql.finalise( stmt_tmp_extract_pathway );

      // sql.query("CREATE INDEX tmp.ind1set ON locset(chr,bp1,bp2); ");
      // sql.query("CREATE INDEX tmp.ind2set ON locset(group_id,name); ");
    }
  
  
  
  //
  // Build a list of var_id's for variants that are required
  //
  // i.e. require:  if a var_id isn't on the list, exclude
  //      exclude:  if a var_id is on the list, exclude
  
  // CHECK -- does this work for files (i.e. they do not flag mask.requires() ?? )
  
  
  if ( mask.requires() ) 
    {

      // Requirements: you must meet *all* requirements
      
      sql.query("CREATE TABLE tmp.require ( var_id INTEGER PRIMARY KEY );");
      
      sql.query("CREATE TABLE tmp.require_loc  ( var_id INTEGER PRIMARY KEY );");
      sql.query("CREATE TABLE tmp.require_var  ( var_id INTEGER PRIMARY KEY );");
      sql.query("CREATE TABLE tmp.require_file ( var_id INTEGER PRIMARY KEY );");
      sql.query("CREATE TABLE tmp.require_reg  ( var_id INTEGER PRIMARY KEY );");
      sql.query("CREATE TABLE tmp.require_id   ( var_id INTEGER PRIMARY KEY );");
      
      int combines = 0;
      int first = 0;
      
      if ( mask.rloc() ) 
	{
	  
	  sql.query( " INSERT OR IGNORE INTO tmp.require_loc (var_id) "
		     " SELECT v.var_id FROM variants AS v , "
		     " ( SELECT chr,bp1,bp2 FROM locdb.loci WHERE group_id IN ( "
		     + mask.loc_require_string() + 
		     " ) ) AS l "
		     "  ON v.chr == l.chr AND v.bp1 BETWEEN l.bp1 AND l.bp2 ;" );
	  
	  ++combines;
	  first = 1;
	  
	}
      
      if ( mask.rvar() ) 
	{
	  
	  sql.query( " INSERT OR IGNORE INTO tmp.require_var (var_id) "
		     " SELECT var_id FROM set_data WHERE set_id IN ( "
		     + mask.var_require_string() + " ); " );

	  
	  ++combines;
	  if ( combines == 1 ) first = 2;
	  
	}
      
      if ( mask.rreg() )
	{
	  
	  sql.query( " INSERT OR IGNORE INTO tmp.require_reg (var_id) "
		     "  SELECT v.var_id FROM variants AS v INNER JOIN tmp.regreq AS zreg "
		     "     ON v.chr == zreg.chr "
		     "    AND v.bp1 BETWEEN zreg.bp1 AND zreg.bp2 ;");		
	  
	  ++combines;
	  if ( combines == 1 ) first = 3;
	  
	}
      
      
      if ( mask.rid() )
	{	  
	  sql.query( " INSERT OR IGNORE INTO tmp.require_id (var_id) "
		     "  SELECT v.var_id FROM variants AS v INNER JOIN tmp.vidreq AS z ON v.name == z.name ;" );	  
	  ++combines;
	  if ( combines == 1 ) first = 4;
	  
	}

      
      // TODO: Should be ~okay, but at some point check whethe this
      // can be done more efficiently in the main SELECT (i.e. and
      // then remove all file-specific filters from this part).  Have
      // a feeling it needs to stay here but can't recall why right
      // now.
      
      // TODO: remind oneself how tmp_reg() is meant to work
      
      if ( mask.files() ) 
	{
	  sql.query( " INSERT OR IGNORE INTO tmp.require_file (var_id) "
		     " SELECT v.var_id FROM variants AS v "
		     " WHERE v.file_id IN ( " + mask.files_include_string() + " ) ; ") ;
	  
	  ++combines;
	  if ( combines == 1 ) first = 5;
	  
	}
      
      
      // Variants that meet all requirements
      
      if ( combines == 1 ) 
	{
	  
	  if ( first == 1 ) 
	    sql.query( " INSERT OR IGNORE INTO tmp.require (var_id) "
		       " SELECT var_id FROM tmp.require_loc; " );
	  else if ( first == 2 )
	    sql.query( " INSERT OR IGNORE INTO tmp.require (var_id) "
		       " SELECT var_id FROM tmp.require_var; " );
	  else if ( first == 3 ) 
	    sql.query( " INSERT OR IGNORE INTO tmp.require (var_id) "
		       " SELECT var_id FROM tmp.require_reg; " );
	  else if ( first == 4 ) 
	    sql.query( " INSERT OR IGNORE INTO tmp.require (var_id) "
		       " SELECT var_id FROM tmp.require_id; " );	  
	  else if ( first == 5 ) 
	    sql.query( " INSERT OR IGNORE INTO tmp.require (var_id) "
		       " SELECT var_id FROM tmp.require_file; " );
	}
      else
	{
	  std::string q = " INSERT OR IGNORE INTO tmp.require (var_id) ";
	  if ( first == 1 )      q+= " SELECT f.var_id FROM tmp.require_loc   AS f ";
	  else if ( first == 2 ) q+= " SELECT f.var_id FROM tmp.require_var   AS f ";
	  else if ( first == 3 ) q+= " SELECT f.var_id FROM tmp.require_reg   AS f ";
	  else if ( first == 4 ) q+= " SELECT f.var_id FROM tmp.require_id    AS f ";
	  else if ( first == 5 ) q+= " SELECT f.var_id FROM tmp.require_file  AS f ";
	  
	  if ( mask.rvar() && first < 2 ) q += " INNER JOIN tmp.require_var ";
	  if ( mask.rreg() && first < 3 ) q += " INNER JOIN tmp.require_reg ";
	  if ( mask.rid() && first < 4 )  q += " INNER JOIN tmp.require_id ";
	  if ( mask.files() && first < 5 ) q += " INNER JOIN tmp.require_file ";
	  
	  q += " USING ( var_id ) ; ";
	  
	  sql.query( q );
	  
	}
      
    }
  

  //
  // Excludes
  //
  
  if ( mask.excludes() ) 
    {
      
      sql.query("CREATE TABLE tmp.exclude ( var_id INTEGER PRIMARY KEY );");
      
      if ( mask.xloc() ) 
	{
	  sql.query( " INSERT OR IGNORE INTO tmp.exclude (var_id) "
		     " SELECT v.var_id FROM variants AS v , "
		     " ( SELECT chr,bp1,bp2 FROM locdb.loci WHERE group_id IN ( "
		     + mask.loc_exclude_string() + 
		     " ) ) AS l "
		     "  ON v.chr == l.chr AND v.bp1 BETWEEN l.bp1 AND l.bp2 ;" );
	}
      
      if ( mask.xvar() ) 
	{
	  sql.query( " INSERT OR IGNORE INTO tmp.exclude (var_id) "
		     " SELECT var_id FROM set_data WHERE set_id IN ( "
		     + mask.var_exclude_string() + " ); ");
	}
      
      if ( mask.xreg() )
	{
	  sql.query( " INSERT OR IGNORE INTO tmp.exclude (var_id) "
		     "  SELECT v.var_id FROM variants AS v INNER JOIN tmp.regex AS zreg "
		     "     ON v.chr == zreg.chr "
		     "    AND v.bp1 BETWEEN zreg.bp1 AND zreg.bp2 ;");		
	}
      
      
      if ( mask.xfiles() ) 
	{
	  sql.query( " INSERT OR IGNORE INTO tmp.exclude (var_id) "
		     " SELECT var_id FROM variants AS v "
		     " WHERE v.file_id IN ( " + mask.files_exclude_string() + " ) ; ") ;
	}
      
    }
  


  //
  // Grouping variable
  //


  if ( mask.named_grouping() ) 
    {
      
      sql.query("CREATE TABLE tmp.grp ( grp VARCHAR(20) NOT NULL );");

      if ( mask.group_loc() ) 
	{
	  
	  if ( mask.loc_exceptions() )
	    {
	      std::string q = " INSERT OR IGNORE INTO tmp.grp (grp) "
		" SELECT DISTINCT t.name FROM tmp.loc AS t "
		" WHERE t.group_id == " + Helper::int2str( mask.group_set() ) + " ; " ;
	      sql.query( q ) ;
	    }
	  else
	    {
	      sql.query( " INSERT OR IGNORE INTO tmp.grp (grp) "
			 " SELECT DISTINCT name FROM locdb.loci "
			 " WHERE group_id == " + Helper::int2str( mask.group_set() ) + " ; ") ;
	    }
	}
      else if ( mask.group_loc_set() )
	{
	  sql.query( " INSERT OR IGNORE INTO tmp.grp (grp) "
		     " SELECT DISTINCT name FROM tmp.locset AS t "
		     " WHERE t.group_id == " + Helper::int2str( mask.group_set() ) + " ; ") ;
	}
      else if ( mask.group_var_set() )
	{
	  sql.query( " INSERT OR IGNORE INTO tmp.grp (grp) "
		     " SELECT name FROM sets WHERE set_id IN ( SELECT DISTINCT set_id FROM superset_data AS s "
		     " WHERE s.superset_id == " + Helper::int2str( mask.group_set() ) + " ) ;" );
	  
	}
    }
  
}
