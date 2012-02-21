
#include "index.h"
#include "cgi.h"
#include "autohtml.h"

#include <set>
#include <iostream>


void show_gene_list( const std::string & loc_set , GStore & g );
void inserter( std::set<std::string> & s , const std::string & t );
void f_display( Variant & , void * );
std::string pp(const std::string & str , const int len = 22 );

struct f_aux 
{ 
  std::set<std::string> meta_var;
  std::set<std::string> meta_geno;
  std::set<std::string> meta_indiv;
  bool show_geno;
  bool show_alt_only;
  int vcounter;
  std::vector<std::string> lines;  

  f_aux() 
  { 
    show_geno = true;
    show_alt_only = true;
    vcounter = 0;
    lines.clear();
  }

};



int main()
{
  
  
  //
  // Core parameters
  //

  //
  // Project path
  //

  std::string project_path = "";


  //
  // Password
  //

  std::string pwd = "(if required)";
  

  //
  // Use defaults Default locus-set ("refseq", "symbol")
  //

  std::string loc_set = PLINKSeq::DEFAULT_LOC_GROUP();
  std::string gene_symbol = PLINKSeq::DEFAULT_GENE_SYMBOL();
  
  
  //
  // Get CGI variables from POST
  //

  //  p = project
  //  v = var-tags
  //  g = genes
  //  i = individuals
  //  r = regions
  //  m = meta-information
  //  q = 'v' or 'i', 's' (variant or individual or set query)

  std::string project = "";
  std::string query_mode = "";
  std::set<std::string> vtags_sel;
  std::set<std::string> genes_sel;
  std::set<std::string> indiv_sel;
  std::set<std::string> meta_sel;
  std::set<Region>      regs_sel;
  std::string genestring = "";

  // some global params
  f_aux aux;

  char **cgivars = getcgivars() ;
  
  for (int i=0; cgivars[i]; i+= 2)
    {
      
      std::string str = cgivars[i];
      
      if      ( str == "vmode" ) query_mode = "v";
      else if ( str == "smode" ) query_mode = "s";
      else if ( str == "imode" ) query_mode = "i";
      else if ( str == "gmode" ) query_mode = "g";
      else if ( str == "showgeno" ) aux.show_geno = true;
      else if ( str == "showaltonly" ) aux.show_alt_only = true;
      else
	{

	  // expecting a key=value pair
	  std::string param = cgivars[i+1];
	  
	  if      ( str == "p"          ) project = param;
	  else if ( str == "genestring" ) genestring = param;

	  else if ( str == "vtags" ) inserter( vtags_sel , param );	  
	  else if ( str == "indiv" ) inserter( indiv_sel , param );
	  else if ( str == "meta"  ) inserter( meta_sel , param );

	}
 
    }

  for (int i=0; cgivars[i]; i++) free(cgivars[i]) ;
  free(cgivars) ;
  

  //
  // Create main page
  //


  HTMLHelper::header( "IndEx" );  
  HTMLHelper::add_css1();
  HTMLHelper::init_combocheckbox();  
  HTMLHelper::init_filterlist();
  HTMLHelper::end_header();
  
  std::cout << "<h2><b>Ind</b>ividual <b>ex</b>ome browser</h2>";


  if ( query_mode == "" )  
    HTMLHelper::end_page( "invalid query mode" );
  

  //
  // Attach project 
  //

  GStore g;  
  g.set_project( project ) ;
  
  if ( ! g.vardb.attached() ) HTMLHelper::end_page( "could not attach project [" + project +"]" );
  if ( ! g.locdb.attached() ) HTMLHelper::end_page( "no locus database attached to project" );


  //
  // Populate entire filter sets (i.e. indiv, genes, var-sets and meta-tags)
  //


  // All possible values
  
  std::set<std::string> indiv;
  std::set<std::string> genes;
  std::set<std::string> vtags;
  std::set<std::string> meta;
  std::map<std::string,std::string> meta_desc;

  // 1) Individuals (from all files)
  
  std::map<int,std::string> files = g.vardb.fetch_files();
  std::map<int,std::string>::iterator i = files.begin();
  while ( i != files.end() )
    {    
      std::vector<std::string> iv = g.vardb.fetch_individuals( i->first );
      for (int j=0; j< iv.size(); j++) indiv.insert( iv[j] );
      ++i;
    }
  
  // 2) Variant tags
  
  std::vector<std::string> sets = g.vardb.get_sets();
  std::vector<std::string> supersets = g.vardb.get_supersets();  
  for (int s = 0 ; s < sets.size() ; s++) vtags.insert( sets[s] );
  for (int s = 0 ; s < supersets.size() ; s++) vtags.insert( supersets[s] );
  

  // 3) Genes
  
  // (no universal list, but parse genestring here --> genes_sel, regs_sel 
  
  std::vector<std::string> tok = Helper::whitespace( genestring );
  
  for (int i=0; i<tok.size(); i++)
    {      
      bool okay = true;      
      Region r(tok[i],okay);
      if ( okay ) regs_sel.insert( r );	  
      else
	{
	  std::string genename = tok[i];
	  Helper::str2upper( genename );
	  genes_sel.insert( genename );
	}
    }

    
  // 4) Meta-fields: Indiv, Variant , Genotype 

  i = files.begin();
  while ( i != files.end() )
    {
      std::vector<std::map<std::string,std::string> > m = g.vardb.fetch_metatypes( i->first );  
      for (int j=0; j<m.size(); j++)
	{
	  std::string name = m[j]["NAME"];
	  std::string grp = m[j]["GRP"];
	  std::string label = name + " | " + grp + " | " + pp( m[j]["DESC"] ) ;
	  meta.insert(name);
	  meta_desc[name] = label;
	  
	  if      ( grp == "Variant" && meta_sel.find(name) != meta_sel.end() ) aux.meta_var.insert(name);
	  else if ( grp == "Genotype" && meta_sel.find(name) != meta_sel.end()) aux.meta_geno.insert(name);
	  else if ( grp == "Individual" && meta_sel.find(name) != meta_sel.end() ) aux.meta_indiv.insert(name);	  
	  
       }
     ++i;
   }
   

  //
  // Four main selection tables
  //

  //                                                  all       current-selection
  //  -----------------------------------------------------------------------------
  //  1) family/individual list (combocheckbox)       indiv     indiv_sel
  //  2) var-tags and supersets (combocheckbox)       vtags     vtags_sel
  //  3) region filter (textbox)                      n/a       regs_sel
  //  4) gene filter (filterlist)                     genes     genes_sel
  //  5) meta-information display (combocheckbox)     meta      meta_sel
  
  
  
  //
  // Create main form
  //

  HTMLHelper::form( "myform" , "./index.cgi" );
  
  //
  // Attach basic information
  //

  HTMLHelper::add_hidden_value( "p" , project );


  //
  // Across the top of the form, have the primary buttons: 
  //
  
  std::cout << "<table border=0 width=100%>";
  
  std::cout << " <tr><td colspan=\"4\" width=100% align=\"center\">"
	    << " <input type=\"submit\" name=\"vmode\" value=\"FILTER variants\"> "
	    << " <input type=\"submit\" name=\"imode\" value=\"List INDIVIDUALs\"> "
	    << " <input type=\"submit\" name=\"smode\" value=\"Display variant SETS\"> "
	    << " <input type=\"submit\" name=\"gmode\" value=\"Display all GENEs\"> ";

  std::cout << "<br>";

  std::string c1 = aux.show_geno ? "CHECKED" : "";
  std::cout << "<input type=\"checkbox\" name=\"showgeno\" " << c1 << "> Show individual genotypes ";

  c1 = aux.show_alt_only ? "CHECKED" : "";
  std::cout << "<input type=\"checkbox\" name=\"showaltonly\" " << c1 << "> Only list non-reference sites ";
  
  std::cout << " <br> ";
  
  std::cout << "</td></tr>";


  // 
  // Main boxes for selections
  //

  std::cout << "<tr>";

  // Var-lists
  std::cout << "<td width=25%>";
  HTMLHelper::add_combocheckbox( "vtags" , vtags , vtags_sel );
  std::cout << "</td>";  

  // Individuals
  std::cout << "<td width=25%>";
  HTMLHelper::add_combocheckbox( "indiv" , indiv , indiv_sel );
  std::cout << "</td>";
    
  // Gene/region list
  std::cout << "<td width=25%>";

  std::cout << "<textarea ";
  std::cout << "name=\"genestring\" rows=\"10\" cols=\"20\">"
	    << genestring
	    << "</textarea></p>";
  std::cout << "</td>";
  
  // Meta-tags
  std::cout << "<td width=25%>";
  HTMLHelper::add_combocheckbox( "meta" , meta_desc , meta_sel );
  std::cout << "</td>"; 
  
  std::cout << "</tr>";


  //
  // DONE 
  //

  std::cout << "</table>";
  HTMLHelper::end_form();

  
  
  //
  // Now display whatever it is we need to display
  //


  //
  // Gene-list mode
  //

  if ( query_mode == "g" )
    {
      show_gene_list( loc_set , g );
      HTMLHelper::end_page();
    }


  //
  // Set query mode
  //

  if ( query_mode == "s" ) 
    {

      std::vector<std::string> sets = g.vardb.get_sets();
      for (int s = 0 ; s < sets.size(); s++)
	std::cout << sets[s] << " " << g.vardb.get_set_size( sets[s] ) << "<br>";
      std::cout << "<hr>";
      
      std::vector<std::string> ss = g.vardb.get_supersets();
      for (int s = 0 ; s < ss.size(); s++)
	{
	  std::cout << "<b>" << ss[s] << "</b> : ";
	  std::vector<std::string> sets = g.vardb.get_sets( ss[s] );
	  for (int t = 0 ; t < sets.size(); t++)
	    std::cout << " " << sets[t] ;
	  std::cout << "<br>";
	}
      std::cout << "<hr>";
	    
      HTMLHelper::end_page();
    }
  
  
  //
  // Variant list
  //

  if ( query_mode == "v" )
    {

      // Create mask

      Mask mask;
      
      // filter individuals
      std::set<std::string>::iterator ii = indiv_sel.begin();
      while ( ii != indiv_sel.end() )
	{
	  mask.include_indiv( *ii );
	  ++ii;
	}
      
      bool anything_to_show = false;
      
      // filter var-sets
      std::set<std::string>::iterator vi = vtags_sel.begin();
      while ( vi != vtags_sel.end() )
	{
	  anything_to_show = true;
	  mask.require_var( *vi );
	  ++vi;
	}

      if ( ! anything_to_show ) 
	HTMLHelper::end_page( "no variants selected" );
      
      // speed-up if no individual information is required
      if ( indiv_sel.size() == 0 ) 
 	{
 	  mask.load_genotype_data( false );      
 	  mask.load_genotype_meta( false );
 	}
      

      // register mask

      g.indmap.populate( g.vardb, g.phmap, mask );

      // append meta-information 
      
      g.vardb.iterate( f_display , &aux , mask );

      // display
      
      std::cout << "<br><em>Found " << aux.vcounter << " variants that match criteria</em><br>&nbsp;<br>";

      std::cout << "<table border=1><tr align=\"center\"><th>Variant</th><th>REF/ALT</th>";

      // meta-information
      std::set<std::string>::iterator mi = meta_sel.begin();
      while ( mi != meta_sel.end() )
	{
	  std::cout << "<th>" << *mi << "</th>";
	  mi++;
	}

      // genotypes
      if ( aux.show_geno ) 
	{
	  const int n = g.indmap.size();
	  for (int i=0;i<n;i++) std::cout << "<th>" << g.indmap(i)->id() << "</th>";
	}

      std::cout << "</tr>";
    
      // table rows
      for (int v=0;v<aux.lines.size();v++)
	std::cout << aux.lines[v] ;
      
      std::cout << "</table>";
      
    }

  
  //
  // Display all phenotype information on selected individual(s)
  //
  
  if ( query_mode == "i" ) 
    {
      std::set<std::string>::iterator ii = indiv_sel.begin();
      while ( ii != indiv_sel.end() )
	{
	  std::cout << "Individual " << *ii << "<br>";
	  std::cout << "<hr>";
	  ++ii;
	}
    }


  // All done
  HTMLHelper::end_page();
    

}


void show_gene_list( const std::string & loc_set , GStore & g )
{
  std::set<std::string> genes = g.locdb.fetch_names( loc_set );
  std::set<std::string>::iterator ig= genes.begin();
  while ( ig != genes.end() ) 
    {
      std::cout << "gene " << *ig << "<br>";
      ++ig;
    }     
}

void inserter( std::set<std::string> & s , const std::string & str )
{
  std::vector<std::string> t = Helper::char_split( str , ',' );
  for (int i=0; i<t.size(); i++) s.insert(t[i]);
}

std::string pp(const std::string & str , const int len )
{
  if ( str.size() < len ) return str;
  return "<abbr title=\"" + str + "\">" + str.substr(0,len) + "...</abbr>";
}

void f_display( Variant & v , void * p )
{

  // main variant display function  
  
  f_aux * aux = (f_aux*)p;


  //
  // Only display variants with 1 or more alternate allele?
  //

  if ( aux->show_alt_only )
    {
      int c, c_tot;
      bool altmin = v.n_minor_allele( &c , &c_tot );
      if ( altmin && c == 0 ) return;
      if ( c == c_tot && !altmin ) return;
    }

  //
  // Generate display line for table
  //

  std::stringstream ss;

  ss << "<tr align=\"center\"><td>" << v << "</td><td>" << v.consensus.reference() << "/" << v.consensus.alternate() << "</td>";

  // meta-information
  std::set<std::string>::iterator mi = aux->meta_var.begin();
  while ( mi != aux->meta_var.end() )
    {
      meta_name_t n = *mi;
      if ( v.meta.has_field( *mi ) )
	ss << "<td>" << v.meta.print( n ) << "</td>";
      else if ( v.consensus.meta.has_field( *mi ) )
	ss << "<td>" << v.consensus.meta.print( n ) << "</td>";
      else 
	ss << "<td>.</td>";
      ++mi;
    }

  // genotypes
  if ( aux->show_geno )
    {
      const int n = v.size();
      for (int i=0; i<n; i++)
	{

	  Genotype & g = v(i);

	  // color-coding
	  if ( g.null() ) ss << "<td align=\"center\" bgcolor=\"gray\">";
	  else if ( g.reference() ) ss << "<td align=\"center\" bgcolor=\"lightgreen\">";
	  else if ( g.heterozygote() ) ss << "<td align=\"center\" bgcolor=\"red\">";
	  else if ( g.alternate_homozygote() ) ss << "<td align=\"center\" bgcolor=\"yellow\">";

	  // actual genotype
	  ss << ( v.flat() ? v.geno_label( v(i) ) : v.label(i) ) ;
	  
	  // genotype meta-information?

	  if ( v.flat() )
	    {
	      std::set<std::string>::iterator k = aux->meta_geno.begin();
	      while ( k != aux->meta_geno.end() )
		{
		  bool prt = false;
		  if ( g.meta.has_field( *k ) )
		    {
		      meta_name_t t = *k;
		      std::string value = g.meta.print( t );
		      if ( value == "" || value == " " ) value = ".";
		      if ( prt ) ss << ";"; else { ss << "<br>"; prt = true; }
		      ss << *k << "=" << value;
		    }
		  ++k;
		}
	    }  
	  else
	    HTMLHelper::end_page( "index does not currently handle projects in which individuals feature more than once" );

// 	  else
// 	    ss << "<td>" << var.gmeta_label(i) << "</td>";
	  
	  // done
	  ss << "</td>";
	}

    }

  ss << "</tr>";

  //
  // All done
  //

  aux->vcounter++;  
  aux->lines.push_back( ss.str() );

}
