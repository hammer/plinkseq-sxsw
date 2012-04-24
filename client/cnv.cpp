#include "cnv.h"
#include "util.h"

#include <iostream>
using namespace std;

extern GStore g;
extern Pseq::Util::Options args;

void f_cnv_denovo_scan(Variant& v, void* p) {

	cout << v.coordinate() << endl;

	const int n = v.size();

	for (int i = 0; i < n; ++i) {
		Individual* pat = v.ind(i)->pat();
		Individual* mat = v.ind(i)->mat();

		// only consider individuals where both parents are present
		if (pat == NULL || mat == NULL)
			continue;

		// there will be a better way to get the actual individual slot
		// but use for now...

		int patn = g.indmap.ind_n(pat->id());
		int matn = g.indmap.ind_n(mat->id());

		Genotype& gChild = v(i);
		Genotype& gPatern = v(patn);
		Genotype& gMatern = v(matn);

		string childID = v.ind(i)->id();
		cout << "child= " << childID << endl;

		double childSQ = gChild.meta.get_double("SQ")[0];

		cout << "childSQ= " << childSQ << endl;

		double paternNQ = gPatern.meta.get_double("NQ")[0];

		cout << "paternNQ= " << paternNQ << endl;

		double maternNQ = gMatern.meta.get_double("NQ")[0];

		cout << "maternNQ= " << maternNQ << endl;

		if (childSQ >= 60 && paternNQ >= 60 && maternNQ >= 60) {
			cout << "child= " << childID << " del_CNV= " << v.coordinate() << endl;
		}

		bool denovo = false;
	}
}

bool Pseq::VarDB::cnv_denovo_scan(Mask& mask) {
	cerr << "Starting CNV de novo scan..." << endl;

	// did we have some special values
	const int n = g.indmap.size();

	// Attach parents
	for (int i = 0; i < n; ++i) {
		Individual* person = g.indmap(i);
		g.inddb.fetch(person);
		Individual* p = g.indmap.ind(person->father());
		Individual* m = g.indmap.ind(person->mother());
		if (p) person->pat(p);
		if (m) person->mat(m);
	}

	g.vardb.iterate(f_cnv_denovo_scan, (void*) NULL, mask);

	return true;
}
