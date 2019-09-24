// GARFIELD - GWAS analysis of regulatory or functional information enrichment with LD correction.
// Copyright (C) 2014 Wellcome Trust Sanger Institute / EMBL - European Bioinformatics Institute
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


using namespace std;
string chrom;
string output = "NULL";

// define type of snpid
// snpid_type must support:
// -stream operators >>,<<
// -comparison <,>,==,!=
// implement parse_snpid specilization
// input files must be sorted with respect to snpid_type
typedef unsigned int snpid_type;
// typedef string snpid_type;
template<class T> T parse_snpid(const string& a);
template<> string parse_snpid<string>(const string& a) {return a;}
template<> unsigned int parse_snpid<unsigned int>(const string& a) {return atoi(a.c_str());}
vector<int> exclude;
int iter=0;

// timing 

clock_t tic() {
	return clock();
}

clock_t toc(clock_t t, const string& a = string()) {
	cerr << "Time elapsed: " << 1.0 * (clock() - t) / CLOCKS_PER_SEC << "s";
	if (!a.empty()) cerr << " (" << a << ")";
	cerr << endl;
	return clock();
}


// open file checks
void safefileopen(ifstream& ifs, const string& fn) {
	ifs.open(fn.c_str());
	if (!ifs.is_open()) {
		ostringstream oss;
		oss << "Could not open file \"" << fn << "\"!" << endl;
		throw runtime_error(oss.str());
	}
}


// read utilities

// need pair >> for reading table fields
namespace std {
	template<class T1, class T2> istream & operator>>(istream & is, pair<T1,T2> & p) {return is >> p.first >> p.second;}
	template<class T1, class T2> ostream & operator<<(ostream & os, const pair<T1,T2> & p) {return os << p.first << p.second;}
}

// helper class for reading a tags field
struct Tags : vector<snpid_type> {
	friend istream & operator>>(istream & is, Tags & tags) {
		tags.clear();
		string a;
		if (is >> a) {
			// copy(istream_iterator<snpid_type>(istringstream(a)), istream_iterator<snpid_type>(), tags);
			istringstream iss(a);
			string x;
			while (getline(iss, x, ',')) tags.push_back(parse_snpid<snpid_type>(x));
			// while (getline(iss, x, ',')) {
			// 	istringstream iss(x);
			// 	snpid_type y;
			// 	if (iss >> y) tags.push_back(y);
			// }
		}
		return is;
	}
};

// helper class for reading an annotation field
struct Annot : vector<bool> {
	friend istream & operator>>(istream & is, Annot & annot) {
		annot.clear();
		string a;
		if (is >> a) {
			annot.resize(a.length());
			for (Annot::size_type i=0; i<a.length(); ++i) annot[i]=a[i]=='1';
	
		}
		return is;
	}
};


// join utilities

// join iterators based on snpid
// assume input iterators sorted by snpid
// iterate through [first1,last1) and [first2,last2) simultaneously
// for every snpid match, insert pair<SnpIdIterator1,SnpIdIterator2> into out
template<class JoinOutputIterator, class SnpIdIterator1, class SnpIdIterator2>
void join(
JoinOutputIterator out,
SnpIdIterator1 first1, SnpIdIterator1 last1,
SnpIdIterator2 first2, SnpIdIterator2 last2)
{
	// assume input iterators sorted by id
	SnpIdIterator1 it1 = first1;
	SnpIdIterator2 it2 = first2;
	while (it1 != last1 && it2 != last2) {
		if (it1->snpid() < it2->snpid()) ++it1;
		else if (it2->snpid() < it1->snpid()) ++it2;
		else {*out++ = make_pair(it1++, it2++);}
	}
}

// count iterator steps, available as member variable i
template<class Iterator>
struct Enumerator : public iterator<input_iterator_tag, typename Iterator::value_type>  {
	size_t i;
	Iterator it;
	Enumerator() : i(0), it() {}
	Enumerator(const Iterator & it) : i(0), it(it) {}
	bool operator==(const Enumerator& a) const {return this->it==a.it;}
	bool operator!=(const Enumerator& a) const {return this->it!=a.it;}
	Enumerator& operator++() {++this->i; ++this->it; return *this;}
	Enumerator operator++(int) {Enumerator tmp(*this); this->operator++(); return tmp;}
	const typename Iterator::value_type& operator*() const {return *this->it;}
	const typename Iterator::value_type* operator->() const {return &*this->it;}
};

// snpid pval container
struct snpid_pval : public pair<snpid_type, double>	{
	snpid_pval() : pair<snpid_type, double>() {}
	snpid_pval(const snpid_type& snpid, const double& pval)
	{this->first=snpid; this->second=pval;}
	snpid_type snpid() const {return this->first;}
	void snpid(const snpid_type& a) {this->first=a;}
	double pval() const {return this->second;}
	void pval(const double& a) {this->second=a;}
	operator snpid_type() const {return this->first;}
};


// snpid containers
// represent data tuples, expose snpid

// snpid pval tags container
struct snpid_pval_tags : public pair<snpid_type, pair<double, Tags> >	{
	snpid_pval_tags(const snpid_pval& a) {this->snpid(a.snpid()); this->pval(a.pval());} 
	// snpid_pval_tags ( const snpid_type& snpid, const double& pval, const Tags& tags )
	// {this->snpid()=snpid; this->pval()=pval; this->tags()=tags;}
	snpid_type snpid() const {return this->first;}
	void snpid(const snpid_type& a) {this->first=a;}
	double pval() const {return this->second.first;}
	void pval(const double& a) {this->second.first=a;}
	const Tags& tags() const {return this->second.second;}
	void tags(const Tags& a) {this->second.second = a;}
	operator snpid_type() const {return this->first;}
	static bool compare_by_pval(const snpid_pval_tags& a, const snpid_pval_tags& b)
	{return a.pval() < b.pval();}
};

// snpid tags container
struct snpid_tags : public pair<snpid_type, Tags> {
	snpid_type snpid() const {return this->first;}
	const Tags& tags() const {return this->second;}
	void tags(const Tags& a) {this->second = a;}
	operator snpid_type() const {return this->first;}
};

// snpid vectorindex container
struct snpid_index : public pair<snpid_type, size_t> {
	snpid_type snpid() const {return this->first;}
	void snpid(const snpid_type& a) {this->first=a;}
	size_t index() const {return this->second;}
	void index(const size_t& a) {this->second=a;}
	operator snpid_type() const {return this->first;}
};

// snpid annotation container
struct snpid_annot : public pair<snpid_type, Annot>	{
	snpid_type snpid() const {return this->first;}
	const Annot& annot() const {return this->second;}
	void annot(const Annot& a) {this->second = a;}
	operator snpid_type() const {return this->first;}
};

// snpid maf tssd container
struct snpid_maf_tssd : public pair<snpid_type, pair<double,int> > {
	snpid_type snpid() const {return this->first;}
	void snpid(const snpid_type& a) {this->first=a;}
	double maf() const {return this->second.first;}
	int tssd() const {return this->second.second;}
	operator snpid_type() const {return this->first;}
};


// JoinOutputIterator assigning tags to vector<snpid_pval_tags>
struct TagsOutputIterator1 : public iterator<output_iterator_tag,void,void,void,void> {
	TagsOutputIterator1& operator=(const pair<const vector<snpid_pval_tags>::iterator&, const istream_iterator<snpid_tags>&>& a)
	{a.first->tags(a.second->tags()); return *this;}
	TagsOutputIterator1& operator++() {return *this;}
	TagsOutputIterator1& operator++(int) {return *this;}
	TagsOutputIterator1& operator*() {return *this;}
};

// JoinOutputIterator assigning tags to vector<Tags>
struct TagsOutputIterator2 : public iterator<output_iterator_tag,void,void,void,void> {
	vector<Tags>* v;
	TagsOutputIterator2(vector<Tags>& v) : v(&v) {}
	TagsOutputIterator2& operator=(const pair<const Enumerator<vector<snpid_pval>::iterator>&, const istream_iterator<snpid_tags>&>& a)
	{(*v)[a.first.i] = a.second->tags(); return *this;}
	TagsOutputIterator2& operator++() {return *this;}
	TagsOutputIterator2& operator++(int) {return *this;}
	TagsOutputIterator2& operator*() {return *this;}
};

// JoinOutputIterator storing annotation
struct TagRefOutputIterator : public iterator<output_iterator_tag,void,void,void,void> {
	vector<Annot>* v;
	TagRefOutputIterator(vector<Annot>& v) : v(&v) {}
	TagRefOutputIterator& operator=(const pair<const vector<snpid_index>::iterator&, const istream_iterator<snpid_annot>&>& a)
	{
		Annot& annot1 = (*v)[a.first->index()];
		Annot annot2 = a.second->annot();
		Annot::size_type i = 0;
		//for (; i<annot1.size(); ++i) annot1[i] = annot1[i] | annot2[i];
		for (; i<annot1.size(); ++i)
		{	
			iter = 0;
			for (int j=0; j<exclude.size(); ++j)
			{
				if (i==exclude[j]) iter++;
			}			
			if (iter==0) 
				annot1[i] = annot1[i] | annot2[i];
		}
		annot1.insert(annot1.end(), annot2.begin()+i, annot2.end());
		return *this;
	}
	TagRefOutputIterator& operator++() {return *this;}
	TagRefOutputIterator& operator++(int) {return *this;}
	TagRefOutputIterator& operator*() {return *this;}
};

// JoinOutputIterator printing data table
struct PrintOutputIterator : public iterator<output_iterator_tag,void,void,void,void> {
	vector<Tags>* vtags;
	vector<Annot>* vannot;
	Annot::size_type nannot;
	PrintOutputIterator(vector<Tags>& vtags, vector<Annot>& vannot,string output)	: vtags(&vtags), vannot(&vannot), nannot(0) {}
	PrintOutputIterator& operator=(const pair<const Enumerator<vector<snpid_pval>::iterator>&, const istream_iterator<snpid_maf_tssd>&>& a)
	{	
		ofstream out(output.c_str(), std::ios_base::app);
		Annot& annot = (*vannot)[a.first.i];
		if (annot.size() != 0) {
			if (nannot == 0) nannot = annot.size();
			if (annot.size() != nannot) throw runtime_error("Error: Annotation lengths do not match!");
			out << chrom <<" "<<a.first->snpid() << " " << a.first->pval()
			<< " " << (*vtags)[a.first.i].size() << " " << a.second->maf()
			<< " " << a.second->tssd() << " ";
			copy(annot.begin(), annot.end(), ostream_iterator<bool>(out));
			out << endl;
		}
		return *this;
	}
	PrintOutputIterator& operator++() {return *this;}
	PrintOutputIterator& operator++(int) {return *this;}
	PrintOutputIterator& operator*() {return *this;}
};




int main(int argc, char* argv[]) {
	clock_t tstart = tic();
	clock_t t = tic();

	// parse options

	if (argc < 6) {
		cerr << "Usage: prep prunetags clumptags maftssd gwas annot" << endl;
		return(1);
	}
	
	// store command line arguments
	string fn_tags_prune = "NULL"; // id t1,t2,..
	string fn_tags_clump = "NULL"; // id t1,t2,..
	string fn_maftssd = "NULL"; // id maf tssd
	string fn_pval0 = "NULL", fn_annot0 = "NULL";
	// file must be sorted by snpid with respect to snpid_type

	chrom = "NA";
	string indices = "0";
	//string output = "NULL";
        
    int i;
    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "-ptags") == 0) {
            fn_tags_prune = argv[i + 1]; i++;
        } else if (strcmp(argv[i], "-ctags") == 0) {
            fn_tags_clump = argv[i + 1]; i++;
        } else if (strcmp(argv[i], "-maftss") == 0) {
            fn_maftssd = argv[i + 1]; i++;
        } else if (strcmp(argv[i], "-pval") == 0) {
            fn_pval0 = argv[i + 1]; i++;
        } else if (strcmp(argv[i], "-ann") == 0) {
            fn_annot0 = argv[i + 1]; i++;
        } else if (strcmp(argv[i], "-chr") == 0) {
            chrom = argv[i + 1]; i++;
        } else if (strcmp(argv[i], "-excl") == 0) {
        	indices = argv[i + 1]; i++;
        } else if (strcmp(argv[i], "-o") == 0) {
            output = argv[i+1];
        }
    }

    if (fn_tags_prune.c_str()=="NULL") {
        cerr << "Error: No pruning LD tags file specified! Use -ptags option." << endl;
        return(1);
    }    
    if (fn_tags_clump.c_str()=="NULL") {
        cerr << "Error: No annotation LD tag file specified! Use -ctags option." << endl;
        return(1);
    }
    if (fn_maftssd.c_str()=="NULL") {
        cerr << "Error: No MAF, TSS distance file specified! Use -maftss option." << endl;
        return(1);
    }
    if (fn_pval0.c_str()=="NULL") {
        cerr << "Error: No pvalue file specified! Use -pval option." << endl;
        return(1);
    }
    if (fn_annot0.c_str()=="NULL") {
        cerr << "Error: No annotation file specified! Use -ann option." << endl;
        return(1);
    }
    if (output.c_str()=="NULL") {
        cerr << "Error: No output file specified! Use -o option." << endl;
        return(1);
    }
    
    //ofstream out(output.c_str());
    //if (!out.is_open()) {
    //    cerr << "Could not open output file!" << endl;
    //    return(1);
    //}
    string fn_pval(fn_pval0); // id pval
	string fn_annot(fn_annot0); // id a1a2...
	stringstream sss(indices);
	string xx;
	while (getline(sss, xx, ',')) {
		exclude.push_back(atoi(xx.c_str()));
	}
	// prune

	vector<snpid_pval> vsnpidpval;
	{
		vector<snpid_pval_tags> vsnpidpvaltags;
		{
            clock_t t = tic();
			{
				//read pval file
				ifstream ifs;
				safefileopen(ifs, fn_pval);
				{
					typedef istream_iterator<snpid_pval> snpidpvalit;
					copy(snpidpvalit(ifs), snpidpvalit(), back_inserter(vsnpidpvaltags));
				}
			}
            t = toc(t, "read pvals");

			{
				ifstream ifs;
				safefileopen(ifs, fn_tags_prune);
				{
					typedef istream_iterator<snpid_tags> it;
					join(TagsOutputIterator1(), vsnpidpvaltags.begin(), vsnpidpvaltags.end(), it(ifs), it());
				}
			}
            t = toc(t, "join tags");
		}

		// go through snps with tags from smallest pvalue to highest pvalue
		// throw away snps that are tags of already seen snps
		sort(vsnpidpvaltags.begin(), vsnpidpvaltags.end(), snpid_pval_tags::compare_by_pval);
		set<snpid_type> sp;
		for (vector<snpid_pval_tags>::iterator it = vsnpidpvaltags.begin(); it != vsnpidpvaltags.end(); ++it) {
			if (sp.find(it->snpid()) == sp.end()) {
				vsnpidpval.push_back(snpid_pval(it->snpid(), it->pval()));
				sp.insert(it->tags().begin(), it->tags().end());
			}
		}

		sort(vsnpidpval.begin(), vsnpidpval.end());

		t = toc(t, "pruning");
	} // --> vsnpidpval contains pruned snp set with pvalues


	// clump

	vector<Tags> vtags;
	vector<Annot> vannot;
	{
		// read tags
		{
			ifstream ifs;
			safefileopen(ifs, fn_tags_clump);
			vtags = vector<Tags>(vsnpidpval.size());
			{
				typedef Enumerator<vector<snpid_pval>::iterator> it1;
				typedef istream_iterator<snpid_tags> it2;
				join(TagsOutputIterator2(vtags),
					it1(vsnpidpval.begin()), it1(vsnpidpval.end()), it2(ifs), it2());
			}
		}

		// expand snp list by tags
		vector<snpid_index> vtr;
		for (vector<snpid_pval>::size_type i=0; i<vsnpidpval.size(); ++i) {
			snpid_index tr;
			tr.snpid(vsnpidpval[i]);
			tr.index(i);
			vtr.push_back(tr);
			for (Tags::iterator itt = vtags[i].begin(); itt != vtags[i].end(); ++itt) {
				tr.snpid(*itt);
				vtr.push_back(tr);
			}
		}
		sort(vtr.begin(), vtr.end());

		// join enrichment annotation
		{
			ifstream ifs;
			safefileopen(ifs, fn_annot);
			vannot = vector<Annot>(vsnpidpval.size());
			{
				typedef istream_iterator<snpid_annot> it;
				join(TagRefOutputIterator(vannot), vtr.begin(), vtr.end(), it(ifs), it());
			}
		}

		t = toc(t, "clumping");
	} // --> vannot contains snps with markers


	// join maf,tssd and print
	{
		ifstream ifs;
		safefileopen(ifs, fn_maftssd);
		{
			typedef Enumerator<vector<snpid_pval>::iterator> it1;
			typedef istream_iterator<snpid_maf_tssd> it2;
			join(PrintOutputIterator(vtags, vannot, output),
				it1(vsnpidpval.begin()), it1(vsnpidpval.end()), it2(ifs), it2());
		}
	}

	toc(tstart);
}
