/*
 * Copyright 2007 Eugene Davydov
 *
 *
 * This file is part of the GERPv2.1 package.
 *
 * GERPv2.1 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * GERPv2.1 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <sstream>
#include <cstring>
#include "etree.h"
#include "emodel.h"
#include "Mseq.h"
#include "MIter.h"
#include "Seq.h"
#include "Vec.h"

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

string tree_fname = "", align_fname = "";
bool align_format_maf = true;
bool skip_prefix = true;
const string prog_name = "gpp-gerpcol";
double total_rate = -1, scale_by = -1, num_tol = 0.001;
string suffix = ".rates";
string refseq = "";
bool proj_ref = false;
double trtv = 2.0;
bool v = false, help = false;


double f(ETree &t, EModel &m, double x) {
  t.computeUp(m,x);
  return -t.computeNorm(m);
}


/* This routine, getBestK, is based on the routine(s) brent from the book
 * Numerical Recipes in C (Cambridge University Press), Copyright (C)
 * 1987-1992 by Numerical Recipes Software.  Used by permission.  Use of
 * this routine other than as an integral part of GERPv2.1 requires an
 * additional license from Numerical Recipes Software.  Further distribution
 * in any form is prohibited.
 */
double getBestK(ETree &t, EModel &m, double ax, double bx, double cx) {
  const int ITMAX  = 100;
  const double CGOLD = 0.3819660f;
  const double ZEPS = 1.0e-10f;

  int iter;

  double tol = num_tol;
  double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;
  
  a = ax;
  b = cx;
  x=w=v=bx;
  fw=fv=fx=f(t,m,x);

  double fa = f(t,m,a);
  double fb = f(t,m,b);
  if (fa < fx) return a;
  if (fb < fx) return b;

  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      return x;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >=
	  q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=f(t,m,u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      v=w;w=x;x=u;
      fv=fw;fw=fx; fx=fu;
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v=u;
	fv=fu;
      }
    }
  }
  return -1;
}

void estimateFreqs(vector<string> &species, double fq[]) {
  ifstream ifs(align_fname.c_str());
  string buffer;

  unsigned int cA = 0;
  unsigned int cC = 0;
  unsigned int cG = 0;
  unsigned int cT = 0;

  while (true) {
    getline(ifs, buffer);
    if (ifs.eof()) break;

    if (buffer.size() > 0 && buffer[0] == 's') {
      stringstream line(buffer);
      string seq;
      unsigned int tn;
      char tch;
      line >> tch >> seq >> tn >> tn >> tch >> tn;
      string s = seq.substr(0, seq.find_first_of('.'));
      if (find(species.begin(), species.end(), s) == species.end()) {
	species.push_back(s);
      }
      line >> seq;

      for (unsigned int i = 0; i < seq.size(); i++) {
	switch(seq[i]) {
	case 'a': case 'A':
	  cA++; break;
    	case 'c': case 'C':
	  cC++; break;
	case 'g': case 'G':
	  cG++; break;
	case 't': case 'T':
	  cT++; break;
	default: break;
	}
      }

    }
  }

  unsigned int cN = cA + cC + cG + cT;
  fq[0] = cT / (double)cN;
  fq[1] = cC / (double)cN;
  fq[2] = cA / (double)cN;
  fq[3] = cG / (double)cN;

  ifs.close();
}

void processMAF(ofstream &fout, ETree &src, double fq[]) {
  vector<string> species;
  if (fq[0] + fq[1] + fq[2] + fq[3] < 0.9)
    estimateFreqs(species, fq);
  EModel m(fq, trtv);
  
  if (v) cout << "Nucleotide frequencies:  A = " << m.f[2] << ", C = " << m.f[1] << ", G = " << m.f[3] << ", T = " << m.f[0] << endl;
  
  vector<string> tsp;
  for (unsigned int j = 0; j < src.foliage.size(); j++) {
    tsp.push_back(src.foliage[j]->sqname);
  }
  sort(species.begin(), species.end());
  sort(tsp.begin(), tsp.end());
  for (unsigned int i = 0, j = 0; i < species.size() || j < tsp.size(); ) {
    if (i < species.size() && (j == tsp.size() || species[i] < tsp[j])) {
      cout << "Alignment species " << species[i] << " not found in tree and therefore ignored." << endl;
      i++;
    } else if (j < tsp.size() && (i == species.size() || tsp[j] < species[i])) {
      cout << "Tree species " << tsp[j] << " not present in alignment and therefore ignored." << endl;
      j++;
    } else {
      i++; j++;
    }
  }
  
  ifstream ifs(align_fname.c_str());
  string buffer;
  
  unsigned int cpos = 0, fskip = 0;
  while (true) {
    getline(ifs, buffer);
    if (ifs.eof()) break;
    if (buffer[0] == 'a') {
      // alignment block
      map<string, string> h;
      unsigned int L = 0;
      unsigned int refstart = 0, reflen = 0;
      
      while (true) {
	getline(ifs, buffer);
	if (ifs.eof() || buffer.size() < 1) break;
	// sequence
	if (buffer[0] == 's') {
	  stringstream line(buffer);
	  string seq;
	  unsigned int tn;
	  char tch;
	  line >> tch >> seq;
	  string spc = seq.substr(0, seq.find_first_of('.'));
	  //cerr << spc << endl;
	  if (spc == refseq) {
	    line >> refstart >> reflen;
	  } else {
	    line >> tn >> tn;
	  }
	  line >> tch >> tn >> seq;
	  h[spc] = seq;
	  if (L == 0) L = seq.size();
	  else if (L != seq.size()) cerr << "alignment block mismatch" << endl;
	}
      }
      
      // end alignment block, do rate computations
      if (refstart == 0 && reflen == 0) cerr << "reference sequence not found" << endl;
      
      // avoids bunch of 0 0 rows at beginning
      if (skip_prefix && cpos == 0)
	fskip = cpos = refstart;

      // fast forward
      while (cpos < refstart) {
	fout << 0.0 << '\t' << 0.0 << endl;
	cpos++;
      }
       
      for (unsigned int i = 0; i < L; i++) {
	if (string("ACGTacgt").find(h[refseq][i]) == string::npos) continue;
	if (proj_ref) h[refseq][i] = '-';
	
	ETree tr(src);
	
	for (unsigned int j = 0; j < tr.foliage.size(); j++) {
	  if (h.find(tr.foliage[j]->sqname) != h.end()) {
	    switch(h[tr.foliage[j]->sqname][i]) {
	    case 't': case 'T': tr.foliage[j]->leafVal = 0; break;
	    case 'c': case 'C': tr.foliage[j]->leafVal = 1; break;
	    case 'a': case 'A': tr.foliage[j]->leafVal = 2; break;
	    case 'g': case 'G': tr.foliage[j]->leafVal = 3; break;
	    default: break;
	    }
	  }
	}
	
	cpos++;
	double rexp = 0, robs = 0;
	if (tr.prep()) {
	  double r = tr.getNeutralRate();
	  tr.scaleBy(1.0 / r);
	  
	  double x = getBestK(tr, m, 0.0, 3.0 * r, 3.1 * r);
	  if (x > 3.0 * r) x = 3.0 * r;
	  robs = x;
	  rexp = r;
	}
	
	fout << rexp << '\t' << (rexp - robs) << endl;
      }
    }
  }
  cout << "Processed alignment of " << cpos - fskip<< " positions." << endl;
}

void processMFA(ofstream &fout, ETree &src) {
  Mseq s(align_fname);
  EModel m(s, trtv);
  if (v) cout << "Nucleotide frequencies:  A = " << m.f[2] << ", C = " << m.f[1] << ", G = " << m.f[3] << ", T = " << m.f[0] << endl;
  
  map<string, Seq*> h;
  
  for (unsigned int i = 0; i < s.getSize(); i++) {
    h[s[i].getTitle()] = &(s[i]);
    
    if (v) {
      bool found = false;
      for (unsigned int j = 0; j < src.foliage.size(); j++) {
	if (s[i].getTitle() == src.foliage[j]->sqname) {
	  found = true;
	  break;
	}
      }
      if (!found) cout << "Alignment species " << s[i].getTitle() << " not found in tree and therefore ignored." << endl;
    }
  }
  
  if (v) {
    for (unsigned int j = 0; j < src.foliage.size(); j++) 
      if (!h[src.foliage[j]->sqname])
	cout << "Tree species " << src.foliage[j]->sqname << " not present in this alignment and therefore ignored." << endl;
  }
  
  unsigned int L = s.getLength();
  if (v) {
    cout << "Processing alignment of " << L << " positions, ";
    cout << "maximum neutral rate is " << src.getNeutralRate() << endl;
  }
  
  for (unsigned int i = 0; i < L; i++) {
    if (string("ACGTacgt").find(h[refseq]->getLetter(i+1)) == string::npos)
      continue;
    if (proj_ref) h[refseq]->changeLetter('-', i+1);

    ETree tr(src);
    for (unsigned int j = 0; j < tr.foliage.size(); j++) {
      Seq* s = h[tr.foliage[j]->sqname];
      if (s) {
	switch(s->getLetter(i+1)) {
	case 't': case 'T': tr.foliage[j]->leafVal = 0; break;
	case 'c': case 'C': tr.foliage[j]->leafVal = 1; break;
	case 'a': case 'A': tr.foliage[j]->leafVal = 2; break;
	case 'g': case 'G': tr.foliage[j]->leafVal = 3; break;
	default: break;
	}
      }
    }
    
    double rexp = 0, robs = 0;
    if (tr.prep()) {
      double r = tr.getNeutralRate();
      tr.scaleBy(1.0 / r);
      
      double x = getBestK(tr, m, 0.0, 3.0 * r, 3.1 * r);
      if (x > 3.0 * r) x = 3.0 * r;
      robs = x;
      rexp = r;
    }
    
    fout << rexp << '\t' << (rexp - robs) << endl;
  }
}

int main (int argc, char *argv[]) {
  double fq[4] = { 0.0 };
  //ifstream nuc_f;

  for ( ; argc > 0; argc--, argv++) {
    if (argv[0][0] != '-') continue;

    switch (argv[0][1]) {
    case 'h':
      cout << prog_name << " options:" << endl << endl;
      cout << " -h \t print help menu" << endl;
      cout << " -v \t verbose mode" << endl;
      cout << " -t <tree filename>" << endl;
      cout << "    \t evolutionary tree" << endl;
      cout << " -f <filename>" << endl;
      cout << "    \t alignment filename" << endl;
      cout << " -a \t alignment in mfa format [default = false]" << endl;
      cout << " -e <reference seq>" << endl;
      cout << "    \t name of reference sequence" << endl;
      cout << " -j \t project out reference sequence" << endl;
      cout << " -r <ratio>" << endl;
      cout << "    \t Tr/Tv ratio [default = 2.0]" << endl;
      //      cout << " -u <filename>" << endl;
      //cout << "    \t nucleotide frequencies [default = compute from alignment]" << endl;
      cout << " -p <precision>" << endl;
      cout << "    \t tolerance for rate estimation [default = 0.001]" << endl;
      cout << " -z force start at position 0 [default = false]" << endl;
      cout << " -n <rate>" << endl;
      cout << "    \t tree neutral rate [default = compute from tree]" << endl;
      cout << " -s <factor>" << endl;
      cout << "    \t tree scaling factor [default = 1.0]" << endl;
      cout << " -x <suffix>" << endl;
      cout << "    \t suffix for naming output files [default = \".rates\"]" << endl << endl;
      cout << "Please see README.txt for more information." << endl << endl;
      help = true;
      break;
    case 'v':
      v = true;
      break;
    case 'z':
      skip_prefix = false;
      break;
    case 'r':
      trtv = strtod(argv[1], NULL);
      break;
      //    case 'u':
      //nuc_f.open(argv[1]);
      //nuc_f >> fq[2] >> fq[1] >> fq[3] >> fq[0];
      //nuc_f.close();
      //break;
    case 's':
      scale_by = strtod(argv[1], NULL);
      break;
    case 't':
      tree_fname = argv[1];
      break;
    case 'a':
      align_format_maf = false;
      break;
    case 'f':
      align_fname = argv[1];
      break;
    case 'e':
      refseq = argv[1];
      break;
    case 'p':
      num_tol = strtod(argv[1], NULL);
      break;
    case 'j':
      proj_ref = true;
      break;
    case 'n':
      total_rate = strtod(argv[1], NULL);
      break;
    case 'x':
      suffix = argv[1];
      break;
    default:
      cerr << "Ignoring unrecognized option " << argv[0] << endl;
      cerr << "Use the -h option to see the help menu." << endl;
      break;
    }
  }

  if (align_fname.length() < 1 || tree_fname.length() < 1 || refseq.length() < 1) {
    if (help) exit(0);
    cerr << "ERROR:  one or more required arguments missing." << endl;
    cerr << "Use 'gerpcol -h' to see the help menu." << endl;
    exit(1);
  }

  if (v) {
    cout << prog_name << endl;
    system("date");
  }

  ETree ttmp(tree_fname);
  ETree src(ttmp);
  if (v) cout << "Neutral rate computed from tree file = " << src.getNeutralRate() << endl;
  if (scale_by > 0) src.scaleBy(scale_by);
  if (total_rate > 0) src.scaleBy(total_rate / src.getNeutralRate());
  if (v) cout << "Neutral rate after rescaling = " << src.getNeutralRate() << endl << endl;

  if (v) cout << "Processing " << align_fname << ", output will be written to " << align_fname + suffix << endl;
  ofstream fout;
  fout.open((align_fname + suffix).c_str(), ios::out);
  fout.precision(3);
  fout.unsetf(ios::scientific);

  if (align_format_maf) {
    processMAF(fout, src, fq);
  } else {
    processMFA(fout, src);
  }

  fout.close();
  if (v) cout << "Finished processing " << align_fname << endl << endl;

  return 0;
}
