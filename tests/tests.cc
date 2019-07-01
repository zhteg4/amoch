// ============================================================================
// tests.cc -- Unit tests
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#include <cstdio>
#include <cmath>
#include <cstring>
// Include all files to verify dependencies
#include "config.h"
#include "defs.h"
#include "ivec.h"
#include "os.h"
#include "param.h"
#include "timer.h"
#include "rng.h"
#include "triple.h"
#include "logfile.h"
#include "types.h"
#include "monomer.h"
#include "builder.h"
#include "cbmc.h"
#include "chaintype.h"
#include "cell.h"
#include "system.h"
#include "chain.h"
#include "energy.h"
#include "lj.h"
#include "fene.h"
#include "restart.h"
#include "output.h"
#include "exclusion.h"
#include "slab.h"

using std::vector;
using AMOCH::Real;
using AMOCH::Param;
using AMOCH::RNG;
using AMOCH::Triple;
using AMOCH::PartType;
using AMOCH::BondType;
using AMOCH::match_part_type;
using AMOCH::match_bond_type;
using AMOCH::Ivec;
using AMOCH::Monomer;
using AMOCH::Chaintype;
using AMOCH::System;
using AMOCH::Energy;
using AMOCH::Chain;
using AMOCH::Cell;
using AMOCH::Exclusion;

static int nfunc = 0;
static int ntest = 0;
static int nsucc = 0;
static int nfail = 0;
static Real tol = 1e-9;
static Real orig_tol;

#define TESTFILE "testfile"

#define LOWER_TOL(t)  orig_tol = tol; tol = t
#define RESET_TOL     tol = orig_tol

#define ADD_TEST   ntest++
#define SUCCESS    nsucc++
#define FAILURE    nfail++;  printf("X(%d)", __LINE__); fflush(stdout)

#define TEST_EQUAL(a,b)  \
   ADD_TEST; if ((a) == (b)) {SUCCESS;} else {FAILURE;}

#define TEST_NOT_EQUAL(a,b) \
   ADD_TEST; if ((a) != (b)) {SUCCESS;} else {FAILURE;}

#define TEST_CLOSE(a,b) \
   ADD_TEST; if (fabs((a)-(b)) < tol) {SUCCESS;} else {FAILURE;}

#define TEST_TRUE(expr) \
   ADD_TEST; if ((expr)) {SUCCESS;} else {FAILURE;}

// ============================================================================
// Return true if str == contents of path, else return false
static bool
compare_str_file(const char *str,
                 const char *path)
{
   size_t len = strlen(str);
   char *buf = new char[len+1];
   FILE *f = fopen(path, "r");
   unsigned i;

   for (i = 0; i < len; i++)
      buf[i] = fgetc(f);
   buf[i] = '\0';
   fclose(f);

   /*
   printf("string = [%s]\n", str);
   printf("  file = [%s]\n", buf);
   */

   bool retval = (0 == strcmp(str, buf));
   delete[] buf;
   return retval;
}

// ============================================================================
static void
test_defs()
{
   TEST_CLOSE(M_PI, DEG2RAD(180.0));
   TEST_CLOSE(180.0, RAD2DEG(M_PI));
   TEST_CLOSE(0.59616123, kB*300.0);  // kT at room temp, kcal/mol
}

// ============================================================================
static void
test_os()
{
#ifdef __unix
   std::string path = AMOCH::join_path("mydir", "myfile");

   TEST_EQUAL("mydir/myfile", path);
#endif
   TEST_EQUAL(false, AMOCH::check_file(TESTFILE));

   FILE *f = fopen(TESTFILE, "w");
   fclose(f);

   TEST_EQUAL(true, AMOCH::check_file(TESTFILE));

   TEST_EQUAL(true, AMOCH::check_dir("src"));
   TEST_EQUAL(false, AMOCH::check_dir("bug"));

   // XXX Not tested: AMOCH::create_dir()
}

// ============================================================================
extern bool 
parse_input(Param *proot,
            const char *path);

static void
test_parser()
{
   FILE *f = fopen(TESTFILE, "w");
   Param proot("root");
   const Param *p;

   fprintf(f, "#\n#\n# test input file { [, }]\n\n\n");
   fprintf(f, "key1 value1\n");
   fprintf(f, "key2\t\t  \n value2\n");
   fprintf(f, "key3 [3a,        3b,\n3c,3d \t]\n");
   fprintf(f, "key4 value4 {key5 value5\n");
   fprintf(f, "   key6 \nvalue6 \n{\n\n\n");
   fprintf(f, "   key7 [7a, 7b, 7c]}}\n\n");
   fprintf(f, "key2 [value8a, value8b]\n");
   fprintf(f, "key2 value9\n");
   fprintf(f, "key10 ' value10a | value10b :'\n");
   fprintf(f, "key11 \"value11a \t value11b\"\n");
   fclose(f);

   TEST_EQUAL(true, parse_input(&proot, TESTFILE));

   TEST_EQUAL("root", proot.key());

   p = proot.find("key1");
   TEST_TRUE(p != NULL); 
   TEST_EQUAL(1, p->nvalues());
   TEST_EQUAL("value1", p->value());
   TEST_EQUAL(NULL, p->find_next("key1"));

   p = proot.find("key2");
   TEST_TRUE(p != NULL);
   TEST_EQUAL(1, p->nvalues());
   TEST_EQUAL("value2", p->value());
   p = p->find_next("key2");
   TEST_TRUE(p != NULL);
   TEST_EQUAL(2, p->nvalues());
   TEST_EQUAL("value8b", p->value(1));
   p = p->find_next("key2");
   TEST_TRUE(p != NULL);
   TEST_EQUAL(1, p->nvalues());
   TEST_EQUAL("value9", p->value());
   TEST_EQUAL(NULL, p->find_next("key2"));

   p = proot.find("key3");
   TEST_EQUAL(4, p->nvalues());
   TEST_EQUAL("3b", p->value(1));
   TEST_EQUAL("3d", p->value(3));

   p = proot.find("key4");
   TEST_EQUAL("value4", p->value());
   const Param *ps = p->find("key6");
   TEST_EQUAL("value6", ps->value());
   p = ps->find("key7");
   TEST_EQUAL(3, p->nvalues());
   TEST_EQUAL("7c", p->value(2));

   p = proot.find("key10");
   TEST_EQUAL(" value10a | value10b :", p->value());
   p = proot.find("key11");
   TEST_EQUAL("value11a \t value11b", p->value());
}

// ============================================================================
static void
test_timer()
{
   AMOCH::Timer t;

   t.start();
   for (int i = 0; i < 1000000; i++)
      ;
   TEST_TRUE(t.elapsed() > 0.0);
}

// ============================================================================
static void
test_rng()
{
   RNG rng(192837465);

   // Reproducible sequence, with checkpoint
   TEST_CLOSE(0.78037583990954, rng.yield());
   TEST_CLOSE(0.50407616258599, rng.yield());
   TEST_CLOSE(0.58125772280619, rng.yield());

   FILE *f = fopen(TESTFILE, "w");

   rng.write(f);
   fclose(f);

   RNG rng2(12345);

   f = fopen(TESTFILE, "r");
   rng2.read(f);
   fclose(f);

   TEST_CLOSE(0.50521647231653, rng2.yield());
   TEST_CLOSE(0.90742635866627, rng2.yield());
   TEST_CLOSE(0.98193687410094, rng2.yield());
   TEST_CLOSE(0.21377842850052, rng2.yield());
   
   // Uniform distribution
   const int nsamples = 15000000;
   const int nbins = 100;
   int count[nbins];
   int i;
   Real target = (Real)(nsamples/nbins);

   for (i = 0; i < nbins; i++)
      count[i] = 0;
   for (i = 0; i < nsamples; i++)
      count[(int)((Real)nbins*rng.yield())]++;
   for (i = 0; i < nbins; i++) {
      if (fabs(((Real)count[i]-target)/target) > 1e-2)
         break;
   }
   TEST_EQUAL(i, nbins);

   // Normal distribution
   int nvals = 5000000;
   const Real mean = 2.463;
   const Real stddev = 1.147;
   Real sample_mean, sample_var;
   Real *vals = new Real[nvals];
   Real z;

   rng.set_normal(mean, stddev);
   sample_mean = 0.0;
   for (i = 0; i < nvals; i++)
   {
      vals[i] = rng.yield_normal();
      sample_mean += vals[i];
   }
   sample_mean /= (Real)nvals;
   LOWER_TOL(1e-3);
   TEST_CLOSE(mean, sample_mean);
   sample_var = 0.0;
   for (i = 0; i < nvals; i++) {
      z = vals[i] - sample_mean;
      sample_var += z*z;
   }
   sample_var /= (Real)(nvals-1);
   TEST_CLOSE(stddev, sqrt(sample_var));
   RESET_TOL;
   delete[] vals;
   
   // Select by weight
   vector<Real> w(5);

   w[0] = 0.1;
   w[1] = 0.06;
   w[2] = 0.44;
   w[3] = 0.11;
   w[4] = 0.29;

   int wc[5] = {0,0,0,0,0};

   for (i = 0; i < nsamples; i++)
      wc[rng.select_by_weight(w)]++;
   for (i = 0; i < 5; i++) {
      target = w[i] * (Real)nsamples;
      if (fabs(((Real)wc[i]-target)/target) > 1e-2)
         break;
   }
   TEST_EQUAL(5, i);
}

// ============================================================================
static void
test_triple()
{
   Triple a(2.3, -0.9883, 1.5e4);
   Triple b(991.3, 4.4e-8, 1.0);
   Triple c;

   c = a+b;
   TEST_CLOSE(993.6, c.x);
   TEST_CLOSE(-0.988299956, c.y);
   TEST_CLOSE(15001.0, c.z);

   c = a-b;
   TEST_CLOSE(-989.0, c.x);
   TEST_CLOSE(-0.988300044, c.y);
   TEST_CLOSE(14999.0, c.z);

   c = a * 1.3;
   TEST_CLOSE(2.99, c.x);
   TEST_CLOSE(-1.28479, c.y);
   TEST_CLOSE(1.95e4, c.z);
   
   TEST_CLOSE(380250010.590785344, c.sqlength());
   TEST_CLOSE(19500.000271559, c.length());
   c.normalize();
   TEST_CLOSE(1.0, c.sqlength());
   TEST_CLOSE(1.0, c.length());

   TEST_CLOSE(17279.9899999565148, dot(a, b));
   TEST_CLOSE(17279.9899999565148, dot(b, a));

   c = cross(a, b);
   TEST_CLOSE(-0.98896, c.x);
   TEST_CLOSE(14869497.7, c.y);
   TEST_CLOSE(979.7017901012, c.z);

   c.zero();
   c += Triple(1,1,1);
   TEST_CLOSE(1.0, c.x);
   TEST_CLOSE(1.0, c.y);
   TEST_CLOSE(1.0, c.z);
}

// ============================================================================
static void
test_logfile()
{
   AMOCH::set_logfile(stdout);
   TEST_EQUAL(stdout, AMOCH::get_logfile());
   AMOCH::open_logfile(TESTFILE);
   AMOCH::log_message("This is a test: %d\n", 42);
   AMOCH::close_logfile();
   TEST_TRUE(compare_str_file("This is a test: 42\n", TESTFILE));
}

// ============================================================================
static void
test_types()
{
   vector<PartType> part_types(0);

   part_types.push_back(PartType(6, "C", 4));        // C_3
   part_types.push_back(PartType(1, "H", 1, true));  // H__Hb
   part_types.push_back(PartType(6, "C", 3));        // C_2
   part_types.push_back(PartType(6, "C", 3, false, true));  // C_R

   TEST_EQUAL(6, part_types[0].element());
   TEST_EQUAL("C", part_types[3].symbol());
   TEST_EQUAL("H", part_types[1].symbol());
   TEST_EQUAL(3, part_types[2].nbonds());
   TEST_EQUAL("C", part_types[3].symbol());
   TEST_EQUAL(true, part_types[1].hbond());
   TEST_EQUAL(true, part_types[3].resonant());

   TEST_EQUAL(0, match_part_type(part_types, 6, 4));
   TEST_EQUAL(2, match_part_type(part_types, 6, 3));
   TEST_EQUAL(3, match_part_type(part_types, 6, 3, false, true));
   TEST_EQUAL(1, match_part_type(part_types, 1, 1, true));
   TEST_EQUAL(-1, match_part_type(part_types, 1, 1));

   vector<BondType> bond_types(0);

   bond_types.push_back(BondType(0, 1, 1));
   bond_types.push_back(BondType(1, 2, 1));
   bond_types.push_back(BondType(3, 3, 2));
   TEST_EQUAL(1, bond_types[1].ptype1());
   TEST_EQUAL(3, bond_types[2].ptype2());
   TEST_EQUAL(2, bond_types[2].order());

   TEST_EQUAL(0, match_bond_type(bond_types, 1, 0, 1));
   TEST_EQUAL(-1, match_bond_type(bond_types, 1, 0));
   TEST_EQUAL(-1, match_bond_type(bond_types, 2, 1, 2));
   TEST_EQUAL(1, match_bond_type(bond_types, 2, 1, 1));
   TEST_EQUAL(2, match_bond_type(bond_types, 3, 3, 2));
   TEST_EQUAL(-1, match_bond_type(bond_types, 3, 3));
}

// ============================================================================
static void
test_monomer()
{
   Monomer m("m1");
   int i;

   TEST_EQUAL(0, strcmp("m1", m.name().c_str()));
   TEST_EQUAL(0, m.natoms());

   TEST_EQUAL(true, m.read("tests/pe.pdb"));
   TEST_EQUAL(8, m.natoms());

   m.identify_backbone(2, 7);
   TEST_EQUAL(4, m.nbackbone());

   vector<PartType> part_types(0);
   Param pmon("monomer");

   TEST_EQUAL(true, m.setup(pmon, part_types));
   TEST_EQUAL(2, part_types.size());
   TEST_CLOSE(30.06904, m.mass());
   TEST_CLOSE(3.116199768949, m.length());  // with head, tail
   TEST_CLOSE(26.03728, m.bbmass());
   TEST_EQUAL(1, part_types[0].element());  // H_
   TEST_EQUAL(1, part_types[0].nbonds());
   TEST_EQUAL(false, part_types[0].hbond());
   TEST_EQUAL(false, part_types[0].resonant());
   TEST_EQUAL(6, part_types[1].element());  // C_3
   TEST_EQUAL(4, part_types[1].nbonds());
   TEST_EQUAL(false, part_types[1].hbond());
   TEST_EQUAL(false, part_types[1].resonant());

   m.set_torsion(2, 0.0);
   m.set_torsion(3, 129.34);

   // Atom types
   int atypes[8] = {0, 1, 1, 0, 0, 0, 0, 0};

   for (i = 0; i < 8; i++) {
      if (atypes[i] != m.type(i))
         break;
   }
   TEST_EQUAL(8, i);
   
   vector<BondType> bond_types(0);

   m.find_bonds(bond_types);

   TEST_EQUAL(4, m.max_bonds_per_atom());
   TEST_EQUAL(7, m.nbonds());
   TEST_EQUAL(2, bond_types.size());
   TEST_CLOSE(1.534215760576, m.bond_length(2));

   // zmatrix; all bond lengths, torsion angles
   const char zmstr[] =
      "  H\n"
      "  C   1  1.10\n"
      "  C   2  1.53   1   111.76\n"
      "  H   3  1.11   2   111.08   1   129.34\n"
      "  H   2  1.12   1   109.59   3  -121.26\n"
      "  H   2  1.12   1   109.56   3   121.19\n"
      "  H   3  1.11   2   111.18   1    60.07\n"
      "  H   3  1.11   2   111.15   1   -60.16\n";

   FILE *f = fopen(TESTFILE, "w");
   m.write_zmatrix(f, part_types);
   fclose(f);
   TEST_TRUE(compare_str_file(zmstr, TESTFILE));

   // Update zmatrix to support bb torsion rotations
   m.adjust_internal_coords(part_types);
   const char zmadjstr[] =
      "  H\n"
      "  C   1  1.10\n"
      "  C   2  1.53   1   111.76\n"
      "  H   3  1.11   2   111.08   1   129.34\n"
      "  H   2  1.12   1   109.59   3  -121.26\n"
      "  H   2  1.12   1   109.56   3   121.19\n"
      "  H   3  1.11   4   107.68   2  -121.96\n"
      "  H   3  1.11   4   107.68   2   121.92\n";
   f = fopen(TESTFILE, "w");
   m.write_zmatrix(f, part_types);
   fclose(f);
   TEST_TRUE(compare_str_file(zmadjstr, TESTFILE));

   TEST_EQUAL(0, match_bond_type(bond_types, 1, 1, 1));  // C_3 - C_3 
   TEST_EQUAL(1, match_bond_type(bond_types, 0, 1, 1));  // C_3 - H_
   TEST_EQUAL(1, match_bond_type(bond_types, 1, 0, 1));  // C_3 - H_

   int b1[7] = {1, 1, 1, 1, 2, 2, 2};
   int b2[7] = {2, 0, 4, 5, 6, 7, 3};
   int bt[7] = {0, 1, 1, 1, 1, 1, 1};
   Ivec b(3);

   for (i = 0; i < 7; i++) {
      b = m.bond(i);
      if ((b[0] != b1[i]) || (b[1] != b2[i]))
         break;
      if (b[2] != bt[i])
         break;
   }
   TEST_EQUAL(7, i);

   // Bonded interactions
   TEST_EQUAL(3, m.bond_separation(0, 3));
   TEST_EQUAL(2, m.bond_separation(4, 5));
   TEST_EQUAL(3, m.bond_separation(4, 7));
   TEST_EQUAL(3, m.bond_separation(0, 6));
   TEST_EQUAL(1, m.bond_separation(1, 2));
   TEST_EQUAL(1, m.bond_separation(2, 7));
   TEST_EQUAL(3, m.bond_separation(3, 5));
   TEST_EQUAL(0, m.bond_separation(2, 2));
   
   // Rotate torsions, update positions
   m.set_anchor(Triple(5.518, 4.320, 4.106));
   m.position[0].x = 6.281;
   m.position[0].y = 3.823;
   m.position[0].z = 5.336;
   m.position[1].x = 6.966;
   m.position[1].y = 2.484;
   m.position[1].z = 5.055;
   m.set_torsion(2, 90.0);
   TEST_CLOSE(90.0, m.torsion_angle(2));
   m.set_torsion(3, 180.0);
   TEST_CLOSE(180.0, m.torsion_angle(3));

   m.update_positions();
   const char *xyzpath = TESTFILE".xyz";
   m.write(xyzpath);
   const char xyzstr[] = 
      "8\n\n"
      "H          6.28100        3.82300        5.33600\n"
      "C          6.96600        2.48400        5.05500\n"
      "C          8.38586        2.67561        4.50626\n"
      "H          8.87372        1.69433        4.30582\n"
      "H          6.36256        1.90249        4.31890\n"
      "H          7.01335        1.88391        5.99491\n"
      "H          8.37431        3.24939        3.55075\n"
      "H          9.02559        3.22985        5.23099\n";
   TEST_TRUE(compare_str_file(xyzstr, xyzpath));

   Triple pos[6] = {
      Triple(8.381, 2.674, 4.507),
      Triple(8.867, 1.696, 4.308),
      Triple(6.359, 1.900, 4.317),
      Triple(7.013, 1.882, 5.999),
      Triple(8.369, 3.244, 3.555),
      Triple(9.019, 3.227, 5.228)
   };

   LOWER_TOL(1e-2);
   for (i = 0; i < 6; i++) {
      TEST_CLOSE(m.position[i+2].x, pos[i].x);
      TEST_CLOSE(m.position[i+2].y, pos[i].y);
      TEST_CLOSE(m.position[i+2].z, pos[i].z);
   }
   RESET_TOL;

   // specify_backbone() -- same monomer, same backbone
   // set torsion 3 in input file, rather than set_torsion()
   Monomer m2("m2");
   vector<int> bb(4, -1);

   TEST_EQUAL(true, m2.read("tests/pe.pdb"));
   bb[0] = 3;
   bb[1] = 1;
   bb[2] = 2;
   bb[3] = 8;
   TEST_EQUAL(0, m2.nbackbone());
   m2.specify_backbone(bb);
   TEST_EQUAL(4, m2.nbackbone());
   part_types.clear();

   f = fopen(TESTFILE, "w");
   fprintf(f, "set torsions {index 4\nangle 129.34}");
   fclose(f);
   (void)parse_input(&pmon, TESTFILE);
   TEST_EQUAL(true, m2.setup(pmon, part_types));

   f = fopen(TESTFILE, "w");
   m2.write_zmatrix(f, part_types);
   fclose(f);
   TEST_TRUE(compare_str_file(zmstr, TESTFILE));
}

// ============================================================================
static void
test_chaintype()
{
   Chaintype ct("chtype1");
   FILE *f = fopen(TESTFILE, "w");

   fprintf(f, "monomers [m1, m2]\n");
   fprintf(f, "monomer_weights [0.63, 0.37]\n");
   fprintf(f, "term m3\n");
   fclose(f);

   Param pchtype("chtype1");

   (void)parse_input(&pchtype, TESTFILE);

   vector<Monomer> monomers(0);

   monomers.push_back(Monomer("m1"));
   monomers.push_back(Monomer("m3"));
   monomers.push_back(Monomer("m2"));

   TEST_EQUAL(true, ct.setup(pchtype, monomers));
   TEST_EQUAL("chtype1", ct.name());
   TEST_EQUAL(2, ct.nmonomers());
   TEST_EQUAL(0, ct.monomer_index(0));  // m1 is index 0 in monomers[]
   TEST_EQUAL(2, ct.monomer_index(1));  // m2 is index 2 in monomers[]
   TEST_EQUAL(1, ct.term_monomer_index());
   TEST_EQUAL(true, ct.select_by_weight());
   TEST_CLOSE(0.63, ct.monomer_weight(0));
   TEST_CLOSE(0.37, ct.monomer_weight(1));

   f = fopen(TESTFILE, "w");
   fprintf(f, "monomers [m2, m1, m3]\n");
   fclose(f);

   Chaintype ct2("chtype2");
   Param pchtype2("chtype2");

   (void)parse_input(&pchtype2, TESTFILE);
   TEST_EQUAL(true, ct2.setup(pchtype2, monomers));
   TEST_EQUAL(3, ct2.nmonomers());
   TEST_EQUAL(2, ct2.monomer_index(0));  // m2 is index 2 in monomers[]
   TEST_EQUAL(0, ct2.monomer_index(1));  // m1 is index 0 in monomers[]
   TEST_EQUAL(1, ct2.monomer_index(2));  // m3 is index 1 in monomers[]
   TEST_EQUAL(-1, ct2.term_monomer_index());
   TEST_EQUAL(false, ct2.select_by_weight());

   // XXX Not tested: ct.monomer_index_by_weight() 
}

// ============================================================================
static void
test_cell()
{
   Cell c(2,2);
   Triple p, tp;

   TEST_EQUAL(0, c.nparticles());
   p.x = 7.7;
   tp.x = -2.3;

   TEST_EQUAL(0, c.add_particle(2, 1, 6, p, tp));
   p.x = 5.4;
   tp.x = -4.6;
   TEST_EQUAL(1, c.add_particle(4, 3, 5, p, tp));
   TEST_EQUAL(2, c.nparticles());
   TEST_EQUAL(2, c.chain(0));
   TEST_EQUAL(4, c.chain(1));
   TEST_EQUAL(1, c.index(0));
   TEST_EQUAL(3, c.index(1));
   TEST_EQUAL(6, c.ptype(0));
   TEST_EQUAL(5, c.ptype(1));
   TEST_CLOSE(7.7, c.pos[0].x);
   TEST_CLOSE(-2.3, c.tpos[0].x);
   TEST_CLOSE(5.4, c.pos[1].x);
   TEST_CLOSE(-4.6, c.tpos[1].x);

   TEST_EQUAL(0, c.nbonds(0));
   TEST_EQUAL(0, c.nbonds(1));

   c.add_bond(0, 1, 2);
   c.add_bond(1, 0, 2);

   TEST_EQUAL(1, c.nbonds(0));
   TEST_EQUAL(1, c.nbonds(1));

   vector<Ivec> bv;

   bv = c.bonds(0);
   TEST_EQUAL(2, bv.size());
   TEST_EQUAL(1, bv[0][0]);
   TEST_EQUAL(2, bv[0][1]);
   bv = c.bonds(1);
   TEST_EQUAL(2, bv.size());
   TEST_EQUAL(0, bv[0][0]);
   TEST_EQUAL(2, bv[0][1]);

   TEST_EQUAL(true, c.bonded(0, 1));
   TEST_EQUAL(true, c.bonded(1, 0));
   TEST_EQUAL(false, c.bonded(0, 2));
   TEST_EQUAL(false, c.bonded(1, 2));

   TEST_EQUAL(0, Cell::nresize);
   TEST_EQUAL(2, c.add_particle(7, 2, 1, p, tp));
   TEST_EQUAL(1, Cell::nresize);
   c.add_bond(1, 3, 3);
   c.add_bond(1, 4, 4);
   TEST_EQUAL(2, Cell::nresize);
}

// ============================================================================
static void
test_chain()
{
   Chain c(2, 3, 100);

   TEST_EQUAL(2, c.type());
   TEST_EQUAL(0, c.nparticles());
   TEST_EQUAL(0, c.nmonomers());
   TEST_EQUAL(3, c.max_monomers());

   c.add_monomer(6);
   c.add_monomer(5);

   c.record_torsion(37.23);
   c.record_torsion(181.9);
   c.record_torsion(224.1);

   c.add_particle_cell(99, 17, 31);

   c.tail_index = 17;

   c.record_particle(0, Triple(0,0,0));
   c.record_particle(2, Triple(1,2,1));

   FILE *f = fopen(TESTFILE, "w");

   c.write(f);
   fclose(f);

   Chain c2(1, 3, 100);
   c2.tail_index = 4;

   f = fopen(TESTFILE, "r");
   c2.read(f);
   fclose(f);

   TEST_EQUAL(2, c2.type());
   TEST_EQUAL(17, c2.tail_index);

   c2.add_monomer(6);
   c2.record_torsion(181.2);
   c2.record_torsion(331.1);

   TEST_EQUAL(17, c2.nparticles());
   TEST_EQUAL(3, c2.nmonomers());
   TEST_EQUAL(1, c2.torsion_count(37));
   TEST_EQUAL(2, c2.torsion_count(181));
   TEST_EQUAL(1, c2.torsion_count(224));
   TEST_EQUAL(1, c2.torsion_count(331));
   TEST_EQUAL(17, c2.cell(99));
   TEST_EQUAL(31, c2.slot(99));

   c2.record_particle(0, Triple(0,0,0));
   c2.record_particle(2, Triple(1,2,1));
   c2.record_particle(6, Triple(3,4,2));
   c2.record_particle(11, Triple(5,8,4));

   Triple p(4,6,3);
   c2.record_particle(15, p);
   c2.finalize(p);

   TEST_CLOSE(5.385164807134504, c2.length(0));
   TEST_CLOSE(10.2469507659596,  c2.length(1));
   TEST_CLOSE(7.810249675906654, c2.length(2));

   TEST_CLOSE(1.0/6.0, c2.com(0).x);
   TEST_CLOSE(1.0/3.0, c2.com(0).y);
   TEST_CLOSE(1.0/6.0, c2.com(0).z);
   TEST_CLOSE(0.6, c2.com(1).x);
   TEST_CLOSE(0.8, c2.com(1).y);
   TEST_CLOSE(0.4, c2.com(1).z);
   TEST_CLOSE(1.5, c2.com(2).x);
   TEST_CLOSE(7.0/3.0, c2.com(2).y);
   TEST_CLOSE(7.0/6.0, c2.com(2).z);
}

// ============================================================================
static void
test_system()
{
   Triple min(1.0, 1.0, 1.0);
   Triple max(7.0, 7.0, 7.0);
   Real cell_size = 2.0;
   System s(min, max, cell_size, 128, 4);  // 3 x 3 x 3 Cells

   TEST_CLOSE(1.0, s.min().x);
   TEST_CLOSE(7.0, s.max().y);
   TEST_CLOSE(6.0, s.size().y);
   TEST_EQUAL(27, s.ncells());
   TEST_EQUAL(0, s.cells[26].nparticles());

   // cell_index()
   Triple p(min);

   TEST_EQUAL(0, s.cell_index(p));

   p.x = 5.1;
   p.y = 3.1;
   TEST_EQUAL(5, s.cell_index(p));

   p.x = 3.5;
   p.y = 2.2;
   p.z = 4.1;
   TEST_EQUAL(10, s.cell_index(p));

   p.x = 2.5;
   p.y = 6.2;
   TEST_EQUAL(15, s.cell_index(p));

   p.x = 6.9;
   p.y = 4.0;
   p.z = 5.8;
   TEST_EQUAL(23, s.cell_index(p));

   p.y = 6.0;
   TEST_EQUAL(26, s.cell_index(p));

   // find_neighbor_cells()
   vector<int> nbrs(26, -1);
   int i;

   s.find_neighbor_cells(0, nbrs);
   int nbrs0[26] = {      1,  4,  3,  5,  2,  8,  6,  7,
                     18, 19, 22, 21, 23, 20, 26, 24, 25,
                      9, 10, 13, 12, 14, 11, 17, 15, 16};
   for (i = 0; i < 26; i++) {
      if (nbrs[i] != nbrs0[i])
         break;
   }
   TEST_EQUAL(26, i);

   s.find_neighbor_cells(13, nbrs);
   int nbrs13[26] = {    14, 17, 16, 15, 12,  9, 10, 11,
                      4,  5,  8,  7,  6,  3,  0,  1,  2,
                     22, 23, 26, 25, 24, 21, 18, 19, 20};
   for (i = 0; i < 26; i++) {
      if (nbrs[i] != nbrs13[i])
         break;
   }
   TEST_EQUAL(26, i);

   s.find_neighbor_cells(26, nbrs);
   int nbrs26[26] = {    24, 18, 20, 19, 25, 22, 23, 21,
                     17, 15,  9, 11, 10, 16, 13, 14, 12,
                      8,  6,  0,  2,  1,  7,  4,  5,  3};
   for (i = 0; i < 26; i++) {
      if (nbrs[i] != nbrs26[i])
         break;
   }
   TEST_EQUAL(26, i);

   // single bin -- no neighbors!
   System s1(Triple(0,0,0), Triple(1,1,1), 1.0, 16, 2);

   s1.find_neighbor_cells(0, nbrs);
   for (i = 0; i < 26; i++) {
      if (nbrs[i] != -1)
         break;
   }
   TEST_EQUAL(26, i);

   // wrap()
   p.x = 0.9;
   p.y = -8.3;
   p.z = 7.0;
   s.wrap(p);
   TEST_CLOSE(6.9, p.x);
   TEST_CLOSE(3.7, p.y);
   TEST_CLOSE(1.0, p.z);

   // min_sqdist()
   Triple p1(1.5,1.5,1.5);
   Triple p2(6.5,6.5,6.5);

   TEST_CLOSE(3.0, s.min_sqdist(p1, p2));
   TEST_CLOSE(3.0, s.min_sqdist(p2, p1));

   // bond_separation():  0 - 1 - 2 - 3 - 4
   Chain c(0, 2, 16);

   c.add_particle_cell(0, 0, 0);  // particle 0, cell 0, slot 0
   c.add_particle_cell(1, 0, 1);  // particle 1, cell 0, slot 1
   s.cells[0].add_bond(0, 1, 0);  // particle in slot 0 bonded to particle 1
   s.cells[0].add_bond(1, 0, 0);  // particle in slot 1 bonded to particle 0

   c.add_particle_cell(2, 1, 0);  // particle 2, cell 1, slot 0
   s.cells[0].add_bond(1, 2, 0);  // particle in slot 1 bonded to particle 2
   s.cells[1].add_bond(0, 1, 0);  // particle in slot 0 bonded to particle 1

   c.add_particle_cell(3, 2, 0);  // particle 3, cell 2, slot 0
   s.cells[2].add_bond(0, 2, 0);  // particle in slot 0 bonded to particle 2
   s.cells[1].add_bond(0, 3, 0);  // particle in slot 1 bonded to particle 3
   c.add_particle_cell(4, 2, 1);  // particle 4, cell 2, slot 1
   s.cells[2].add_bond(1, 3, 0);  // particle in slot 1 bonded to particle 3
   s.cells[2].add_bond(0, 4, 0);  // particle in slot 0 bonded to particle 0

   TEST_EQUAL(true, s.bond_separation(c, 0, 4, 4));
   TEST_EQUAL(true, s.bond_separation(c, 4, 0, 4));
   TEST_EQUAL(false, s.bond_separation(c, 4, 0, 3));
   TEST_EQUAL(true, s.bond_separation(c, 1, 3, 2));
   TEST_EQUAL(true, s.bond_separation(c, 3, 1, 2));
   TEST_EQUAL(false, s.bond_separation(c, 3, 1, 1));
   TEST_EQUAL(true, s.bond_separation(c, 3, 2, 1));
   TEST_EQUAL(true, s.bond_separation(c, 2, 3, 1));
   TEST_EQUAL(true, s.bond_separation(c, 4, 2, 4));
   TEST_EQUAL(true, s.bond_separation(c, 2, 4, 4));
   TEST_EQUAL(false, s.bond_separation(c, 0, 2, 1));
   TEST_EQUAL(false, s.bond_separation(c, 2, 0, 1));
}

// ============================================================================
static void
test_lj()
{
   FILE *f = fopen(TESTFILE, "w");

   fprintf(f, "params C{epsilon 0.08559\nsigma 3.125691}\n");
   fprintf(f, "params H{epsilon 0.01368\nsigma 2.561778}\n");
   fclose(f);

   Energy *energy = new AMOCH::LJ(2, 1.2);
   vector<PartType> part_types(0);

   part_types.push_back(PartType(1, "H"));
   part_types.push_back(PartType(6, "C"));

   Param pen("energy");

   TEST_EQUAL(true, parse_input(&pen, TESTFILE));
   TEST_EQUAL(true, energy->setup(pen, part_types));

   energy->scale();

   TEST_CLOSE(1.44, energy->pair_range_sq());
   // C-H at 1.8 Ang
   TEST_CLOSE(31.1867196147, energy->hard_energy(0, 1, 3.24));
   TEST_CLOSE(31.1867196147, energy->hard_energy(1, 0, 3.24));
   TEST_CLOSE(29.1206634449, energy->pair_energy(1, 0, 3.24));
   TEST_CLOSE(29.1206634449, energy->pair_energy(0, 1, 3.24));
   TEST_CLOSE(1.0, energy->bond_energy(0, 1, 3.24)+1.0);
   TEST_CLOSE(1.0, energy->bond_energy(1, 0, 3.24)+1.0);

   delete energy;
}

// ============================================================================
static void
test_fene()
{
   FILE *f = fopen(TESTFILE, "w");

   fprintf(f, "params C{epsilon 1\nsigma 1\nk 30\nR0 4}\n");
   fclose(f);

   Energy *energy = new AMOCH::FENE(1, 1.2);
   vector<PartType> part_types(0);

   part_types.push_back(PartType(6, "C"));

   Param pen("energy");

   TEST_EQUAL(true, parse_input(&pen, TESTFILE));
   TEST_EQUAL(true, energy->setup(pen, part_types));

   energy->scale();

   TEST_CLOSE(1.44, energy->pair_range_sq());
   // C beads at 1.1 Ang
   TEST_CLOSE(0.0166275506, energy->pair_energy(0, 0, 1.21));
   TEST_CLOSE(1.2745232708, energy->hard_energy(0, 0, 1.21));
   TEST_CLOSE(18.8729869241, energy->bond_energy(0, 0, 1.21));
}

// ============================================================================
static void
test_restart()
{
   const Real a = 8.72838;
   System s(Triple(0,0,0), Triple(a,a,a), a, 64, 4);
   vector<Chain> chains(2, Chain(0, 2, 17));
   int i, j, n;
   FILE *f;

   chains[0].add_monomer(17);
   chains[1].add_monomer(17);

   f = fopen("tests/pmma.dump", "r");
   TEST_EQUAL(true, AMOCH::read_lammps_dump(f, s, chains));
   fclose(f);
   TEST_EQUAL(34, s.cells[0].nparticles());
   TEST_CLOSE(8.76167, s.cells[0].tpos[0].z);
   TEST_CLOSE(0.03329, s.cells[0].pos[0].z);
   TEST_CLOSE(s.cells[0].tpos[33].x, s.cells[0].pos[33].x);
   TEST_EQUAL(0, s.cells[0].chain(16));
   TEST_EQUAL(16, s.cells[0].index(16));
   TEST_EQUAL(1, s.cells[0].chain(17));
   TEST_EQUAL(0, s.cells[0].index(17));
   TEST_EQUAL(3, s.cells[0].ptype(11));
   TEST_EQUAL(4, s.cells[0].ptype(29));
   for (i = 0; i < 2; i++) {
      for (j = 0; j < 17; j++) {
         n = 17*i + j;
         if (chains[i].cell(j) != 0)
            break;
         if (chains[i].slot(j) != n)
            break;
      }
      TEST_EQUAL(17, j);
   }

   f = fopen("tests/bonds.dat", "r");
   TEST_EQUAL(true, AMOCH::read_bonds(f, s, chains));
   fclose(f);
   TEST_EQUAL(1, s.cells[0].nbonds(33));

   vector<Ivec> bonds(0, Ivec(0));

   bonds = s.cells[0].bonds(33);
   TEST_EQUAL(13, bonds[0][0]);
   TEST_EQUAL(0, bonds[0][1]);

   TEST_EQUAL(4, s.cells[0].nbonds(30));
   bonds = s.cells[0].bonds(30);
   TEST_EQUAL(12, bonds[0][0]);
   TEST_EQUAL(5,  bonds[0][1]);
   TEST_EQUAL(14, bonds[1][0]);
   TEST_EQUAL(0,  bonds[1][1]);
   TEST_EQUAL(15, bonds[2][0]);
   TEST_EQUAL(0,  bonds[2][1]);
   TEST_EQUAL(16, bonds[3][0]);
   TEST_EQUAL(0,  bonds[3][1]);
}

// ============================================================================
static void
test_cbmc()
{
   AMOCH::Builder *builder = new AMOCH::CBMC;

   Monomer m("mon1");
   Param pbuild("builder");
   Param pmon("monomer");
   FILE *f = fopen(TESTFILE, "w");
   fprintf(f, "build 1\n");
   fprintf(f, "head 3\ntail 8\n");
   fprintf(f, "configs 2\n");
   fclose(f);
   (void)parse_input(&pbuild, TESTFILE);
   (void)m.read("tests/pe.pdb");
   builder->setup_monomer(m, pbuild);
   TEST_EQUAL(4, m.nbackbone());

   vector<PartType> part_types(0);
   vector<BondType> bond_types(0);

   (void)m.setup(pmon, part_types);
   m.find_bonds(bond_types);
   m.adjust_internal_coords(part_types);
   builder->update_types(m, m, part_types, bond_types);
   vector<Monomer> monomers(1, m);

   TEST_EQUAL(2, part_types.size());
   TEST_EQUAL(2, bond_types.size());
   TEST_EQUAL(true, builder->setup(pbuild, part_types, bond_types));
   TEST_EQUAL(false, builder->backbone_only());

   System sys(Triple(0,0,0), Triple(10,10,10), 5, 32, 4);

   Chain c(0, 2, 16);
   vector<Chain> chains(1, c);

   Chaintype ct("chtype");
   Param pcht("chtype2");
   f = fopen(TESTFILE, "w");
   fprintf(f, "monomers mon1");
   fclose(f);
   (void)parse_input(&pcht, TESTFILE);
   (void)ct.setup(pcht, monomers);
   vector<Chaintype> chtypes(1, ct);

   f = fopen(TESTFILE, "w");
   fprintf(f, "params C{epsilon 0.08559\nsigma 3.125691}\n");
   fprintf(f, "params H{epsilon 0.01368\nsigma 2.561778}\n");
   fclose(f);
   Energy *energy = new AMOCH::LJ(2, 1.2);
   Param pen("energy");
   (void)parse_input(&pen, TESTFILE);
   (void)energy->setup(pen, part_types);
   energy->scale();

   vector<const Exclusion *> exclusions(0, NULL);
   RNG rng(12345);
   vector<int> monselect(1, 0);

   builder->build(sys, chains, monomers, monselect, chtypes, exclusions, 
                  *energy, rng);
   TEST_EQUAL(1, chains[0].nmonomers());
   TEST_EQUAL(7, chains[0].nparticles());
   TEST_EQUAL(1, monselect[0]);
   int n = 0;
   for (int i = 0; i < sys.ncells(); i++)
      n += sys.cells[i].nparticles();
   TEST_EQUAL(7, n);

   f = fopen(TESTFILE, "w");
   builder->write(f);
   fclose(f);
   delete builder;

   builder = new AMOCH::CBMC;
   TEST_EQUAL(true, builder->setup(pbuild, part_types, bond_types));
   builder->build(sys, chains, monomers, monselect, chtypes, exclusions, 
                  *energy, rng);
   TEST_EQUAL(2, chains[0].nmonomers());
   TEST_EQUAL(14, chains[0].nparticles());
   TEST_EQUAL(2, monselect[0]);
   n = 0;
   for (int i = 0; i < sys.ncells(); i++)
      n += sys.cells[i].nparticles();
   TEST_EQUAL(14, n);

   delete builder;
   delete energy;
}

// ============================================================================
static void
test_output()
{
   System sys(Triple(0,0,0), Triple(6,6,6), 6.0, 8, 4);
   Chain c(0, 1, 8);
   vector<Chain> chains(1, c);
   int i;
   int ptypes[8] = {0, 1, 1, 0, 0, 0, 0, 0};
   Triple pos[8] = {
      Triple(2.167, 0.027, 0.725),
      Triple(2.976, 0.719, 0.983),
      Triple(4.244, -0.026, 1.420),
      Triple(5.059, 0.688, 1.679),
      Triple(2.643, 1.393, 1.807),
      Triple(3.213, 1.362, 0.102),
      Triple(4.621, -0.690, 0.608),
      Triple(4.050, -0.658, 2.317)
   };
   Triple p;
   vector<PartType> part_types(0);

   part_types.push_back(PartType(1, "H", 1));
   part_types.push_back(PartType(6, "C", 4));

   for (i = 0; i < 8; i++) {
      p = pos[i];
      sys.wrap(p);
      sys.cells[0].add_particle(0, i, ptypes[i], p, pos[i]);
      chains[0].add_particle_cell(i, 0, i);
   }
   sys.cells[0].add_bond(0, 1, 0);
   sys.cells[0].add_bond(1, 0, 0);
   sys.cells[0].add_bond(1, 2, 1);
   sys.cells[0].add_bond(2, 1, 1);
   sys.cells[0].add_bond(1, 4, 0);
   sys.cells[0].add_bond(4, 1, 0);
   sys.cells[0].add_bond(1, 5, 0);
   sys.cells[0].add_bond(5, 1, 0);
   sys.cells[0].add_bond(2, 3, 0);
   sys.cells[0].add_bond(3, 2, 0);
   sys.cells[0].add_bond(2, 6, 0);
   sys.cells[0].add_bond(6, 2, 0);
   sys.cells[0].add_bond(2, 7, 0);
   sys.cells[0].add_bond(7, 2, 0);

   chains[0].add_monomer(8);

   const char wxyzstr[] = 
      "8\n"
      "PE out test wrapped\n"
      "   H          2.167          0.027          0.725\n"
      "   C          2.976          0.719          0.983\n"
      "   C          4.244          5.974          1.420\n"
      "   H          5.059          0.688          1.679\n"
      "   H          2.643          1.393          1.807\n"
      "   H          3.213          1.362          0.102\n"
      "   H          4.621          5.310          0.608\n"
      "   H          4.050          5.342          2.317\n";
   FILE *f = fopen(TESTFILE, "w");
   AMOCH::write_xyz(f, "PE out test wrapped", sys, chains, part_types, true);
   fclose(f);
   TEST_EQUAL(true, compare_str_file(wxyzstr, TESTFILE));

   const char uxyzstr[] = 
      "8\n"
      "PE out test unwrapped\n"
      "   H          2.167          0.027          0.725\n"
      "   C          2.976          0.719          0.983\n"
      "   C          4.244         -0.026          1.420\n"
      "   H          5.059          0.688          1.679\n"
      "   H          2.643          1.393          1.807\n"
      "   H          3.213          1.362          0.102\n"
      "   H          4.621         -0.690          0.608\n"
      "   H          4.050         -0.658          2.317\n";
   f = fopen(TESTFILE, "w");
   AMOCH::write_xyz(f, "PE out test unwrapped", sys, chains, part_types, false);
   fclose(f);
   TEST_EQUAL(true, compare_str_file(uxyzstr, TESTFILE));
      
   const char wpdbstr[] = 
"REMARK   1\n"
"REMARK   1 PDB out test wrapped\n"
"ATOM      1    H                 2.167   0.027   0.725                       H  \n"
"ATOM      2    C                 2.976   0.719   0.983                       C  \n"
"ATOM      3    C                 4.244   5.974   1.420                       C  \n"
"ATOM      4    H                 5.059   0.688   1.679                       H  \n"
"ATOM      5    H                 2.643   1.393   1.807                       H  \n"
"ATOM      6    H                 3.213   1.362   0.102                       H  \n"
"ATOM      7    H                 4.621   5.310   0.608                       H  \n"
"ATOM      8    H                 4.050   5.342   2.317                       H  \n"
"TER       9\n"
"CONECT    1    2\n"
"CONECT    2    1    3    5    6\n"
"CONECT    3    2    4    7    8\n"
"CONECT    4    3\n"
"CONECT    5    2\n"
"CONECT    6    2\n"
"CONECT    7    3\n"
"CONECT    8    3\n"
"END\n";
   f = fopen(TESTFILE, "w");
   AMOCH::write_pdb(f, "PDB out test wrapped", sys, chains, part_types, true);
   fclose(f);
   TEST_EQUAL(true, compare_str_file(wpdbstr, TESTFILE));

   const char updbstr[] = 
"REMARK   1\n"
"REMARK   1 PDB out test unwrapped\n"
"ATOM      1    H                 2.167   0.027   0.725                       H  \n"
"ATOM      2    C                 2.976   0.719   0.983                       C  \n"
"ATOM      3    C                 4.244  -0.026   1.420                       C  \n"
"ATOM      4    H                 5.059   0.688   1.679                       H  \n"
"ATOM      5    H                 2.643   1.393   1.807                       H  \n"
"ATOM      6    H                 3.213   1.362   0.102                       H  \n"
"ATOM      7    H                 4.621  -0.690   0.608                       H  \n"
"ATOM      8    H                 4.050  -0.658   2.317                       H  \n"
"TER       9\n"
"CONECT    1    2\n"
"CONECT    2    1    3    5    6\n"
"CONECT    3    2    4    7    8\n"
"CONECT    4    3\n"
"CONECT    5    2\n"
"CONECT    6    2\n"
"CONECT    7    3\n"
"CONECT    8    3\n"
"END\n";
   f = fopen(TESTFILE, "w");
   AMOCH::write_pdb(f, "PDB out test unwrapped", sys, chains, part_types,false);
   fclose(f);
   TEST_EQUAL(true, compare_str_file(updbstr, TESTFILE));
      
   const char astr[] = 
      "8 atoms\n\nATOMS\n\n"
      "1  1  1  0.000000  2.167000000000  0.027000000000  0.725000000000\n"
      "2  1  2  0.000000  2.976000000000  0.719000000000  0.983000000000\n"
      "3  1  2  0.000000  4.244000000000  5.974000000000  1.420000000000\n"
      "4  1  1  0.000000  5.059000000000  0.688000000000  1.679000000000\n"
      "5  1  1  0.000000  2.643000000000  1.393000000000  1.807000000000\n"
      "6  1  1  0.000000  3.213000000000  1.362000000000  0.102000000000\n"
      "7  1  1  0.000000  4.621000000000  5.310000000000  0.608000000000\n"
      "8  1  1  0.000000  4.050000000000  5.342000000000  2.317000000000\n";
   f = fopen(TESTFILE, "w");
   AMOCH::write_atoms_dat(f, sys, chains);
   fclose(f);
   TEST_EQUAL(true, compare_str_file(astr, TESTFILE));

   const char atstr[] = 
      "1  H  1  0  0\n"
      "2  C  4  0  0\n";
   f = fopen(TESTFILE, "w");
   AMOCH::write_atom_type_dat_in(f, part_types);
   fclose(f);
   TEST_EQUAL(true, compare_str_file(atstr, TESTFILE));

   const char bstr[] = 
      "BONDS\n\n"
      "1  1  1  2\n"
      "2  2  2  3\n"
      "3  1  2  5\n"
      "4  1  2  6\n"
      "5  1  3  4\n"
      "6  1  3  7\n"
      "7  1  3  8\n";
   f = fopen(TESTFILE, "w");
   AMOCH::write_bonds_dat(f, sys, chains);
   fclose(f);
   TEST_EQUAL(true, compare_str_file(bstr, TESTFILE));

   vector<BondType> bond_types(0);
   bond_types.push_back(BondType(0, 1, 1));
   bond_types.push_back(BondType(1, 1, 1));

   const char btstr[] = 
      "1  1  2  1\n"
      "2  2  2  1\n";
   f = fopen(TESTFILE, "w");
   AMOCH::write_bond_type_dat_in(f, bond_types);
   fclose(f);
   TEST_EQUAL(true, compare_str_file(btstr, TESTFILE));

   // XXX not tested: AMOCH::show_status()
}

// ============================================================================
static void
test_slab()
{
   FILE *f = fopen(TESTFILE, "w");
   fprintf(f, "min [0,0,0]\nmax[5,5,5]\n");
   fclose(f);

   Param pbuild("exclude");
   (void)parse_input(&pbuild, TESTFILE);

   Exclusion *ex = new AMOCH::Slab;

   TEST_EQUAL(true, ex->setup(pbuild));
   TEST_CLOSE(125.0, ex->volume());
   TEST_EQUAL(false, ex->inside(Triple(0,0,0)));
   TEST_EQUAL(false, ex->inside(Triple(5,5,5)));
   TEST_EQUAL(true, ex->inside(Triple(1,2,4)));
   TEST_EQUAL(false, ex->inside(Triple(1, -2, 4)));
   TEST_EQUAL(false, ex->inside(Triple(1,6,4)));
   TEST_CLOSE(3.0, ex->min_sqdist(Triple(-1,-1,-1)));
   TEST_CLOSE(3.0, ex->min_sqdist(Triple(-1,6,-1)));
   TEST_CLOSE(1.0, ex->min_sqdist(Triple(2,3,6)));
   TEST_CLOSE(2.0, ex->min_sqdist(Triple(6,1,6)));
   delete ex;
}

// ============================================================================
// Unit tests
// ============================================================================
int
main()
{
   int prev_succ, prev_cases;

#define TEST_FUNCTION(func) \
   printf("\n%18.18s: ", #func); fflush(stdout); nfunc++;  \
   prev_succ = nsucc; prev_cases = ntest; func(); \
   if ((nsucc-prev_succ) == (ntest-prev_cases)) printf("\u2713")

   printf("\n==================\n");
   printf("Running unit tests\n");
   printf("==================\n");

   TEST_FUNCTION(test_defs);
   TEST_FUNCTION(test_os);
   TEST_FUNCTION(test_parser);
   TEST_FUNCTION(test_timer);
   TEST_FUNCTION(test_rng);
   TEST_FUNCTION(test_triple);
   TEST_FUNCTION(test_logfile);
   //AMOCH::set_logfile(stdout);
   TEST_FUNCTION(test_types);
   TEST_FUNCTION(test_monomer);
   TEST_FUNCTION(test_chaintype);
   TEST_FUNCTION(test_cell);
   TEST_FUNCTION(test_chain);
   TEST_FUNCTION(test_system);
   TEST_FUNCTION(test_lj);     // Energy
   TEST_FUNCTION(test_fene);   // Energy
   TEST_FUNCTION(test_restart);
   TEST_FUNCTION(test_cbmc);   // Builder
   TEST_FUNCTION(test_output);
   TEST_FUNCTION(test_slab);
   // TODO Sphere, Cylinder

   // Nothing to test: config.h, ivec.h, version.h
   // Not tested: main.cc

   printf("\n\n%d test functions with %d test cases\n", nfunc, ntest);
   printf("   --> %d passed\n", nsucc);
   printf("   --> %d failed\n\n", nfail);
   return 0;
}

#if 0
// ============================================================================
static void
test_exclusion()
{
   // Sphere
   FILE *f = fopen(TESTFILE, "w");

   fprintf(f, "center [1,1,1]\nradius 2\n");
   fclose(f);

   Parser *parser = new Parser;

   TEST_EQUAL(true, parser->parse(TESTFILE));

   Exclusion *ex = new AMOCH::Sphere;

   TEST_EQUAL(true, ex->setup(*parser));
   TEST_CLOSE(33.51032163829112, ex->volume());

   TEST_EQUAL(true, ex->inside(Triple(0,0,0)));
   TEST_EQUAL(true, ex->inside(Triple(2,2,2)));
   TEST_EQUAL(false, ex->inside(Triple(3,1,1)));
   TEST_EQUAL(false, ex->inside(Triple(1,-1,1)));
   TEST_EQUAL(false, ex->inside(Triple(1,1,-4)));

   Triple comp;

   TEST_CLOSE(24.28718707889797, ex->min_sqdist(Triple(5,5,5), comp));
   TEST_CLOSE(-2.845299461620748, comp.x);
   TEST_CLOSE(-2.845299461620748, comp.y);
   TEST_CLOSE(-2.845299461620748, comp.z);

   TEST_CLOSE(1.0, ex->min_sqdist(Triple(1,1,-2), comp));
   TEST_CLOSE(0.0, comp.x);
   TEST_CLOSE(0.0, comp.y);
   TEST_CLOSE(1.0, comp.z);

   delete ex;
   delete parser;

   // Cylinder
   f = fopen(TESTFILE, "w");
   fprintf(f, "point [1,0,0]\naxis[1,0,0]\nradius 1\nlength 9\n");
   fclose(f);
   parser = new Parser;
   TEST_EQUAL(true, parser->parse(TESTFILE));

   ex = new Cylinder;

   TEST_EQUAL(true, ex->setup(*parser));
   TEST_CLOSE(28.274333882308138, ex->volume());

   TEST_EQUAL(false, ex->inside(Triple(5.0, 2.0, 0.5)));
   TEST_EQUAL(true, ex->inside(Triple(5.0, 0.0, 0.5)));
   TEST_EQUAL(false, ex->inside(Triple(5.0, 1.0, 0.0)));
   TEST_EQUAL(false, ex->inside(Triple(8.0, 0.0, -1.0)));
   TEST_EQUAL(false, ex->inside(Triple(0.5, 0.5, 0.5)));
   TEST_EQUAL(false, ex->inside(Triple(10.5, 0.5, 0.5)));
   TEST_EQUAL(true, ex->inside(Triple(8.9, 0.5, 0.5)));
   TEST_EQUAL(false, ex->inside(Triple(1.0, 0.0, 0.0)));

   TEST_CLOSE(9.0, ex->min_sqdist(Triple(4,4,0), comp));
   TEST_CLOSE(0.0, comp.x);
   TEST_CLOSE(-3.0, comp.y);
   TEST_CLOSE(0.0, comp.z);
   TEST_CLOSE(16.0, ex->min_sqdist(Triple(4,0,5), comp));
   TEST_CLOSE(0.0, comp.x);
   TEST_CLOSE(0.0, comp.y);
   TEST_CLOSE(-4.0, comp.z);
   TEST_CLOSE(3.3431457505076, ex->min_sqdist(Triple(4,-2,2), comp));
   // nearest point on cylinder: (4, 1/sqrt(2), 1/sqrt(2))
   TEST_CLOSE(0.0, comp.x);
   TEST_CLOSE(1.29289321881345, comp.y);
   TEST_CLOSE(-1.29289321881345, comp.z);
   delete ex;
   delete parser;

   f = fopen(TESTFILE, "w");
   fprintf(f, "point [0,1,0]\naxis[1,1,0]\nradius 5\nlength 2\n");
   fclose(f);
   parser = new Parser;
   TEST_EQUAL(true, parser->parse(TESTFILE));

   ex = new Cylinder;

   TEST_EQUAL(true, ex->setup(*parser));
   TEST_EQUAL(false, ex->inside(Triple(0,1,0)));
   TEST_EQUAL(true, ex->inside(Triple(0.1, 1.1, 0.0)));
   TEST_EQUAL(true, ex->inside(Triple(-0.1, 2.0, 0.0)));
   TEST_EQUAL(false, ex->inside(Triple(4,4,0)));
   TEST_EQUAL(0, ex->inside(Triple(-0.5, 0.5, 0.0)));
   delete ex;
   delete parser;
   
}
#endif
