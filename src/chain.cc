// ============================================================================
// chain.cc -- AMOCH::Chain methods
// ----------------------------------------------------------------------------
// Copyright
// ----------------------------------------------------------------------------
// LICENSE
// ============================================================================

#include "chain.h" 

using AMOCH::Chain;
using AMOCH::Triple;
using AMOCH::Real;

// ============================================================================
// Constructor
// ============================================================================
Chain::Chain(int type,
             int max_monomers,
             int max_particles) :
   _type(type),
   _max_monomers(max_monomers),
   _max_particles(max_particles),
   _nmonomers(0),
   _nparticles(0),
   _monomer_start(max_monomers, -1),
   _torsion_count(360, 0),
   _length(max_monomers, 0.0),
   _com(max_monomers, Triple()),
   _m0(max_monomers, Triple()),
   _cell(max_particles, -1),
   _slot(max_particles, -1),
   tail_index(-1)
{}

// ============================================================================
// Add cell, slot info for particle n
// ============================================================================
void 
Chain::add_particle_cell(int n,
                         int cell,
                         int slot)
{
   _cell[n] = cell;
   _slot[n] = slot;
}

// ============================================================================
// Add a monomer consisting of nparticles to the Chain
// ============================================================================
void 
Chain::add_monomer(int nparticles)
{
   _monomer_start[_nmonomers++] = _nparticles;
   _nparticles += nparticles;
}

// ============================================================================
// Record a final torsion angle
// ============================================================================
void 
Chain::record_torsion(Real angle)
{
   while (angle < 0.0)
      angle += 360.0;
   while (angle >= 360.0)
      angle -= 360.0;
   _torsion_count[(int)angle]++;
}

// ============================================================================
// Record a final particle position, after all building is complete
// ============================================================================
void 
Chain::record_particle(int index,
                       const Triple& pos)
{
   int j;

   for (j = 0; j < _max_monomers-1; j++) {
      if ((index >= _monomer_start[j]) && (index < _monomer_start[j+1]))
         break;
   }
   // particle is in monomer j
   _com[j] += pos;
   if (index == _monomer_start[j])  // first particle in monomer j
      _m0[j] = pos;
}

// ============================================================================
// Calculate final chain statistics
// ============================================================================
void 
Chain::finalize(const Triple& tail_position)
{
   int j, N;
   Triple d;
   Triple head = _m0[0];

   for (j = 0; j < _max_monomers-1; j++) {
      N = _monomer_start[j+1] - _monomer_start[j];
      if (N > 1)
         _com[j] *= 1.0/(Real)N;
      d = _m0[j+1] - head;
      _length[j] = d.length();
   }
   // j == _max_monomers-1
   N = _nparticles - _monomer_start[j];
   if (N > 1)
      _com[j] *= 1.0/(Real)N;
   d = tail_position - head;
   _length[j] = d.length();
}

// ============================================================================
// Write internal state to the output FILE f
// ============================================================================
void 
Chain::write(FILE *f) const
{
   fwrite((void *)&_type, sizeof(int), 1, f);
   fwrite((void *)&_max_monomers, sizeof(int), 1, f);
   fwrite((void *)&_max_particles, sizeof(int), 1, f);
   fwrite((void *)&_nmonomers, sizeof(int), 1, f);
   fwrite((void *)&_nparticles, sizeof(int), 1, f);
   fwrite((void *)&_monomer_start[0], sizeof(int), _max_monomers, f);
   fwrite((void *)&_torsion_count[0], sizeof(int), 360, f);
   fwrite((void *)&_length[0], sizeof(Real), _max_monomers, f);
   // Do NOT write _com: any scenario that requires reading Chains will also
   // require reading particles (e.g. LAMMPS dump file); any final stats will
   // require calling record_particle() FOR EVERY PARTICLE -- _com sums will
   // be corrupted with duplicates and rescaled incorrectly in finalize()
   fwrite((void *)&_m0[0], sizeof(Triple), _max_monomers, f);
   fwrite((void *)&_cell[0], sizeof(int), _max_particles, f);
   fwrite((void *)&_slot[0], sizeof(int), _max_particles, f);
   fwrite((void *)&tail_index, sizeof(int), 1, f);
}

// ============================================================================
// Read internal state from the input FILE f
// ============================================================================
void 
Chain::read(FILE *f)
{
   fread((void *)&_type, sizeof(int), 1, f);
   fread((void *)&_max_monomers, sizeof(int), 1, f);
   fread((void *)&_max_particles, sizeof(int), 1, f);
   fread((void *)&_nmonomers, sizeof(int), 1, f);
   fread((void *)&_nparticles, sizeof(int), 1, f);
   fread((void *)&_monomer_start[0], sizeof(int), _max_monomers, f);
   fread((void *)&_torsion_count[0], sizeof(int), 360, f);
   fread((void *)&_length[0], sizeof(Real), _max_monomers, f);
   // Do not read _com; see comments above in write()
   fread((void *)&_m0[0], sizeof(Triple), _max_monomers, f);
   fread((void *)&_cell[0], sizeof(int), _max_particles, f);
   fread((void *)&_slot[0], sizeof(int), _max_particles, f);
   fread((void *)&tail_index, sizeof(int), 1, f);
}

