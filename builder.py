#!/usr/bin/env python3
import numpy as np
import argparse

##### constants for types of atoms, masses and bond types
TYPES = {
	'GEM2': 1,
	'GEM4': 2, 
}
MASSES = {
	'GEM2': 1.0,
	'GEM4': 1.0,
}
BONDS = {
	'BOND': 1, # is this equilibrium bond length, or is it something else?
}

##### helper functions
class Lammps_universe():
  '''
  class in which the lammps universe is built and dumped
  '''
  def __init__(self, box_L = 10.0, out_fn = 'data.lmp'):
    self.atoms = []
    self.bonds = []
    self.Nmol = 0
    self.box_L = box_L
    self.out_fn = out_fn

  def add_dimer(self, xyz, bond_length = 1.0):
    '''
    add a dimer to position .xyz with uniformly random orientations, log the atom and bond info
    '''
    self.Nmol += 1
    orientation = np.random.normal(size = 3)
    orientation /= np.linalg.norm(orientation)

    for nind, name in enumerate(['GEM2', 'GEM4']):
      atom_id = len(self.atoms) + 1
      self.atoms.append({
        'type' : TYPES[name],
        'id'   : atom_id,
        'molid': self.Nmol,
        'xyz'  : xyz + (nind - 0.5) * bond_length * orientation,
      })

    bond_id = len(self.bonds) + 1
    self.bonds.append({
      'ids': (atom_id, atom_id - 1),
      'id_bond': bond_id,
      'type': BONDS['BOND'],
    })

  def initialize_random_liquid(self, N_molecules = 100):
    '''
    initialize whole system of dimers at uniformly random positions
    '''
    for n in range(N_molecules):
      center_of_mass = self.box_L * np.random.random(3)
      self.add_dimer(center_of_mass)

  def dump_data(self, comment = ''):
    '''
    dump the current contents of the universe into the data file
    '''
    with open(self.out_fn, 'w') as out_f:
      out_f.write('LAMMPS GEM dimers {:}\n'.format(comment))
      out_f.write('{:d} atoms\n{:d} bonds\n{:d} atom types\n{:d} bond types\n\n'.format(len(self.atoms), len(self.bonds), len(TYPES), len(BONDS)))
      for k in ['x', 'y', 'z']:
        out_f.write('{:.4f} {:.4f} {:}lo {:}hi\n'.format(0, self.box_L, k, k))
      out_f.write('\nMasses\n\n')
      for t in TYPES:
        out_f.write('{:d} {:.6f}\n'.format(TYPES[t], MASSES[t]))
      out_f.write('\nAtoms\n\n')
      for atom in self.atoms:
        out_f.write('{:d} {:d} {:d} {:.6f} {:.6f} {:.6f}\n'.format(atom['id'], atom['molid'], atom['type'], *atom['xyz']))
      out_f.write('\nBonds\n\n')
      for bond in self.bonds:
        out_f.write('{:d} {:d} {:d} {:d}\n'.format(bond['id_bond'], bond['type'], *bond['ids']))
    out_f.close()

##### parsing the input
parser = argparse.ArgumentParser()
parser.add_argument(
            "--N_molecules",
            required = True,
            type = int,
            help = 'number of molecules',
            )
parser.add_argument(
            "--density",
            required = True,
            type = float,
            help = 'density of molecules',
            )
parser.add_argument(
            "--output",
            type = str,
            required = True,
            help = 'input filename',
            )
args = parser.parse_args()

##### generate and dump the system
box_L = (args.N_molecules / args.density)**(1/3)
ls = Lammps_universe(box_L = box_L, out_fn = args.output)
ls.initialize_random_liquid(N_molecules = args.N_molecules)
ls.dump_data()

exit()
