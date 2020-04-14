""" Nuclear Data Scripts

This module is a collection of scripts to convert various nucler data formats

"""
from cream.endf import ENDFTape
from cream.atom_relax import AtomRelax, Z_EL
import click


def read_endf_files(files):
    """ Read a collection of ENDF files from paths and convert to ENDFTapes

    Args:
        files (list of str): Paths to ENDF files

    Result:
        List of open ENDFTapes
    """
    return [ENDFTape(path) for path in list(files)]


def read_relaxation(tapes):
    """ Read atomic relaxation data from a list of tapes

    Args:
        tapes (list of ENDFTape): Search for atomic relaxation data

    Result:
        Dictionary of {MAT: AtomRelax} pairs

    Raises:
        ValueError: If in collection of materials multiple data for the same
            atom (MAT) is provided.
    """
    relax = dict()
    for t in tapes:
        for MAT, _, _ in t.get(MF=28):
            if MAT in relax.keys():
                msg = ('Atomic Relaxation data for {} is repeated. First was'
                       ' in {}. Duplicate data is found in {}'
                       .format(Z_EL[MAT/100], relax[MAT].source, t))
                raise ValueError(msg)
            atom = AtomRelax()
            atom.from_endf(MAT/100, t)
            relax[MAT] = atom
    return relax


def serpent_relaxation_lib(files, E_cut=1.0e3):
    """ Return Serpent-formatted atomic relaxation data library

    Args:
        files (list of str): Paths to ENDF files with relaxation data
        E_cut (float, optional): Minimum Photon Energy to prune needless
            transition data [eV]. Default=1.0e3

    Result:
        (tuple): tuple containing:
            (str) Relaxation data library
            (str) Electronic ground state Library
    """
    atoms = read_relaxation(read_endf_files(files))

    # Prune and check
    [a.prune_energy(E_cut) for a in atoms.values()]
    [a.check() for a in atoms.values()]

    # Build relaxation library
    relaxLib = str()
    for key, val in sorted(atoms.items()):
        relaxLib += val.print_trans()

    # Build ground state library
    groundLib = 'Nelement {}\n'.format(len(atoms))
    for key, val in sorted(atoms.items()):
        groundLib += val.print_ground()
    return relaxLib, groundLib


@click.group()
def data():
    "Manage nuclear data. "
    pass


@data.command()
@click.argument('input',
                type=click.Path(exists=True),
                nargs=-1)
@click.option('-o1', '--output-trans',
              required=True,
              type=click.Path(),
              help='Transition data output file',
              prompt='Transition data output')
@click.option('-o2', '--output-ground',
              required=True,
              type=click.Path(),
              help='Electronic ground state data output file',
              prompt='Ground state data output')
@click.option('-f', '--format',
              type=click.Choice(['SCONE', 'Serpent'], case_sensitive=False),
              default='SCONE',
              show_default=True,
              help='Format of output library')
def relaxation(input, output_trans, output_ground, format):
    """Generate atomic relaxation data library from ENDF.

    [INPUT] is a list of ENDF files you want to include in the library

    """

    # Gracefull exit if no files were given
    if len(input) == 0:
        print('Failed. No files were given')
        return

    # Generate libraries as strings
    if format.lower() == 'scone':
        print('bla')
    elif format.lower() == 'serpent':
        res = serpent_relaxation_lib(input)

    # Print to files
    with open(output_trans, 'w') as f:
        f.write(res[0])
    with open(output_ground, 'w') as f:
        f.write(res[1])

    # Print end
    print('Success!')
