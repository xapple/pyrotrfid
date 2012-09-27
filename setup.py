from distutils.core import setup

setup(
        name             = 'pyrotrfid',
        version          = '0.9.1',
        description      = 'Digital TRFLP with pyrosequencing reads',
        long_description = open('README.txt').read(),
        license          = 'GNU General Public License 3.0',
        url              = 'http://xapple.github.com/pyrotrfid/',
        author           = 'Lucas Sinclair',
        author_email     = 'lucas.sinclair@me.com',
        classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
        install_requires = ['scipy', 'matplotlib', 'biopython', 'pysam'],
        packages         = ['pyrotrfid'],
        scripts          = ['pyrotrfid/pyrotrfid'],
    )