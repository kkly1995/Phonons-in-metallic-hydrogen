import numpy as np

from nexus import settings, generate_physical_system
from nexus import generate_pwscf, job, generate_pw2qmcpack
from nexus import generate_qmcpack, QmcpackInput, run_project
from nexus import loop, linear, vmc, dmc
from qmcpack_input import force

# geometry
lattice_parameter = 2.32 # in bohr
c_over_a = 2.475

# supercell
lattice_vectors = lattice_parameter*np.eye(3)
lattice_vectors[0,0] *= 4
lattice_vectors[1,1] *= 4
lattice_vectors[2,2] *= 2
lattice_vectors[2,2] *= c_over_a
atomic_positions = np.loadtxt('scaled_positions')
atomic_positions = np.matmul(atomic_positions, lattice_vectors)

## just pick one displacement
#d = 0

settings(
        pseudo_dir  = '/u/sciteam/ly1/pseudo',
        results     = '',
        sleep       = 3,
        generate_only = 1,
        machine     = 'bluewaters_xe'
        )

for d in range(70):
    system = generate_physical_system(
            units       = 'B',
            axes        = lattice_vectors,
            elem        = ['H']*128,
            pos         = atomic_positions + np.loadtxt('displacements/%s.dat' % d),
            )

    scf = generate_pwscf(
            system      = system,
            identifier  = 'scf',
            path        = 'LDA/scf%s' % d,
            job         = job(cores=32, app='pw.x'),
            input_type  = 'generic',
            calculation = 'scf',
            ecutwfc     = 300,
            occupations = 'smearing',
            smearing    = 'fd',
            degauss     = 0.0063338,
            pseudos     = ['H.qmc.lda.upf'],
            nosym       = True,
            tprnfor     = True,
            wf_collect  = True,
            kgrid       = (7, 7, 7),
            kshift      = (1, 1, 1)
            )

    p2q = generate_pw2qmcpack(
            identifier  = 'p2q',
            path        = 'LDA/scf%s' % d,
            job         = job(cores=1, app='pw2qmcpack.x'),
            write_psir  = False,
            dependencies= (scf, 'orbitals'),
            )

    # load jastrow
    wf = QmcpackInput('jastrows.xml')
    #wf = QmcpackInput('jastrows_no_kspace.xml')
    jastrows = wf.qmcsystem.wavefunction.jastrows

    #opt = generate_qmcpack(
    #        system      = system,
    #        identifier  = 'opt',
    #        path        = 'opt',
    #        job         = job(cores=32, app='qmcpack'),
    #        input_type  = 'basic',
    #        pseudos     = [],
    #        corrections = [],
    #        jastrows    = jastrows,
    #        twistnum    = 0,
    #        calculations= [loop(max = 5,
    #            qmc     = linear(
    #                timestep    = 0.5,
    #                samples     = 10000,
    #                warmupsteps = 100,
    #                blocks      = 100,
    #                substeps    = 5,
    #                usedrift    = True,
    #                minmethod   = 'OneShiftOnly',
    #                )
    #            )],
    #        dependencies= (p2q, 'orbitals'),
    #        )

    qmc = generate_qmcpack(
            system      = system,
            identifier  = 'qmc',
            path        = 'LDA/qmc%s' % d,
            job         = job(nodes=172, cores=5504, app='qmcpack'),
            input_type  = 'basic',
            pseudos     = [],
            corrections = [],
            jastrows    = jastrows,
            lr_handler  = 'ewald',
            lr_dim_cutoff   = 20,
            estimators  = [force(
                mode    = 'cep',
                addionion   = True,
                rcut    = 1.0,
                nbasis  = 4,
                weightexp   = 2,
                )],
            calculations=[
                vmc(
                    warmupsteps = 100,
                    steps   = 10,
                    blocks  = 800,
                    substeps= 5,
                    timestep= 0.5,
                    usedrift= True,
                    samples = 512,
                ),
                dmc(
                    warmupsteps = 100,
                    steps   = 10,
                    blocks  = 800,
                    timestep= 0.03,
                    )
                ],
            dependencies= [
                (p2q, 'orbitals')
                ]
            )

run_project()
