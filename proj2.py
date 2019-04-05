""" main script for proj2 """
import numpy as np
import matplotlib.pyplot as plt
from meshes import Mesh

from test_FVSolver import Solver
from solver import fluxes
import pickle
# task 1


if False:
    bump = Mesh('meshes/bump0_kfid.gri')
    # test freestream
    for order in [1, 2]:
        FVSolver = Solver(bump, order=order)
        FVSolver.solve(maxIter=1000, freeStreamTest=True)
        FVSolver.plotResiduals()
        print(FVSolver.Rmax[0])
    plt.legend()
    plt.show()

    # run real test
    for order in [1, 2]:
        FVSolver = Solver(bump, order=order)
        FVSolver.solve()
        FVSolver.plotResiduals()
    plt.legend()
    plt.show()

# side project investigation
if False:
    # make the converganece plots
    lvs = 5
    orders = 2

    err_cd = np.zeros((2, lvs))
    err_cl = np.zeros((2, lvs))
    Es = np.zeros((2, lvs))
    dof = np.zeros(lvs)

    for order in [1, 2]:
        for recon_bc_flux in [False, True]:
            for recon_p in [False, True]:
                for lv in range(lvs):
                    bump = Mesh('meshes/bump' + str(lv) + '_kfid.gri', check=False)
                    dof[lv] = bump.nElem

                    FVSolver = Solver(bump, order=order)
                    FVSolver.recon_p = recon_p
                    FVSolver.recon_bc_flux = recon_bc_flux
                    FVSolver.solve()
                    FVSolver.postprocess()
                    # FVSolver.writeSolution('sol_bump_lv' + str(lv) + '_O' + str(order))
                    # fileObject = open('FVSolver_lv'+str(lv) + '_O' + str(order)+'.p', 'wb')
                    # pickle.dump(FVSolver, fileObject)
                    err_cd[order-1, lv] = np.abs(FVSolver.cd - 2.94278e-6)
                    err_cl[order-1, lv] = np.abs(FVSolver.cl - 1.537095)
                    Es[order-1, lv] = FVSolver.Es

                label = 'recon bc: ' + str(recon_bc_flux) + ', recon P: ' + str(recon_p) + ', :'

                plt.figure(1)
                conv_rate = -1*np.polyfit(np.log10(np.sqrt(dof)), np.log10(err_cl[1]), 1)[0]
                plt.loglog(np.sqrt(dof), err_cl[order-1, :], '-o', label=label+"%.2f" % conv_rate)
                print('cl order '+str(order) + label, conv_rate)

                plt.figure(2)
                conv_rate = -1*np.polyfit(np.log10(np.sqrt(dof)), np.log10(err_cd[1]), 1)[0]

                plt.loglog(np.sqrt(dof), err_cd[order-1, :], '-o', label=label+"%.2f" % conv_rate)
                print('cd order '+str(order) + label, conv_rate)

                plt.figure(3)
                conv_rate = -1*np.polyfit(np.log10(np.sqrt(dof)), np.log10(Es[1]), 1)[0]
                plt.loglog(np.sqrt(dof), Es[order-1, :], '-o', label=label+"%.2f" % conv_rate)
                print('Es order '+str(order) + label, conv_rate)

    plt.figure(1)

    plt.ylabel(r'Error')
    plt.xlabel(r'$\sqrt{dof}$')
    plt.title('Convergence of Cl')

    plt.legend()

    plt.figure(2)

    plt.ylabel(r'Error')
    plt.xlabel(r'$\sqrt{dof}$')
    plt.title('Convergence of Cd')
    plt.legend()

    plt.figure(3)

    plt.ylabel(r'Error')
    plt.xlabel(r'$\sqrt{dof}$')
    plt.title('Convergence of Es')

    plt.legend()
    plt.show()

if False:
    # make the converganece plots
    lvs = 5
    orders = 2

    err_cd = np.zeros((orders, lvs))
    err_cl = np.zeros((orders, lvs))
    Es = np.zeros((orders, lvs))
    dof = np.zeros(lvs)

    for order in [1, 2]:
        for lv in range(lvs):
            bump = Mesh('meshes/bump' + str(lv) + '_kfid.gri', check=False)
            dof[lv] = bump.nElem

            FVSolver = Solver(bump, order=order)
            FVSolver.solve()
            FVSolver.postprocess()
            FVSolver.writeSolution('sol_bump_lv' + str(lv) + '_O' + str(order))
            fileObject = open('FVSolver_lv'+str(lv) + '_O' + str(order)+'.p', 'wb')
            pickle.dump(FVSolver, fileObject)
            err_cd[order-1, lv] = np.abs(FVSolver.cd - 2.94278e-6)
            err_cl[order-1, lv] = np.abs(FVSolver.cl - 1.537095)
            Es[order-1, lv] = FVSolver.Es

        label = 'Order:' + str(order) + ' Conv. Rate:'

        plt.figure(1, figsize=(8, 6))
        conv_rate = -1*np.polyfit(np.log10(np.sqrt(dof)), np.log10(err_cl[order-1]), 1)[0]
        plt.loglog(np.sqrt(dof), err_cl[order-1, :], '-o', label=label+"%.2f" % conv_rate)
        print('cl order '+str(order) + label, conv_rate)

        plt.figure(2, figsize=(8, 6))
        conv_rate = -1*np.polyfit(np.log10(np.sqrt(dof)), np.log10(err_cd[order-1]), 1)[0]

        plt.loglog(np.sqrt(dof), err_cd[order-1, :], '-o', label=label+"%.2f" % conv_rate)
        print('cd order '+str(order) + label, conv_rate)

        plt.figure(3, figsize=(8, 6))
        conv_rate = -1*np.polyfit(np.log10(np.sqrt(dof)), np.log10(Es[order-1]), 1)[0]
        plt.loglog(np.sqrt(dof), Es[order-1, :], '-o', label=label+"%.2f" % conv_rate)
        print('Es order '+str(order) + label, conv_rate)

    plt.figure(1)

    plt.ylabel(r'Error', fontsize=14)
    plt.xlabel(r'$\sqrt{dof}$', fontsize=14)
    plt.title('Convergence of Cl', fontsize=14)
    plt.legend()
    plt.savefig('figures/conv_cl', bbox_inches='tight')

    plt.figure(2)

    plt.ylabel(r'Error', fontsize=14)
    plt.xlabel(r'$\sqrt{dof}$', fontsize=14)
    plt.title('Convergence of Cd', fontsize=14)
    plt.legend()
    plt.savefig('figures/conv_cd', bbox_inches='tight')

    plt.figure(3)

    plt.ylabel(r'Error', fontsize=14)
    plt.xlabel(r'$\sqrt{dof}$', fontsize=14)
    plt.title('Convergence of Es', fontsize=14)

    plt.legend()
    plt.savefig('figures/conv_Es', bbox_inches='tight')
    plt.show()


if True:

    for order in [1, 2]:
        bump = Mesh('meshes/bump2_kfid.gri', check=False)

        FVSolver = Solver(bump, order=order)
        FVSolver.solve()
        FVSolver.postprocess()
        FVSolver.plotCP()
    plt.legend()
    plt.gca().invert_yaxis()

    plt.show()
    # import ipdb
    # ipdb.set_trace()
