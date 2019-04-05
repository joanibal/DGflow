""" main script for proj2 """
import numpy as np
import matplotlib.pyplot as plt
from meshes import Mesh

from cfdsolvers import DGSolver
import pickle
# task 1


def bumpShape(x):
    return 0.0625*np.exp(-25*x**2)


if False:
    # bump = Mesh('meshes/bump0_kfid.gri')
    # test freestream

    bump = Mesh('meshes/bump0_kfid.gri', wallGeomFunc=bumpShape)
    for order in [0, 1, 2]:
        CFDSolver = DGSolver(bump, order=order)
        CFDSolver.testFreestream()
        CFDSolver.plotResiduals()
        print(CFDSolver.Rmax[0])
    plt.legend()
    plt.show()

    # # run real test
    for order in [0, 1, 2]:
        CFDSolver = DGSolver(bump, order=order)
        print(1.94*2.02**-order)
        CFDSolver.solve(cfl=1.94*2.02**-order)
        CFDSolver.plotResiduals()

    plt.legend()
    plt.show()


if False:
    # make the converganece plots
    # lvs = 2
    orders = 3

    # err_cd = np.zeros((orders, lvs))
    # err_cl = np.zeros((orders, lvs))
    # Es = np.zeros((orders, lvs))
    # wallTime = np.zeros((orders, lvs))
    # dof = np.zeros(lvs)

    data = {}

    # bump = Mesh('meshes/bump0_kfid.gri', wallGeomFunc=bumpShape)

    for order in range(orders):
        nlvs = 5 - order
        data[order] = {
            'err_cd': np.zeros(nlvs),
            'err_cl': np.zeros(nlvs),
            'Es': np.zeros(nlvs),
            'wallTime': np.zeros(nlvs),
            'dof': np.zeros(nlvs)
        }
        bump = Mesh('meshes/bump0_kfid.gri', wallGeomFunc=bumpShape)

        for lv in range(nlvs):
            if lv > 0:
                bump.refine()

            print('===========', order, lv, '====================')
            CFDSolver = DGSolver(bump, order=order)
            data[order]['dof'][lv] = CFDSolver.dof
            CFDSolver.solve(maxIter=100000, cfl=((1.94 - lv*0.25)*2.02**-order))

            CFDSolver.postprocess()
            CFDSolver.writeSolution('sol_bump_lv' + str(lv) + '_Ord' + str(order))
            data[order]['err_cd'][lv] = np.abs(CFDSolver.cd - 2.94278e-6)
            data[order]['err_cl'][lv] = np.abs(CFDSolver.cl - 1.537095)
            data[order]['Es'][lv] = CFDSolver.Es
            data[order]['wallTime'][lv] = CFDSolver.wallTime

        label = 'Order:' + str(order) + ' Conv. Rate:'

        plt.figure(1, figsize=(8, 6))
        conv_rate = -1*np.polyfit(np.log10(np.sqrt(data[order]['dof'])), np.log10(data[order]['err_cl']), 1)[0]
        plt.loglog(np.sqrt(data[order]['dof']), data[order]['err_cl'], '-o', label=label+"%.2f" % conv_rate)
        print('cl order '+str(order) + label, conv_rate)

        plt.figure(2, figsize=(8, 6))
        conv_rate = -1*np.polyfit(np.log10(np.sqrt(data[order]['dof'])), np.log10(data[order]['err_cd']), 1)[0]

        plt.loglog(np.sqrt(data[order]['dof']), data[order]['err_cd'], '-o', label=label+"%.2f" % conv_rate)
        print('cd order '+str(order) + label, conv_rate)

        plt.figure(3, figsize=(8, 6))
        conv_rate = -1*np.polyfit(np.log10(np.sqrt(data[order]['dof'])), np.log10(data[order]['Es']), 1)[0]
        plt.loglog(np.sqrt(data[order]['dof']), data[order]['Es'], '-o', label=label+"%.2f" % conv_rate)
        print('Es order '+str(order) + label, conv_rate)

        plt.figure(4, figsize=(8, 6))
        plt.loglog(data[order]['wallTime'], data[order]['Es'], '-o', label='Order:' + str(order))

    plt.figure(1)

    plt.ylabel('Error', fontsize=14)
    plt.xlabel(r'$\sqrt{dof}$', fontsize=14)
    plt.title('Convergence of Cl', fontsize=14)
    plt.legend()
    plt.savefig('figures/conv_cl', bbox_inches='tight')

    plt.figure(2)

    plt.ylabel('Error', fontsize=14)
    plt.xlabel(r'$\sqrt{dof}$', fontsize=14)
    plt.title('Convergence of Cd', fontsize=14)
    plt.legend()
    plt.savefig('figures/conv_cd', bbox_inches='tight')

    plt.figure(3)

    plt.ylabel('Error', fontsize=14)
    plt.xlabel(r'$\sqrt{dof}$', fontsize=14)
    plt.title('Convergence of Es', fontsize=14)

    plt.legend()
    plt.savefig('figures/conv_Es', bbox_inches='tight')

    plt.figure(4)

    plt.ylabel('Es Error', fontsize=14)
    plt.xlabel('Walltime [s]', fontsize=14)
    plt.title('Runtime cost for accuracy', fontsize=14)

    plt.legend()
    plt.savefig('figures/walltime', bbox_inches='tight')

    plt.show()

    fileObject = open('convergence_data.p', 'wb')
    pickle.dump(data, fileObject)
    fileObject.close()


if True:
    # compare the performance of DG with FVM
    fileObject = open('./data/FV_convergence_data.p', 'rb')
    dataFV = pickle.load(fileObject)
    fileObject.close()
    fileObject = open('./data/DG_convergence_data.p', 'rb')
    dataDG = pickle.load(fileObject)
    fileObject.close()

    # plot DG dat
    for order in dataDG.keys():

        label = 'DG p=' + str(order)

        plt.figure(1, figsize=(8, 6))
        plt.loglog(np.sqrt(dataDG[order]['dof']), dataDG[order]['err_cl'], '-o', label=label)

        plt.figure(2, figsize=(8, 6))
        plt.loglog(np.sqrt(dataDG[order]['dof']), dataDG[order]['err_cd'], '-o', label=label)

        plt.figure(3, figsize=(8, 6))
        plt.loglog(np.sqrt(dataDG[order]['dof']), dataDG[order]['Es'], '-o', label=label)

        plt.figure(4, figsize=(8, 6))
        plt.loglog(dataDG[order]['wallTime'], dataDG[order]['Es'], '-o', label=label)

    # plot FV data
    for order in dataFV.keys():

        label = 'FV Order:' + str(order)

        plt.figure(1, figsize=(8, 6))
        plt.loglog(np.sqrt(dataFV[order]['dof']), dataFV[order]['err_cl'], '--s', label=label)

        plt.figure(2, figsize=(8, 6))
        plt.loglog(np.sqrt(dataFV[order]['dof']), dataFV[order]['err_cd'], '--s', label=label)

        plt.figure(3, figsize=(8, 6))
        plt.loglog(np.sqrt(dataFV[order]['dof']), dataFV[order]['Es'], '--s', label=label)

        plt.figure(4, figsize=(8, 6))
        plt.loglog(dataFV[order]['wallTime'], dataFV[order]['Es'], '--s', label=label)

    plt.figure(1)

    plt.ylabel('Error', fontsize=14)
    plt.xlabel(r'$\sqrt{dof}$', fontsize=14)
    plt.title('Comparison of Cl Convergence', fontsize=14)
    plt.legend()
    plt.savefig('figures/comp_cl', bbox_inches='tight')

    plt.figure(2)

    plt.ylabel('Error', fontsize=14)
    plt.xlabel(r'$\sqrt{dof}$', fontsize=14)
    plt.title('Comparison of Cd Convergence', fontsize=14)
    plt.legend()
    plt.savefig('figures/comp_cd', bbox_inches='tight')

    plt.figure(3)

    plt.ylabel('Error', fontsize=14)
    plt.xlabel(r'$\sqrt{dof}$', fontsize=14)
    plt.title('Comparison of Es Convergence', fontsize=14)

    plt.legend()
    plt.savefig('figures/comp_Es', bbox_inches='tight')

    plt.figure(4)

    plt.ylabel('Es Error', fontsize=14)
    plt.xlabel('Walltime [s]', fontsize=14)
    plt.title('Comparison of Runtime Cost for Accuracy', fontsize=14)

    plt.legend()
    plt.savefig('figures/comp_walltime', bbox_inches='tight')

    plt.show()


if False:

    bump = Mesh('meshes/bump0_kfid.gri', wallGeomFunc=bumpShape)
    bump.refine()
    bump.refine()
    for order in [0, 1, 2]:
        CFDSolver = DGSolver(bump, order=order)
        CFDSolver.solve(maxIter=100000, cfl=((1.56)*2.02**-order))

        # FVSolver.solve()
        CFDSolver.postprocess()
        CFDSolver.plotCP()
    plt.legend()
    plt.gca().invert_yaxis()

    plt.show()
    # import ipdb
    # ipdb.set_trace()
