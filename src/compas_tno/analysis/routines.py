import time
from compas_tno.shapes import Shape
from compas_tno.diagrams import FormDiagram
import compas_tno
import os
import math


def limit_analysis_GSF(analysis, thk, thk_reduction, thk_refined=None, limit_equal=0.01, printout=True, save_forms=False):
    """Routine to compute the succesive max/min optimisation optimisation.

    Parameters
    ----------
    analysis : :class:`~compas_tno.analysis.Analysis`
        The Analysis object
    thk : float
        The initial thickness of the structure
    thk_reduction : float
        Step of reduction of the thickness
    thk_refined : float, optional
        If require a smaller step close to convergence, by default None
    limit_equal : float, optional
        Tolerance for the stopping criteria where the min/max solutions are considered equal, by default 0.01
    printout : bool, optional
        Whether or not printing in the screen the solutions obtained, by default True
    save_forms : bool, optional
        Whether or not solutions should be saved, by default False

    Returns
    -------
    [list, list], [list, list]
        Lists with the thicknesses obtained in min/max computation and the soulutions normalised of min/max thrust.
    """

    solutions_min = []  # empty lists to keep track of the solutions for min thrust
    solutions_max = []  # empty lists to keep track of the solutions for max thrust
    thicknesses_min = []
    thicknesses_max = []
    t0 = thk
    form0 = analysis.form.copy()
    thk_reduction0 = thk_reduction
    data_diagram = analysis.form.parameters
    data_shape = analysis.shape.datashape
    ro = analysis.shape.ro
    last_min = 0
    last_max = 100
    last_thk_min = t0
    objectives = ['min', 'max']
    exitflag = 0
    address = None

    print('Limit Analysis - GSF: For ({0}) with diagram ({1})'.format(data_shape['type'], data_diagram['type']))

    while (last_max - last_min) > limit_equal:

        if exitflag == 0:
            pass
        else:
            objectives = objectives[::-1]
            print('Warning: Did not Achieve precision required. Stopping Anyway.\n')
            break

        for analysis.optimiser.settings['objective'] in objectives:

            analysis.form = form0
            exitflag = 0
            thk = t0
            thk_reduction = thk_reduction0
            count = 0

            print('\n----- Starting the [', analysis.optimiser.settings['objective'], '] problem for intial thk:', thk)
            print('THK  |   Solved  |   Opt.Val |   Opt/W   |   THK red.  |   Setup time  |   Run time')

            while exitflag == 0 and thk > 0:

                # analysis.form = form0  # added for general case
                data_shape['thk'] = thk
                time0 = time.time()
                analysis.shape = Shape.from_library(data_shape)
                analysis.shape.ro = ro
                swt = analysis.shape.compute_selfweight()
                analysis.apply_selfweight()
                analysis.apply_envelope()
                analysis.apply_reaction_bounds()
                analysis.set_up_optimiser()
                setup_time = time.time() - time0
                analysis.run()
                run_time = time.time() - time0 - setup_time
                exitflag = analysis.optimiser.exitflag  # get info if optimisation was succeded ot not
                fopt = analysis.optimiser.fopt  # objective function optimum value
                fopt_over_weight = fopt/swt  # divide by selfweight

                if exitflag == 0:
                    if analysis.optimiser.settings['objective'] == 'min':
                        solutions_min.append(fopt_over_weight)
                        thicknesses_min.append(thk)
                        last_min = fopt_over_weight
                        last_thk_min = thk
                    else:
                        solutions_max.append(fopt_over_weight)
                        thicknesses_max.append(thk)
                        last_max = abs(fopt_over_weight)
                    if save_forms:
                        address = save_forms + '_' + analysis.optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'
                    else:
                        address = os.path.join(compas_tno.get('/temp/'), 'form_' + analysis.optimiser.settings['objective'] + '.json')
                    analysis.form.to_json(address)
                    print('{0}  |   True  |   {1:.1f} |   {2:.6f} |   {3}   | {4:.2f}s  |   {5:.2f}s'.format(thk, fopt, fopt_over_weight, thk_reduction, setup_time, run_time))
                    thk = round(thk - thk_reduction, 5)
                else:
                    if not address:
                        if analysis.optimiser.settings['objective'] == objectives[0]:
                            print('Failed in Initial Thickness and objective ({0})'.format(analysis.optimiser.settings['objective']))
                            objectives = objectives[::-1]
                    elif analysis.optimiser.settings['objective'] == 'max' and last_thk_min <= thk:
                        print('{0}  |   False  |   XXXX |   XXXX |   Load next [min]  | {1:.2f}s   |  {2:.2f}s'.format(thk, setup_time, run_time))
                        # print('Try loading next [min] optimisation - 1')
                        next_min = None
                        for i in range(len(thicknesses_min)):
                            if thicknesses_min[i] < thk or (thicknesses_min[i] == thk and count == 0):
                                next_min = thicknesses_min[i]
                                break
                        if next_min:
                            address = save_forms + '_' + 'min' + '_thk_' + str(100*next_min) + '.json'
                            analysis.form = FormDiagram.from_json(address)
                            thk = next_min
                            exitflag = 0
                        else:
                            print('Tried loading all [min] optimisation')
                            print('---------- End of process -----------', '\n')
                    elif thk_reduction > thk_refined:  # or (analysis.optimiser.settings['objective'] == 'max' and last_thk_min < last_thk_max and thk_reduction > thk_refined/2):
                        analysis.form = FormDiagram.from_json(address)  # reload last solved
                        thk_ = thk
                        thk = round(thk + thk_reduction, 5)
                        thk_reduction = thk_reduction / 2
                        print('{0}  |   False  |   XXXX |   XXXX |   {1}    | {2:.2f}s   |  {3:.2f}s'.format(thk_, thk_reduction, setup_time, run_time))
                        thk = round(thk - thk_reduction, 5)
                        exitflag = 0
                    else:
                        print('---------- End of process -----------', '\n')
                count += 1

    if printout:
        print('------------------ SUMMARY ------------------ ')
        print('ANALYSIS FOR THICKNESS t0:', t0)
        print('thicknesses min ({0}):'.format(len(thicknesses_min)))
        print(thicknesses_min)
        print('thicknesses max ({0}):'.format(len(thicknesses_max)))
        print(thicknesses_max)
        print('min solutions ({0}):'.format(len(solutions_min)))
        print(solutions_min)
        print('max solutions ({0}):'.format(len(solutions_max)))
        print(solutions_max)

    return [thicknesses_min, thicknesses_max],  [solutions_min, solutions_max]


def thk_minmax_GSF(analysis, thk_max, thk_step=0.05, printout=True, save_forms=None, skip_minthk=False):
    """Routine to compute the succesive max/min optimisations starting from the minimum thickness.

    Parameters
    ----------
    analysis : :class:`~compas_tno.analysis.Analysis`
        The Analysis object
    thk_max : float
        Maximum thickness, in which the algorithm should stop
    thk_step : float, optional
        The increasing steps of the algorithm, by default 0.05
    printout : bool, optional
        If prints should be displayed in the screen, by default True
    save_forms : bool, optional
        If forms should be saved as JSON, by default None
    skip_minthk : bool, optional
        If the minimum thickness computation can be skiped (it assumes the minimum thickness as it is), by default False

    Returns
    -------
    [list, list], [list, list]
        Lists with the thicknesses obtained in min/max computation and the soulutions normalised of min/max thrust.
    """

    solutions_min = []  # empty lists to keep track of the solutions for min thrust
    solutions_max = []  # empty lists to keep track of the solutions for max thrust
    thicknesses_min = []
    thicknesses_max = []
    objectives = ['min', 'max']
    ro = analysis.shape.ro
    data_shape = analysis.shape.datashape

    # Find extreme (min thickness) solution:

    print('\n----- Starting the min thk optimisation for starting thk: {0:.4f}'.format(analysis.shape.data['thk']))

    if not skip_minthk:
        time0 = time.time()
        analysis.optimiser.settings['variables'] = ['ind', 'zb', 't']
        analysis.optimiser.settings['objective'] = 't'
        analysis.apply_selfweight()
        analysis.apply_envelope()
        analysis.apply_reaction_bounds()
        analysis.set_up_optimiser()
        setup_time = time.time() - time0
        analysis.run()
        run_time = time.time() - time0 - setup_time
        exitflag = analysis.optimiser.exitflag
    else:
        exitflag = 0
        setup_time = 0.0
        run_time = 0.0

    if exitflag == 0:
        thk_min = analysis.form.attributes['thk']
        # compute_reactions(analysis.form)
        swt = analysis.form.lumped_swt()
        T = analysis.form.thrust()
        T_over_swt = T/swt
        thk_increase = int(thk_step*100)*math.ceil(thk_min*100.0/int(thk_step*100))/100.0 - thk_min
        print('Min THK  |   Solved  |   Thrust |   T/W   |  THK incr. |  Setup time  |   Run time')
        print('{0:.5f}  |   True  |   {1:.1f} |   {2:.6f} |   {3:.2f}   | {4:.2f}s  |   {5:.2f}s'.format(thk_min, T, T_over_swt, thk_increase, setup_time, run_time))

        # STORE

        address0 = os.path.join(compas_tno.get('/temp/'), 'form0.json')
        analysis.form.to_json(address0)

        solutions_min.append(T_over_swt)
        solutions_max.append(-1 * T_over_swt)
        thicknesses_min.append(thk_min)
        thicknesses_max.append(thk_min)

        # SAVE
        if save_forms:
            address_min = save_forms + '_' + 'min' + '_thk_' + str(100*thk_min) + '.json'
            address_max = save_forms + '_' + 'max' + '_thk_' + str(100*thk_min) + '.json'
            analysis.form.to_json(address_min)
            analysis.form.to_json(address_max)

        thk = round(thk_min + thk_increase, 5)
        thk_increase = thk_step
        thk0 = thk

        if thk_min == thk_max:
            print('Warning: Minimum THK Optimisation found optimum at start. Try rescalling the problem.')
            return
    else:
        print('Error: Minimum THK Optimisation did not find a solution: Try scalling the optimisation, or starting in a different thickness value')
        return

    # Start inverse loop for min/max:

    analysis.optimiser.settings['variables'] = ['ind', 'zb']
    for analysis.optimiser.settings['objective'] in objectives:

        analysis.form = FormDiagram.from_json(address0)

        exitflag = 0
        thk = thk0
        count = 0
        first_fail = True

        print('\n----- Starting the inverse [', analysis.optimiser.settings['objective'], '] problem for minimum thk:', thk)
        print('THK  |   Solved  |   Opt.Val |   Opt/W   |   THK incr.  |   Setup time  |   Run time')

        while exitflag == 0 and thk <= thk_max:

            data_shape['thk'] = thk

            time0 = time.time()
            analysis.shape = Shape.from_library(data_shape)
            analysis.shape.ro = ro

            swt = analysis.shape.compute_selfweight()

            lumped_swt = analysis.form.lumped_swt()

            analysis.form.scale_form(swt/lumped_swt)
            lumped_swt = analysis.form.lumped_swt()

            analysis.apply_envelope()
            analysis.apply_reaction_bounds()
            analysis.set_up_optimiser()
            setup_time = time.time() - time0
            analysis.run()
            run_time = time.time() - time0 - setup_time
            exitflag = analysis.optimiser.exitflag  # get info if optimisation was succeded ot not
            fopt = analysis.optimiser.fopt  # objective function optimum value
            fopt_over_weight = fopt/lumped_swt  # divide by selfweight

            if exitflag == 0:
                if analysis.optimiser.settings['objective'] == 'min':
                    solutions_min.append(fopt_over_weight)
                    thicknesses_min.append(thk)
                else:
                    solutions_max.append(fopt_over_weight)
                    thicknesses_max.append(thk)
                if save_forms:
                    address = save_forms + '_' + analysis.optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'
                else:
                    address = os.path.join(compas_tno.get('/temp/'), 'form_' + analysis.optimiser.settings['objective'] + '.json')
                analysis.form.to_json(address)
                print('{0:.5f}  |   True  |   {1:.1f} |   {2:.6f} |   {3}   | {4:.2f}s  |   {5:.2f}s'.format(thk, fopt, fopt_over_weight, thk_increase, setup_time, run_time))
                thk = thk + thk_increase
                first_fail = True
            else:
                print('Failed Optimisation for [{0}] with thk: {1}'.format(analysis.optimiser.settings['objective'], thk))
                other_objective = 'min' if analysis.optimiser.settings['objective'] == 'max' else 'max'
                if thk < thk_max:
                    if not first_fail:
                        thk = thk + thk_increase
                        first_fail = True
                        print('Second Fail with [{0}] opt, will load [{1}] opt in incr. thk: {2}'.format(analysis.optimiser.settings['objective'], other_objective, thk))
                    else:
                        print('First Fail with [{0}] opt, will load [{1}] opt in thk: {2}'.format(analysis.optimiser.settings['objective'], other_objective, thk))
                        first_fail = False
                    address_other_obj = save_forms + '_' + other_objective + '_thk_' + str(100*thk) + '.json'
                    analysis.form = FormDiagram.from_json(address_other_obj)
                    exitflag = 0
            count += 1
        print('---------- End of process -----------', '\n')

    thicknesses_min.reverse()
    thicknesses_max.reverse()
    solutions_min.reverse()
    solutions_max.reverse()

    if printout:
        print('------------------ SUMMARY ------------------ ')
        print('ANALYSIS FOUND MIN THK:', thk_min)
        print('thicknesses min ({0}):'.format(len(thicknesses_min)))
        print(thicknesses_min)
        print('thicknesses max ({0}):'.format(len(thicknesses_max)))
        print(thicknesses_max)
        print('min solutions ({0}):'.format(len(solutions_min)))
        print(solutions_min)
        print('max solutions ({0}):'.format(len(solutions_max)))
        print(solutions_max)

    return [thicknesses_min, thicknesses_max],  [solutions_min, solutions_max]


def max_n_minmax_GSF(analysis, n_step=0.01, printout=True, save_forms=False):
    """Routine to compute max/min thrust optimisation problems for non-analytic surfaces, considering the nomal vectors 'n'.

    Parameters
    ----------
    analysis : :class:`~compas_tno.analysis.Analysis`
        The Analysis object
    n_step : float, optional
        The step to decrease the normal vector, by default 0.01
    printout : bool, optional
        If prints are added to the screen, by default True
    save_forms : bool, optional
        If forms are saved as JSON, by default False

    Returns
    -------
    [list, list], [list, list]
        Lists with the thicknesses obtained in min/max computation and the soulutions normalised of min/max thrust.
    """

    solutions_min = []  # empty lists to keep track of the solutions for min thrust
    solutions_max = []  # empty lists to keep track of the solutions for max thrust
    thicknesses_min = []
    thicknesses_max = []
    objectives = ['min', 'max']
    ro = analysis.shape.ro
    data_shape = analysis.shape.datashape
    thk0 = data_shape['thk']
    t = data_shape['t']
    initial_intrados = analysis.shape.intrados.copy()
    initial_extrados = analysis.shape.extrados.copy()
    middle = analysis.shape.middle.copy()

    # Find extreme (min thickness) solution:

    print('\n----- Starting the min thk optimisation for approx. thk: {0:.4f}'.format(thk0))

    time0 = time.time()
    analysis.optimiser.settings['variables'] = ['ind', 'zb', 'n']
    analysis.optimiser.settings['objective'] = 'n'
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.apply_reaction_bounds()
    analysis.set_up_optimiser()
    setup_time = time.time() - time0
    analysis.run()
    run_time = time.time() - time0 - setup_time

    exitflag = analysis.optimiser.exitflag

    if exitflag == 0:
        n = - 1 * analysis.optimiser.fopt
        thk = thk0 - 2 * n
        thk_min = thk
        T = analysis.form.thrust()
        swt = analysis.shape.compute_selfweight()
        T_over_swt = T/swt
        n_reduction = n - int(n_step*100)*math.floor(n*100.0/int(n_step*100))/100.0
        print('n (offset) | Min THK  |   Solved  |   Thrust |   T/W   |  n decr. |  Setup time  |   Run time')
        print('{0:.5f}  | {1:.5f}  |   True  |   {2:.1f} |   {3:.6f} |   {4:.4f}   | {5:.2f}s  |   {6:.2f}s'.format(n, thk_min, T, T_over_swt,
                                                                                                                    n_reduction, setup_time, run_time))

        # STORE

        address0 = os.path.join(compas_tno.get('/temp/'), 'form0.json')
        analysis.form.to_json(address0)

        solutions_min.append(T_over_swt)
        solutions_max.append(-1 * T_over_swt)
        thicknesses_min.append(thk)
        thicknesses_max.append(thk)

        # SAVE
        if save_forms:
            address_min = save_forms + '_' + 'min' + '_thk_' + str(100*thk) + '.json'
            address_max = save_forms + '_' + 'max' + '_thk_' + str(100*thk) + '.json'
            analysis.form.to_json(address_min)
            analysis.form.to_json(address_max)

        if n < 0:
            print('Warning: Reduction of cross section (n), or offset is negative:', n_reduction)
            return

        n0 = round(n - n_reduction, 5)
        print('from here we ewill use this n:', n0)
        n_reduction = n_step

        print('Changing solver to IPOPT')
        analysis.optimiser.settings['solver'] == 'IPOPT'

    else:
        print('Error: Minimum THK Optimisation did not find a solution: Try entering a different geometry')
        return

    # Start inverse loop for min/max:

    analysis.optimiser.settings['variables'] = ['ind', 'zb']
    for analysis.optimiser.settings['objective'] in objectives:

        analysis.form = FormDiagram.from_json(address0)

        exitflag = 0
        n = n0
        thk = thk0 - 2 * n
        # print('resulting in this start thickness', thk)
        count = 0
        first_fail = True

        print('\n----- Starting the inverse [', analysis.optimiser.settings['objective'], '] problem for minimum thk:', thk)
        print('n (offset) | THK  |   Solved  |   Opt.Val |   Opt/W   |   THK incr.  |   Setup time  |   Run time')

        while exitflag == 0 and thk <= thk0:

            time0 = time.time()
            intrados = initial_intrados.offset_mesh(n=n, direction='up')
            extrados = initial_extrados.offset_mesh(n=n, direction='down')
            analysis.shape = Shape.from_meshes_and_formdiagram(analysis.form, intrados, extrados, middle=middle, data={'type': 'general', 't': t, 'thk': thk})
            analysis.shape.ro = ro
            swt = analysis.shape.compute_selfweight()

            lumped_swt = analysis.form.lumped_swt()
            analysis.apply_selfweight()
            analysis.form.scale_form(swt/lumped_swt)

            analysis.apply_envelope()
            analysis.apply_reaction_bounds()
            analysis.set_up_optimiser()
            setup_time = time.time() - time0
            analysis.run()
            run_time = time.time() - time0 - setup_time
            exitflag = analysis.optimiser.exitflag  # get info if optimisation was succeded ot not
            fopt = analysis.optimiser.fopt  # objective function optimum value
            fopt_over_weight = fopt/swt  # divide by selfweight

            if exitflag == 0:
                if analysis.optimiser.settings['objective'] == 'min':
                    solutions_min.append(fopt_over_weight)
                    thicknesses_min.append(thk)
                else:
                    solutions_max.append(fopt_over_weight)
                    thicknesses_max.append(thk)
                if save_forms:
                    address = save_forms + '_' + analysis.optimiser.settings['objective'] + '_thk_' + str(100*thk) + '.json'
                else:
                    address = os.path.join(compas_tno.get('/temp/'), 'form_' + analysis.optimiser.settings['objective'] + '.json')
                analysis.form.to_json(address)
                print('{0:.5f}  | {1:.5f}  |   True  |   {2:.1f} |   {3:.6f} |   {4}   | {5:.2f}s  |   {6:.2f}s'.format(
                    n, thk, fopt, fopt_over_weight, n_reduction, setup_time, run_time))
                n = n - n_reduction
                thk = thk0 - 2*n
                first_fail = True
            else:
                print('Failed Optimisation for [{0}] with thk: {1}'.format(analysis.optimiser.settings['objective'], thk))
                other_objective = 'min' if analysis.optimiser.settings['objective'] == 'max' else 'max'
                if thk < thk0:
                    if not first_fail:
                        n = n - n_reduction
                        thk = thk0 - 2*n
                        first_fail = True
                        print('Second Fail with [{0}] opt, will load [{1}] opt in incr. thk: {2}'.format(analysis.optimiser.settings['objective'], other_objective, thk))
                    else:
                        print('First Fail with [{0}] opt, will load [{1}] opt in thk: {2}'.format(analysis.optimiser.settings['objective'], other_objective, thk))
                        first_fail = False
                    address_other_obj = save_forms + '_' + other_objective + '_thk_' + str(100*thk) + '.json'
                    analysis.form = FormDiagram.from_json(address_other_obj)
                    exitflag = 0
            count += 1
        print('---------- End of process -----------', '\n')

    if printout:
        print('------------------ SUMMARY ------------------ ')
        print('ANALYSIS FOUND MIN THK:', thk_min)
        print('thicknesses min ({0}):'.format(len(thicknesses_min)))
        print(thicknesses_min)
        print('thicknesses max ({0}):'.format(len(thicknesses_max)))
        print(thicknesses_max)
        print('min solutions ({0}):'.format(len(solutions_min)))
        print(solutions_min)
        print('max solutions ({0}):'.format(len(solutions_max)))
        print(solutions_max)

    return [thicknesses_min, thicknesses_max],  [solutions_min, solutions_max]
