import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess

#This script extracts tensor flow validation outputs for an ensemble
#   of neural network potentials, including atomic energy RMSE and
#   force component RMSE. At the end, it prints out both to allow
#   choice of NNP minimizing test set RMSE.

mother_directory = '/pscratch/sd/p/pedrogm/olsalz-project/validate/alex_md_structures'
val_interv = 1000

def get_last_step(folder_name):
    
    directory = mother_directory + folder_name + '/test_mol_set'

    %cd {directory}
    command = r"ls ./ | sed 's/[^0-9]*\([0-9]\+\).*/\1/g' | sort -n | tail -1"
    shell_output_last_step = subprocess.check_output(command, shell=True)
    last_step =int(shell_output_last_step.decode("utf-8").strip())

    return last_step

def collect_validation(val_interv, folder_name):

    directory = mother_directory + folder_name + '/test_mol_set'
    %cd {directory}
    last_step = get_last_step(folder_name)
    print(last_step)

    steps = [i for i in range(0, last_step + 1, val_interv)]

    e_rmse = []
    force_rmse = []
    force_norm_rmse = []

    for i in steps:

        df = pd.read_csv(str(i) + '.dat', delim_whitespace=True)
        number = np.array(df["n_atoms"])
        e_ref = np.array(df["e_ref"])
        e_nn = np.array(df["e_nn"])

        e_ref_p_atom = np.divide(e_ref, number)
        e_nn_p_atom = np.divide(e_nn, number)
        e_diff = e_ref_p_atom - e_nn_p_atom
        e_rmse.append(np.sqrt(1/np.size(e_diff)*np.sum(np.square(e_diff))))

        df_forces = pd.read_csv(str(i) + '_forces.dat', delim_whitespace=True)
        fx_error = abs(np.asarray(df_forces["fx_nn"]) -  np.asarray(df_forces["fx_ref"]))
        fy_error = abs(np.asarray(df_forces["fy_nn"]) -  np.asarray(df_forces["fy_ref"]))
        fz_error = abs(np.asarray(df_forces["fz_nn"]) -  np.asarray(df_forces["fz_ref"]))
        force_rmse_value = np.sqrt(
            1/(3*np.size(fx_error)) * (
                np.sum(np.square(fx_error)) +
                np.sum(np.square(fy_error)) +
                np.sum(np.square(fz_error))
            )
        )
        f_norm_nn = np.sqrt(
            np.square(np.asarray(df_forces["fx_nn"])) +
            np.square(np.asarray(df_forces["fy_nn"])) +
            np.square(np.asarray(df_forces["fz_nn"]))
        )
        f_norm_ref = np.sqrt(
            np.square(np.asarray(df_forces["fx_ref"])) +
            np.square(np.asarray(df_forces["fy_ref"])) +
            np.square(np.asarray(df_forces["fz_ref"]))
        )
        force_norm_rmse_value = np.sqrt(1/np.size(f_norm_nn)*np.sum(np.square(f_norm_nn - f_norm_ref)))

        force_rmse.append(force_rmse_value)
        force_norm_rmse.append(force_norm_rmse_value)

    return steps, e_rmse, force_rmse, force_norm_rmse

#Extract data
steps_1, e_1, f_1, f_1_norm = collect_validation(val_interv, '/train1')
steps_2, e_2, f_2, f_2_norm = collect_validation(val_interv, '/train2')
steps_3, e_3, f_3, f_3_norm = collect_validation(val_interv, '/train3')
steps_4, e_4, f_4, f_4_norm = collect_validation(val_interv, '/train4')
steps_5, e_5, f_5, f_5_norm = collect_validation(val_interv, '/train5')

#Find the best performing neural network
print("The NN_1 with min atomic E RMSE is after " + str(steps_1[e_1.index(min(e_1))]) + " steps at " + str(np.round(1000*(np.min(e_1)), 3)) + " meV with atomic F RMSE of " + str(np.round(1000*f_1[e_1.index(min(e_1))], 2)) + " meV/Å")
print("The NN_1 with min atomic E RMSE is after " + str(steps_2[e_2.index(min(e_2))]) + " steps at " + str(np.round(1000*(np.min(e_2)), 3)) + " meV with atomic F RMSE of " + str(np.round(1000*f_2[e_2.index(min(e_2))], 5)) + " meV/Å")
print("The NN_1 with min atomic E RMSE is after " + str(steps_3[e_3.index(min(e_3))]) + " steps at " + str(np.round(1000*(np.min(e_3)), 3)) + " meV with atomic F RMSE of " + str(np.round(1000*f_3[e_3.index(min(e_3))], 5)) + " meV/Å")
print("The NN_1 with min atomic E RMSE is after " + str(steps_4[e_4.index(min(e_4))]) + " steps at " + str(np.round(1000*(np.min(e_4)), 3)) + " meV with atomic F RMSE of " + str(np.round(1000*f_4[e_4.index(min(e_4))], 5)) + " meV/Å")
print("The NN_1 with min atomic E RMSE is after " + str(steps_5[e_5.index(min(e_5))]) + " steps at " + str(np.round(1000*(np.min(e_5)), 3)) + " meV with atomic F RMSE of " + str(np.round(1000*f_5[e_5.index(min(e_5))], 5)) + " meV/Å")
print(" ")
print("The NN_1 with min atomic F RMSE is after " + str(steps_1[f_1.index(min(f_1))]) + " steps at " + str(round(1000*min(f_1), 2)) + " meV/Å")
print("The NN_2 with min atomic F RMSE is after " + str(steps_2[f_2.index(min(f_2))]) + " steps at " + str(round(1000*min(f_2), 2)) + " meV/Å")
print("The NN_3 with min atomic F RMSE is after " + str(steps_3[f_3.index(min(f_3))]) + " steps at " + str(round(1000*min(f_3), 2)) + " meV/Å")
print("The NN_4 with min atomic F RMSE is after " + str(steps_4[f_4.index(min(f_4))]) + " steps at " + str(round(1000*min(f_4), 2)) + " meV/Å")
print("The NN_5 with min atomic F RMSE is after " + str(steps_5[f_5.index(min(f_5))]) + " steps at " + str(round(1000*min(f_5), 2)) + " meV/Å")
