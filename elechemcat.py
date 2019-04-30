from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
from scipy import stats
from math import log
import pandas as pd
import numpy as np
import sys
import os

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 11

font = FontProperties()
font.set_family('serif')
font.set_name('Times New Roman')
font.set_size('large')

err_names = ["", " ", None, "NaN"]


def ecc_analysis(inp_file):
    '''
    Main funciton caller. redirect input file to parser and perform analysis
    based on input of the file.

    **Parameters**
        inp_file        : *string*

    **Returns**

        None:           No return, but either will show or save plots from sub-funcs.
    '''
    main_dir = inp_file.split("/")[1]
    GLOBALS, XRD, CD, CV, ECSA, REACTION, KL = read_inp_file(inp_file)

    plots_path = GLOBALS["plots_path"][0]

    def process_CV_LSV(CV, CV_xy_range):
        '''
        Data processor for CV and LSV data including reactions Electrochemical
        Surface Area (ECSA), etc.

        **Parameters**
            CV              : *dict*
            CV_xy_range     : *list, float*

        **Returns**

            None:           No return, but either will show or save plots
        '''
        x_shifts_list = []
        y_divisors_list = []
        if len(CV["cv_files_list"]) == len(CV["x_shifts_list"]):
            for i in CV["x_shifts_list"]:
                if i == "REF_ELEC_RHE":
                    x_shifts_list.append(float(GLOBALS["REF_ELEC_RHE"][0]))
                else:
                    x_shifts_list.append(float(i))
            for i in CV["y_divisors_list"]:
                if i == "GEO_AREA_FACTOR":
                    y_divisors_list.append(float(GLOBALS["GEO_AREA_FACTOR"][0]))
                else:
                    y_divisors_list.append(float(i))
        elif (len(CV["cv_files_list"]) > 1 and len(CV["x_shifts_list"]) == 1):
            if CV["x_shifts_list"][0] == "REF_ELEC_RHE":
                x_shifts_list = [float(GLOBALS["REF_ELEC_RHE"][0])] * len(CV["cv_files_list"])
            else:
                x_shifts_list = [float(CV["x_shifts_list"][0])] * len(CV["cv_files_list"])
            if CV["y_divisors_list"][0] == "GEO_AREA_FACTOR":
                y_divisors_list = [float(GLOBALS["GEO_AREA_FACTOR"][0])] * len(CV["cv_files_list"])
            else:
                y_divisors_list = [float(CV["y_divisors_list"][0])] * len(CV["cv_files_list"])

        if len(CV) > 12:
            # ECSA
            cv_fig, cv_ax = generate_echem_plot(
                GLOBALS["working_path"][0],
                CV["cv_files_list"],
                CV["SAMPLE_NAME"],
                CV["x_column_header"],
                CV["y_column_header"],
                CV["control_column_header"],
                CV["control_type"][0],
                CV["outfile"],
                x_shifts_list,
                y_divisors_list,
                CV["num_of_scans"],
                CV["intg_from_pot"],
                CV["intg_to_pot"],
                CV["bg_intg_shift"],
                CV["ScanRate"],
                CV["StrChrg"])
        else:
            # CV and LSV
            cv_fig, cv_ax = generate_echem_plot(
                GLOBALS["working_path"][0],
                CV["cv_files_list"],
                CV["cv_legends_list"],
                CV["x_column_header"],
                CV["y_column_header"],
                CV["control_column_header"],
                CV["control_type"][0],
                CV["outfile"],
                x_shifts_list,
                y_divisors_list)

        plt.xlabel(CV["x_label"][0], fontproperties=font)
        plt.ylabel(CV["y_label"][0], fontproperties=font)

        x_axis = [float(i) for i in GLOBALS[CV_xy_range][:4]]
        plt.xlim(x_axis[0], x_axis[1])
        cv_ax.xaxis.set_major_locator(plt.MultipleLocator(x_axis[2]))
        cv_ax.xaxis.set_minor_locator(plt.MultipleLocator(x_axis[3]))

        y_axis = [float(i) for i in GLOBALS[CV_xy_range][4:]]
        plt.ylim(y_axis[0], y_axis[1])
        cv_ax.yaxis.set_major_locator(plt.MultipleLocator(y_axis[2]))
        cv_ax.yaxis.set_minor_locator(plt.MultipleLocator(y_axis[3]))

        cv_ax.grid(True, alpha=0.25, linestyle=':')

        if CV["control_type"][0] != "gradient":
            cv_ax.legend(frameon=False)

        if GLOBALS["auto_open_status"][0] == "True":
            plt.show()
        else:
            out_file = os.path.join(main_dir, plots_path, CV["outfile"][0] + '.png')
            plt.savefig(out_file, format='png', dpi=1000, bbox_inches='tight')

    def generate_echem_plot(
        source_dir,
        cv_files_list,
        cv_legends_list,
        x_column_header,
        y_column_header,
        control_column_header,
        control_type,
        out_label,
        x_shift=[0],
        y_divisor=[1],
        num_of_scans=0,
        intg_from_pot=0,
        intg_to_pot=0,
        bg_intg_shift=0,
        ScanRate=0,
        StrChrg=0,
    ):
        '''
        General plot processor and generator

        **Parameters**
            source_dir              : *string*
            cv_files_list           : *list, string*
            cv_legends_list         : *list, string*
            x_column_header         : *string*
            y_column_header         : *string*
            control_column_header   : *string*
            control_type            : *string*
            out_label               : *string*
            x_shift=[0]             : *optional*, *list, float*
            y_divisor=[1]           : *optional*, *list, float*
            num_of_scans=0          : *optional*, *float*
            intg_from_pot=0         : *optional*, *float*
            intg_to_pot=0           : *optional*, *float*
            bg_intg_shift=0         : *optional*, *float*
            ScanRate=0              : *optional*, *float*
            StrChrg=0               : *optional*, *float*

        **Returns**

            fig, ax: matplotlib figure and axes
        '''
        fig, ax = plt.subplots()
        cv_main_colors = ['r', 'b', 'k', 'g', 'm', 'c']
        i = 0
        for cv_i in cv_files_list:
            if cv_i in err_names:
                i = i + 1
                continue
            try:
                cv_i_data = pd.read_csv(os.path.join(main_dir, source_dir, cv_i), delimiter=";")
                cv_x = cv_i_data[x_column_header].values[:, 0] + x_shift[i]
                cv_y = cv_i_data[y_column_header].values[:, 0] / y_divisor[i]
                cv_c = cv_i_data[control_column_header].values[:, 0]
                cv_cs = set(cv_c.tolist())
                cv_cs_count = len(cv_cs)
                df = pd.DataFrame({
                    'y': cv_y,
                    'x': cv_x,
                    'control': cv_c})
                for control, control_df in df.groupby('control'):
                    if control_type == "gradient":
                        if control == 1:
                            continue
                        opc = (float(control + i * cv_cs_count) / float(len(cv_files_list) * cv_cs_count + 1))
                        m_size = 0.5
                        m_color = "darkblue"
                        m_label = None
                        plt.plot(control_df.x, control_df.y, '-', alpha=opc, linewidth=m_size, label=m_label, color=m_color)
                    elif control_type == "integration":
                        CO_from_intg = abs(cv_i_data.loc[cv_i_data['Scan'] == 1, 'Potential applied (V)'] - (float(intg_from_pot[0]) - x_shift[i])).idxmin()
                        CO_to_intg = abs(cv_i_data.loc[cv_i_data['Scan'] == 1, 'Potential applied (V)'] - (float(intg_to_pot[0]) - x_shift[i])).idxmin()
                        CO_pots = cv_x[CO_from_intg:CO_to_intg]
                        CO_currs = cv_y[CO_from_intg:CO_to_intg]
                        CO_baseline_values = np.zeros(len(CO_currs))
                        CO_currs = CO_currs.tolist()
                        CO_baseline_values = CO_baseline_values.tolist()

                        BG_from_intg = CO_from_intg + int(bg_intg_shift[0])
                        BG_to_intg = CO_to_intg + int(bg_intg_shift[0])
                        BG_pots = cv_x[BG_from_intg:BG_to_intg]
                        BG_currs = cv_y[BG_from_intg:BG_to_intg]
                        BG_baseline_values = np.zeros(len(BG_currs))
                        BG_currs = BG_currs.tolist()
                        BG_baseline_values = BG_baseline_values.tolist()

                        CO_peak = np.trapz(CO_currs, CO_pots) - np.trapz(CO_baseline_values, CO_pots)
                        BG_paek = np.trapz(BG_currs, BG_pots) - np.trapz(BG_baseline_values, BG_pots)

                        CO_Stripping = CO_peak - BG_paek
                        Charge = CO_Stripping / float(ScanRate[0])
                        ecsa_calc = Charge * 1000000 / float(StrChrg[0])
                        RFactor = ecsa_calc / float(GLOBALS["GEO_AREA"][0])

                        plt.plot(cv_x, cv_y, '-', alpha=1, linewidth=1, label='CV', color='k')
                        plt.fill_between(CO_pots, BG_currs, CO_currs, alpha=0.25, linewidth=0.7, label='CO stripping', color='g')
                        plt.text(0.3, 0.0004, 'ECSA = %.1f $\mathregular{cm^2}$\nRoughness Factor = %.1f' % (ecsa_calc, RFactor), fontsize=10)
                        break
                    elif (control_type == "linear" or int(control_type) == control):
                        opc = 1
                        m_size = 1
                        m_color = cv_main_colors[i]
                        m_label = cv_legends_list[i]
                        plt.plot(control_df.x, control_df.y, '-', alpha=opc, linewidth=m_size, label=m_label, color=m_color)
                    else:
                        continue
                i = i + 1
            except IOError:
                print("file: \"" + str(cv_i) + "\" does not exist.")
        return fig, ax

    def process_generate_KL(KL):
        '''
        Data processor and plotter for Koutecky-Levich eqaution. The function prcesses
        steady state electrochemical measurements to generate KL plot, Reaction number
        (number of electrons involved in the reaction), and Tafal plot

        **Parameters**

            KL:         *dict*
                        KL dictionary provided by the parser

        **Returns**

            None:       No return, but either will show or save plots
        '''
        KL_traces = []
        KL_main_colors = ['r', 'b', 'k', 'g', 'm', 'c']
        fig, ax = plt.subplots()

        BG_data = pd.read_csv(os.path.join(main_dir, KL["source_path"][0], "BG.txt"), delimiter=";")
        BG_curr_values = np.asarray(BG_data['WE(1).Current (A)']) / float(GLOBALS["GEO_AREA"][0]) * 1000

        for i in range(int(KL["num_rpms"][0])):
            file_name_i = "LSV(" + str(i + 1) + ").txt"
            data_i = pd.read_csv(os.path.join(main_dir, KL["source_path"][0], file_name_i), delimiter=";")
            rotation_i = np.asarray(data_i[KL["rotation_header"][0]])[0] * 1000
            pot_values_i = np.asarray(data_i[KL["LSV_x_column_header"][0]]) + float(GLOBALS["REF_ELEC_RHE"][0])
            curr_values_i = np.asarray(data_i[KL["y_column_header"][0]]) / float(GLOBALS["GEO_AREA"][0]) * 1000

            LSV_pots = pot_values_i
            LSV_currs = curr_values_i - BG_curr_values
            KL_pots = []
            KL_currs = []

            for j in range(60):
                file_name_j = "KL(" + str((i * 60) + j + 1) + ").txt"
                data_j = pd.read_csv(os.path.join(main_dir, KL["source_path"][0], file_name_j), delimiter=";")
                pot_values_j = np.asarray(data_j[KL["KL_x_column_header"][0]])
                curr_values_j = np.asarray(data_j[KL["y_column_header"][0]])
                curr_avg_j = np.mean(curr_values_j[:]) / float(GLOBALS["GEO_AREA"][0]) * 1000
                KL_pots.append(pot_values_j[0] + float(GLOBALS["REF_ELEC_RHE"][0]))
                KL_currs.append(curr_avg_j)

            KL_traces.append([rotation_i, KL_pots, KL_currs])
            plt.plot(LSV_pots, LSV_currs, '-', alpha=1, linewidth=0.5, label='LSV', color=KL_main_colors[i])
            plt.plot(KL_pots, KL_currs, 'o', alpha=1, markersize=3, label='KL', color=KL_main_colors[i])

        plt.xlabel(CV["x_label"][0], fontproperties=font)
        plt.ylabel(CV["y_label"][0], fontproperties=font)

        x_axis = [float(i) for i in GLOBALS["KL_plot_xy_range"][:4]]
        plt.xlim(x_axis[0], x_axis[1])
        ax.xaxis.set_major_locator(plt.MultipleLocator(x_axis[2]))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(x_axis[3]))

        y_axis = [float(i) for i in GLOBALS["KL_plot_xy_range"][4:]]
        plt.ylim(y_axis[0], y_axis[1])
        ax.yaxis.set_major_locator(plt.MultipleLocator(y_axis[2]))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(y_axis[3]))

        ax.grid(True, alpha=0.25, linestyle=':')

        if GLOBALS["auto_open_status"][0] == "True":
            plt.show()
        else:
            out_file = os.path.join(main_dir, plots_path, KL["outfile"][0] + '-SSvsLSV.png')
            plt.savefig(out_file, format='png', dpi=1000, bbox_inches='tight')

        ###############################################################################################
        ###############################################################################################

        fig, ax = plt.subplots()

        rpms = [KL_traces[j][0] for j in range(5)]
        KL_a_rotations = [1.0 / (i ** (0.5)) for i in rpms]
        KL_n = []

        const0 = 0.620  # A/cm2/(rad/s)^1/2
        const = const0 * (np.pi / 30.0) ** (1.0 / 2.0)  # A / cm2 / (rpm) ^ 1 / 3
        Fc = 96485.33  # A s / mol
        Co = 1.21 / 1000000  # mol / cm3
        Do = 1.93 / 100000  # cm2 / s
        v = 1.09 / 100  # cm2 / s
        Bo = const * Fc * Co * ((Do) ** (2.0 / 3.0)) * ((v) ** (-1.0 / 6.0))  # A/cm2/(rpm)^1/2

        for i in range(60):
            currs = []
            for j in range(5):
                currs.append(-1 / (KL_traces[j][2][i] / 1000))

            slope, intercept, r_value, p_value, std_err = stats.linregress(KL_a_rotations, currs)

            B = float(1.0 / float(slope))  # B
            jK = float(1.0 / intercept)  # jK
            n = float(1.0 / (slope * Bo))  # n
            k = jK / n / Fc / Co  # k

            KL_n.append(n)

        plt.plot(KL_pots, KL_n, 'o', alpha=1, markersize=3.5, label='Reaction Number', color='b')

        plt.xlabel(CV["x_label"][0], fontproperties=font)
        plt.ylabel("Reaction Number, n", fontproperties=font)

        x_axis = [float(i) for i in GLOBALS["KL_plot_xy_range"][:4]]
        plt.xlim(x_axis[0], x_axis[1])
        ax.xaxis.set_major_locator(plt.MultipleLocator(x_axis[2]))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(x_axis[3]))

        ax.grid(True, alpha=0.25, linestyle=':')

        if GLOBALS["auto_open_status"][0] == "True":
            plt.show()
        else:
            out_file = os.path.join(main_dir, plots_path, KL["outfile"][0] + '-RxnNum.png')
            plt.savefig(out_file, format='png', dpi=1000, bbox_inches='tight')

        ###############################################################################################
        ###############################################################################################

        fig, ax = plt.subplots()

        for i in range(5):
            if int(KL_traces[i][0]) != 1600:
                continue
            logcurr = [log(abs(x), 10) for x in KL_traces[i][2]]
            overpot = [(float(KL["RXN_pot"][0]) - pot) * 1000 for pot in KL_traces[i][1]]

        op_1 = [y for y in [y for y in overpot if y <= 280] if y >= 229]
        op_2 = [s for s in [s for s in overpot if s <= 500] if s > 380]

        lc_1 = logcurr[:]
        lc_2 = logcurr[:]

        lc_1 = lc_1[overpot.index(op_1[0]):overpot.index(op_1[-1]) + 1]
        lc_2 = lc_2[overpot.index(op_2[0]):overpot.index(op_2[-1]) + 1]

        # calculating  and plotting trendlines
        z1 = np.polyfit(lc_1, op_1, 1)
        p1 = np.poly1d(z1)
        x1 = np.linspace(lc_1[-1] + 1.5, lc_1[0] - 1.5, 1000)
        plt.plot(lc_1, p1(lc_1), '-', alpha=1, linewidth=1, color='r')
        plt.plot(x1, p1(x1), ':', alpha=0.5, linewidth=1, color='r')
        l1 = np.array((lc_1[-1] * 0.98, p1(lc_1[-1]) * 0.98))
        plt.text(l1[0], l1[1], "%.2f mV/dec" % z1[0], fontsize=11, rotation=-15, color='r')

        z2 = np.polyfit(lc_2, op_2, 1)
        p2 = np.poly1d(z2)
        x2 = np.linspace(lc_2[-1], lc_2[0], 1000)
        plt.plot(lc_2, p2(lc_2), '-', alpha=1, linewidth=1, color='g')
        plt.plot(x2, p2(x2), ':', alpha=0.5, linewidth=1, color='g')
        l2 = np.array((lc_2[-1] * 1.05, p2(lc_2[-1])))
        plt.text(l2[0], l2[1], "%.2f mV/dec" % z2[0], fontsize=11, rotation=-88, color='g')

        plt.plot(logcurr, overpot, 'o', alpha=1, markersize=4, label='Steady State Reaction', color='b')
        plt.xlabel("log[io]", fontproperties=font)
        plt.ylabel("Overpotential (mV)", fontproperties=font)

        plt.xlim(-2.0, 1.5)
        plt.ylim(525, 200)

        ax.grid(True, alpha=0.25, linestyle=':')

        if GLOBALS["auto_open_status"][0] == "True":
            plt.show()
        else:
            out_file = os.path.join(main_dir, plots_path, KL["outfile"][0] + '-Tafal.png')
            plt.savefig(out_file, format='png', dpi=1000, bbox_inches='tight')

    if len(CV) > 0:
        print("Processing and Plotting CV ....")
        process_CV_LSV(CV, "CV_plot_xy_range")

    if len(CD) > 0:
        print("Processing and Plotting CD ....")
        CD_i = {}
        for i in range(int(CD["NUMBER_OF_SAMPLES"][0])):
            for cd_label in CD:
                if (cd_label != "NUMBER_OF_SAMPLES" and cd_label != "CD_NUM"):
                    CD_i[cd_label] = CD[cd_label][i]
            process_CV_LSV(CD_i, "CD_plot_xy_range")

    if len(REACTION) > 0:
        print("Processing and Plotting Reaction ....")
        process_CV_LSV(REACTION, "reaction_plot_xy_range")

    if len(ECSA) > 0:
        print("Processing and Plotting ECSA ....")
        ECSA_i = {}
        for i in range(int(ECSA["NUMBER_OF_SAMPLES"][0])):
            for ecsa_label in ECSA:
                if (ecsa_label != "NUMBER_OF_SAMPLES" and ecsa_label != "ECSA_NUM"):
                    ECSA_i[ecsa_label] = ECSA[ecsa_label][i]
            process_CV_LSV(ECSA_i, "ECSA_plot_xy_range")

    if len(KL) > 0:
        print("Processing and Plotting KL ....")
        process_generate_KL(KL)

    print("Done.")


def read_inp_file(inp_file):
    '''
    Input files loader and parser. The proper and a full exmaple of the input file
        can be found in the Example folder

    **Parameters**

        inp_file:   *string*
                    path to input file

    **Returns**

        valid:  *dict*, *dict*, *dict*, *dict*, *dict*, *dict*, *dict*
                dictonaries of the infomration parsed from the input file in
                GLOBALS, XRD, CD, CV, ECSA, REACTION, KL
    '''
    GLOBALS, XRD, CD, CV, ECSA, REACTION, KL = {}, {}, {}, {}, {}, {}, {}

    with open(inp_file) as input_file:
        GLOBALS_switch = False
        XRD_switch = False
        CD_switch = False
        CD_iter_switch = False
        CV_switch = False
        ECSA_switch = False
        ECSA_iter_switch = False
        REACTION_switch = False
        KL_switch = False
        for line in input_file:
            if (line.strip() == "#" or line.strip() == ""):
                continue
            elif line.strip() == "&GLOBALS":
                GLOBALS_switch = True
            elif line.strip() == "&XRD":
                XRD_switch = True
            elif line.strip() == "&CD":
                CD_switch = True
            elif line.strip() == "&CV":
                CV_switch = True
            elif line.strip() == "&ECSA":
                ECSA_switch = True
            elif line.strip() == "&REACTION":
                REACTION_switch = True
            elif line.strip() == "&KL":
                KL_switch = True
            elif line.strip() == "///":
                GLOBALS_switch = False
                XRD_switch = False
                CD_switch = False
                CV_switch = False
                ECSA_switch = False
                REACTION_switch = False
                KL_switch = False
            elif GLOBALS_switch:
                GLOBALS[line.strip().split()[0]] = [i.split(", ") for i in line.strip().split("= ")[1:]][0]
            elif XRD_switch:
                XRD[line.strip().split()[0]] = [i.split(", ") for i in line.strip().split("= ")[1:]][0]
            elif CD_switch:
                if (str(line.strip().split()[0]) == "CD_NUM" and int(line.strip().split("= ")[1:][0]) > 1):
                    CD_iter_switch = True
                    continue
                if CD_iter_switch:
                    CD[line.strip().split()[0]] = [CD[line.strip().split()[0]], [i.split(", ") for i in line.strip().split("= ")[1:]][0]]
                else:
                    CD[line.strip().split()[0]] = [i.split(", ") for i in line.strip().split("= ")[1:]][0]
            elif CV_switch:
                CV[line.strip().split()[0]] = [i.split(", ") for i in line.strip().split("= ")[1:]][0]
            elif ECSA_switch:
                if (str(line.strip().split()[0]) == "ECSA_NUM" and int(line.strip().split("= ")[1:][0]) > 1):
                    ECSA_iter_switch = True
                    continue
                if ECSA_iter_switch:
                    ECSA[line.strip().split()[0]] = [ECSA[line.strip().split()[0]], [i.split(", ") for i in line.strip().split("= ")[1:]][0]]
                else:
                    ECSA[line.strip().split()[0]] = [i.split(", ") for i in line.strip().split("= ")[1:]][0]
            elif REACTION_switch:
                REACTION[line.strip().split()[0]] = [i.split(", ") for i in line.strip().split("= ")[1:]][0]
            elif KL_switch:
                KL[line.strip().split()[0]] = [i.split(", ") for i in line.strip().split("= ")[1:]][0]
    return GLOBALS, XRD, CD, CV, ECSA, REACTION, KL


systemArg = sys.argv
if len(systemArg) < 2:
    print("ERROR: no input file is specified. ex: (elechecat.py input.inp) ")
elif len(systemArg) > 2:
    print("CAUTION: Only one input file will be processed: %s" % systemArg[1])
    ecc_analysis(systemArg[1])
else:
    ecc_analysis(systemArg[1])

if __name__ == "__main__":
    # ecc_analysis(os.path.join("Exmaple/input.inp"))
    pass
