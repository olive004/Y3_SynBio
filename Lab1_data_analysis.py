# Synbio foundations Lab 1 Data Analysis: Promoter Strength


import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from collections.abc import Iterable        # for list fxns
import seaborn as sns                       # matplotlib sux for heatmaps
# matplotlib inline
plt.style.use('Solarize_Light2')
# plt.style.use('seaborn-pastel')


def flatten(lis):
     for item in lis:
         if isinstance(item, Iterable) and not isinstance(item, str):
             for x in flatten(item):
                 yield x
         else:        
             yield item


def choose_cols(data, label_specification):         # https://stackoverflow.com/questions/17485747/how-to-convert-a-nested-list-into-a-one-dimensional-list-in-python
    """Only include columns in data that are values 
    in the dictionary """
    cols = flatten(list(label_specification.values()))
    return data[cols]

def get_red_mean(data, label_specification):
    """ Get mean between experimentally redundant columns that as values
    correspond to a key in the label_specification """
    if isinstance(data, pd.Series):
        data_m = pd.Series(index=label_specification.keys())
    else:
        data_m = pd.DataFrame(index=data.index, columns=label_specification.keys())

    for key, val in label_specification.items():
        if isinstance(data, pd.DataFrame):      # Differentiate between Series, which only have 1 axis not 2
            if val[0] in data.columns:
                data_m[key] = data[val].mean(axis=1, skipna=True)   
            else: continue
        else: 
            if val[0] in data.index.tolist():
                data_m[key] = data[val].mean(skipna=True)
            else: continue
    return data_m

def get_red_std(data, label_specification):
    """ Get std between experimentally redundant columns that as values
    correspond to a key in the label_specification """
    if isinstance(data, pd.Series):
        data_std = pd.Series(index=label_specification.keys(), dtype=np.dtype('Float64') )
    else:
        data_std = pd.DataFrame(index=data.index, columns=label_specification.keys())

    for key, val in label_specification.items():
        if isinstance(data, pd.DataFrame):      # Differentiate between Series, which only have 1 axis not 2
            if val[0] in data.columns:
                data_std[key] = data[val].std(axis=1, skipna=True)
            else: continue
        else: 
            if val[0] in data.index.tolist():
                data_std[key] = data[val].std(skipna=True)
            else: continue
    return data_std



def clean_negs(df):
    df[df == np.nan] = 0
    df[(df == np.inf) | df == - np.inf] = 0
    df[df < 0] = 0
    return df



def expand_labels(group_labels, label_specification):
    """ Some labels correspond to a group of sub-labels. 
    Takes in array of group labels and uses the dict 
    label_specification to get corresponding sub-labels. """
    sub_labels = []
    label_ind = 0
    for key, val in label_specification.items():
        if key in group_labels:
            for i in val:
                sub_labels.append(i)
        label_ind = label_ind +1
    return sub_labels


def print_summary(series_summary, intro_sentence):
    print(intro_sentence)
    for i in range(0,len(series_summary)):
        print(str(series_summary.index[i]), str(series_summary[i]))
    





if __name__=="__main__":

    # Calibration info from output: ABS (0.5/0.58) = OD
    # TRUE	Reader
    # OD600	ABS600
    # 0.500	0.580
    abs2od_conv_factor = (0.5/0.58) 

    # Load in data
    filename_ABS = 'ABS600.csv'
    filename_GFP = 'GFP.csv'
    promoters = ['A1', 'B2', 'C2', 'D2', 'D2+', 'S','N']
    promoter_cols = {'MINUTES': ['MINUTES'], 'LB': ['C1','C2','C3'], 'A1':['C4','C5','C6'], 'B2': ['C7','C8','C9'], 'C2': ['C10', 'C11'],'D2':['D1','D2','D3'], 'D2+':['D4','D6'], 'S': ['D7','D8','D9']	, 'N':['D10','D11','D12']}
    data_ABS = pd.read_csv(filename_ABS)
    data_GFP = pd.read_csv(filename_GFP)
    data_ABS = choose_cols(data_ABS, promoter_cols)
    full_time = data_ABS['MINUTES']
    data_GFP = choose_cols(data_GFP, promoter_cols)     # only using columns in promoter_cols dict  

    # convert ABS to OD
    data_ABS = data_ABS.mul(abs2od_conv_factor)
    data_ABS['MINUTES'] = full_time
    
    # Data processing
    # Get mean out of 3 samples
    data_ABS_m = get_red_mean(data_ABS, promoter_cols)
    data_ABS_m['MINUTES'] = data_ABS['MINUTES']

    # CORRECT background ABS600 (values)
    ABS_correction_sample = data_ABS_m['LB']
    # subtract
    data_ABS_corr_m = data_ABS_m.subtract(ABS_correction_sample, axis='rows')
    # >0
    data_ABS_corr_m[data_ABS_corr_m <= 0] = 0
    # Minutes
    data_ABS_corr_m['MINUTES'] = data_ABS['MINUTES']  


    ########################
    # Plot curves of ABS600   
    # savefig(fname, dpi=None, facecolor='w', edgecolor='w',
    #     orientation='portrait', papertype=None, format=None,
    #     transparent=False, bbox_inches=None, pad_inches=0.1,
    #     frameon=None, metadata=None)
    time_ABS = data_ABS['MINUTES']  
    plt.figure()
    for i in range(1, (len(data_ABS_corr_m.columns))):
        plt.plot(time_ABS, data_ABS_corr_m.iloc[:,i])
        curr_col = data_ABS_corr_m.columns[i]

    plt.legend(data_ABS_corr_m.columns[1:])
    plt.title('ABS growth curve')
    plt.xlabel('time (min)')
    plt.ylabel('OD Absorbance (AU)')
    plt.savefig('lab1_ABS.eps', format='eps')
    plt.show()
    ########################



    # LOG PHASE: chosen interval minutes            # Some samples had weird growth values; 120-180 and 480-540 (C2,D2+)
    log_phase_range = [120,180]         # Minutes of max growth (has to be accurate)
    log_phase_range_bad = [480, 540]    # 645] was used before lol

    # DF indices corresponding to growth range
    # --> Gets iloc index of dataframe by making boolean df and using argwhere() to get the indices of truths
    log_ph_ok_bool = (data_ABS_corr_m['MINUTES']==log_phase_range[0]) | (data_ABS_corr_m['MINUTES']==log_phase_range[1])
    log_ph_bad_bool = (data_ABS_corr_m['MINUTES']==log_phase_range_bad[0]) | (data_ABS_corr_m['MINUTES']==log_phase_range_bad[1])
    log_ph_ok_ind = np.argwhere(log_ph_ok_bool.values).flatten()
    log_ph_bad_ind = np.argwhere(log_ph_bad_bool.values).flatten()
    
    # Columns good/bad
    ok_cols = ['LB','A1','B2','D2','S','N']    # Have to find manually from experimental data
    bad_cols = ['LB','C2','D2+']
    # Wells that correspond to each sample label
    ok_cols_expanded = expand_labels(ok_cols, promoter_cols) 

    # Get ABS log phase samples (i am sorry for this)
    #OLD
    ABS_log_phase = data_ABS_corr_m.iloc[min(log_ph_bad_ind[0],log_ph_ok_ind[0]):(max(log_ph_bad_ind[1], log_ph_ok_ind[1])+1), :]
    ABS_log_phase_ok = data_ABS_corr_m.iloc[log_ph_ok_ind[0]:log_ph_ok_ind[1]+1, :]
    ABS_log_phase_bad = data_ABS_corr_m.iloc[log_ph_bad_ind[0]:log_ph_bad_ind[1]+1, :]
    #NEW
    ABS_log_phase = data_ABS_corr_m[ok_cols[:]].iloc[log_ph_ok_ind[0]:log_ph_ok_ind[1]+1, :]
    # Insert the lagged promoters with the other promoters using .values
    ABS_log_phase.insert(3, bad_cols[1], ABS_log_phase_bad[bad_cols[1]].values)     # C2 values
    ABS_log_phase.insert(5, bad_cols[2], ABS_log_phase_bad[bad_cols[2]].values)     # D2+ values

    # &&&&&&&&&&&&&&&&&&&&&&&&&&& PLOTTING ABS &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Plot ABS log phase
    plt.figure()
    time_log_ABS = ABS_log_phase_ok['MINUTES']
    for i in range(0, (len(ABS_log_phase.columns))):
        plt.plot(time_log_ABS, ABS_log_phase.iloc[:,i])
        curr_col = ABS_log_phase.columns[i]

    plt.legend(ABS_log_phase.columns[:])
    plt.title('ABS log phase all')
    plt.xlabel('time (min)')
    plt.ylabel('OD Absorbance (AU)')
    plt.savefig('lab1_logABS_all.eps', format='eps')
    plt.show()

    # Plot ABS log phase: GOOD PROMOTERS
    plt.figure()
    time_log_ABS = ABS_log_phase_ok['MINUTES']
    for i in ok_cols:
        plt.plot(time_log_ABS, ABS_log_phase_ok[i])
    plt.legend(ABS_log_phase_ok[ok_cols])
    plt.title('ABS log phase (good samples)')
    plt.xlabel('time (min)')
    plt.ylabel('OD Absorbance (AU)')
    plt.savefig('lab1_logABS_good.eps', format='eps')
    plt.show()

    # Plot ABS log phase: BAD PROMOTERS
    plt.figure()
    time_log_ABS = ABS_log_phase_bad['MINUTES']
    for i in bad_cols:
        plt.plot(time_log_ABS, ABS_log_phase_bad[i])
    plt.legend(ABS_log_phase_bad[bad_cols])
    plt.title('ABS log phase (lagged samples)')
    plt.xlabel('time (min)')
    plt.ylabel('OD Absorbance (AU)')
    plt.savefig('lab1_logABS_bad.eps', format='eps')
    plt.show()
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&










    ########################
    # 2. GFP
    # Correct for sample measuring bg noise
    GFP_correction_sample = data_GFP['C3'] 
    # Skip 'MINUTES' with 1: then add 'MINUTES' back
    data_GFP_corr = data_GFP.iloc[:, 1:].subtract(GFP_correction_sample, axis='rows')
    data_GFP_corr.insert(0, 'MINUTES', data_GFP['MINUTES'])

    data_GFP_corr = clean_negs(data_GFP_corr)

    # Select log phase
    log_ph_ok_bool = (data_GFP_corr['MINUTES']==log_phase_range[0]) | (data_GFP_corr['MINUTES']==log_phase_range[1])
    log_ph_bad_bool = (data_GFP_corr['MINUTES']==log_phase_range_bad[0]) | (data_GFP_corr['MINUTES']==log_phase_range_bad[1])
    log_ph_ok_ind = np.argwhere(log_ph_ok_bool.values).flatten()
    log_ph_bad_ind = np.argwhere(log_ph_bad_bool.values).flatten()

    # Columns for dGFP calculations
    ok_cols_gfp = ['LB','A1','B2','D2','S','N']    # Have to find manually from experimental data
    bad_cols_gfp = ['C2','D2+']
    # Wells that correspond to each sample label; there's 3 wells per promoter
    ok_cols_gfp = expand_labels(ok_cols_gfp, promoter_cols) 
    bad_cols_gfp = expand_labels(bad_cols_gfp, promoter_cols) 

    # Separate vectors for good and bad promoters since we're looking at different times for htem
    # GFP_log_phase = data_GFP_corr.iloc[min(log_ph_bad_ind[0],log_ph_ok_ind[0]):(max(log_ph_bad_ind[1], log_ph_ok_ind[1])+1), :]
    # Get good promoters first so the time access is correct, then add in bad promoters from different time interval
    GFP_log_phase = data_GFP_corr[ok_cols_gfp[:]].iloc[log_ph_ok_ind[0]:log_ph_ok_ind[1]+1, :]
    # GFP_log_phase_ok = data_GFP_corr[ok_cols_gfp[:]].iloc[log_ph_ok_ind[0]:log_ph_ok_ind[1]+1, :]
    GFP_log_phase_bad = data_GFP_corr[bad_cols_gfp].iloc[log_ph_bad_ind[0]:log_ph_bad_ind[1]+1, :]
    # Insert the lagged promoters with the other promoters using .values
    for i in range(0,len(promoter_cols["C2"])):      # C2
        pos_c2 = (3*len(promoter_cols["C2"]))        # pos of C2 in expanded columns
        GFP_log_phase.insert(pos_c2+i, bad_cols_gfp[i], GFP_log_phase_bad[bad_cols_gfp[i]].values) 
    for i in range(len(promoter_cols["C2"]),len(promoter_cols["C2"])+len(promoter_cols["D2+"])):      # D2+
        pos_d2p = (4*len(promoter_cols["D2+"]))
        GFP_log_phase.insert(pos_d2p+i, bad_cols_gfp[i], GFP_log_phase_bad[bad_cols_gfp[i]].values) 

    # Get dGFP (change in GFP) 
    # Again, do everything for good and bad samples
    dGFP = GFP_log_phase.iloc[-1, :] - GFP_log_phase.iloc[0, :]
    # dGFP_ok = GFP_log_phase_ok.iloc[-1, :] - GFP_log_phase_ok.iloc[0, :]
    # dGFP_bad = GFP_log_phase_bad.iloc[-1, :] - GFP_log_phase_bad.iloc[0, :]

    # Change in GFP per cell
    #idx
    half_ind = min(log_ph_bad_ind[0],log_ph_ok_ind[0]) + int(np.floor(len(GFP_log_phase)/2)-1)
    # half_ind_ok = log_ph_ok_ind[0] + int(np.floor(len(GFP_log_phase_ok)/2)-1)
    half_ind_bad = log_ph_bad_ind[0] + int(np.floor(len(GFP_log_phase_bad)/2)-1)
    #S
    dGFP_cell = dGFP[ok_cols_gfp] / data_ABS_corr_m['S'].iloc[half_ind]
    dGFP_cell_bad = dGFP[bad_cols_gfp] / data_ABS_corr_m['S'].iloc[half_ind_bad]
    dGFP_cell = dGFP_cell.append(dGFP_cell_bad)
    #old
    # dGFP_cell_ok = dGFP_ok / data_ABS_corr_m['S'].iloc[half_ind_ok]
    # dGFP_cell_bad = dGFP_bad / data_ABS_corr_m['S'].iloc[half_ind_bad]

    # clean up data
    dGFP_cell = clean_negs(dGFP_cell)
    # dGFP_cell_ok = clean_negs(dGFP_cell_ok)
    # dGFP_cell_bad = clean_negs(dGFP_cell_bad)


    # Get mean & std of dGFP per cell
    mean_dGFP_cell = get_red_mean(dGFP_cell, promoter_cols)
    std_dGFP_cell = get_red_std(dGFP_cell, promoter_cols)


    # Plot Bar chart of promoter stats
    plt.bar(promoters, mean_dGFP_cell[2:], yerr=std_dGFP_cell[2:], ecolor='green', capsize=10)
    plt.title('Mean dGFP by promoter')
    plt.xlabel('Promoter')
    plt.ylabel('dGFP (RFU)')
    plt.savefig('lab2_dGFP.eps', format='eps')
    plt.show()


    # RPU

    # Get Relative Promoter Units (RPU) by dividing by standard (S) promoter
    # Subtract measure from promoter N first
    rpu_dGFP_m = mean_dGFP_cell.subtract(mean_dGFP_cell['N'])
    rpu_dGFP_m = clean_negs(rpu_dGFP_m)
    rpu_dGFP_m = rpu_dGFP_m / rpu_dGFP_m['S']
    # For scaling std error's 
    rpu_scaling_factor = 1 / mean_dGFP_cell['S']

    # Standard Error 
    # is the standard deviation divided by the square-root 
    # of the number of data points (n). In this 
    # experiment you did n=3.
    n=3
    std_dGFP_cell = std_dGFP_cell.subtract(std_dGFP_cell['N']) 
    std_err_GFP = (std_dGFP_cell * rpu_scaling_factor) / (np.sqrt(n))

    # Plot Bar chart
    plt.bar(promoters, rpu_dGFP_m[2:], yerr=std_err_GFP[promoters], ecolor='red', capsize=10) #, yerr=std_dGFP_cell[2:], ecolor='green', capsize=10)
    plt.title('Relative promoter strengths')
    plt.xlabel('Promoter')
    plt.ylabel('RPU')
    plt.savefig('lab2_dGFP_rpu.eps', format='eps')
    plt.show()

    print_summary(mean_dGFP_cell[2:], "The mean dGFP strength in RFU for each promoter is \n")
    print_summary(rpu_dGFP_m, "The mean dGFP strength in RPU for each promoter is \n")


    

    # CLASS COMPARISON
    # data (first column ours rest unsorted)
    class_data_raw = {"Team 1": rpu_dGFP_m[promoters[0:5]].values,
        "Team 2": [np.nan,  0.511,      np.nan,     np.nan,     np.nan],
        "Team 3": [np.nan,  00.59716,   np.nan,     0.038395,   0.032766],
        "Team 4": [np.nan,  0.417,      0.0044,     0.017,      -0.0033],
        "Team 5": [np.nan,  np.nan,     -0.0035,    np.nan,     np.nan],
        "Team 10": [np.nan, 0.458,      np.nan,     0.004,      0.008],
        "Team 11": [2.53,   np.nan,     0.009,      np.nan,     np.nan],
        "Team 12": [2.018,  np.nan,     np.nan,     np.nan,     np.nan],
        "Team 13": [np.nan, np.nan,     -0.006,     0.021,      0.006],
        "Team 14": [2.688,  np.nan,     np.nan,     np.nan,     np.nan]}
    class_data = pd.DataFrame(class_data_raw, index=["A1", "B2", "C2", "D2", "D2+"], columns=["Team 1","Team 2","Team 3","Team 4","Team 5", "Team 10","Team 11","Team 12","Team 13","Team 14"])

    # stats
    mean_class = class_data.mean(axis=1)
    std_class = class_data.std(axis=1)

    # plot class
    plt.bar(mean_class.index, mean_class, yerr=std_class, ecolor='black', capsize=10, zorder=1) #, yerr=std_dGFP_cell[2:], ecolor='green', capsize=10)
    plt.scatter(mean_class.index, class_data_raw["Team 1"], c="red", zorder=2)
    plt.title('Comparing RPU to class')
    plt.xlabel('Promoter')
    plt.ylabel('RPU')
    plt.savefig('lab3_class_comp.eps', format='eps')
    plt.show()

    print_summary(mean_class, "The mean of each promoter for the class is \n")
    print_summary(std_class, "The standard deviation of each promoter for the class is \n")


    class_A2 = {"2": 0.874, "4": 0.7929, "5":0.9157,"6":0.607, "9":0.88}

    print("A2 mean: ", str(np.mean(list(class_A2.values()))))



    # # Get volume for the ng/ul concentration of promoter we got
    # # conv ng/ul to g/L
    # concentrations = [194.0,77.6,29.1,27.3]
    # concentrations = concentrations*(0.001)





    # CLASS 4: Heatmap of final results
    plt.style.use('Solarize_Light2')

    data_rfp_4 = pd.read_csv('Lab4_red.csv')
    data_gfp_4 = pd.read_csv('Lab4_green.csv')
    data_od_4 = pd.read_csv('Lab4_OD.csv')

    iptg = [0, 0.5,1,2,3,4,5,6]
    ahl = [0,	0.05,	0.1,	0.125,	0.15]

    # Seaborn automatically interprets these for annotation
    data_rfp_4.columns = ahl
    data_rfp_4.index = iptg
    data_gfp_4.columns = ahl
    data_gfp_4.index = iptg
    data_od_4.columns = ahl
    data_od_4.index = iptg

    data_rfp_4 = data_rfp_4.div(data_od_4)
    data_gfp_4 = data_gfp_4.div(data_od_4)

    # rfp
    fig, ax = plt.subplots()
    ax = sns.heatmap(data_rfp_4, annot=True, fmt=".6g")
    ax.set_title("Heatmap red fluorescence")
    ax.set_xlabel("IPTG (mM)")
    ax.set_ylabel("AHL (uM)")
    plt.savefig('lab4_rfp.eps', format='eps')
    plt.show()

    # gfp
    fig, ax = plt.subplots()
    ax = sns.heatmap(data_gfp_4, annot=True, fmt=".6g")
    ax.set_title("Heatmap green fluorescence")
    ax.set_xlabel("IPTG (mM)")
    ax.set_ylabel("AHL (uM)")
    plt.savefig('lab4_gfp.eps', format='eps')
    plt.show()

    # od
    fig, ax = plt.subplots()
    ax = sns.heatmap(data_od_4, annot=True, fmt=".4g")
    ax.set_title("Heatmap optical density (OD)")
    ax.set_xlabel("IPTG (mM)")
    ax.set_ylabel("AHL (uM)")
    plt.savefig('lab4_od.eps', format='eps')
    plt.show()

