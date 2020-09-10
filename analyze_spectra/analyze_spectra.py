import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fname = sys.argv[-1]
spectra = pd.read_csv(fname, index_col = 0) # the first column
                                            # beomes the index. 
                                            # Just 0 indexed for 
                                            # data1.csv but this may 
                                            # not always be the case

snrs = [np.mean(spectra.iloc[:, i]) / np.std(spectra.iloc[:, i]) \
        for i in range(spectra.shape[1])]   # estimated SNR for 
                                            # each frequency bin

print("File {} has {} observations and {} frequency channels."\
        .format(fname, spectra.shape[0], spectra.shape[1]))
print("This script will make two plots of power (one in the time "\
        +"domain and one in the frequency domain) and one plot of "\
        +"SNR vs. frequency channel.")
print("You will be asked to input the frequency channel indices to "\
        +"plot in the time domain, and the time indices to plot in "\
        +"the frequency domain.")

def get_indices_from_input(indices_type, index_max):
    """ Prompt user to input a range or list of indices within the 
        range [0, index_max) to be returned as an object of type list.
    Args:
        indices_type (str): Phrase indicating the variable (axis) of 
                            the indices to be gathered.
        index_max (int): The maximum index (non-inclusive) to be 
                         accepted.
    Returns:
        list: The list of indices chosen by the user.
    """
    while(True):
        while(True):
            answer = input(("To input a range of {} indices, enter "\
                    +"'range', to input a list of {} indices enter "\
                    +"'list': ").format(indices_type, indices_type))
            if (answer in ['range', 'list']):
                break
            else:
                print("Invalid input, please try again.")
           
        indices = None  # variable to be returned
        if (answer == 'range'):
            while(True):
                indices = input(("Enter the range of {} indices in "\
                        +"start,stop,step format where start is "\
                        +"included and stop is not (the range must "\
                        +"be a subset of the range 0,{},1): ")\
                        .format(indices_type, index_max))
                indices = indices.split(',')
                if (len(indices) != 3):
                    print("Invalid format, please try again.")
                    continue
                try:
                    indices = [int(index) for index in indices]
                except(ValueError):
                    print("One or more of the values you entered "\
                            +"was not a number, please try again.")
                else:
                    indices = list(range(indices[0], indices[1], \
                            indices[2]))
                    break
        else: # answer == 'list'
            while(True):
                indices = input(("Enter the list of {} indices in "\
                        +"ascending order using format: 1,2,3,4,5,.."\
                        +". (all indices must be in list(range({})))"\
                        +": ").format(indices_type, index_max))
                indices = indices.split(',')
                try:
                    indices = [int(index) for index in indices]
                    break
                except(ValueError):
                    print("One or more of the values you entered "\
                            +"was not a number, please try again.")
        if (sorted(indices) != indices):
            print("The list of indices you've selected is not in "\
                    +"ascending order, please try again.")
            continue
        if (indices[0] < 0 or indices[-1] > index_max):
            print("Your selected range of indices is not within the "\
                    +"required bounds, please try again")
            continue
        print(("You've selected the {} indices given by the "\
                +"following list: ").format(indices_type))
        print(indices)
        answer = input("To confirm your selection, enter 'confirm', "\
                +"enter anything else to try again: ")
        if (answer == 'confirm'):
            return indices

freq_indices = get_indices_from_input("frequency channel", \
        spectra.shape[1])
time_indices = get_indices_from_input("time", spectra.shape[0])

fig = plt.figure(figsize = (12, 8))
ax1 = fig.add_subplot(2, 1, 1)  # time domain plot
ax2 = fig.add_subplot(2, 1, 2)  # frequency domain plot
time_axis = list(spectra.index)
freq_axis = list(range(spectra.shape[1]))
for index in freq_indices:
    ax1.plot(time_axis, spectra.iloc[:, index], label = str(index))
for index in time_indices:
    ax2.plot(freq_axis, spectra.iloc[index, :], label = str(index))
ax1.set_xlabel('Time [unknown]')
ax1.set_ylabel('Power [unknown]')
ax2.set_xlabel('Frequency Channel')
ax2.set_ylabel('Power [unknown]')
# Legend specifications chosen to fit every column (row) in data1.csv
# for the time domain (freq. domain) plot Thus, it may need to be 
# modified for a different file
ax1.legend(loc = 'upper right', title = 'freq. ch.', \
        fontsize = 'xx-small', ncol = 5, bbox_to_anchor = (1.35, 1))
ax2.legend(loc = 'upper right', title = 'time', \
        fontsize = 'xx-small', ncol = 5, bbox_to_anchor = (1.35, 1))
fig.tight_layout(pad = 3.0)
fig2 = plt.figure(figsize = (9, 9))
plt.plot(freq_axis, snrs)
plt.title('SNR vs. Frequency Channel')
plt.xlabel('Frequency Channel')
plt.ylabel('Estimated Signal-to-Noise Ratio')
plt.show()
            
