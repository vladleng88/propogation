import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def graphPlot(ettaRef, dSdX, potentialRef, withemRef, T_BOOM, P_BOOM, T_DOP, FF_DOP, ttt, FFF):
    fig, axs = plt.subplots(3, 2)

    axs[0, 0].plot(list(ettaRef.values()), list(dSdX.values()), linewidth=2.0, color='red')
    axs[0, 0].set_xlabel('etta')  # Add an x-label to the axes.
    axs[0, 0].set_ylabel('dS/dX')  # Add a y-label to the axes.

    axs[0, 1].plot(list(ettaRef.values()), list(potentialRef.values()), linewidth=2.0, color='red')
    axs[0, 1].set_xlabel('etta')  # Add an x-label to the axes.
    axs[0, 1].set_ylabel('Potential')  # Add a y-label to the axes.

    axs[2, 0].plot(list(ettaRef.values()), list(withemRef.values()), linewidth=2.0, color='red')
    axs[2, 0].set_xlabel('etta')  # Add an x-label to the axes.
    axs[2, 0].set_ylabel('Withem')  # Add a y-label to the axes.

    axs[1, 0].plot(list(T_DOP.values()), list(FF_DOP.values()), linewidth=2.0, color='red')
    axs[1, 0].plot(list(ttt.values()), list(FFF.values()), linewidth=2.0, color='blue')
    axs[1, 0].set_xlabel('t_wave')  # Add an x-label to the axes.
    axs[1, 0].set_ylabel('Ф_wave')  # Add a y-label to the axes.
    axs[1, 0].set_title("FFF_WAVE")  # Add a title to the axes.

    axs[1, 1].plot(list(T_BOOM.values()), list(P_BOOM.values()), linewidth=2.0, color='red')
    axs[1, 1].set_xlabel('t_boom')  # Add an x-label to the axes.
    axs[1, 1].set_ylabel('p_boom')  # Add a y-label to the axes.

    plt.show()