import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.constants import labels_Earth, Z_Earth, A_Earth
from utils.plot_utils import plot_Ye, save_fig

folder_plots = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "plots"
)

print(folder_plots)
fig4 = plot_Ye(Z_Earth, A_Earth, labels_Earth)
save_fig(fig4, folder_plots, "/Figure4.png")

