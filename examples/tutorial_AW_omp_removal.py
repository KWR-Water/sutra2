# set working directory to main directory of package
import codecs
import imp

from matplotlib.pyplot import ylabel
from set_cwd_to_project_root import project_root
from pathlib import Path

from pandas import RangeIndex, read_csv
import pandas as pd
from pandas._testing import assert_frame_equal
import numpy as np
import datetime as dt
import os
import ast  # abstract syntax trees
import sys
import copy
# path = os.getcwd()  # path of working directory
from pathlib import Path

from zmq import zmq_version_info

# Plotting modules
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors

# sutra2 modules
import sutra2.Analytical_Well as AW
import sutra2.ModPath_Well as mpw
import sutra2.Transport_Removal as TR
