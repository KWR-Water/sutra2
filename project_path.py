"""  Convenience file to easily get the full path of the project. """
import sys
from pathlib import Path
import os

module_path = Path(__file__).parent.parent
if module_path not in sys.path:
    sys.path.append(str(module_path))

