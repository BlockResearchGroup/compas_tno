import os

__author__ = ["Ricardo Maia Avelino"]
__copyright__ = "Block Research Group - ETH Zurich"
__license__ = "MIT License"
__email__ = "maia@arch.ethz.ch"
__version__ = "0.3.0"


HERE = os.path.dirname(__file__)

HOME = os.path.abspath(os.path.join(HERE, "../../"))
DATA = os.path.abspath(os.path.join(HOME, "data"))
DOCS = os.path.abspath(os.path.join(HOME, "docs"))
TEMP = os.path.abspath(os.path.join(HOME, "temp"))


__all__ = ["HOME", "DATA", "DOCS", "TEMP"]
