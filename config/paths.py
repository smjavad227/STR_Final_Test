import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.paths import *

import os
DIABETES_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIABETES_DIR = os.path.join(DIABETES_DIR, 'diabetes_extension')
DIABETES_FIGURES_DIR = os.path.join(DIABETES_DIR, 'results', 'figures')
DIABETES_TABLES_DIR = os.path.join(DIABETES_DIR, 'results', 'tables')
OUTPUT_DIR = os.path.join(DIABETES_DIR, 'output')
