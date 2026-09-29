import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import synthetic_methods as synthetic  # noqa: E402

# The combination package needs the datasets folder when it is imported
if 'INTOGEN_DATASETS' not in os.environ:
    folder = tempfile.mkdtemp(prefix='intogen_datasets_')
    synthetic.make_datasets(folder)
    os.environ['INTOGEN_DATASETS'] = folder
