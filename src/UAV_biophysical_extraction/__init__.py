
try:
    from ._version import version as __version__  # written by setuptools_scm (optional)
except Exception:
    __version__ = "0.0.0"

# Export primary classes; update import paths if filenames differ.
from .processing.RF_processing import RFProcessor
from .processing.NN_processing import NNProcessor
from .processing.splits_processing import DatasetProcessor
from .utils.metrics import error_metrics

__all__ = ["RFProcessor", "NNProcessor", "DatasetProcessor", "error_metrics", "__version__"]