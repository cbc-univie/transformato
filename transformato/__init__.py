"""
transformato
Workflow to set up a relative free energy calculation of ligands with a common core scaffold
"""

# Add imports here
from .utils import load_config_yaml
from .system import SystemStructure
from .state import IntermediateStateFactory
from .mutate import ProposeMutationRoute

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
