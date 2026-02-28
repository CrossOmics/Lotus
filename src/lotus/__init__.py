# Eagerly import every submodule and auto-wrap their public functions
# with the lightweight ``tracked`` decorator so that lineage tracking
# is active for *any* Lotus function the user calls.
from lotus import (
    annotation_analysis,
    clustering,
    core_analysis,
    io,
    preprocessing,
    visualization,
    lineagetracker
)
# from lotus.lineagetracker import track_module

# for _mod in (
#     annotation_analysis,
#     clustering,
#     core_analysis,
#     io,
#     preprocessing,
#     visualization,
# ):
#     track_module(_mod)

# del _mod  # clean up module-level namespace

__all__ = [
    'annotation_analysis',
    'clustering',
    'core_analysis',
    'io',
    'lineagetracker',
    'preprocessing',
    'visualization',
]
