"""blade — python format library + gated stage runner for the
Wolfram-free Blade C pipeline (redg1/fitrel/dumppoints/ssolve/recmod/dynamicrr).

See README.md for scope and honest caveats.
"""

from .formats import (BladeFormatError, Database, DegreeInfo, DegreesFile,
                      EvalFile, EvalList, FitTable, Kinematics, Manifest,
                      OutSolution, PointsFile, RecCoeff, RecMono, RecPart,
                      RRRes, Scheme, SystemFile, TmpTemplate,
                      assemble_rational_functions, big_uint_primes, classify,
                      flags_size, prime_id, round_trip)
from .pipeline import (BladeGateError, PipelineConfig, StageResult, run_all,
                       run_dumppoints, run_dynamicrr, run_fitrel, run_recmod,
                       run_redg1, run_ssolve, write_system_and_evallists)
from .ratrec import reconstruct, ratfun_evaluator, poly_eval_mod
from .semibl import (Orderings, SemiBLJob, SemiBLTable,
                     gather_sample_points, run_permutation_tests,
                     sort_degrees_by_complexity)

__version__ = "1.0.0"
