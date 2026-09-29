"""Contain native solver stderr without changing exceptions or numerical inputs.

POSIX fd redirection is process-wide. Solver calls using this helper are
serialized; campaigns already evaluate sequentially. Unrelated threads must
not concurrently write diagnostics to fd 2 during a solve.
"""
from contextlib import contextmanager
import ctypes
import os
import sys
import tempfile
import threading

_LOCK = threading.RLock()


def _flush():
    sys.stderr.flush()
    ctypes.CDLL(None).fflush(None)


@contextmanager
def capture_native_stderr():
    """Capture actual fd 2 writes, restoring it even on BaseException."""
    captured = {'text': ''}
    with _LOCK, tempfile.TemporaryFile(mode='w+b') as stream:
        _flush()
        saved = os.dup(2)
        try:
            os.dup2(stream.fileno(), 2)
            yield captured
        finally:
            try: _flush()
            finally:
                os.dup2(saved, 2)
                os.close(saved)
            stream.seek(0)
            captured['text'] = stream.read().decode('utf-8', errors='replace')


def quiet_solve_ivp(*args, **kwargs):
    from scipy.integrate import solve_ivp
    captured = {}
    try:
        with capture_native_stderr() as captured:
            result = solve_ivp(*args, **kwargs)
    except BaseException as error:
        # Preserve the original exception object, type and message.
        if captured.get('text'): error.native_solver_stderr = captured['text']
        raise
    if captured['text']: result.native_solver_stderr = captured['text']
    return result
