'''Expression‑based objectives for Dada‑Engine research.

This module defines :class:`ExpressionObjective`, which allows users to specify
arbitrary arithmetic expressions in the study TOML.  The expression may use any
symbol present in the result dictionary returned by the thermodynamic evaluator
plus the derived ``ratio`` variable (micro‑tube length‑to‑diameter).  The caret
operator ``^`` is treated as exponentiation (mapped to ``**``).
'''  

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping, Any

from asteval import Interpreter  # safe sandboxed evaluator

from .geometry_helpers import exchanger_length_to_diameter


@dataclass(frozen=True, slots=True)
class ExpressionObjective:
    """Represent a user‑defined arithmetic objective.

    Parameters
    ----------
    name:
        Identifier used in the TOML file (e.g. ``perf``).
    expr:
        Arithmetic expression as a string.  Supported operators are the usual
        Python arithmetic ones; ``^`` is automatically converted to ``**``.
    direction:
        Either ``"minimize"`` or ``"maximize"``; defaults to ``"minimize"``.
    """

    name: str
    expr: str
    direction: Literal["minimize", "maximize"] = "minimize"

    def evaluate(self, results: Mapping[str, float], geometry: Any) -> float:
        """Safely evaluate the expression.

        ``results`` is a mapping of result symbols (e.g. ``pressure_drop``) to
        numeric values produced by a campaign run. ``geometry`` is the machine
        geometry object required to compute the ``ratio`` helper.
        """
        # Prepare the symbol table – copy the results and add the derived ratio.
        symtable: dict[str, float] = dict(results)
        try:
            symtable["ratio"] = exchanger_length_to_diameter(geometry)
        except Exception as exc:
            raise ValueError(f"Failed to compute ratio for expression objective: {exc}") from exc

        # Convert caret exponentiation to Python's ``**``.
        safe_expr = self.expr.replace("^", "**")

        # ``minimal=True`` disables builtins and only the provided symbols are
        # available, keeping the evaluation sandboxed.
        aeval = Interpreter(symtable=symtable, minimal=True, use_numpy=False)
        value = aeval(safe_expr)
        if aeval.error:
            msgs = "; ".join(str(e.get_error()) for e in aeval.error)
            raise ValueError(f"Error evaluating expression '{self.expr}': {msgs}")
        if not isinstance(value, (int, float)):
            raise ValueError(
                f"Expression did not evaluate to a numeric type: {value!r}"
            )
        return float(value)
