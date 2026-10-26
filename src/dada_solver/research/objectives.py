'''Expression‑based objectives for Dada‑Engine research.

This module defines :class:`ExpressionObjective`, which allows users to specify
arbitrary arithmetic expressions in the study TOML.  The expression may use any
symbol present in the result dictionary returned by the thermodynamic evaluator
plus the derived ``ratio`` variable (micro‑tube length‑to‑diameter).  The caret
operator ``^`` is treated as exponentiation (mapped to ``**``).
'''

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal, Mapping, Any

from asteval import Interpreter  # safe sandboxed evaluator

from dada_solver.sizing.objectives import ObjectiveValue
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

    def evaluate(
        self,
        target: Mapping[str, float] | Any,
        geometry: Any = None,
        **context: Any,
    ) -> float | ObjectiveValue:
        """Evaluate the arithmetic objective.

        Supports two invocation patterns:
        1. Standalone / dictionary evaluation:
           ``evaluate(results_mapping, geometry=None) -> float``
        2. Solver / campaign evaluation:
           ``evaluate(evaluation, **context) -> ObjectiveValue``
        """
        safe_expr = self.expr.replace("^", "**")

        # Case 1: Standalone evaluation with Mapping
        if isinstance(target, Mapping):
            symtable: dict[str, float] = dict(target)
            if "ratio" not in symtable and geometry is not None:
                try:
                    symtable["ratio"] = exchanger_length_to_diameter(geometry)
                except Exception as exc:
                    raise ValueError(f"Failed to compute ratio for expression objective: {exc}") from exc

            aeval = Interpreter(usersyms=symtable, use_numpy=False)
            value = aeval(safe_expr)
            if aeval.error:
                msgs = "; ".join(str(e.get_error()) for e in aeval.error)
                raise ValueError(f"Error evaluating expression '{self.expr}': {msgs}")
            if not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(
                    f"Expression did not evaluate to a finite numeric type: {value!r}"
                )
            return float(value)

        # Case 2: Sizing / campaign evaluation with DesignEvaluation (or SimpleNamespace)
        evaluation = target
        usable = bool(getattr(evaluation, "usable", False))
        perf = getattr(evaluation, "performance", None)
        if not usable or perf is None:
            return ObjectiveValue(self.name, None, available=False)

        symtable = {
            "power": perf.gas_power,
            "gas_power": perf.gas_power,
            "indicated_power_w": perf.gas_power,
            "heat_input_w": perf.heat_in_power,
            "heat_in_power": perf.heat_in_power,
            "heat_out_w": perf.heat_out_power,
            "heat_out_power": perf.heat_out_power,
            "indicated_thermal_efficiency": perf.thermal_efficiency,
            "thermal_efficiency": perf.thermal_efficiency,
            "mechanical_input_power": perf.mechanical_input_power,
            "indicated_mechanical_input_power_w": perf.mechanical_input_power,
        }
        if perf.motor_power is not None:
            symtable["motor_power"] = perf.motor_power
            symtable["motor_power_w"] = perf.motor_power
        if perf.cooling_power is not None:
            symtable["cooling_power"] = perf.cooling_power
            symtable["cooling_power_w"] = perf.cooling_power
        if perf.cooling_cop is not None:
            symtable["cooling_cop"] = perf.cooling_cop
            symtable["cop"] = perf.cooling_cop

        # Microtube productivity context if provided
        hi_count = context.get("heat_in_microtube_count")
        ho_count = context.get("heat_out_microtube_count")
        if hi_count is not None and ho_count is not None and (hi_count + ho_count) > 0:
            total_tubes = hi_count + ho_count
            if perf.cooling_power is not None:
                cooling_prod = perf.cooling_power / total_tubes
                symtable["cooling_power_per_total_microtube_w"] = cooling_prod
                if perf.cooling_cop is not None:
                    symtable["cooling_cop_times_power_per_total_microtube_w"] = perf.cooling_cop * cooling_prod
            if perf.motor_power is not None:
                motor_prod = perf.motor_power / total_tubes
                symtable["motor_power_per_total_microtube_w"] = motor_prod
                if perf.thermal_efficiency is not None:
                    symtable["thermal_efficiency_times_power_per_total_microtube_w"] = perf.thermal_efficiency * motor_prod

        # Exchanger geometry ratio
        geom = getattr(evaluation, "design", None) or getattr(evaluation, "model", None) or geometry
        if geom is not None and "ratio" not in symtable:
            try:
                symtable["ratio"] = exchanger_length_to_diameter(geom)
            except Exception:
                pass

        aeval = Interpreter(usersyms=symtable, use_numpy=False)
        value = aeval(safe_expr)
        if aeval.error or not isinstance(value, (int, float)) or not math.isfinite(value):
            return ObjectiveValue(self.name, None, available=False)

        score = -float(value) if self.direction == "maximize" else float(value)
        return ObjectiveValue(self.name, score, available=True)
