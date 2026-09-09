"""Geometry-based heat-exchanger screening and sizing tools."""

from dada_solver.exchangers.models import (
    ExchangerOperatingPoint,
    ExchangerPerformance,
    ParallelRectangularChannels,
    ThermalResistanceModel,
    TransportProperties,
    evaluate_parallel_channels,
)
from dada_solver.exchangers.sizing import (
    ExchangerRequirements,
    ExchangerSizingResult,
    screen_parallel_channel_designs,
)
from dada_solver.exchangers.cycle_adapter import (
    ExchangerSide,
    conservative_port_flow_operating_point,
)
from dada_solver.exchangers.optimization import (
    ContinuousGeometryBounds,
    ExchangerObjective,
    ExchangerOptimizationResult,
    optimize_parallel_channel_geometry,
)
from dada_solver.exchangers.coupling import (
    CoupledExchangerSizer,
    CoupledExchangerSizingResult,
    CoupledSizingSettings,
    CoupledSizingStatus,
    ExchangerGeometrySearch,
)
from dada_solver.exchangers.hydraulic_validity import (
    HydraulicQuasiSteadyValidity,
    HydraulicValidityThresholds,
    assess_hydraulic_quasi_steady_validity,
)
from dada_solver.exchangers.hydraulic_comparison import (
    HydraulicComparisonPoint,
    compare_hydraulic_models,
)
from dada_solver.exchangers.two_sided import (
    ChannelSidePerformance,
    IncompressibleFluidProperties,
    LiquidLoopOperatingPoint,
    TwoSidedExchangerPerformance,
    combine_channel_sides,
    channel_side_from_gas_performance,
    evaluate_incompressible_channel_side,
    evaluate_two_sided_exchanger,
)

__all__ = [
    "ExchangerOperatingPoint",
    "ExchangerObjective",
    "ExchangerOptimizationResult",
    "ExchangerGeometrySearch",
    "ExchangerPerformance",
    "ExchangerRequirements",
    "ExchangerSizingResult",
    "ExchangerSide",
    "ParallelRectangularChannels",
    "ContinuousGeometryBounds",
    "CoupledExchangerSizer",
    "CoupledExchangerSizingResult",
    "CoupledSizingSettings",
    "CoupledSizingStatus",
    "HydraulicQuasiSteadyValidity",
    "HydraulicComparisonPoint",
    "HydraulicValidityThresholds",
    "ChannelSidePerformance",
    "IncompressibleFluidProperties",
    "LiquidLoopOperatingPoint",
    "TwoSidedExchangerPerformance",
    "ThermalResistanceModel",
    "TransportProperties",
    "evaluate_parallel_channels",
    "conservative_port_flow_operating_point",
    "optimize_parallel_channel_geometry",
    "assess_hydraulic_quasi_steady_validity",
    "compare_hydraulic_models",
    "combine_channel_sides",
    "channel_side_from_gas_performance",
    "evaluate_incompressible_channel_side",
    "evaluate_two_sided_exchanger",
    "screen_parallel_channel_designs",
]
