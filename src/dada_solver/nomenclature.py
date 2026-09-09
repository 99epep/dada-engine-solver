"""Study names for display; legacy Python/TOML identifiers remain readable."""


def study_name(identifier: str) -> str:
    return {
        "C": "H_i",
        "H": "H_o",
        "small_to_cold": "S_to_H_i",
        "large_to_hot": "L_to_H_o",
        "cold_to_large": "H_i_to_L",
        "hot_to_small": "H_o_to_S",
        "cold_exchanger_infeasible": "H_i_exchanger_infeasible",
        "hot_exchanger_infeasible": "H_o_exchanger_infeasible",
    }.get(identifier, identifier)
