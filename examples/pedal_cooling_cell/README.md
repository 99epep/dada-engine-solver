# Pedal-powered cooling cell: three mechanical approaches

This example explores a **low-tech, human-powered refrigeration cell**: a closed air circuit, a nominal **1 bar air-fill objective**, and a mechanism intended to be driven by a pedal crank. The aim is to investigate what cooling performance is attainable without a compressor requiring a high-pressure refrigerant charge. These are **simulation results**, not measurements from a built machine.

The three studies test a mechanical hypothesis: **more articulated linkages can better approximate a favorable piston-motion law**. Each of the two pistons is driven by its own mechanism. In the animated reports, the small- and large-cylinder linkages have deliberately been given different geometries to illustrate the range of possible implementations, rather than to suggest that both sides must look alike.

## Explore the results

Open the interactive, self-contained HTML reports first; each includes animated linkages and the thermodynamic traces:

- [Slider-crank report](slider/slider_coolcell.html)
- [Four-bar report](fourbar/fourbar_coolcell.html)
- [Six-bar report](sixbar/sixbar_coolcell.html)

Each directory also contains the following, using the same `*_coolcell` prefix:

- `.json`: evaluated result, including exact COP, powers, validity checks and diagnostic data.
- `.toml`: Research study definition, with design parameters, search ranges and constraints.
- `.basis.json`: portable starting machine/thermodynamic configuration used by the study.
- `.small.mechanism.json` and `.large.mechanism.json`: the two physical linkage definitions.
- `.html`: a standalone report combining the results, plots and mechanism animations.

## Common thermodynamic setup

The DADA cell circulates air between **two cylinders of different swept volumes** through hot and cold microtube heat exchangers. Passive valves select the flow path. Compression and expansion are interleaved with heat-transfer phases, aiming toward a Brayton/Joule-like refrigeration cycle. The heat exchangers have finite gas-side transfer and pressure losses, finite wall thermal storage and prescribed external fluid streams.

All three studies use **air**, a **−15 °C cold source** and a **25 °C hot source** (a 40 K temperature lift), at **0.933 Hz** (about 56 rpm). Their saved configurations specify a fixed inventory of approximately **20.7 g of air**, rather than imposing 1 bar throughout the cycle: internal pressures necessarily vary during operation. Both external streams are modeled as water-like liquids with prescribed flow and heat-transfer conductance.

At these source temperatures, the **Carnot refrigeration COP** is approximately **6.45**. The comparisons below report cooling output and *indicated* mechanical input: bearing, gear and pedal losses, as well as external pump/fan power, are **not included**.

## Results, in increasing order of cooling COP

| Piston mechanisms | Cooling COP | Fraction of Carnot | Cooling power | Indicated input |
| --- | ---: | ---: | ---: | ---: |
| [Slider-crank](slider/slider_coolcell.html) | **1.53** | **23.7%** | **212 W** | **139 W** |
| [Four-bar](fourbar/fourbar_coolcell.html) | **1.57** | **24.3%** | **227 W** | **145 W** |
| [Six-bar](sixbar/sixbar_coolcell.html) | **1.60** | **24.8%** | **237 W** | **148 W** |

### 1. Slider-crank — the simplest reference

Two independently phased slider-crank mechanisms provide a mechanically straightforward baseline. With about **23 L total swept volume**, a **0.76 small/large swept-volume ratio**, and approximately **217 cold-side / 223 hot-side microtubes**, this design reaches a COP of **1.53**.

### 2. Four-bar — more control over piston motion

Adding four-bar linkages changes the piston displacement and dwell profile over each revolution. The jointly retuned thermal design uses roughly **25 L swept volume**, a **0.76 swept-volume ratio**, and **223 / 219 microtubes** (cold / hot). Its COP rises to **1.57**.

### 3. Six-bar — the most flexible motion

The six-bar mechanisms provide still more freedom to shape the two piston trajectories. With approximately **25 L swept volume**, a **0.77 swept-volume ratio**, and **225 / 218 microtubes** (cold / hot), the best of these three configurations achieves a COP of **1.60**.

The exchanger tube lengths and diameters were also adjusted slightly between optimizations. Thus, the progression supports the value of **more flexible mechanical synthesis**, but it is **not a controlled proof** that the linkage complexity alone explains every COP gain, nor a proof of global optimality. The precise dimensions, values and simulation conditions are preserved in the accompanying files.

**Scope:** The modeled external streams have generous prescribed capacities and their circulation power is excluded. Actual pedal effort, mechanical reliability, exchanger construction and achievable real-world performance still require engineering validation.
