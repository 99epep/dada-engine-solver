# Motor exchangers: duty envelope and initial literature comparison

## Current gas-model source audit (2026-09-16)

The implementation ledger is [MICROTUBE_GAS_MODEL.md](MICROTUBE_GAS_MODEL.md);
measured-versus-predicted errors are in [EXCHANGER_VALIDATION.md](EXCHANGER_VALIDATION.md).
The older design survey below is retained as historical context, not as a
statement of the current reservoir temperatures or current gas-model capability.

| Source | Evidence and role | Calibration decision |
| --- | --- | --- |
| Ewart, Perrier, Graur et al., CFM 2007, [HAL](https://hal.science/hal-03361781v1), local `HAL_microtubes.pdf` | Steady isothermal N2/Ar/He in fused silica. Eq.3/5 and Tables 1–3 verify the second-order slip implementation and its VHS convention. | Published fits reproduced; metal accommodation remains unknown. |
| Aubert (1999), local `50015150.pdf`, chapter III §3.4.3.2, PDF pp.181–186; annex A13–A22 | Experimental pressure transmission through silica microtubes, including parallel tubes. Gain agrees better than phase. | Unsteady diagnostic motivation only; no h(t) measurement or thermal multiplier. |
| [Doty et al. 1991](https://dotynmr.com/download/pubs/1991_HTE_Doty_HeatExchanger.pdf), Eq.7, Tables 1–2 | Closest supplied integral gas/gas exchanger data. N2/He flow, pressure, four temperatures, UA, effectiveness and tube loss remain in the original CSVs. | Primary integral benchmark; discrepancies retained, no opaque calibration. |
| [Yang, Hong, Morini & Asako 2014](https://cris.unibo.it/handle/11585/332323), DOI 10.1016/j.ijheatmasstransfer.2014.07.017, pp.732–740 | Experimental/numerical heated gas in stainless tubes, 750/510/170 micrometres, Re 3000–12000. | Low-compressibility conventional baseline and compressibility warning; no enhancement formula inferred from an abstract. |
| [Turner 2003](https://digitalcommons.uri.edu/oa_diss/1625/), DOI 10.23860/diss-1625, abstract | Experimental etched silicon/glass microchannels, hydraulic diameters 5–96 micrometres; friction and heat transfer depend on Ma/Kn. | Different geometry/boundary; no tube Nu calibration. This is Turner, not Aubert. |
| [Gavasane et al.](https://doi.org/10.1080/10407782.2024.2351199), abstract, pp.7410–7428 | **Numerical DSMC**, Kn 0.05–5, pressure ratios 1.2–4 and walls 323–723 K. | Domain warning, not experimental data. No DSMC in cycle integration. |
| Lelea–Nishio–Takano, supplied [distilled-water paper](https://www.academia.edu/18316510/The_experimental_research_on_microtube_heat_transfer_and_fluid_flow_of_distilled_water) | Liquid-water experiments; entry/conjugate-conduction methodology. | `out_of_scope_liquid_calibration`; no gas coefficient inferred. Full linked text was not used to extract equations. |
| Qi et al., [S001793100700542X](https://www.sciencedirect.com/science/article/abs/pii/S001793100700542X), DOI 10.1016/j.ijheatmasstransfer.2007.08.017 | Liquid nitrogen flow boiling, Part II. | **`out_of_scope_two_phase`**. |
| [Wang & Dai 2014](https://doi.org/10.1007/s12209-014-2138-5), abstract, pp.21–26 | Micro-sized exchanger experiments, Re 500–1800. Available abstract does not establish a transferable DADA gas model. | No coefficients extracted from inaccessible full text. |
| [Qpedia December 2007](https://www.qats.com/Qpedia-Article/Qpedia_Dec07_Flow_and_heat_transfer_correlations_for_mini_and_macrochannels) | Engineering survey, not a new gas experiment. | Background only. |
| [Thermopedia 1121](https://www.thermopedia.com/content/1121/) | Access was blocked during this audit. | No equation extracted or calibration claimed. |
| [S1738573323004746](https://www.sciencedirect.com/science/article/pii/S1738573323004746), [S2214157X24005896](https://www.sciencedirect.com/science/article/pii/S2214157X24005896) | Full source contents were not retrievable. | `unverified_source_not_used`; no invented identification or coefficients. |
| Supplied ResearchGate records: [experimental uncertainty thesis](https://www.researchgate.net/publication/45351722_Analyse_d%27ecoulements_liquides_ou_gazeux_en_micro-conduites_decouplage_des_incertitudes_experimentales), [single-phase review](https://www.researchgate.net/publication/279392413_Single-phase_fluid_flow_and_heat_transfer_in_microtubes) | Methodological/review leads; no quantitative equation extracted from these records. | Not used as a calibration dataset. |

Production equations are attributed individually in the model ledger, including
transport tables and their temperature envelopes. Bibliographic presence here
is not evidence that a source validates the DADA pulse waveform or a metal/gas
surface pair.

## Historical architecture survey

Initial research: 2026-09-08. Updated user requirements: 2026-09-09.
Motor efficiency guides sizing under explicit compactness and useful-power
requirements. Ideal chronology remains the target; isothermality is deferred
and has no veto in the core validity assessment. The human-powered cooling
cell's low-tech fabrication objective does not constrain the motor study.

The current target is **2–10 Hz**, with reservoirs at **298.15 K (25 C)** and
**448.15 K (175 C)**. Prefer a documented commercial exchanger and adapt the
machine around it. If suitable commercial information cannot be found, the user
authorizes extrapolation from Doty's measurements, with the assumptions and
uncertainty kept explicit. See [MOTOR_DEMONSTRATOR.md](MOTOR_DEMONSTRATOR.md).

## 1. Historical controlled-cycle duty

The figures below belong to the original 0.1515 Hz, 600/300 K example. They
are a reproducible baseline, **not the requirements of the new motor**.

```console
PYTHONPATH=src python3 -m dada_solver.exchangers.duty \
  examples/motor_controlled_example.toml \
  --output-prefix outputs/motor_exchanger_duty
```

The JSON contains summaries. The CSV retains simultaneous states and flows
at both ports of each exchanger, using a dry ideal gas with air-like constants.
This is neither an optimized motor nor a final hardware specification.

The period is **6.600 s**, corresponding to **0.1515 Hz** or **9.09 rpm**.
The intended flow is **intermittent and unidirectional**. Cycle frequency is
not the frequency of a zero-mean sinusoidal flow. Pulse shape, active duration
and rise/fall times must accompany comparisons with published experiments.

| Calculated quantity | H_i, heat input | H_o, heat rejection |
|---|---:|---:|
| Reservoir temperature | 600 K | 300 K |
| Lumped exchanger gas temperature | 576.87–599.94 K | 298.53–323.66 K |
| Absolute gas pressure | 1.239–3.251 bar | 1.239–3.251 bar |
| Mean heat received by gas | +1079.93 W | −887.56 W |
| Peak forward inlet-port flow | 20.22 g/s | 19.63 g/s |
| Peak forward outlet-port flow | 16.97 g/s | 20.65 g/s |
| Mean outlet-port flow | 5.896 g/s | 5.896 g/s |
| Total active outlet-flow duration | 4.290 s | 3.292 s |
| Active outlet-flow fraction | 65.0% | 49.9% |
| Assumed UA, to replace with hardware data | 125 W/K | 125 W/K |
| Assumed gas volume | 0.660 L | 0.660 L |

Here, active flow exceeds **1% of that port's peak forward flow**. This is an
explicit, configurable reporting threshold, not a physical valve threshold.
Durations are accumulated over a cycle and need not describe one uninterrupted
pulse. Extrema need not be simultaneous: use the CSV for joint calculations.
Adjacent-cylinder temperatures are included to distinguish incoming gas from
the lumped exchanger temperature.

### Discrepancy from the intended unidirectional flow

The H_i→L and H_o→S check valves never produce negative flow. The existing
S↔H_i and L↔H_o internal links allow local reflux. Integrated reverse mass is
approximately **0.333%** and **2.134%** of forward mass at the H_i and H_o inlet
ports, respectively. These are model results, not measurements of a real motor.
Signed histories preserve them rather than hiding them behind absolute values.

No additional check valves are silently introduced. Literature selection still
targets intermittent unidirectional flow; the discrepancies help assess the
distance from ideal chronology and the limits of instantaneous steady maps.

## 2. Plate-fin and parallel-microtube candidates

| Architecture | Evidence found | Decisive motor checks |
|---|---|---|
| Aluminum plate-fin | Industrial hardware, manufacturer tests, several fin types | Joint temperature/pressure rating, headers, pressure loss and gas volume |
| Stainless-steel or suitable-alloy plate-fin | Manufacturer documentation for hot processes | Thermal/hydraulic maps for a product at the required scale |
| Parallel microtubes | Gas prototype measurements of UA, effectiveness and pressure drop | Tube/header flow distribution, temperature and flow scaling |

**No architecture has been declared the winner.** Compare recalculated motor
efficiency at the same energy boundary, accounting jointly for heat transfer,
pressure loss and gas volume. Exchanger effectiveness is not motor efficiency.
Auxiliaries must be included to move from gas work to useful system output.

### Primary sources and scope

1. **ALPEMA — aluminum plate-fin applications.** The association states a
   temperature range up to **204 C**. The historical 600 K reservoir exceeded
   that range. The new 175 C reservoir is within it, but this does not certify
   every product, its seals, or its combined pressure/temperature rating.
   [Source](https://alpema.net/applications.html).

2. **Sumitomo Precision Products — Stainless Steel Plate-Fin Heat Exchanger.**
   The manufacturer presents stainless steel for hot processes and also names
   nickel-rich and titanium-based alloys. This establishes a fabrication route,
   not a UA/pressure-drop map for DADA.
   [Manufacturer sheet](https://www.spp.co.jp/netsu/__assets/en/pdf/sumitomo_downloads_02_eng.pdf).

3. **Chart — Brazed Aluminum Heat Exchangers**, printed pages 8–10. The
   manufacturer describes custom designs, transient-analysis capability and
   experimentally characterized fins. The brochure itself is not a numerical
   performance table for the solver.
   [Brochure](https://files.chartindustries.com/Brazed-Aluminum-Heat-Exchangers.pdf).

4. **Doty et al. (1991), The Microtube Strip Heat Exchanger**, Heat Transfer
   Engineering 12(3), printed page 38, Tables 1 and 2. Nitrogen and helium
   prototype tests provide both thermal and hydraulic measurements. Three
   nitrogen points at 322 kPa are transcribed in SI in
   `examples/data/doty_1991_nitrogen_reference.csv`, with reported uncertainties.
   The first measured UA is 3.7 W/K versus 7.1 W/K calculated in the paper.
   These are prototype results, not proof of present commercial availability
   or performance under interrupted flow.
   [Paper hosted by the industrial author](https://dotynmr.com/download/pubs/1991_HTE_Doty_HeatExchanger.pdf).

5. **Zohir et al. (2006), An experimental investigation of heat transfer to
   pulsating pipe air flow with different amplitudes.** Air tests cover Re
   750–12,320 and 1–10 Hz with imposed wall heat flux. Frequency now overlaps
   the proposed motor envelope; waveform, Reynolds number and boundary
   conditions still require comparison. Verify the waveform and absence of
   reflux in the full experiment before transferring a correlation. No universal
   enhancement factor is applied.
   [Authors' institution](https://pure.kfupm.edu.sa/en/publications/an-experimental-investigation-of-heat-transfer-to-pulsating-pipe-/).

6. **Bhowmik and Lee (2009), Analysis of heat transfer and pressure drop
   characteristics in an offset strip fin heat exchanger.** This steady
   numerical study uses water; it is not gas-pulse validation. It helps identify
   geometry-dependent fin correlations. Offset-strip fins differ from the
   smooth rectangular channels in the current module.
   [Paper](https://www.sciencedirect.com/science/article/pii/S0735193308002443).

## 3. Measurements and authorized extrapolation

The Doty data remain measured prototype references, with lower flow rates,
different temperatures and a specific geometry. T1–T4 retain the source's
sensor names. Even the three equal-pressure rows have different temperatures,
so they are not a controlled UA-versus-flow sweep at fixed boundary conditions.

The user now authorizes extrapolation if adequate commercial data are missing.
That supersedes the earlier no-extrapolation research boundary. The estimate
must retain source conditions, specify gas-property and geometric scaling,
identify the changes in Reynolds number and pulsed operation, and report
sensitivity scenarios. Extrapolated values must never be labeled measurements
or manufacturer ratings. The current 125 W/K assumptions are not automatically
replaced merely because a reference table exists.

Before coupling an exchanger estimate, record its geometry, passages and header
volume, gas and external-fluid conditions, thermal transfer and pressure-drop
range, measurement boundary, uncertainty, and pulse assumptions. Further
transient modeling remains deferred until its benefit and cost are explained.
