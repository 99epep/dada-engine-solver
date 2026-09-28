"""Geometry diagnostics using production joint states; all screens are sampled.

Triangle-clearance helpers preserve the Stage 2F convention.
"""
import math
import numpy as np

def segment_distance_origin_batch(A: np.ndarray, B: np.ndarray):
    AB = B - A
    den = np.sum(AB*AB, axis=1)
    t = -np.sum(A*AB, axis=1) / np.maximum(den, 1e-30)
    t = np.clip(t, 0.0, 1.0)
    Q = A + t[:, None]*AB
    return np.linalg.norm(Q, axis=1)


def signed_triangle_clearance_origin(E: np.ndarray, F: np.ndarray, H: np.ndarray):
    d = np.minimum.reduce([
        segment_distance_origin_batch(E, F),
        segment_distance_origin_batch(F, H),
        segment_distance_origin_batch(H, E),
    ])

    def cross_to_origin(A, B):
        AB = B - A
        AO = -A
        return AB[:, 0]*AO[:, 1] - AB[:, 1]*AO[:, 0]

    s1 = cross_to_origin(E, F)
    s2 = cross_to_origin(F, H)
    s3 = cross_to_origin(H, E)
    eps = 1e-12

    inside = (
        ((s1 >= -eps) & (s2 >= -eps) & (s3 >= -eps))
        | ((s1 <= eps) & (s2 <= eps) & (s3 <= eps))
    )
    signed = np.where(inside, -d, d)
    return float(np.min(signed)), int(np.sum(inside))



def zero_crossings(velocity):
    """Cyclic sign changes after removing numerically zero samples, without smoothing."""
    values=np.asarray(velocity)
    scale=float(np.max(np.abs(values)))
    nonzero=values[np.abs(values)>max(scale*1e-12,1e-15)]
    signs=np.sign(nonzero)
    return int(np.count_nonzero(signs!=np.roll(signs,1))) if len(signs) else 0


def six_bar_metrics(mechanism, samples=1440):
    states=[mechanism.joint_state(float(t)) for t in np.linspace(0,2*math.pi,samples,endpoint=False)]
    joints={name:np.array([s['joints'][name] for s in states]) for name in states[0]['joints']}
    stroke=mechanism.stroke_over_crank
    transverse=np.array([s['transverse'] for s in states])
    clearance,inside=signed_triangle_clearance_origin(joints['E'],joints['F'],joints['H'])
    return dict(stroke_over_crank=stroke,
        minimum_primary_transmission_sine=min(s['primary_transmission_sine'] for s in states),
        minimum_secondary_transmission_sine=min(s['secondary_transmission_sine'] for s in states),
        minimum_rod_axis_cosine=min(s['rod_axis_cosine'] for s in states),
        EH_over_crank=mechanism.link_ef*math.hypot(mechanism.h_along_over_ef,mechanism.h_normal_over_ef),
        H_axis_lateral_rms_over_stroke=float(np.std(transverse))/stroke,
        H_axis_lateral_span_over_stroke=float(np.ptp(transverse))/stroke,
        crank_axis_to_EFH_clearance_over_crank=clearance,
        crank_axis_inside_EFH_frames=inside,
        zero_crossing_count=zero_crossings([s['derivative'] for s in states]))


def four_bar_metrics(model, side, samples=1440, envelope_frame_angle_rad=0.):
    assembly=getattr(model,side+'_assembly')
    normalized=getattr(model,'_'+side)
    states=[assembly.evaluate(*model._crank(float(t))) for t in np.linspace(0,2*math.pi,samples,endpoint=False)]
    positions=[]
    for state in states:
        slider=assembly.slider
        p=(slider.axis_origin_x+state.coordinate*math.cos(slider.axis_angle),
           slider.axis_origin_y+state.coordinate*math.sin(slider.axis_angle))
        positions.extend(((0.,0.),(assembly.loop.rocker_pivot_x,assembly.loop.rocker_pivot_y),
                          state.crank_pin,state.coupler_joint,state.output_point,p))
    c,s=math.cos(envelope_frame_angle_rad),math.sin(envelope_frame_angle_rad)
    positions=np.asarray(positions) @ np.array([[c,s],[-s,c]])
    envelope=float(np.linalg.norm(np.ptp(positions,axis=0)))
    return dict(stroke_over_crank=normalized.stroke/model.crank_radius,
        minimum_primary_transmission_sine=min(abs(s.four_bar_cross_product) for s in states)/(assembly.loop.coupler_length*assembly.loop.rocker_length),
        minimum_rod_axis_cosine=min(s.connecting_rod_transverse_margin for s in states)/assembly.slider.connecting_rod_length,
        stroke_over_envelope=normalized.stroke/envelope,
        zero_crossing_count=zero_crossings([s.coordinate_derivative for s in states]))
