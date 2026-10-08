"""Whole-machine report animation using production closures and replayed gas states.

Independent visible-stroke scales are presentation only, not a physical shared
shaft assembly. Cylinder placement follows the validated machine drawing.
"""
import io
import math
import numpy as np

FRAME_COUNT = 66
FRAME_DURATION_MS = 85
EXCHANGER_GAP = 9.6
MACHINE_SCALE = 0.90
CYLINDER_LINE_WIDTHS = dict(small=6.5,large=7.)
VERTICAL_MARGIN_PX = 20
# Presentation calibration: 100-dpi reference drawing, 1176 horizontal pixels
# over 28.8730100164066 drawing units. Keep assembly strokes in drawing units
# so viewport zoom scales them together with conduits and exchanger glyphs.
REFERENCE_PIXELS_PER_UNIT = 40.73007972953834
REFERENCE_DPI = 100.
TRANSFER_ROW_OFFSET = .155 * 1.25


def temperature_rgb(T, Tmin, Tmax):
    a = 0.5 if Tmax <= Tmin else (float(T) - Tmin) / (Tmax - Tmin)
    a = min(1.0, max(0.0, a))
    return (a, 0.0, 1.0 - a)


def candidate_mechanisms(study, physical):
    """Resolve exact candidate coordinates through the current Research adapter."""
    from .schema import compile_study
    from .artifacts import MechanismArtifact
    from .families import PHYSICAL_FAMILIES, parameter_specs
    definition = compile_study(study)
    values = dict(definition.fixed_parameters, **physical)
    artifacts = {}
    for side in ('small', 'large'):
        settings = study.settings[side]
        if settings['family'] not in PHYSICAL_FAMILIES:
            return None
        settings = {k:v for k,v in settings.items() if k not in ('artifact', 'sha256')}
        geometry = {name:values[f'kinematics.{"shared" if study.data['kinematics']['coupling']=='shared_crank' and name in ('phase_rad','crank_direction') else side}.{name}']
                    for name in parameter_specs(settings, side)}
        source = study.artifacts.get(side)
        if source is not None and source.scientific['geometry'] == geometry:
            artifacts[side] = source
        else:
            artifacts[side] = MechanismArtifact.create(settings['family'], geometry, settings=settings,
                constraints=source.scientific['constraints'] if source else ())
    design = definition.adapter.build(values)
    return artifacts, design.configuration


def _oriented(model, angle, side):
    state = model.state(float(angle))
    geometry = model.geometry
    if model.family == 'six_bar':
        axis = geometry.slider_axis_angle
        origin = np.array((-math.sin(axis), math.cos(axis))) * geometry.slider_axis_offset
    elif model.family == 'four_bar':
        slider = geometry.slider
        axis = slider.axis_angle
        origin = np.array((slider.axis_origin_x, slider.axis_origin_y))
    else:
        axis = 0.
        origin = np.array((0., geometry.offset_over_crank))
    rotation = (0. if side == 'small' else math.pi) - axis
    c, s = math.cos(rotation), math.sin(rotation)
    matrix = np.array(((c,-s),(s,c)))
    dy = (matrix @ origin)[1]
    joints = {name:matrix @ np.asarray(point) - (0.,dy)
              for name,point in state['joints'].items()}
    # Gas is on the inner-head side: increasing volume means outward travel.
    # Six-bar's production normalization already decreases with coordinate.
    if model.family != 'six_bar' and model.artifact.scientific['geometry']['volume_increases_with_coordinate']:
        joints = {name:point*np.array((-1.,1.)) for name,point in joints.items()}
    return joints, state['links']


def _linkage_hits_cylinder(envelope, links, shift, head, outer, width, gap):
    """Test production bar segments against the three cylinder walls.

    The outward cylinder end is open. P is the piston articulation, not a
    separate obstacle: a linkage fitting through that opening needs no rod.
    Expand wall rectangles by the drawing clearance before clipping segments.
    """
    segments = np.asarray([[state[a],state[b]] for state in envelope for a,b in links],float)
    segments[:,:,0] += shift
    low,high = sorted((head,outer))
    walls = ((low-gap,high+gap,width/2-gap,width/2+gap),
             (low-gap,high+gap,-width/2-gap,-width/2+gap),
             (head-gap,head+gap,-width/2-gap,width/2+gap))
    start = segments[:,0,:]; direction = segments[:,1,:]-start
    for x0,x1,y0,y1 in walls:
        enter = np.zeros(len(start)); leave = np.ones(len(start))
        possible = np.ones(len(start),dtype=bool)
        for axis,lo,hi in ((0,x0,x1),(1,y0,y1)):
            parallel = np.abs(direction[:,axis]) < 1e-14
            possible &= ~parallel | ((start[:,axis]>=lo)&(start[:,axis]<=hi))
            divisor = np.where(parallel,1.,direction[:,axis])
            t0=(lo-start[:,axis])/divisor; t1=(hi-start[:,axis])/divisor
            enter=np.maximum(enter,np.where(parallel,-np.inf,np.minimum(t0,t1)))
            leave=np.minimum(leave,np.where(parallel,np.inf,np.maximum(t0,t1)))
        if np.any(possible & (enter<=leave)):
            return True
    return False


def _place(samples, side, head, width, gap, *, envelope, visible_stroke, links):
    """Keep the piston in its chamber; translate the rigid linkage outward.

    A collision-free linkage can attach directly; otherwise the whole-cycle
    envelope determines an offset placing the linkage outside the cylinder.
    P remains the production rod articulation; piston_positions are the rigidly
    connected piston faces. Their constant separation is graphical only.
    """
    outward = -1 if side == 'small' else 1
    inner_piston = max(s['P'][0] for s in envelope) if side == 'small' else min(s['P'][0] for s in envelope)
    shift = head + outward*.035*width - inner_piston
    placed = [{k:v+np.array((shift,0.)) for k,v in state.items()} for state in samples]
    piston_positions = np.array([state['P'][0] for state in placed])
    outer = head + outward*1.10*visible_stroke
    # Include P itself and every other joint: the entire connecting rod and
    # linkage must stay outside, not merely its moving upstream endpoint.
    envelope_x = np.array([point[0]+shift for state in envelope for point in state.values()])
    delta = min(0.,outer-gap-envelope_x.max()) if side == 'small' else max(0.,outer+gap-envelope_x.min())
    # No mandatory extension: retain direct piston articulation if the full
    # cycle fits through the open cylinder end without touching its walls.
    wall_head = head - outward*_stroke_world(CYLINDER_LINE_WIDTHS[side])
    if not _linkage_hits_cylinder(envelope,links,shift,wall_head,outer,width,gap):
        delta = 0.
    for state in placed:
        for name in state:
            state[name] = state[name]+np.array((delta,0.))
    return placed, dict(x_inner=head,x_outer=outer,width=width,
                        piston_positions=piston_positions,rod_length=abs(delta),clearance=gap)


def machine_geometry(artifacts, study_angles, *, width_px=1200):
    """Production joints, fixed cylinder layout, and whole-cycle rod clearance."""
    from .mechanism_view import MechanismModel
    raw, envelopes, reference, links = {}, {}, {}, {}
    # Preserve the existing drawing's cylinder dimensions independently of the
    # new animation cadence. Dense closure samples determine linkage clearance.
    reference_angles = np.linspace(0.,2*math.pi,48,endpoint=False)
    envelope_angles = np.linspace(0.,2*math.pi,1440,endpoint=False)
    for side in ('small','large'):
        model = MechanismModel(artifacts[side])
        if model.primary:
            raise ValueError('A complete physical piston mechanism is required.')
        frames = [_oriented(model,angle,side) for angle in study_angles]
        raw[side] = [frame[0] for frame in frames]
        links[side] = frames[0][1]
        reference[side] = [_oriented(model,angle,side)[0] for angle in reference_angles]
        extrema = (model.law.model.stationary_points(model.law.side) if model.family=='four_bar'
                   else model.geometry.stationary_points())
        angles = np.concatenate((envelope_angles,[p['angle_rad'] for p in extrema],np.asarray(study_angles)))
        envelopes[side] = [_oriented(model,angle,side)[0] for angle in angles]
    strokes = {side:float(np.ptp([state['P'][0] for state in frames])) for side,frames in reference.items()}
    if min(strokes.values()) <= 1e-9:
        raise ValueError('Mechanism has no visible piston stroke.')
    common = .90 * (1.5 * .5 * sum(strokes.values()))
    widths = dict(small=common,large=common/.84)
    heads = dict(small=-MACHINE_SCALE*EXCHANGER_GAP/4,large=MACHINE_SCALE*EXCHANGER_GAP/4)
    local = {side:[{k:v*common/strokes[side] for k,v in state.items()} for state in frames]
             for side,frames in raw.items()}
    envelopes = {side:[{k:v*common/strokes[side] for k,v in state.items()} for state in frames]
                 for side,frames in envelopes.items()}
    def place(gap):
        return {side:_place(local[side],side,heads[side],widths[side],gap,
                           envelope=envelopes[side],visible_stroke=common,links=links[side]) for side in local}
    probes = place(0.)
    cloud = np.array([p for frames,_ in probes.values() for state in frames for p in state.values()])
    xmin = min(cloud[:,0].min(),*(c['x_outer'] for _,c in probes.values()))
    xmax = max(cloud[:,0].max(),*(c['x_outer'] for _,c in probes.values()))
    # Six pixels of visual margin plus the production drawing's joint radius.
    gap = 6.*(xmax-xmin)/max(1.,width_px-12.)+.085
    placed = place(gap)
    # Frame the dense envelope too: framing must not depend on which frames
    # happen to be exported, nor crop the new outboard mechanisms and rods.
    clouds = []
    for side,(frames,cyl) in placed.items():
        delta = frames[0]['P']-local[side][0]['P']
        clouds.extend(point+delta for state in envelopes[side] for point in state.values())
    cloud = np.array(clouds)
    xmin = min(cloud[:,0].min(),*(c['x_outer'] for _,c in placed.values()))
    xmax = max(cloud[:,0].max(),*(c['x_outer'] for _,c in placed.values()))
    ymin = min(cloud[:,1].min()-.31*.40*widths['large'],-widths['large']/2)
    ymax = max(cloud[:,1].max()+.18,widths['large']/2)
    from .transfer_block import BLOCK_WIDTH,HX_H,TRI_H
    block_scale = (heads['large']-heads['small'])/BLOCK_WIDTH
    row_extent = .31*widths['large']+max(HX_H,TRI_H)*block_scale
    ymin,ymax = min(ymin,-row_extent),max(ymax,row_extent)
    pad = .03*max(xmax-xmin,ymax-ymin)
    return dict(sides=placed,links=links,heads=heads,
                limits=(xmin-pad,xmax+pad,ymin-pad,ymax+pad))


def layout_for_configuration(configuration):
    """Order is Ho then Hi; operating direction does not determine placement."""
    return ''.join('U' if value == 'upstream' else 'D' for value in
                   (configuration.heat_out_valve_placement, configuration.heat_in_valve_placement))


def frame_series(cycle, frames):
    if type(frames) is not int or frames < 2:
        raise ValueError('At least two animation frames are required.')
    required = {'S_temperature_K','L_temperature_K','Ho_temperature_K',
                'Hi_temperature_K','Ho_to_S_valve_open','Hi_to_L_valve_open'}
    if not required.issubset(cycle['series']):
        raise ValueError('Machine replay is missing required temperatures or valve states.')
    angles = np.asarray(cycle['angle_rad'],float)
    if len(angles) < 2 or np.any(np.diff(angles)<0):
        raise ValueError('Replay cycle angles must be ordered.')
    query = np.linspace(0.,2*math.pi,frames,endpoint=False)
    result = {}
    for name,values in cycle['series'].items():
        values = np.asarray(values,float)
        if values.shape != angles.shape or not np.all(np.isfinite(values)):
            raise ValueError(f'Invalid machine replay series: {name}.')
        if name.endswith('_valve_open') and np.any((values < 0) | (values > 1)):
            raise ValueError(f'Invalid valve opening: {name}.')
        if name.endswith('_valve_open') and np.all(np.isin(values,[0.,1.])):
            # Model states are binary: retain transitions, never invent valve lift.
            result[name] = values[np.clip(np.searchsorted(angles,query,side='right')-1,0,len(values)-1)]
        else:
            result[name] = np.interp(query,angles,values)
    return query, result


def _stroke_world(reference_points):
    return reference_points * REFERENCE_DPI / 72. / REFERENCE_PIXELS_PER_UNIT


def _stroke_points(ax, reference_points):
    """Convert the calibrated assembly stroke to this viewport's point size."""
    pixels = abs(ax.transData.transform((1.,0.))[0]-ax.transData.transform((0.,0.))[0])
    return _stroke_world(reference_points) * pixels * 72. / ax.figure.dpi


def _cylinder_head(ax, cylinder, side):
    """Keep the calibrated cylinder face fixed independently of viewport zoom."""
    return cylinder['x_inner'] + (1 if side=='small' else -1)*_stroke_world(CYLINDER_LINE_WIDTHS[side])


def _transfer_endpoints(ax, geometry):
    """Reach each wall's outer edge: half a stroke overlap, no overshoot."""
    return {side:_cylinder_head(ax,cyl,side) + (-1 if side=='small' else 1)*
            .5*_stroke_world(CYLINDER_LINE_WIDTHS[side])
            for side,(_,cyl) in geometry['sides'].items()}


def draw_machine_frame(ax, geometry, index, data, layout, Tmin, Tmax):
    from matplotlib.patches import Circle, Rectangle
    from .transfer_block import PlacedAxes, LAYOUTS, draw_row
    ax.clear()
    ax.set_facecolor('#f0f0f0')
    xmin,xmax,ymin,ymax = geometry['limits']
    ax.set_xlim(xmin,xmax); ax.set_ylim(ymin,ymax)
    ax.set_aspect('equal',adjustable='box'); ax.apply_aspect()
    colors = ['#1f77b4','#ff7f0e','#2ca02c','#9467bd','#d62728','#444444']
    for side,temperature in (('small','S_temperature_K'),('large','L_temperature_K')):
        frames,cyl = geometry['sides'][side]
        state = frames[index]
        width,head,outer,piston = cyl['width'],cyl['x_inner'],cyl['x_outer'],cyl['piston_positions'][index]
        # Move only the cylinder head inward by its full rendered stroke
        # thickness. Keeping the outer end fixed extends it by that amount;
        # pistons and linkages retain their exact positions.
        linewidth = _stroke_points(ax,CYLINDER_LINE_WIDTHS[side])
        head = _cylinder_head(ax,cyl,side)
        x0,x1 = sorted((head,piston))
        ax.add_patch(Rectangle((x0,-width/2),x1-x0,width,
                              facecolor=temperature_rgb(data[temperature][index],Tmin,Tmax),
                              edgecolor='none',alpha=.92,zorder=1))
        ax.plot([piston,piston],[-width/2,width/2],color='#555555',lw=_stroke_points(ax,9.5 if side=='small' else 10.),zorder=4,solid_capstyle='butt')
        for xs,ys in (([head,head],[-width/2,width/2]),([head,outer],[width/2]*2),([head,outer],[-width/2]*2)):
            ax.plot(xs,ys,color='#111111',lw=linewidth,zorder=6,solid_capstyle='round')
        if cyl['rod_length'] > 0.:
            ax.plot([piston,state['P'][0]],[0.,0.],color='#555555',lw=_stroke_points(ax,4.),zorder=4,solid_capstyle='butt')
        ax.add_patch(Circle(state['P'],.075,facecolor='white',edgecolor='black',lw=1.1,zorder=8))
        if 'F' in state:
            _draw_six_bar(ax,state,width,lw=2.55 if side=='small' else 2.65)
        else:
            A,B = state['A'],state['B']
            ax.add_patch(Circle(A,float(np.linalg.norm(B-A)),fill=False,edgecolor='black',lw=.85,zorder=0))
            for i,(a,b) in enumerate(geometry['links'][side]):
                ax.plot(*zip(state[a],state[b]),color=colors[min(i,len(colors)-1)],lw=2.55,solid_capstyle='round',zorder=4)
            for name,point in state.items():
                if name != 'P':
                    ax.add_patch(Circle(point,.075,facecolor='white',edgecolor='black',lw=1.1,zorder=8))
            for name in ('A','D','G'):
                if name in state:
                    x,y = state[name]; scale=.40*width
                    ax.plot([x,x],[y-.08*scale,y-.22*scale],color='black',lw=1.8)
                    ax.plot([x-.18*scale,x+.18*scale],[y-.24*scale]*2,color='black',lw=1.8)
                    for k in range(5):
                        xi=x-.18*scale+k*(.36*scale/4)
                        ax.plot([xi-.04*scale,xi+.04*scale],[y-.31*scale,y-.24*scale],color='black',lw=1.2)
    large_width = geometry['sides']['large'][1]['width']
    endpoints = _transfer_endpoints(ax,geometry)
    for branch,direction,y in (('Ho','left',TRANSFER_ROW_OFFSET*large_width),('Hi','right',-TRANSFER_ROW_OFFSET*large_width)):
        placed = PlacedAxes(ax,endpoints['small'],endpoints['large'],y,zorder_offset=6)
        opening = data[branch+('_to_S_valve_open' if branch=='Ho' else '_to_L_valve_open')][index]
        draw_row(placed,y=0.,name=branch,**LAYOUTS[layout][branch],valve_direction=direction,
                 valve_open=float(opening),
                 left_color=temperature_rgb(data['S_temperature_K'][index],Tmin,Tmax),
                 right_color=temperature_rgb(data['L_temperature_K'][index],Tmin,Tmax),
                 hx_color=temperature_rgb(data[branch+'_temperature_K'][index],Tmin,Tmax))
    xmin,xmax,ymin,ymax = geometry['limits']
    ax.set_xlim(xmin,xmax); ax.set_ylim(ymin,ymax)
    ax.set_aspect('equal',adjustable='box'); ax.axis('off')


def render_machine_webp(artifacts, cycle, *, layout, motor_operation=False, frames=FRAME_COUNT):
    """Encode a deterministic fixed-frame whole-machine animation, without graphs."""
    from PIL import Image, features
    from PIL import _webp
    if not features.check('webp') or not hasattr(_webp, 'WebPAnimEncoder'):
        raise ValueError('Pillow animated WebP support is unavailable.')
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from .transfer_block import LAYOUTS
    if layout not in LAYOUTS:
        raise ValueError(f'Unknown transfer layout {layout!r}.')
    query,data = frame_series(cycle,frames)
    geometry = machine_geometry(artifacts,-query if motor_operation else query)
    temperatures = np.concatenate([np.asarray(cycle['series'][name]) for name in
                                   ('S_temperature_K','L_temperature_K','Ho_temperature_K','Hi_temperature_K')])
    Tmin,Tmax = float(temperatures.min()),float(temperatures.max())
    fig = Figure(figsize=(12,5),dpi=100,facecolor='#f0f0f0')
    canvas = FigureCanvasAgg(fig); ax = fig.add_axes((.01,.01,.98,.98))
    images = []
    top,bottom = 500,0
    try:
        for i in range(frames):
            draw_machine_frame(ax,geometry,i,data,layout,Tmin,Tmax)
            canvas.draw()
            pixels = np.asarray(canvas.buffer_rgba())
            rows = np.flatnonzero(np.any(pixels[:,:,:3]!=240,axis=(1,2)))
            if rows.size:
                top,bottom = min(top,int(rows[0])),max(bottom,int(rows[-1])+1)
            images.append(Image.fromarray(pixels.copy()).convert('RGB'))
        # Crop the same vertical window across all frames. This preserves every
        # other element's pixel scale and adds exactly 20 raw pixels of margin.
        height = bottom-top+2*VERTICAL_MARGIN_PX
        for i,image in enumerate(images):
            cropped = Image.new('RGB',(1200,height),'#f0f0f0')
            cropped.paste(image,(0,VERTICAL_MARGIN_PX-top))
            image.close(); images[i] = cropped
        output = io.BytesIO()
        images[0].save(output,format='WEBP',save_all=True,append_images=images[1:],
                       duration=FRAME_DURATION_MS,loop=0,quality=82,method=4)
        return output.getvalue(), dict(layout=layout,frames=frames,frame_duration_ms=FRAME_DURATION_MS,cycle_duration_ms=frames*FRAME_DURATION_MS,width=1200,height=height,vertical_margin_px=VERTICAL_MARGIN_PX,
            cylinder_inner_faces={side:_cylinder_head(ax,cyl,side) for side,(_,cyl) in geometry['sides'].items()},
            transfer_block_endpoints=_transfer_endpoints(ax,geometry),temperature_range_K=[Tmin,Tmax],
            method='One-cycle replay from saved periodic state; production linkage closures. Independent visible-stroke scaling is illustrative, not a shared-shaft assembly.')
    finally:
        for image in images: image.close()
        fig.clear()


_COLORS = {
    "crank": "#1f77b4",
    "coupler": "#ff7f0e",
    "rocker": "#2ca02c",
    "secondary_dyad": "#9467bd",
    "secondary_plate": "#d62728",
    "rod": "#444444",
    "ground": "#111111",
    "joint_fill": "#ffffff",
    "piston": "#555555",
    "cylinder": "#111111",
    "panel_bg": "#f0f0f0",
    "S": "#1f77b4",
    "L": "#ff7f0e",
    "Hi": "#d62728",
    "Ho": "#2ca02c",
}

def draw_ground_symbol(ax, p, scale, lw=2.0):
    x, y = p
    ax.plot([x, x], [y - 0.08 * scale, y - 0.22 * scale], color=_COLORS["ground"], lw=lw, solid_capstyle="round")
    base_y = y - 0.24 * scale
    half = 0.18 * scale
    ax.plot([x - half, x + half], [base_y, base_y], color=_COLORS["ground"], lw=lw, solid_capstyle="round")
    for k in range(5):
        xi = x - half + k * (2 * half / 4)
        ax.plot([xi - 0.04 * scale, xi + 0.04 * scale], [base_y - 0.07 * scale, base_y], color=_COLORS["ground"], lw=lw * 0.7, solid_capstyle="round")


def draw_joint(ax, p, r=0.075, lw=1.1):
    from matplotlib.patches import Circle
    ax.add_patch(Circle((p[0], p[1]), r, facecolor=_COLORS["joint_fill"], edgecolor="black", lw=lw, zorder=8))


def _draw_six_bar(ax, s, width, lw=2.5):
    from matplotlib.patches import Circle
    A, B, C, D, E = s["A"], s["B"], s["C"], s["D"], s["E"]
    F, G, H, P = s["F"], s["G"], s["H"], s["P"]

    crank_radius = float(np.linalg.norm(B - A))
    ax.add_patch(Circle((A[0], A[1]), crank_radius, fill=False, edgecolor="black", lw=0.85, zorder=0))

    ax.plot(*zip(A, B), color=_COLORS["crank"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(B, C), color=_COLORS["coupler"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(C, D), color=_COLORS["rocker"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(B, E), color=_COLORS["coupler"], lw=lw * 0.92, solid_capstyle="round")
    ax.plot(*zip(C, E), color=_COLORS["coupler"], lw=lw * 0.92, solid_capstyle="round")

    ax.plot(*zip(G, F), color=_COLORS["secondary_dyad"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(E, F), color=_COLORS["secondary_plate"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(E, H), color=_COLORS["secondary_plate"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(F, H), color=_COLORS["secondary_plate"], lw=lw, solid_capstyle="round")
    ax.plot(*zip(H, P), color=_COLORS["rod"], lw=lw, solid_capstyle="round")

    gs = 0.40 * width
    draw_ground_symbol(ax, A, gs, lw=lw * 0.72)
    draw_ground_symbol(ax, D, gs, lw=lw * 0.72)
    draw_ground_symbol(ax, G, gs, lw=lw * 0.72)

    jr = 0.075 if width < 7 else 0.085
    for pt in (A, B, C, D, E, F, G, H):
        draw_joint(ax, pt, r=jr, lw=lw * 0.42)


