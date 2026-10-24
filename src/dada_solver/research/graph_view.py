"""Browse and export existing Research report curves in one Matplotlib figure."""
from pathlib import Path
from .plot_data import PLOT_NAMES

GRAPH_NAMES=tuple(name for name in PLOT_NAMES if name!='mechanisms')


def draw_graph(ax, records, name):
    """Draw report samples directly, without resampling or scientific evaluation."""
    ax.clear();messages=[];domains=set();units=set();count=0
    for record in records:
        graph=record.get('plots',{}).get(name,{})
        if graph.get('angle') is not None:
            domains.add(graph.get('angle_domain','solver_cycle_angle'))
            units.add(graph.get('unit',graph.get('volume_unit','')))
    if len(domains)>1 or len(units)>1:
        raise ValueError('Cannot overlay curves with different angle conventions or units.')
    for record in records:
        graph=record.get('plots',{}).get(name,{})
        label=record['candidate_id'][:12]
        if graph.get('unavailable') or not graph.get('angle'):
            messages.append(label+': '+graph.get('unavailable','No recorded curve available.'));continue
        series=graph.get('series')
        if series is None and name=='volumes':
            series=[dict(label=side,values=graph[side]) for side in ('small','large')]
        for row in series or []:
            ax.plot(graph['angle'],row['values'],label=label+' / '+row['label'])
            count+=1
        if not series: messages.append(label+': No series available.')
    ax.set_title('Indicated gas-pressure torques' if name=='torques' else name.replace('_',' ').title())
    ax.set_xlabel('Cycle angle [deg]');ax.set_ylabel(next(iter(units),''))
    ax.grid(True,alpha=.25)
    if count: ax.legend(fontsize='small')
    if messages:
        ax.text(.02,.98,'\n'.join(messages),transform=ax.transAxes,va='top',fontsize='small',wrap=True)
    return messages


class GraphBrowser:
    """One window, cached curves and native Matplotlib view history."""
    def __init__(self,data,*,names=GRAPH_NAMES,save_directory=None,svg_vector_conduits=True):
        import matplotlib.pyplot as plt
        from matplotlib.axes import Axes
        browser=self
        class GraphAxes(Axes):
            def _get_view(self):
                return browser.index,super()._get_view()

            def _set_view(self,view):
                index,limits=view
                if browser.index!=index:
                    browser.index=index
                    browser.show_current(redraw=False)
                super()._set_view(limits)
        self.records=[row for row in data.get('selected',data['records']) if 'plots' in row]
        self.names=tuple(name for name in names for _ in (self.records if name=='mechanisms' else (None,)))
        self.machine_records={}
        self.machines={}
        record_index=0
        for index,name in enumerate(self.names):
            if name=='mechanisms':
                self.machine_records[index]=self.records[record_index%len(self.records)]
                record_index+=1
        self.frame_index=0;self.paused=False;self.animation=None
        self.svg_vector_conduits=svg_vector_conduits
        if not self.records: raise ValueError('No candidates with plot data; choose an available candidate.')
        if not self.names or any(name not in PLOT_NAMES for name in self.names):
            raise ValueError('Choose report plots, including mechanisms for live machine animation.')
        self.index=0;self.save_directory=Path(save_directory) if save_directory is not None else None
        self.figure=plt.figure(figsize=(11,7))
        self._curve_facecolor=self.figure.get_facecolor()
        self.ax=self.figure.add_subplot(axes_class=GraphAxes)
        self.figure.subplots_adjust(bottom=.15,top=.88)
        self.status=self.figure.text(.5,.035,'',ha='center',fontsize='small')
        self.heading=self.figure.suptitle(data.get('name','Research candidate curves'))
        self.figure._research_graph_browser=self
        if self.save_directory is not None:
            self.save_directory.mkdir(parents=True,exist_ok=True)
            self.figure.canvas.get_default_filename=lambda: str(self.save_directory.resolve()/(self.current_graph+'.png'))
        toolbar=getattr(self.figure.canvas.manager,'toolbar',None)
        if toolbar is not None: self.initialize_toolbar(toolbar)
        else: self.show_current()
        self.figure.canvas.mpl_connect('key_press_event',self.on_key)
        if self.machine_records:
            from matplotlib.animation import FuncAnimation
            from .machine_render import FRAME_COUNT, FRAME_DURATION_MS
            self.animation=FuncAnimation(self.figure,self.update_frame,frames=FRAME_COUNT,
                interval=FRAME_DURATION_MS,blit=False,cache_frame_data=False)
            self.figure.canvas.mpl_connect('draw_event',lambda _:self._sync_timer())
            self.figure.canvas.mpl_connect('close_event',lambda _:self.animation.event_source.stop() if self.animation.event_source is not None else None)
        native_save=self.figure.savefig
        def savefig(filename,*args,**kwargs):
            kind=kwargs.get('format')
            if kind is None and isinstance(filename,(str,Path)): kind=Path(filename).suffix[1:]
            if self.current_graph=='mechanisms' and kind and kind.lower()=='svg':
                from .machine_animation import animated_machine_svg
                payload=animated_machine_svg(self._machine(),vector_conduits=self.svg_vector_conduits)
                if hasattr(filename,'write'): filename.write(payload)
                else: Path(filename).write_bytes(payload)
                return
            return native_save(filename,*args,**kwargs)
        self.figure.savefig=savefig


    def initialize_toolbar(self,toolbar):
        """Seed native Back/Forward history with graph-aware Axes views."""
        self.toolbar=toolbar
        for index in range(len(self.names)):
            self.index=index;self.show_current(redraw=False);toolbar.push_current()
        for _ in range(len(self.names)-1): toolbar.back()
        self.figure.canvas.draw_idle()

    @property
    def current_graph(self): return self.names[self.index]

    def show_current(self,*,redraw=True):
        machine=self.current_graph=='mechanisms'
        self.heading.set_visible(not machine);self.status.set_visible(not machine)
        if self.current_graph=='mechanisms':
            try:
                from .machine_animation import frame_machine_axes
                frame_machine_axes(self.ax,self._machine()['geometry'])
                self._draw_machine()
            except (ValueError,KeyError,RuntimeError,ImportError) as error:
                self.ax.clear();self.ax.axis('off')
                self.ax.text(.02,.98,f'Machine animation unavailable: {error}',transform=self.ax.transAxes,va='top',wrap=True)
        else:
            self.figure.set_facecolor(self._curve_facecolor)
            self.figure.set_size_inches(11,7,forward=True)
            self.figure.subplots_adjust(left=.125,right=.9,bottom=.15,top=.88)
            draw_graph(self.ax,self.records,self.current_graph)
        methods=list(dict.fromkeys(row.get('plots',{}).get(self.current_graph,{}).get('method','') for row in self.records))
        self.status.set_text(f'{self.index+1} / {len(self.names)} — '+ '\n'.join(m for m in methods if m))
        self._sync_timer()
        if redraw: self.figure.canvas.draw_idle()

    def _machine(self):
        from .machine_animation import prepare_animation
        if self.index not in self.machines:
            payload=self.machine_records[self.index]['plots'].get('mechanisms',{})
            if payload.get('unavailable'): raise ValueError(payload['unavailable'])
            if payload.get('kind')!='machine_animation': raise ValueError('Live drawing requires vector report inputs, not an encoded WebP.')
            self.machines[self.index]=prepare_animation(payload)
        return self.machines[self.index]

    def _draw_machine(self):
        from .machine_animation import draw_animation_frame
        draw_animation_frame(self.ax,self._machine(),self.frame_index)

    def _sync_timer(self):
        if self.animation is None or self.animation.event_source is None: return
        if self.paused or self.current_graph!='mechanisms' or self.index not in self.machines: self.animation.event_source.stop()
        else: self.animation.event_source.start()

    def update_frame(self,index):
        self.frame_index=index
        if self.current_graph=='mechanisms' and not self.paused:
            limits=self.ax.get_xlim(),self.ax.get_ylim()
            try:
                self._draw_machine()
                self.ax.set_xlim(limits[0]);self.ax.set_ylim(limits[1])
            except (ValueError,KeyError,RuntimeError,ImportError): self._sync_timer()
        return []

    def next(self):
        self.index=(self.index+1)%len(self.names);self.show_current()
        if getattr(self,'toolbar',None) is not None: self.toolbar.push_current()

    def previous(self):
        self.index=(self.index-1)%len(self.names);self.show_current()
        if getattr(self,'toolbar',None) is not None: self.toolbar.push_current()

    def on_key(self,event):
        # Interactive backends already route arrows/Home through their toolbar.
        if event.key==' ':
            self.paused=not self.paused;self._sync_timer();return
        if event.key=='pagedown': self.next();return
        if event.key=='pageup': self.previous();return
        if getattr(self,'toolbar',None) is not None: return
        if event.key=='right': self.next()
        elif event.key=='left': self.previous()
        elif event.key in ('home','end'):
            self.index=0 if event.key=='home' else len(self.names)-1;self.show_current()

    def save_current(self,path):
        """Export only the graph; refuse to overwrite an existing file."""
        import matplotlib.pyplot as plt
        path=Path(path)
        if path.suffix.lower() not in ('.png','.svg','.pdf'):
            raise ValueError('Choose a .png, .svg or .pdf output.')
        if path.exists(): raise ValueError('Graph output already exists; choose a new filename.')
        if self.current_graph=='mechanisms' and path.suffix.lower()=='.svg':
            from .machine_animation import animated_machine_svg
            payload=animated_machine_svg(self._machine(),vector_conduits=self.svg_vector_conduits)
            path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('xb') as stream: stream.write(payload)
            return path
        figure,axes=plt.subplots(figsize=(11,6))
        try:
            if self.current_graph=='mechanisms':
                from .machine_animation import draw_animation_frame, frame_machine_axes
                frame_machine_axes(axes,self._machine()['geometry'])
                draw_animation_frame(axes,self._machine(),self.frame_index)
            else: draw_graph(axes,self.records,self.current_graph)
            axes.set_xlim(self.ax.get_xlim());axes.set_ylim(self.ax.get_ylim())
            if self.current_graph!='mechanisms': figure.tight_layout()
            path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('xb') as stream:
                figure.savefig(stream,format=path.suffix[1:].lower(),dpi=180)
        finally: plt.close(figure)
        return path



def browse_graphs(data,**kwargs):
    browser=GraphBrowser(data,**kwargs)
    return browser.figure,browser
