"""Deterministic center-first Sobol scheduling in the unchanged global space."""
from .candidate import canonical_json
from .strategy import SobolStrategy
from .parameters import ChoiceParameter


class ScheduledSobol:
    def __init__(self, space, settings, *, index=0):
        self.space, self.settings, self._index = space, settings, index
        self.local = settings['domain'] == 'local_regions_v1'
        self.regions = settings['regions'] if self.local else [dict(id='global',center={p.name:p.initial for p in space.parameters})]
        self.centers = []
        known = {}
        for region in self.regions:
            key = canonical_json(region['center'])
            if key in known:
                self.centers[known[key]]['region_ids'].append(region['id'])
            else:
                known[key] = len(self.centers)
                self.centers.append(dict(physical=region['center'],region_ids=[region['id']]))
        samples = max(0,index-len(self.centers))
        n = len(self.regions)
        self.engines = [SobolStrategy(len(space.parameters),seed=settings['seed'],scramble=settings['scramble'],
            index=samples//n+(i<samples%n)) for i in range(n)]
        self.last_physical = self.last_origin = None

    @property
    def index(self): return self._index

    def next_point(self):
        if self.index < len(self.centers):
            center = self.centers[self.index]
            physical = dict(center['physical'])
            origin = dict(kind='center' if self.local else 'initial',region_id=center['region_ids'][0],
                          region_ids=center['region_ids'],region_sequence_index=None)
        else:
            i = (self.index-len(self.centers)) % len(self.regions)
            region,engine = self.regions[i],self.engines[i]
            origin = dict(kind='sobol',region_id=region['id'],region_sequence_index=engine.index)
            u = engine.next_point()
            if self.local:
                radius = self.settings['radius_fraction']
                center = self.space.encode(region['center'])
                u = tuple(x if isinstance(p,ChoiceParameter) and self.settings.get('choice_scope')=='declared_choices'
                          else max(0.,z-radius)+x*(min(1.,z+radius)-max(0.,z-radius))
                          for p,z,x in zip(self.space.parameters,center,u))
            physical = self.space.decode(u)
        # Canonical physical encoding makes integer duplicates (and identical
        # overlaps across regions) share identity; center floats stay exact.
        self.last_physical, self.last_origin = physical, origin
        self._index += 1
        return self.space.encode(physical)

    def state(self):
        samples = max(0,self.index-len(self.centers))
        return dict(type='scheduled_sobol_v1',index=self.index,
            centers_completed=min(self.index,len(self.centers)),
            next_region=samples%len(self.regions),
            regions=[dict(id=r['id'],index=e.index) for r,e in zip(self.regions,self.engines)])
