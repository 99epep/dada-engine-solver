"""Freeze a CoolProp 8 oracle; optional development dependency, never solver code."""
import argparse
import json
from pathlib import Path
import CoolProp
from CoolProp.CoolProp import PropsSI

TEMPERATURES = dict(Air=[100,125,150,175,199.9,200,250,300,500,1000],
                    Helium=[50,75,100,150,199.9,200,300,500,1000])


def generate():
    if CoolProp.__version__.split('.')[0]!='8':
        raise ValueError('Generate the reference using CoolProp 8.x.')
    rows=[]
    for fluid,temperatures in TEMPERATURES.items():
        for t in temperatures:
            properties={}
            sensitivity={}
            for field,output in [('viscosity_pa_s','V'),('conductivity_w_m_k','L'),('cp_j_kg_k','C')]:
                value=PropsSI(output,'T',t,'Dmass',1e-10,fluid)
                smaller=PropsSI(output,'T',t,'Dmass',1e-11,fluid)
                properties[field]=value
                sensitivity[field]=abs(value/smaller-1)
            rows.append(dict(species=fluid.lower(),temperature_k=t,**properties,
                             density_reduction_relative_change=sensitivity))
    return dict(schema_version=1,oracle='CoolProp',version=CoolProp.__version__,
        git_revision=CoolProp.__gitrevision__,density_kg_m3=1e-10,
        density_check_kg_m3=1e-11,scope='Dilute limit oracle; not a DADA EOS or phase validation.',rows=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('tests/data/coolprop8_dilute_transport.json'))
    args=parser.parse_args()
    args.output.write_text(json.dumps(generate(),indent=2,allow_nan=False)+'\n')
