"""Authoritative numerical primitives shared by Python and optional Numba.

No optional imports or configuration objects here. Public model methods retain
validation and error semantics. Keep the compiled kernel and its transitive
numerical helpers in this file: cache identity hashes the complete source.
"""
import math
import numpy as np

SPECIES = ('air', 'nitrogen', 'argon', 'helium')

def molar_cp(t, oxygen=False):
    if oxygen:
        a,b,c,d,e=((31.32234,-20.23531,57.86644,-36.50624,-.007374) if t<700 else
                   (30.03235,8.772972,-3.988133,.788313,-.741599))
    else:
        a,b,c,d,e=((28.98641,1.853978,-9.647459,16.63537,.000117) if t<500 else
                   (19.50583,19.88705,-8.598535,1.369784,.527601))
    x=t/1000
    return a+b*x+c*x*x+d*x**3+e/x**2

def helium_property(t, conductivity):
    temperatures=(200,225,250,273.15,275,300,325,350,375,400,450,500,600,700,800,900,1000)
    values=((118.5,128.2,137.7,146.2,146.9,155.9,164.6,173.2,181.6,189.9,206.,221.7,251.9,280.9,308.9,336.1,362.6)
            if conductivity else (15.1,16.4,17.6,18.7,18.8,19.9,21.,22.1,23.2,24.3,26.3,28.3,32.2,35.9,39.5,43.,46.4))
    return float(np.interp(t,np.array(temperatures),np.array(values)))*(1e-3 if conductivity else 1e-6)

def gas_constant(species):
    return (287.05,296.803,208.132,2077.1)[species]


def viscosity(t,species):
    if species==3: return helium_property(t,False)
    mu,s=((1.716e-5,111.),(1.663e-5,107.),(2.125e-5,114.))[species]
    return mu*(t/273.)**1.5*(273.+s)/(t+s)


def conductivity(t,species):
    if species==3: return helium_property(t,True)
    k,s=((.0241,194.),(.0242,150.),(.0163,170.))[species]
    return k*(t/273.)**1.5*(273.+s)/(t+s)


def transport_cp(t,species):
    if species>=2: return 2.5*gas_constant(species)
    if species==1: return molar_cp(t)/.0280134
    return (.79*molar_cp(t)+.21*molar_cp(t,True))/(.79*.0280134+.21*.0319988)


def temperature(m,u,cv):
    return u/(m*cv)


def pressure(m,r,t,v):
    return m*r*t/v


def requires_slip(kn):
    return kn>=.001


def mean_free_path(mu,p,r,t):
    return mu/p*math.sqrt(math.pi*r*t/2)


def poiseuille(p1,p2,t,mu,length,diameter,count,r):
    return count*math.pi*diameter**4*(p1-p2)*(p1+p2)/(256*mu*length*r*t)


def hausen(graetz):
    return 3.66+.0668*graetz/(1+.04*graetz**(2/3))


def minor_loss_coefficient(header,rho,area,cda):
    quadratic=header/(4*rho*area**2)
    if cda>0: quadratic+=1/(2*rho*cda**2)
    return quadratic


def laminar_network_flow(dp,linear,quadratic):
    return 2*dp/(linear+math.sqrt(linear**2+4*quadratic*dp))


def orifice_flow(pin,pout,t,r,gamma,area):
    ratio=pout/pin
    if ratio<=(2.0/(gamma+1.0))**(gamma/(gamma-1.0)):
        factor=math.sqrt(gamma/(r*t))
        factor*=(2.0/(gamma+1.0))**((gamma+1.0)/(2.0*(gamma-1.0)))
        return area*pin*factor,True
    radicand=2.0*gamma/(r*t*(gamma-1.0))*(ratio**(2.0/gamma)-ratio**((gamma+1.0)/gamma))
    return area*pin*math.sqrt(max(0.0,radicand)),False


def flow_numbers(flow,p1,p2,t,d,length,area,r,mu,k,cp):
    rho=min(p1,p2)/(r*t)
    u=abs(flow)/(rho*area);re=abs(flow)*d/(area*mu);pr=cp*mu/k
    speed=math.sqrt(cp/(cp-r)*r*t);ma=u/speed
    gz=re*pr*d/length;ratio=max(p1,p2)/min(p1,p2)
    compressibility=2*(ratio-1)/(ratio+1)
    return rho,u,re,pr,speed,ma,gz,ratio,compressibility


def thermal_kind(re,pr):
    if re<2300: return 0
    if 4000<=re<=5e6 and .5<=pr<=2000: return 1
    return 2


def domain_flags(re,pr,ma,kn,drop,length,d,max_ma,max_drop):
    flags=0
    if ma>max_ma: flags|=1
    if requires_slip(kn): flags|=2
    if kn>.1: flags|=4
    if drop>max_drop: flags|=8
    if not .5<=pr<=2000: flags|=16
    kind=thermal_kind(re,pr)
    if kind==0 and length<.05*re*d: flags|=32
    if kind==1 and length<10*d: flags|=64
    if kind==2: flags|=128
    return flags


def film_port_conductance(area,nu,k,d):
    return .5*area*nu*k/d


def series_conductance(conductance,resistance):
    return 1/(1/conductance+resistance)


def wall_heat_rates(effective,inlet,wall_temperature,conductance,gas_temperature):
    air_heat=effective*(inlet-wall_temperature)
    gas_heat=conductance*(wall_temperature-gas_temperature)
    return air_heat,gas_heat,air_heat-gas_heat


def enthalpy(cp,t):
    return cp*t


def accumulate_transfer(mass,energy,source,destination,flow,specific_enthalpy):
    enthalpy_flow=flow*specific_enthalpy
    mass[source]-=flow;mass[destination]+=flow
    energy[source]-=enthalpy_flow;energy[destination]+=enthalpy_flow


def apply_piston_work(energy,pressures,volume_rates):
    energy[0]-=pressures[0]*volume_rates[0]
    energy[1]-=pressures[1]*volume_rates[1]
    return pressures[0]*volume_rates[0]+pressures[1]*volume_rates[1]



def laminar_diagnostics(flow,p1,p2,t,d,length,area,maximum_mach,maximum_drop,thermal_entry,species):
    mu=viscosity(t,species);k=conductivity(t,species);cp=transport_cp(t,species)
    r=gas_constant(species)
    rho,u,re,pr,speed,ma,gz,ratio,drop=flow_numbers(flow,p1,p2,t,d,length,area,r,mu,k,cp)
    kn=mean_free_path(mu,min(p1,p2),r,t)/d
    valid=domain_flags(re,pr,ma,kn,drop,length,d,maximum_mach,maximum_drop)==0 and thermal_kind(re,pr)==0
    nu=hausen(gz) if thermal_entry else 3.66
    return valid,nu


def directed_flow(pin,pout,t,p,r,gamma):
    # p: diameter, length, count, area, multiplier, header K, valve CdA,
    # maximum Mach, maximum pressure drop, thermal-entry flag, Tmin, Tmax.
    if pout>=pin:
        return True,0.
    d,length,count,area,multiplier,header,cda=p[:7]
    if not p[10]<=t<=p[11]:
        return False,0.
    species=int(p[12])
    mu=viscosity(t,species)
    kn=mean_free_path(mu,(pin+pout)/2,gas_constant(species),t)/d
    if requires_slip(kn):
        return False,0.
    rho=(pin+pout)/(2*r*t)
    nominal=poiseuille(pin,pout,t,mu,length/2,d,count,r)/multiplier
    linear=(pin-pout)/nominal
    quadratic=minor_loss_coefficient(header,rho,area,cda)
    dp=pin-pout
    flow=laminar_network_flow(dp,linear,quadratic)
    # Conservative compiled domain: fall back even for unsupported report-only
    # states. Python owns transition/turbulent/slip treatment and failure strings.
    valid,_=laminar_diagnostics(flow,pin,pout,t,d,length,area,p[7],p[8],p[9],species)
    if not valid:
        return False,0.
    cap_area=min(area,cda) if cda>0 else area
    cap,_=orifice_flow(pin,pout,t,r,gamma,cap_area)

    return True,min(flow,cap)


def wall_kernel(values,volumes,volume_rates,gas,links,walls,sources,destinations,one_way,ports):
    result=np.zeros(15)
    t=np.empty(4);pressures=np.empty(4)
    r,cv,cp,gamma,omega=gas
    for j in range(4):
        m,u=values[2*j],values[2*j+1]
        if not math.isfinite(m) or not math.isfinite(u) or m<=0 or u<=0 or not math.isfinite(volumes[j]) or volumes[j]<=0:
            return False,result
        t[j]=temperature(m,u,cv)
        pressures[j]=pressure(m,r,t[j],volumes[j])
    flows=np.empty(4)
    # Same transport accumulation order as ThermodynamicModel.assemble_rates.
    mass=result[:8:2];energy=result[1:8:2]
    for link in range(4):
        a,b=sources[link],destinations[link]
        reverse=pressures[a]<pressures[b]
        if reverse and one_way[link]:
            flows[link]=0.
            continue
        if reverse: a,b=b,a
        ok,flow=directed_flow(pressures[a],pressures[b],t[a],links[link],r,gamma)
        if not ok: return False,result
        flows[link]=-flow if reverse else flow
        accumulate_transfer(mass,energy,a,b,flow,enthalpy(cp,t[a]))
    for side in range(2):
        j=side+2
        # wall: capacity, air inlet, air conductance, film resistance, internal
        # area, diameter, length, flow area, Mach/drop guards, entry, Tmin/Tmax.
        w=walls[side]
        tw=values[8+side]/w[0]
        if not math.isfinite(tw) or tw<=0 or not w[11]<=t[j]<=w[12]:
            return False,result
        species=int(w[13])
        k=conductivity(t[j],species);conductance=0.
        for port in range(2):
            link=ports[side,port]
            a,b=sources[link],destinations[link]
            p1,p2=pressures[a],pressures[b]
            if flows[link]==0: p1,p2=pressures[j],pressures[j]
            ok,nu=laminar_diagnostics(flows[link],p1,p2,t[j],w[5],w[6],w[7],w[8],w[9],w[10],species)
            if not ok: return False,result
            conductance+=film_port_conductance(w[4],nu,k,w[5])
        overall=series_conductance(conductance,w[3])
        air_heat,gas_heat,wall_rate=wall_heat_rates(w[2],w[1],tw,overall,t[j])
        result[2*j+1]+=gas_heat
        result[8+side]=wall_rate
        result[10+side]=air_heat;result[12+side]=gas_heat
    result[14]=apply_piston_work(energy,pressures,volume_rates)
    return True,result/omega
