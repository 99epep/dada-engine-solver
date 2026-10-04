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

# Adapted from CoolProp 8.0.0 (MIT), see THIRD_PARTY_NOTICES.md.
# Only rho -> 0 terms: Arp/McCarty/Friend viscosity, Hands/Arp conductivity.
def helium_property(t, conductivity):
    if conductivity:
        exponent=3.739232544/t-26.20316969/t**2+59.82252246/t**3-49.26397634/t**4
        return 2.7870034e-3*t**.7034007057*math.exp(exponent)
    if t<=100:
        x=math.log(t)
        exponent=-.135311743/x+1.00347841+1.20654649*x-.149564551*x*x+.012520841*x*x*x
        return math.exp(exponent)*1e-7
    return 196*t**.71938*math.exp(12.451/t-295.67/t/t-4.1249)*1e-7

def gas_constant(species):
    return (287.05,296.803,208.132,2077.1)[species]


def viscosity(t,species):
    if species==3: return helium_property(t,False)
    if species==0:
        # Lemmon/Jacobsen collision integral, Air.json (CoolProp MIT).
        x=math.log(t/103.3)
        integral=math.exp(.431-.4623*x+.08406*x**2+.005341*x**3-.00331*x**4)
        return 2.66958e-8*math.sqrt(28.9586*t)/(.36**2*integral)
    mu,s=((1.716e-5,111.),(1.663e-5,107.),(2.125e-5,114.))[species]
    return mu*(t/273.)**1.5*(273.+s)/(t+s)


def conductivity(t,species):
    if species==3: return helium_property(t,True)
    if species==0:
        # Lemmon/Jacobsen eta0_and_poly; EOS reducing temperature, not Tc.
        tau=132.6312/t
        return .001308*viscosity(t,0)*1e6+.001405*tau**(-1.1)-.001036*tau**(-.3)
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


def bennett_mean_nusselt(reynolds, prandtl, diameter_over_length):
    """Author HeatLib constant-T average branch; domain checked by the caller.

    Bennett (2020a), DOI 10.1115/1.4047834. Equations and applicability are
    transcribed in docs/MICROTUBE_GAS_MODEL.md; no external runtime dependency.
    """
    if reynolds == 0: return 3.66
    z=1/(reynolds*prandtl*diameter_over_length)
    offset=(3.66-6.54)/4.35
    exponent=(3.66+41.0)/13.3
    leveque=.40377*(64/z)**(1/3)
    graetz=(leveque**exponent+(3.66-offset)**exponent)**(1/exponent)+offset
    g=1.10*(1+.140/prandtl**(2/3))**(3/4)
    modified=.40377*(5.312/(g*math.sqrt(prandtl*z)*z))**(1/3)
    return graetz/math.tanh(leveque/modified*(1+.565*(prandtl*z)**(1/3)))


def shah_entry_excess(xplus):
    """Cumulative excess over Poiseuille; Shah-London Eq.192 (Fanning)."""
    if xplus == 0: return 0.
    # Algebraically equivalent form, stable at both zero and large xplus.
    if xplus <= 1:
        return (1.25*xplus*xplus+.00021*(13.76*math.sqrt(xplus)-64*xplus))/(xplus*xplus+.00021)
    return (1.25+.00021*(13.76/xplus**1.5-64/xplus))/(1+.00021/xplus**2)


def shah_apparent_darcy(reynolds, xplus):
    return (64+shah_entry_excess(xplus)/xplus)/reynolds


def shah_segment_pressure_loss(flow, diameter, area, mu, rho, x1, x2):
    if flow == 0: return 0.
    reynolds=abs(flow)*diameter/(area*mu)
    excess=shah_entry_excess(x2/(diameter*reynolds))-shah_entry_excess(x1/(diameter*reynolds))
    return excess*flow*flow/(2*rho*area*area)


def turbulent_darcy(reynolds, relative_roughness=0.):
    """Existing Haaland closure, shared with the public exchanger model."""
    return (-1.8*math.log10((relative_roughness/3.7)**1.11+6.9/reynolds))**-2


def turbulent_nusselt(reynolds, prandtl, friction):
    numerator=(friction/8.0)*(reynolds-1000.0)*prandtl
    denominator=1.0+12.7*math.sqrt(friction/8.0)*(prandtl**(2.0/3.0)-1.0)
    return numerator/denominator


def transition_friction(reynolds, slip_factor=1., entrance_darcy=0.):
    low=64/(2300*slip_factor)+entrance_darcy
    return low+(turbulent_darcy(4000)-low)*(reynolds-2300)/1700


def transition_heat(reynolds, prandtl, diameter_over_length, thermal_entry=True):
    low=bennett_mean_nusselt(2300,prandtl,diameter_over_length) if thermal_entry else 3.66
    high=turbulent_nusselt(4000,prandtl,turbulent_darcy(4000))
    return low+(high-low)*(reynolds-2300)/1700


def tube_network_residual(flow, diameter, area, mu, linear, quadratic, multiplier, length, rho, dp, axial_start=0., entrance=False):
    reynolds=flow*diameter/(area*mu)
    if reynolds<2300:
        loss=linear*flow
        if entrance:
            loss+=shah_segment_pressure_loss(flow,diameter,area,mu,rho,axial_start,axial_start+length)
    else:
        if reynolds<4000:
            entrance_darcy=0.
            if entrance:
                # Match the cumulative laminar segment excess at the endpoint;
                # never evaluate Shah at a transitional Reynolds number.
                excess=shah_entry_excess((axial_start+length)/(diameter*2300))-shah_entry_excess(axial_start/(diameter*2300))
                entrance_darcy=excess*diameter/(multiplier*length)
            friction=transition_friction(reynolds,1.,entrance_darcy)
        else:
            friction=turbulent_darcy(reynolds)
        loss=multiplier*friction*length/diameter*flow*flow/(2*rho*area*area)
    return loss+quadratic*flow*flow-dp


def continuum_network_flow(upper, diameter, area, mu, linear, quadratic, multiplier, length, rho, dp, axial_start=0., entrance=False):
    """Bracketed solve of the existing no-slip network, no domain extrapolation.

    The public reference uses scipy.brentq with xtol=1e-15 and its default
    rtol=4*epsilon. Bisection resolves the bracket to floating-point precision; failure
    returns to Python for the authoritative domain exception. Transition matches
    the complete laminar segment endpoint, including entrance excess.
    """
    high=min(upper,5e6*area*mu/diameter)
    low=0.
    # The analytic upper bound can undershoot by one rounding unit at
    # nearly equal pressures, where the entrance correction is negligible.
    upper_residual=tube_network_residual(high,diameter,area,mu,linear,quadratic,multiplier,length,rho,dp,axial_start,entrance)
    if abs(upper_residual)<=8*2.220446049250313e-16*dp:
        return True,high
    if tube_network_residual(high,diameter,area,mu,linear,quadratic,multiplier,length,rho,dp,axial_start,entrance)<0:
        return False,0.
    for _ in range(100):
        middle=(low+high)/2
        residual=tube_network_residual(middle,diameter,area,mu,linear,quadratic,multiplier,length,rho,dp,axial_start,entrance)
        if residual==0 or middle==low or middle==high:
            return True,middle
        if residual<0: low=middle
        else: high=middle
    return False,0.


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
    if 2300<=re<4000 and .5<=pr<=2000: return 3
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
    if kind in (0,3) and pr>500: flags|=16
    if kind==0 and re>0:
        if length/(d*re*pr)<=1e-6: flags|=256
    # Transition interpolates endpoint closures; retain both entry guards.
    if kind==3 and length<.05*2300*d: flags|=32
    if kind in (1,3) and length<10*d: flags|=64
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
    nu=bennett_mean_nusselt(re,pr,d/length) if thermal_entry else 3.66
    return valid,nu


def continuum_diagnostics(flow,p1,p2,t,d,length,area,maximum_mach,maximum_drop,thermal_entry,species):
    """Available no-slip regimes only; invalid states retain Python diagnostics."""
    mu=viscosity(t,species);k=conductivity(t,species);cp=transport_cp(t,species)
    r=gas_constant(species)
    rho,u,re,pr,speed,ma,gz,ratio,drop=flow_numbers(flow,p1,p2,t,d,length,area,r,mu,k,cp)
    kn=mean_free_path(mu,min(p1,p2),r,t)/d
    if domain_flags(re,pr,ma,kn,drop,length,d,maximum_mach,maximum_drop)!=0:
        return False,0.
    kind=thermal_kind(re,pr)
    if kind==0: nu=bennett_mean_nusselt(re,pr,d/length) if thermal_entry else 3.66
    elif kind==3: nu=transition_heat(re,pr,d/length,thermal_entry)
    elif kind==1: nu=turbulent_nusselt(re,pr,turbulent_darcy(re))
    else: return False,0.
    return True,nu


def directed_flow(pin,pout,t,p,r,gamma,reverse=False):
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
    half=int(p[13]) if len(p)>13 else 0
    if reverse: half=1-half
    if flow>0:
        ok,flow=continuum_network_flow(flow,d,area,mu,linear,quadratic,multiplier,length/2,rho,dp,half*length/2,True)
        if not ok: return False,0.
    # All domain failures (including report-only states) still use the public
    # Python path; no guard or authoritative rejection message is bypassed.
    valid,_=continuum_diagnostics(flow,pin,pout,t,d,length,area,p[7],p[8],p[9],species)
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
    return wall_balance_kernel(values, volume_rates, gas, links, walls, sources,
        destinations, one_way, ports, t, pressures, None)


def wall_balance_kernel(values,volume_rates,gas,links,walls,sources,destinations,one_way,ports,t,pressures,enthalpies):
    """Common conservative balances after backend-specific state reconstruction."""
    result=np.zeros(15)
    r,cv,cp,gamma,omega=gas
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
        ok,flow=directed_flow(pressures[a],pressures[b],t[a],links[link],r,gamma,reverse)
        if not ok: return False,result
        flows[link]=-flow if reverse else flow
        accumulate_transfer(mass,energy,a,b,flow,enthalpy(cp,t[a]) if enthalpies is None else enthalpies[a])
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
            ok,nu=continuum_diagnostics(flows[link],p1,p2,t[j],w[5],w[6],w[7],w[8],w[9],w[10],species)
            if not ok: return False,result
            conductance+=film_port_conductance(w[4],nu,k,w[5])
        overall=series_conductance(conductance,w[3])
        air_heat,gas_heat,wall_rate=wall_heat_rates(w[2],w[1],tw,overall,t[j])
        result[2*j+1]+=gas_heat
        result[8+side]=wall_rate
        result[10+side]=air_heat;result[12+side]=gas_heat
    result[14]=apply_piston_work(energy,pressures,volume_rates)
    return True,result/omega
