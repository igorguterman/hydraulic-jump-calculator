import math
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq
from hydraulic_jump.core import Channel, G, calculate_jump, normal_depth
from hydraulic_jump.profiles import Profile, solve_basin
from hydraulic_jump.report import payload, html_report


@pytest.mark.parametrize('fr',[1.001,1.2,1.7,2.5,4.5,9,17,40])
def test_belanger_and_energy_identity(fr):
    y,b=.7,4.
    q=fr*math.sqrt(G*y)*b*y
    r=calculate_jump(q,y,b)
    assert r.y2/y==pytest.approx((math.sqrt(1+8*fr*fr)-1)/2,rel=1e-10)
    assert r.energy_loss==pytest.approx(r.specific_energy1-r.specific_energy2,abs=1e-10)
    assert r.momentum_residual<1e-10
    assert r.froude2<1


@pytest.mark.parametrize('z',[0,.5,1.5,5])
def test_trapezoid_pressure_integral_and_conjugacy(z):
    r=calculate_jump(10,.45,3,section='trapezoidal',side_slope=z)
    def impulse(y):
        pressure=quad(lambda h:h*(3+2*z*(y-h)),0,y)[0]
        return 100/(G*y*(3+z*y))+pressure
    assert impulse(r.y1)==pytest.approx(impulse(r.y2),rel=1e-10)
    c=Channel(3,'trapezoidal',z)
    assert c.fr2(10,r.critical_depth)==pytest.approx(1,rel=1e-10)
    assert r.energy_loss>0
    if z:
        assert r.jump_length is None
    else:
        assert r.y2==pytest.approx(calculate_jump(10,.45,3).y2)


def test_published_fhwa_example():
    # HEC-14 pp.7-40/41, Q=14.2, b=4.3, y1=.375;
    # printed values rounded: V=8.81, Fr=4.6, yc=1.036, y2=2.255.
    r=calculate_jump(14.2,.375,4.3)
    assert r.v1==pytest.approx(8.81,abs=.005)
    assert r.froude1==pytest.approx(4.6,abs=.01)
    assert r.critical_depth==pytest.approx(1.036,abs=.0005)
    assert r.y2==pytest.approx(2.255,abs=.001)


@pytest.mark.parametrize('field,bad',[('q',0),('q',float('inf')),('q',float('nan')),('y1',-1),('width',0),('side_slope',-1),('tailwater',-1)])
def test_reject_invalid(field,bad):
    kwargs=dict(q=10,y1=.45,width=3,side_slope=1.5,tailwater=2)
    kwargs[field]=bad
    with pytest.raises(ValueError):
        calculate_jump(**kwargs)


def test_subcritical_rejected():
    with pytest.raises(ValueError):
        calculate_jump(1,3,4)


def test_normal_and_uniform_profile():
    c=Channel(3,'trapezoidal',1.5)
    yn=normal_depth(10,c,.02,.001)
    assert c.friction(10,yn,.02)==pytest.approx(.001,rel=1e-10)
    p=Profile(10,c,yn,0,50,.02,.001)
    assert p.depth(50)==pytest.approx(yn,rel=1e-9)


def test_horizontal_profile_against_independent_depth_integral():
    c=Channel(3)
    p=Profile(10,c,.45,0,10,.015)
    y=float(p.depth(10))
    def dxdy(y):
        a=3*y
        radius=a/(3+2*y)
        fr2=100/(G*9*y**3)
        sf=.015**2*100/(a*a*radius**(4/3))
        return (fr2-1)/sf
    distance=quad(dxdy,.45,y,epsabs=1e-10)[0]
    assert distance==pytest.approx(10,abs=1e-6)


def test_manufactured_jump_location_and_refinement():
    c=Channel(3)
    # Independent scalar quadrature yields both boundary depths for xj=12.
    def dxdy(y):
        a=3*y
        return (100/(G*9*y**3)-1)/(.015**2*100/(a*a*(a/(3+2*y))**(4/3)))
    y_at=brentq(lambda y:quad(dxdy,.45,y)[0]-12,.45001,.9)
    after=calculate_jump(10,y_at,3).y2
    tw=brentq(lambda y:quad(dxdy,after,y)[0]-18,1.2,after)
    r=solve_basin(10,c,.45,tw,30,.015)
    fine=solve_basin(10,c,.45,tw,30,.015,rtol=1e-10)
    assert r.status=='found'
    assert r.x_jump==pytest.approx(12,abs=1e-5)
    assert abs(r.x_jump-fine.x_jump)<1e-5
    assert r.momentum_residual<1e-9
    assert r.contained


def test_tailwater_change_and_no_jump_status():
    c=Channel(3)
    low=solve_basin(10,c,.45,1.8,30,.015)
    high=solve_basin(10,c,.45,1.9,30,.015)
    assert high.x_jump<low.x_jump
    assert solve_basin(10,c,.45,3.,30,.015).status=='outside'
    assert solve_basin(10,c,.45,1.2,30,.015).status=='outside'


def test_indeterminate_no_friction():
    r=calculate_jump(10,.45,3)
    assert solve_basin(10,Channel(3),.45,r.y2,30,0).status=='indeterminate'


def test_uncontained_and_missing_length():
    c=Channel(3)
    r=solve_basin(10,c,.45,1.8,30,.015,supplied_jump_length=100)
    assert r.contained is False
    assert r.reserve<0
    # Trapezoid with z=0 reproduces rectangle exactly.
    assert solve_basin(10,Channel(3,'trapezoidal',0),.45,1.8,30,.015).x_jump==r.x_jump


def test_profile_stops_at_critical_without_extrapolation():
    p=Profile(10,Channel(3),.45,0,1000,.03)
    assert p.stopped
    assert p.depth(p.hi)<p.yc
    with pytest.raises(ValueError):
        p.depth(1000)


def test_report_escapes_user_content():
    r=calculate_jump(10,.45,3)
    doc=html_report(payload({'length_source':'<script>alert(1)</script>'},r))
    assert '<script>' not in doc
    assert '&lt;script&gt;' in doc
    assert '1.0.0' in doc
