"""SI, prismatic horizontal jump; hydrostatic far sections, alpha=beta=1."""
from dataclasses import dataclass, asdict
from math import sqrt, tanh, isfinite
from scipy.optimize import brentq

G = 9.80665
VERSION = '1.0.0'
HEC14 = 'https://www.fhwa.dot.gov/engineering/hydraulics/pubs/06086/hec14.pdf'


def positive(name, value, zero=False):
    if not isinstance(value, (int, float)) or not isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(f'{name}: требуется конечное положительное число (ноль допускается только в специально отмеченных полях).')


@dataclass(frozen=True)
class Channel:
    width: float
    section: str = 'rectangular'
    side_slope: float = 0.0

    def __post_init__(self):
        positive('Ширина дна', self.width)
        positive('Откос z', self.side_slope, zero=True)
        if self.section not in ('rectangular', 'trapezoidal'):
            raise ValueError('Неизвестная форма сечения.')

    @property
    def z(self):
        return self.side_slope if self.section == 'trapezoidal' else 0.0

    def area(self, y):
        return y * (self.width + self.z*y)

    def top(self, y):
        return self.width + 2*self.z*y

    def perimeter(self, y):
        return self.width + 2*y*sqrt(1+self.z**2)

    def radius(self, y):
        return self.area(y)/self.perimeter(y)

    def momentum(self, q, y):
        return q*q/(G*self.area(y)) + self.width*y*y/2 + self.z*y**3/3

    def fr2(self, q, y):
        return q*q*self.top(y)/(G*self.area(y)**3)

    def energy(self, q, y):
        return y+(q/self.area(y))**2/(2*G)

    def friction(self, q, y, n):
        return (n*q/(self.area(y)*self.radius(y)**(2/3)))**2

    def critical(self, q):
        positive('Q',q)
        scale = (q*q/(G*self.width**2))**(1/3)
        return scale*brentq(lambda t:self.fr2(q,scale*t)-1,1e-12,2,xtol=1e-13)


def normal_depth(q, channel, n, slope):
    positive('Q',q)
    positive('n',n)
    positive('Уклон',slope)
    hi=max(channel.critical(q),1.)
    while channel.friction(q,hi,n)>slope:
        hi*=2
        if hi>1e8:
            raise ValueError('Нормальная глубина вне вычислительного диапазона.')
    return brentq(lambda y:channel.friction(q,y,n)-slope,hi*1e-12,hi)


@dataclass(frozen=True)
class JumpResult:
    q: float
    y1: float
    y2: float
    v1: float
    v2: float
    froude1: float
    froude2: float
    critical_depth: float
    energy_loss: float
    specific_energy1: float
    specific_energy2: float
    jump_length: float | None
    regime: str
    momentum1: float
    momentum2: float
    momentum_residual: float
    dissipation_percent: float
    power_kw: float
    tailwater_ratio: float | None
    warnings: tuple[str,...]

    def to_dict(self):
        return asdict(self)


def calculate_jump(q,y1,width,*,section='rectangular',side_slope=0.,tailwater=None):
    positive('Q',q)
    positive('y1',y1)
    if tailwater is not None:
        positive('Глубина нижнего бьефа',tailwater)
    c=Channel(width,section,side_slope)
    yc=c.critical(q)
    fr1=sqrt(c.fr2(q,y1))
    if fr1<=1+1e-6:
        raise ValueError('Fr1 ≤ 1 (или слишком близко к 1): сверхкритический вход для прыжка отсутствует.')
    m1=c.momentum(q,y1)
    if c.z==0:
        y2=y1*4*fr1**2/(sqrt(1+8*fr1**2)+1)
    else:
        hi=max(2*yc,2*y1)
        while c.momentum(q,hi)<m1:
            hi*=2
        y2=yc*brentq(lambda t:(c.momentum(q,t*yc)-m1)/m1,1,hi/yc,xtol=1e-13)
    m2=c.momentum(q,y2)
    e1,e2=c.energy(q,y1),c.energy(q,y2)
    loss=(y2-y1)**3/(4*y1*y2) if c.z==0 else e1-e2
    if loss < -1e-9*e1 or y2<=yc:
        raise ValueError('Решение не прошло физическую проверку энергии/глубины.')
    labels=[(1.7,'Волнистый переход'),(2.5,'Слабый прыжок'),(4.5,'Колеблющийся прыжок'),(9,'Устойчивый прыжок'),(float('inf'),'Сильный прыжок')]
    regime=next(name for upper,name in labels if fr1<upper)
    warnings=[]
    # HEC-14 eq.7.25: rectangular free jump; application guard is 1.7..17.
    length=220*y1*tanh((fr1-1)/22) if c.z==0 and 1.7<=fr1<=17 else None
    if length is None:
        warnings.append('Длина не вычислена: HEC-14 (7.25) применяется здесь только к прямоугольнику при 1,7 ≤ Fr1 ≤ 17. Для проверки размещения нужна подтверждённая длина.')
    if fr1<1.7:
        warnings.append('Волнистый переход: модель локализованного турбулентного прыжка ограничена.')
    if 2.5<=fr1<4.5:
        warnings.append('Колеблющийся прыжок: возможны волны и неустойчивость положения.')
    ratio=tailwater/y2 if tailwater is not None else None
    if ratio is not None:
        if ratio<.99:
            warnings.append('В заданном сечении нижний бьеф ниже сопряжённой глубины. Для определения положения рассчитайте профили.')
        elif ratio>1.01:
            warnings.append('Нижний бьеф выше сопряжённой глубины. Это сравнение глубин, а не расчёт затопленного прыжка.')
    return JumpResult(q,y1,y2,q/c.area(y1),q/c.area(y2),fr1,sqrt(c.fr2(q,y2)),yc,max(0,loss),e1,e2,length,regime,m1,m2,abs(m2-m1)/m1,max(0,loss)/e1*100,G*q*max(0,loss),ratio,tuple(warnings))
