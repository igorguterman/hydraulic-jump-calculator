"""GVF branches, matched as a point jump on a horizontal basin.

dy/dx=(S0-Sf)/(1-Fr^2); Manning SI; alpha=beta=1.
Finite jump length is a separate footprint check (HEC-14 section 7.4.1).
"""
from dataclasses import dataclass
from math import isfinite
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from .core import Channel, calculate_jump, positive


class Profile:
    def __init__(self,q,channel,y0,x0,x1,n,slope=0,*,rtol=1e-8):
        positive('Q',q)
        positive('Начальная глубина',y0)
        positive('n',n,zero=True)
        if not all(isfinite(v) for v in (x0,x1,slope)) or x0==x1:
            raise ValueError('Недопустимые координаты профиля.')
        self.yc=channel.critical(q)
        branch=-1 if y0<self.yc else 1
        if abs(y0/self.yc-1)<1e-4:
            raise ValueError('Граничная глубина слишком близка к критической для GVF.')
        def ode(x,values):
            # Guard internal trial steps. Terminal events prevent publishing
            # continuations across critical depth or dry-bed singularities.
            u=min(max(float(values[0]),1e-8),1-1e-7) if branch<0 else max(float(values[0]),1+1e-7)
            y=u*self.yc
            return [(slope-channel.friction(q,y,n))/(1-channel.fr2(q,y))/self.yc]
        def critical(x,y):
            return branch*(y[0]-1)-1e-4
        def dry(x,y):
            return y[0]-1e-7
        critical.terminal=dry.terminal=True
        self.solution=solve_ivp(ode,(x0,x1),[y0/self.yc],dense_output=True,rtol=rtol,
                                atol=rtol*.01,max_step=abs(x1-x0)/200,events=(critical,dry))
        if not self.solution.success:
            raise ValueError('Профиль не сошёлся: '+self.solution.message)
        self.lo=min(x0,float(self.solution.t[-1]))
        self.hi=max(x0,float(self.solution.t[-1]))
        self.stopped=self.solution.status==1

    def depth(self,x):
        if np.any(np.asarray(x)<self.lo-1e-9) or np.any(np.asarray(x)>self.hi+1e-9):
            raise ValueError('Экстраполяция профиля запрещена.')
        return self.solution.sol(x)[0]*self.yc


@dataclass
class BasinResult:
    status: str
    message: str
    x_jump: float | None
    jump: object | None
    reserve: float | None
    contained: bool | None
    momentum_residual: float | None
    upstream: Profile
    downstream: Profile
    candidates: tuple[float,...]


def solve_basin(q,channel:Channel,y_in,y_tail,length,n,*,rtol=1e-8,supplied_jump_length=None):
    positive('Длина бассейна',length)
    positive('Глубина нижнего бьефа',y_tail)
    if supplied_jump_length is not None:
        positive('Заданная длина прыжка',supplied_jump_length)
    yc=channel.critical(q)
    if y_in>=yc or y_tail<=yc:
        raise ValueError('Нужны yвх < yкр < yнб: сверхкритический вход и подкритический выход.')
    up=Profile(q,channel,y_in,0,length,n,rtol=rtol)
    down=Profile(q,channel,y_tail,length,0,n,rtol=rtol)
    lo,hi=max(up.lo,down.lo),min(up.hi,down.hi)
    def empty(status,msg,roots=()):
        return BasinResult(status,msg,None,None,None,None,None,up,down,roots)
    if hi<=lo:
        return empty('no_overlap','Ветви не имеют общей области до критической глубины. Прыжок не локализован.')
    def residual(x):
        m1=channel.momentum(q,float(up.depth(x)))
        m2=channel.momentum(q,float(down.depth(x)))
        return (m2-m1)/max(m1,m2)
    xs=np.linspace(lo,hi,401)
    rs=[residual(x) for x in xs]
    if max(abs(r) for r in rs)<1e-9:
        return empty('indeterminate','Положение неопределённо: импульсы совпадают на всём участке, например при n=0.')
    roots=[]
    for i in range(len(xs)-1):
        if abs(rs[i])<1e-10:
            roots.append(float(xs[i]))
        elif rs[i]*rs[i+1]<0:
            roots.append(brentq(residual,float(xs[i]),float(xs[i+1]),xtol=1e-9))
    if abs(rs[-1])<1e-10:
        roots.append(float(xs[-1]))
    roots=tuple(sorted({round(x,8) for x in roots}))
    if len(roots)>1:
        return empty('multiple','Несколько корней. Автоматический выбор не выполнен.',roots)
    if not roots:
        if min(rs)>0:
            msg='Импульс подпорной ветви выше на общей области: возможно смещение к входу/затопление. Свободный прыжок внутри не найден.'
        else:
            msg='Импульс подпорной ветви ниже: свободный прыжок внутри не найден, возможно вынесен за выход.'
        return empty('outside',msg)
    x=roots[0]
    jump=calculate_jump(q,float(up.depth(x)),channel.width,section=channel.section,side_slope=channel.z)
    lj=supplied_jump_length if supplied_jump_length is not None else jump.jump_length
    reserve=length-x-lj if lj is not None else None
    return BasinResult('found','Положение найдено по равенству импульсов. Длина прыжка проверяется отдельно.',x,jump,reserve,
                       reserve>=0 if reserve is not None else None,abs(residual(x)),up,down,roots)


def basin_table(result,length,points=201):
    rows=[]
    for x in np.linspace(0,length,points):
        yu=float(result.upstream.depth(x)) if result.upstream.lo<=x<=result.upstream.hi else None
        yd=float(result.downstream.depth(x)) if result.downstream.lo<=x<=result.downstream.hi else None
        actual=(yu if x<result.x_jump else yd) if result.x_jump is not None else None
        rows.append({'x, м':float(x),'Сверхкритическая ветвь, м':yu,'Подпорная ветвь, м':yd,'Составной профиль, м':actual})
    return rows
