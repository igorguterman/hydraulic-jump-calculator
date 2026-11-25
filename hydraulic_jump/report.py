"""Portable, escaped HTML reports; print to PDF using any browser."""
from datetime import datetime, timezone
from html import escape
import hashlib
import json
import platform
from importlib.metadata import version
from .core import VERSION, HEC14

ASSUMPTIONS = (
    'Установившийся одномерный поток воды; SI; g=9,80665 м/с²; плотность 1000 кг/м³.',
    'Призматическое русло, горизонтальное дно в зоне прыжка; гидростатика в удалённых сечениях; α=β=1.',
    'Трение в самом прыжке не учитывается; на плавно изменяющихся участках используется Маннинг.',
    'Положение — точечный разрыв по равенству импульсов. Конечная длина проверяется отдельно, валик не разрешается.',
    'HEC-14 (7.25) для длины — эмпирика прямоугольного русла; программный диапазон 1,7–17 по Fr. Это не универсальный диапазон стандарта.',
    'Запас длины означает только геометрическое размещение. Не проверены размыв, прочность, кавитация, воздухововлечение, волны и высота стенок.',
    'Гасители, пороги, расширения, затопленный/нестационарный прыжок и наклонное дно требуют отдельных методик.',
)
LABELS = {
 'width':'Ширина дна, м','section':'Сечение','side_slope':'Откос z, горизонталь/вертикаль',
 'tailwater':'Глубина на выходе, м (0 — неизвестна)','length':'Длина бассейна, м','n':'n Маннинга, с/м⅓',
 'supplied_length':'Внешняя длина прыжка, м (0 — автоматически)','length_source':'Источник внешней длины',
 'q':'Расход Q, м³/с','y1':'Глубина до y₁, м','y2':'Сопряжённая глубина y₂, м',
 'v1':'Скорость до, м/с','v2':'Скорость после, м/с','froude1':'Fr₁','froude2':'Fr₂',
 'critical_depth':'Критическая глубина, м','energy_loss':'Потеря удельной энергии, м',
 'specific_energy1':'Удельная энергия до, м','specific_energy2':'Удельная энергия после, м',
 'jump_length':'Длина HEC-14 (7.25), м','regime':'Форма прыжка','momentum1':'M₁, м³','momentum2':'M₂, м³',
 'momentum_residual':'Относительная невязка импульса','dissipation_percent':'Доля рассеянной энергии, %',
 'power_kw':'Рассеиваемая мощность при ρ=1000, кВт','tailwater_ratio':'yнб/y₂','warnings':'Замечания',
}


def result_rows(result):
    return [{'Показатель':LABELS.get(k,k),'Значение':v if isinstance(v,str) else ('—' if v is None else f'{v:.8g}')}
            for k,v in result.to_dict().items() if k!='warnings']


def payload(inputs,result,basin=None,profile=None,convergence=None):
    canonical=json.dumps(inputs,sort_keys=True,ensure_ascii=False,allow_nan=False)
    return {'version':VERSION,'software':{'python':platform.python_version(),'scipy':version('scipy'),'numpy':version('numpy'),'streamlit':version('streamlit')},'calculation_id':hashlib.sha256(canonical.encode()).hexdigest()[:12],
            'created_utc':datetime.now(timezone.utc).isoformat(),'inputs':inputs,
            'jump':result.to_dict(),'basin':basin,'profile':profile,'convergence':convergence,
            'assumptions':ASSUMPTIONS,'sources':[HEC14,
            'https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.0/overview-of-optional-capabilities/mixed-flow-regime-calculations']}


def html_report(data):
    def table(mapping):
        return '<table>'+''.join('<tr><th>'+escape(LABELS.get(k,k))+'</th><td>'+escape(str(v))+'</td></tr>' for k,v in mapping.items())+'</table>'
    profile=data.get('profile') or []
    svg=''
    if profile:
        columns=['Сверхкритическая ветвь, м','Подпорная ветвь, м']
        max_y=max((r[k] or 0) for r in profile for k in columns)*1.1
        max_x=max(r['x, м'] for r in profile)
        svg='<svg viewBox="0 0 800 260" role="img" aria-label="Расчётные ветви глубины">'
        svg+='<path d="M40 10 V230 H780" fill="none" stroke="black"/>'
        for k,color in zip(columns,['#148','orange']):
            points=' '.join(f'{40+730*r["x, м"]/max_x:.2f},{230-210*r[k]/max_y:.2f}' for r in profile if r[k] is not None)
            svg+=f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2"/>'
        svg+='<text x="350" y="255">x, м; глубина от горизонтального дна</text></svg>'
    body='<h1>Расчёт гидравлического прыжка</h1>'
    body+=f'<p>Версия {escape(data["version"])} · ID {escape(data["calculation_id"])} · UTC {escape(data["created_utc"])}</p>'
    body+='<p>Библиотеки: '+escape(str(data['software']))+'</p>'
    body+='<h2>Исходные данные</h2>'+table(data['inputs'])
    body+='<h2>Прыжок</h2>'+table(data['jump'])
    if data.get('basin'):
        body+='<h2>Горизонтальный бассейн</h2>'+table(data['basin'])+svg
    if data.get('convergence'):
        body+='<h2>Численная проверка</h2>'+table(data['convergence'])
    body+='<h2>Формулы</h2><p>A=y(b+zy); T=b+2zy; P=b+2y√(1+z²); R=A/P.<br>M=Q²/(gA)+by²/2+zy³/3; M₁=M₂.<br>Fr²=Q²T/(gA³); E=y+V²/(2g); ΔE=E₁−E₂.<br>Sf=(nQ/(AR^(2/3)))²; dy/dx=(S₀−Sf)/(1−Fr²).<br>L=220y₁ tanh((Fr₁−1)/22), HEC-14 (7.25), только прямоугольник.</p>'
    body+='<h2>Допущения и ограничения</h2><ul>'+''.join('<li>'+escape(s)+'</li>' for s in data['assumptions'])+'</ul>'
    body+='<h2>Источники</h2>'+''.join('<p><a href="'+escape(s,quote=True)+'">'+escape(s)+'</a></p>' for s in data['sources'])
    return '<!doctype html><html lang="ru"><meta charset="utf-8"><title>Расчёт прыжка</title><style>body{font:15px Arial;max-width:950px;margin:32px auto;color:#172331}table{border-collapse:collapse;width:100%;margin-bottom:16px}td,th{border:1px solid #ccd;padding:7px;text-align:left;overflow-wrap:anywhere}th{width:45%}h2{margin-top:28px}svg{width:100%}@media print{body{margin:0}tr{break-inside:avoid}h2{break-after:avoid}}</style>'+body+'</html>'
