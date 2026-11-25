import json
from dataclasses import asdict
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st
from hydraulic_jump.core import Channel, calculate_jump, normal_depth, VERSION, HEC14
from hydraulic_jump.profiles import solve_basin, basin_table
from hydraulic_jump.report import payload, html_report, result_rows, ASSUMPTIONS

st.set_page_config(page_title='Гидравлический прыжок',page_icon='🌊',layout='wide')
st.title('Гидравлический прыжок и горизонтальный бассейн')
st.caption(f'Версия {VERSION} · SI · расчёт по импульсу и одномерным профилям потока')
DEFAULT={'q':10.,'y1':.45,'width':3.,'section':'rectangular','side_slope':1.5,'tailwater':1.8,'length':30.,'n':.015,'supplied_length':0.,'length_source':''}
with st.sidebar:
    st.header('Расчётный проект')
    uploaded=st.file_uploader('Загрузить сохранённый проект JSON',type=['json'])
    if uploaded and st.button('Применить проект'):
        try:
            raw=json.loads(uploaded.getvalue())
            inp=raw.get('inputs',raw)
            if not isinstance(inp,dict) or set(inp)-set(DEFAULT):
                raise ValueError('Неизвестные поля проекта.')
            config={**DEFAULT,**inp}
            for key,default in DEFAULT.items():
                value=config[key]
                if isinstance(default,float) and (type(value) not in (int,float) or not np.isfinite(value) or value<0 or value>1e6):
                    raise ValueError(f'Недопустимое значение {key}.')
                if isinstance(default,str) and not isinstance(value,str):
                    raise ValueError(f'Недопустимый тип {key}.')
            Channel(config['width'],config['section'],config['side_slope'])
            if config['q']<=0 or config['y1']<=0 or config['length']<=0:
                raise ValueError('Q, y1, length должны быть положительными.')
            if any(config[k]<.000001 for k in ('q','y1','width')) or config['length']<.001 or config['n']>1:
                raise ValueError('Значение за пределами полей интерфейса.')
            for key,value in config.items():
                st.session_state[key]=value
            st.success('Исходные данные загружены.')
        except (ValueError,TypeError,OverflowError) as exc:
            st.error(f'Проект не применён: {exc}')
    for key,value in DEFAULT.items():
        st.session_state.setdefault(key,value)
    q=st.number_input('Расход Q, м³/с',min_value=.000001,max_value=1e6,key='q',format='%.6f')
    y1=st.number_input('Входная глубина y₁, м',min_value=.000001,max_value=1e6,key='y1',format='%.6f')
    section=st.selectbox('Сечение',['rectangular','trapezoidal'],key='section',format_func=lambda s:'Прямоугольное' if s=='rectangular' else 'Трапецеидальное')
    width=st.number_input('Ширина дна b, м',min_value=.000001,max_value=1e6,key='width')
    z=st.number_input('Откос z, горизонталь / вертикаль',min_value=0.,max_value=1e6,key='side_slope',disabled=section=='rectangular')
    tailwater=st.number_input('Глубина на выходе yнб, м (0 — неизвестна)',min_value=0.,max_value=1e6,key='tailwater')
    length=st.number_input('Длина горизонтального бассейна B, м',min_value=.001,max_value=1e6,key='length')
    n=st.number_input('Шероховатость Маннинга n, с/м⅓',min_value=0.,max_value=1.,key='n',format='%.5f',step=.001)
    with st.expander('Подтверждённая длина из другого источника'):
        supplied=st.number_input('Длина прыжка L, м (0 — автоматически)',min_value=0.,max_value=1e6,key='supplied_length')
        length_source=st.text_input('Источник заданной длины / обоснование',key='length_source')
    st.caption('Все глубины отсчитываются от одного горизонтального днища. Значение n=0 — идеализация без трения.')
inputs={k:st.session_state[k] for k in DEFAULT}
channel=Channel(width,section,z)
try:
    result=calculate_jump(q,y1,width,section=section,side_slope=z,tailwater=tailwater or None)
except (ValueError,OverflowError,ZeroDivisionError) as exc:
    st.error(str(exc))
    st.info('Проверьте данные. Для прыжка нужен сверхкритический вход. Измените расход, глубину или сечение.')
    st.download_button('Сохранить исходные данные',json.dumps({'inputs':inputs},ensure_ascii=False,indent=2),file_name='jump_project.json')
    st.stop()

tabs=st.tabs(['Прыжок','Бассейн и профили','Сценарии','Методики и отчёт'])
with tabs[0]:
    cols=st.columns(4)
    for col,label,value in zip(cols,['Fr₁','Сопряжённая глубина y₂','Потеря энергии','Рассеяние'],[f'{result.froude1:.4f}',f'{result.y2:.4f} м',f'{result.energy_loss:.4f} м',f'{result.dissipation_percent:.2f} %']):
        col.metric(label,value)
    st.write(f'**{result.regime}** · невязка импульса {result.momentum_residual:.2e}')
    st.dataframe(pd.DataFrame(result_rows(result)),hide_index=True,use_container_width=True)
    for warning in result.warnings:
        st.warning(warning)
    if result.jump_length is not None:
        st.caption('Длина рассчитана по HEC-14 (7.25); это эмпирическая зависимость, а не точное разрешение валика.')
    with st.expander('Графики удельной энергии и импульса'):
        ys=np.geomspace(max(y1*.4,1e-6),result.y2*1.5,250)
        fig=make_subplots(rows=1,cols=2,subplot_titles=('Удельная энергия E, м','Функция импульса M, м³'))
        for i,fn in enumerate([channel.energy,channel.momentum],start=1):
            fig.add_trace(go.Scatter(x=[fn(q,y) for y in ys],y=ys,name='E(y)' if i==1 else 'M(y)'),row=1,col=i)
            fig.add_trace(go.Scatter(x=[fn(q,y1),fn(q,result.y2)],y=[y1,result.y2],mode='markers',name='Сечения 1 / 2'),row=1,col=i)
        fig.update_yaxes(title_text='Глубина, м')
        st.plotly_chart(fig,use_container_width=True)
    with st.expander('Вспомогательный расчёт нормальной глубины подходящего русла'):
        st.caption('Для отдельного равномерного участка с положительным уклоном. Уклон не применяется к горизонтальному бассейну.')
        slope=st.number_input('Уклон подходящего участка S₀, м/м',min_value=.000001,value=.001,format='%.6f')
        if n>0:
            yn=normal_depth(q,channel,n,slope)
            st.write(f'Нормальная глубина: {yn:.6f} м · Fr={channel.fr2(q,yn)**.5:.4f}')
        else:
            st.info('При n=0 нормальная глубина по Маннингу не определяется.')

basin=None
profile_rows=None
basin_summary=None
convergence=None
with tabs[1]:
    st.write('Свободный прыжок в горизонтальном призматическом бассейне без порогов и гасителей. Профили интегрируются от обоих граничных сечений.')
    if tailwater<=0:
        st.info('Задайте глубину на выходе для поиска положения прыжка.')
    elif supplied>0 and not length_source.strip():
        st.error('Для внешней длины укажите источник или обоснование. Проверка бассейна не выполнена.')
    else:
        try:
            basin=solve_basin(q,channel,y1,tailwater,length,n,supplied_jump_length=supplied or None)
            basin_summary={'status':basin.status,'Сообщение':basin.message,'x прыжка от входа, м':basin.x_jump,'Запас B−x−L, м':basin.reserve,'Размещение по длине':basin.contained,'Невязка импульса':basin.momentum_residual,'Длина из внешнего источника, м':supplied or None,'Источник длины':length_source if supplied else 'HEC-14 (7.25), если применимо','Параметры в месте прыжка':basin.jump.to_dict() if basin.jump else None}
            st.info(basin.message)
            if basin.x_jump is not None:
                st.metric('Начало прыжка от входа',f'{basin.x_jump:.4f} м')
                if basin.contained is True:
                    st.success(f'Проверка длины выполнена: запас {basin.reserve:.3f} м. Это проверка размещения, а не всех условий безопасности сооружения.')
                elif basin.contained is False:
                    st.error(f'Прыжок не помещается по длине: дефицит {-basin.reserve:.3f} м.')
                else:
                    st.warning('Положение найдено, но размещение не проверено: нет применимой длины прыжка.')
                st.write('**Параметры непосредственно перед и после найденного прыжка:**')
                st.dataframe(pd.DataFrame(result_rows(basin.jump)),hide_index=True,use_container_width=True)
                for warning in basin.jump.warnings:
                    st.warning(warning)
                fine=solve_basin(q,channel,y1,tailwater,length,n,rtol=1e-10,supplied_jump_length=supplied or None)
                delta=abs(fine.x_jump-basin.x_jump) if fine.x_jump is not None else None
                convergence={'rtol основной':1e-8,'rtol уточнённый':1e-10,'Изменение x, м':delta,'Статус уточнённого расчёта':fine.status}
                st.caption(f'Уточнение точности: изменение координаты {delta:.3g} м.' if delta is not None else 'Уточнённый расчёт не подтвердил положение — требуется исследование.')
            if basin.upstream.stopped or basin.downstream.stopped:
                st.warning('Одна из ветвей остановлена у критической глубины. Экстраполяция за неё не выполняется.')
            profile_rows=basin_table(basin,length)
            frame=pd.DataFrame(profile_rows)
            fig=go.Figure()
            for col in ['Сверхкритическая ветвь, м','Подпорная ветвь, м']:
                fig.add_trace(go.Scatter(x=frame['x, м'],y=frame[col],name=col,line={'dash':'dash'},connectgaps=False))
            if basin.x_jump is not None:
                xj=basin.x_jump
                left=frame[frame['x, м']<xj]
                right=frame[frame['x, м']>xj]
                fig.add_trace(go.Scatter(x=list(left['x, м'])+[xj,xj]+list(right['x, м']),y=list(left['Сверхкритическая ветвь, м'])+[basin.jump.y1,basin.jump.y2]+list(right['Подпорная ветвь, м']),name='Составной профиль (точечный разрыв)',line={'width':4}))
                lj=supplied or basin.jump.jump_length
                if lj is not None:
                    fig.add_vrect(x0=xj,x1=xj+lj,fillcolor='orange',opacity=.12,line_width=0,annotation_text='Проверяемая длина L')
            fig.add_trace(go.Scatter(x=[0,length],y=[0,0],name='Дно',line={'color':'#765'}))
            fig.update_layout(xaxis_title='x от входа, м',yaxis_title='Глубина от дна, м',height=470)
            st.plotly_chart(fig,use_container_width=True)
            st.caption('Штриховые линии — отдельные расчётные ветви, не одновременно существующие поверхности. Вертикаль — точечная модель прыжка. Оранжевая зона — его эмпирическая протяжённость.')
            st.download_button('Профили CSV',frame.to_csv(index=False).encode('utf-8-sig'),file_name='profiles.csv',mime='text/csv')
        except (ValueError,OverflowError,ZeroDivisionError) as exc:
            st.error(f'Расчёт бассейна не выполнен: {exc}')
            basin_summary={'status':'error','Сообщение':str(exc)}

with tabs[2]:
    st.write('В каждой строке задавайте свой расход, входную глубину и нижний бьеф. Геометрия, n и длина бассейна берутся из боковой панели.')
    st.caption('Начальные строки используют одинаковые глубины только как шаблон. Замените их данными режимов сооружения.')
    data=pd.DataFrame([{'Q, м³/с':v,'yвх, м':y1,'yнб, м':tailwater} for v in [q*.75,q,q*1.25]])
    edited=st.data_editor(data,num_rows='dynamic',hide_index=True,key='scenarios',use_container_width=True)
    signature=json.dumps(inputs,sort_keys=True)+edited.to_json()
    if st.button('Рассчитать сценарии',type='primary'):
        rows=[]
        for _,row in edited.iterrows():
            rec=row.to_dict()
            try:
                ri=calculate_jump(float(row.iloc[0]),float(row.iloc[1]),width,section=section,side_slope=z,tailwater=float(row.iloc[2]) or None)
                rec.update({'Fr₁':ri.froude1,'y₂, м':ri.y2,'ΔE, м':ri.energy_loss,'L, м':ri.jump_length,'Статус':'Прыжок рассчитан'})
                if row.iloc[2]>ri.critical_depth:
                    if supplied>0 and not length_source.strip():
                        raise ValueError('Не указан источник внешней длины.')
                    rb=solve_basin(float(row.iloc[0]),channel,float(row.iloc[1]),float(row.iloc[2]),length,n,supplied_jump_length=supplied or None)
                    rec.update({'x, м':rb.x_jump,'Запас, м':rb.reserve,'Статус':rb.status,'Комментарий':rb.message})
            except (ValueError,OverflowError,ZeroDivisionError) as exc:
                rec['Статус']=str(exc)
            rows.append(rec)
        st.session_state['batch_result']=(signature,rows)
    if st.session_state.get('batch_result',('',[]))[0]==signature:
        out=pd.DataFrame(st.session_state['batch_result'][1])
        st.dataframe(out,hide_index=True,use_container_width=True)
        st.download_button('Сценарии CSV',out.to_csv(index=False).encode('utf-8-sig'),file_name='scenarios.csv',mime='text/csv')

with tabs[3]:
    st.subheader('Методики и границы применения')
    for note in ASSUMPTIONS:
        st.write('• '+note)
    st.latex(r'M=\frac{Q^2}{gA}+\frac{by^2}{2}+\frac{zy^3}{3},\quad M_1=M_2')
    st.latex(r'Fr^2=\frac{Q^2T}{gA^3},\quad E=y+\frac{V^2}{2g},\quad \frac{dy}{dx}=\frac{S_0-S_f}{1-Fr^2}')
    st.latex(r'L=220y_1\tanh\left(\frac{Fr_1-1}{22}\right)')
    st.markdown(f'[FHWA HEC-14, 2006]({HEC14}): импульс — §6.2, формула (6.1); длина прямоугольного прыжка — §7.4.1, (7.25).')
    st.markdown('[USACE HEC-RAS: смешанные режимы и поиск прыжка по импульсу](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.0/overview-of-optional-capabilities/mixed-flow-regime-calculations)')
    st.caption('Программные пределы полей — ограничения ввода, не подтверждение физической применимости. Порог сравнения глубин ±1% — только индикатор, не нормативный запас.')
    st.subheader('Сохранение проекта и отчёт')
    report=payload(inputs,result,basin_summary,profile_rows,convergence)
    st.write(f'Идентификатор расчёта: **{report["calculation_id"]}** · версия {VERSION}')
    st.download_button('Сохранить проект и результаты JSON',json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False),file_name='jump_project.json',mime='application/json')
    st.download_button('Скачать расчётный отчёт HTML',html_report(report).encode('utf-8'),file_name='jump_report.html',mime='text/html')
    st.caption('Откройте HTML-отчёт в браузере → Ctrl+P → «Сохранить как PDF». Отчёт содержит исходные данные, результаты, формулы, источники и ограничения.')
