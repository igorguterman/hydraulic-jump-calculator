import io
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from hydraulic_jump import calculate_jump

st.set_page_config(page_title="Гидравлический прыжок", page_icon="🌊", layout="wide")
st.title("Расчёт гидравлического прыжка")
st.caption("Предварительный инженерный расчёт по сохранению количества движения")

with st.sidebar:
    st.header("Исходные данные")
    q = st.number_input("Расход Q, м³/с", min_value=0.0001, value=10.0, step=1.0)
    y1 = st.number_input("Глубина до прыжка y₁, м", min_value=0.0001, value=0.45, step=0.05)
    section = st.selectbox("Сечение", ["rectangular", "trapezoidal"], format_func=lambda x: "Прямоугольное" if x == "rectangular" else "Трапецеидальное")
    width = st.number_input("Ширина дна / канала, м", min_value=0.0001, value=3.0, step=0.1)
    z = st.number_input("Откос z (горизонталь:вертикаль)", min_value=0.0, value=1.5, step=0.1, disabled=section == "rectangular")
    tailwater = st.number_input("Глубина нижнего бьефа, м (0 — не задана)", min_value=0.0, value=0.0, step=0.05)

try:
    result = calculate_jump(q, y1, width, section=section, side_slope=z, tailwater=tailwater or None)
except ValueError as exc:
    st.error(str(exc))
    st.stop()

st.subheader("Результаты")
cols = st.columns(4)
cols[0].metric("Fr₁", f"{result.froude1:.3f}")
cols[1].metric("Сопряжённая глубина y₂", f"{result.y2:.3f} м")
cols[2].metric("Потеря энергии", f"{result.energy_loss:.3f} м")
cols[3].metric("Оценка длины прыжка", f"{result.jump_length:.2f} м")

st.write(f"**Режим:** {result.regime}")
st.write(f"**Рекомендация по бассейну:** {result.basin_recommendation}")
st.table({
    "Показатель": ["Скорость до прыжка", "Скорость после прыжка", "Fr₂", "Критическая глубина", "Энергия до", "Энергия после"],
    "Значение": [f"{result.v1:.3f} м/с", f"{result.v2:.3f} м/с", f"{result.froude2:.3f}", f"{result.critical_depth:.3f} м", f"{result.specific_energy1:.3f} м", f"{result.specific_energy2:.3f} м"],
})

for warning in result.warnings:
    st.warning(warning)

st.info("Результат предназначен для предварительной проверки. Для рабочего проекта нужны проверка по действующему стандарту, расчёт профиля и при сложной геометрии — 2D/CFD или физическая модель.")

st.subheader("Схематический профиль прыжка")
x1 = result.jump_length * 0.25
x2 = result.jump_length * 0.75
fig = go.Figure()
fig.add_trace(go.Scatter(x=[0, x1, x2, result.jump_length], y=[result.y1, result.y1, result.y2, result.y2], mode="lines", name="Свободная поверхность", line={"width": 4}))
fig.add_trace(go.Scatter(x=[0, result.jump_length], y=[0, 0], mode="lines", name="Днище", line={"color": "saddlebrown"}))
fig.update_layout(xaxis_title="Координата, м", yaxis_title="Глубина, м", height=380, margin={"l": 20, "r": 20, "t": 20, "b": 20})
st.plotly_chart(fig, use_container_width=True)

st.subheader("Проверка по диапазону расходов")
q_min, q_max, q_count = st.columns(3)
q_min = q_min.number_input("Q от, м³/с", min_value=0.0001, value=max(q / 2, 0.0001), step=1.0)
q_max = q_max.number_input("Q до, м³/с", min_value=0.0001, value=q * 1.5, step=1.0)
q_count = q_count.number_input("Количество точек", min_value=2, max_value=50, value=5, step=1)
if q_max <= q_min:
    st.warning("Q до должно быть больше Q от.")
else:
    rows = []
    for qi in [q_min + i * (q_max - q_min) / (q_count - 1) for i in range(q_count)]:
        try:
            ri = calculate_jump(qi, y1, width, section=section, side_slope=z, tailwater=tailwater or None)
            rows.append({"Q, м³/с": qi, "Fr₁": ri.froude1, "y₂, м": ri.y2, "Потери энергии, м": ri.energy_loss, "Длина прыжка, м": ri.jump_length})
        except ValueError:
            rows.append({"Q, м³/с": qi, "Fr₁": None, "y₂, м": None, "Потери энергии, м": None, "Длина прыжка, м": None})
    sweep = pd.DataFrame(rows)
    st.dataframe(sweep, hide_index=True, use_container_width=True)

report = f"""РАСЧЁТ ГИДРАВЛИЧЕСКОГО ПРЫЖКА\n\nQ = {result.q:.4f} м³/с\ny₁ = {result.y1:.4f} м\ny₂ = {result.y2:.4f} м\nFr₁ = {result.froude1:.4f}\nFr₂ = {result.froude2:.4f}\nПотеря энергии = {result.energy_loss:.4f} м\nОценка длины прыжка = {result.jump_length:.4f} м\nРежим: {result.regime}\nРекомендация: {result.basin_recommendation}\n"""
st.download_button("Скачать текстовый отчёт", data=io.BytesIO(report.encode("utf-8")), file_name="hydraulic_jump_report.txt", mime="text/plain")
