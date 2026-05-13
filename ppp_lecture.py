"""
PPP/CPP Lecture PDF Generator
Uses matplotlib mathtext for formula rendering + reportlab for document assembly.
Only uses matplotlib-supported math commands (no \bigl, \cot, \begin{cases}, \mathbb).
"""

import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle,
    PageBreak, HRFlowable,
)
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY

W, H = A4
MARGIN = 2.0 * cm


# ─── Formula renderer ──────────────────────────────────────────────────────────

def formula_block(tex, fontsize=14, dpi=200):
    """Render LaTeX via matplotlib mathtext → inline Image in a left-padded table."""
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.patch.set_alpha(0)
    text = fig.text(0, 0, f"${tex}$", fontsize=fontsize, color='black',
                    ha='left', va='baseline', fontfamily='DejaVu Sans')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bb = text.get_window_extent(renderer=renderer)
    pad = 6
    dpi_f = float(dpi)
    fig_w = (bb.width  + 2*pad) / dpi_f
    fig_h = (bb.height + 2*pad) / dpi_f
    if fig_w < 0.1: fig_w = 0.1
    if fig_h < 0.1: fig_h = 0.1
    fig.set_size_inches(fig_w, fig_h)
    text.set_position((pad / (dpi_f * fig_w),
                       pad / (dpi_f * fig_h)))
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi_f,
                bbox_inches='tight', pad_inches=pad/dpi_f,
                transparent=True)
    plt.close(fig)
    buf.seek(0)
    img = Image(buf)
    img.drawWidth  *= 0.72
    img.drawHeight *= 0.72
    indent = 1.2*cm
    avail  = W - 2*MARGIN - indent
    tbl = Table([[img]], colWidths=[avail])
    tbl.setStyle(TableStyle([
        ('ALIGN',        (0,0), (-1,-1), 'LEFT'),
        ('LEFTPADDING',  (0,0), (-1,-1), 0),
        ('RIGHTPADDING', (0,0), (-1,-1), 0),
        ('TOPPADDING',   (0,0), (-1,-1), 4),
        ('BOTTOMPADDING',(0,0), (-1,-1), 4),
    ]))
    return tbl


# ─── Styles ───────────────────────────────────────────────────────────────────

styles = getSampleStyleSheet()

title_style = ParagraphStyle('T', parent=styles['Title'],
    fontSize=22, leading=28, spaceAfter=6,
    textColor=colors.HexColor('#1a1a2e'))

h1_style = ParagraphStyle('H1', parent=styles['Heading1'],
    fontSize=16, leading=20, spaceBefore=14, spaceAfter=6,
    textColor=colors.HexColor('#16213e'))

h2_style = ParagraphStyle('H2', parent=styles['Heading2'],
    fontSize=12, leading=15, spaceBefore=9, spaceAfter=3,
    textColor=colors.HexColor('#0f3460'))

h3_style = ParagraphStyle('H3', parent=styles['Heading3'],
    fontSize=10.5, leading=13, spaceBefore=7, spaceAfter=2,
    textColor=colors.HexColor('#533483'))

body_style = ParagraphStyle('B', parent=styles['Normal'],
    fontSize=10, leading=14, spaceAfter=4, alignment=TA_JUSTIFY)

cap_style = ParagraphStyle('C', parent=styles['Normal'],
    fontSize=9, leading=11, spaceAfter=6,
    textColor=colors.HexColor('#555555'), alignment=TA_CENTER)

def P(text, style=None): return Paragraph(text, style or body_style)
def H1(text):            return Paragraph(text, h1_style)
def H2(text):            return Paragraph(text, h2_style)
def H3(text):            return Paragraph(text, h3_style)
def SP(n=6):             return Spacer(1, n)
def HR():                return HRFlowable(width='100%', thickness=0.5,
                                           color=colors.HexColor('#cccccc'), spaceAfter=4)


def var_table(rows, headers=None, col_widths=None):
    avail = W - 2*MARGIN
    if col_widths is None:
        col_widths = [3.8*cm, avail - 3.8*cm]
    data = ([headers] if headers else []) + list(rows)
    tbl = Table(data, colWidths=col_widths)
    st = [
        ('FONTSIZE',      (0,0), (-1,-1), 9),
        ('LEADING',       (0,0), (-1,-1), 12),
        ('TOPPADDING',    (0,0), (-1,-1), 2),
        ('BOTTOMPADDING', (0,0), (-1,-1), 2),
        ('VALIGN',        (0,0), (-1,-1), 'TOP'),
        ('GRID',          (0,0), (-1,-1), 0.3, colors.HexColor('#cccccc')),
        ('BACKGROUND',    (0,0), (-1,-1), colors.HexColor('#fafafa')),
        ('TEXTCOLOR',     (0,0), (0,-1),  colors.HexColor('#0f3460')),
        ('FONTNAME',      (0,0), (0,-1),  'Helvetica-Bold'),
    ]
    if headers:
        st += [
            ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#16213e')),
            ('TEXTCOLOR',  (0,0), (-1,0), colors.white),
            ('FONTNAME',   (0,0), (-1,0), 'Helvetica-Bold'),
        ]
    tbl.setStyle(TableStyle(st))
    return tbl


def wide_table(rows, headers=None, col_widths=None):
    avail = W - 2*MARGIN
    if col_widths is None:
        n = len((headers or rows[0]))
        col_widths = [avail/n]*n
    data = ([headers] if headers else []) + list(rows)
    tbl = Table(data, colWidths=col_widths)
    st = [
        ('FONTSIZE',      (0,0), (-1,-1), 9),
        ('LEADING',       (0,0), (-1,-1), 12),
        ('TOPPADDING',    (0,0), (-1,-1), 3),
        ('BOTTOMPADDING', (0,0), (-1,-1), 3),
        ('VALIGN',        (0,0), (-1,-1), 'MIDDLE'),
        ('GRID',          (0,0), (-1,-1), 0.3, colors.HexColor('#aaaaaa')),
        ('ROWBACKGROUNDS', (0,1), (-1,-1),
         [colors.HexColor('#f5f8ff'), colors.HexColor('#ffffff')]),
    ]
    if headers:
        st += [
            ('BACKGROUND', (0,0), (-1,0), colors.HexColor('#16213e')),
            ('TEXTCOLOR',  (0,0), (-1,0), colors.white),
            ('FONTNAME',   (0,0), (-1,0), 'Helvetica-Bold'),
            ('ALIGN',      (0,0), (-1,0), 'CENTER'),
        ]
    tbl.setStyle(TableStyle(st))
    return tbl


# ─── Content ──────────────────────────────────────────────────────────────────

def build_story():
    s = []

    # ── TITLE PAGE ────────────────────────────────────────────────────────────
    s += [SP(60)]
    s.append(Paragraph("Точное позиционирование PPP и CPP",
        ParagraphStyle('BT', parent=title_style, fontSize=26, leading=32,
                       alignment=TA_CENTER)))
    s += [SP(8)]
    s.append(Paragraph("Математическая лекция на основе кода RTKLIB",
        ParagraphStyle('BS', parent=body_style, fontSize=13, alignment=TA_CENTER,
                       textColor=colors.HexColor('#555555'))))
    s += [SP(12)]
    s.append(HRFlowable(width='60%', thickness=1.5,
                         color=colors.HexColor('#16213e'), spaceAfter=12))
    s.append(Paragraph("Precise Point Positioning &amp; Collaborative Precise Positioning",
        ParagraphStyle('BE', parent=body_style, fontSize=11, alignment=TA_CENTER,
                       textColor=colors.HexColor('#16213e'))))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 1. ИЗМЕРИТЕЛЬНАЯ МОДЕЛЬ
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("1. Измерительная модель ГНСС"))
    s.append(HR())

    s.append(H2("1.1 Уравнения наблюдений"))
    s.append(P("В основе PPP лежат два типа наблюдений: <b>псевдодальность</b> "
               "(code) и <b>фазовые измерения</b> (carrier phase). "
               "Для спутника s на частоте f_j:"))
    s.append(SP(4))
    s.append(P("<b>Псевдодальность (в метрах):</b>"))
    s.append(formula_block(
        r"P_j^s = \rho^s + c\,\delta t_r - c\,\delta t^s + T^s + I_j^s"
        r" + b_{r,j} - b_j^s + \varepsilon_{P_j}"))
    s.append(P("<b>Фазовое измерение (в метрах):</b>"))
    s.append(formula_block(
        r"L_j^s = \rho^s + c\,\delta t_r - c\,\delta t^s + T^s - I_j^s"
        r" + \lambda_j N_j^s + b_{r,j}^\phi - b_j^{s,\phi} + \varepsilon_{L_j}"))
    s.append(P("Ионосферная задержка I_j^s <b>положительна</b> для P и "
               "<b>отрицательна</b> для L. Неоднозначность присутствует только в фазе."))
    s += [SP(4)]
    s.append(var_table([
        [r'ρˢ',            'Геометрическое расстояние спутник–приёмник, м'],
        ['c',              'Скорость света, 299 792 458 м/с'],
        ['δt_r',           'Смещение часов приёмника, с'],
        ['δtˢ',            'Смещение часов спутника (из точных эфемерид), с'],
        ['Tˢ',             'Тропосферная задержка, м'],
        ['I_j^s',          'Ионосферная задержка на частоте f_j, м'],
        ['λ_j = c / f_j',  'Длина волны на частоте j, м'],
        ['N_j^s',          'Целочисленная фазовая неоднозначность, циклы'],
        ['b_{r,j}, b_j^s', 'Приборные задержки (DCB) приёмника и спутника'],
    ]))

    s += [SP(8)]
    s.append(H2("1.2 Геометрическое расстояние и релятивистские поправки"))
    s.append(formula_block(
        r"\rho^s = \|\mathbf{r}^s(t^s) - \mathbf{r}_r(t_r)\| + \delta\rho_\mathrm{Shapiro}"))
    s.append(P("<b>Поправка Шапиро</b> (задержка в гравитационном поле Земли):"))
    s.append(formula_block(
        r"\delta\rho_\mathrm{Shapiro} = \frac{2\mu}{c^2}"
        r"\ln\frac{r^s + r_r + \rho^s}{r^s + r_r - \rho^s}"))
    s.append(P("μ = GM = 3.986×10¹⁴ м³/с² — гравитационный параметр Земли."))
    s.append(P("<b>Поворот Земли</b> за время распространения сигнала (эффект Саньяка):"))
    s.append(formula_block(
        r"\mathbf{r}^s_\mathrm{corr} = R_z(\Omega_E \cdot \Delta t)\,\mathbf{r}^s(t^s)"
        r",\quad \Delta t = \rho^s/c"))
    s.append(P("Ω_E = 7.2921151467×10⁻⁵ рад/с — угловая скорость вращения Земли."))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 2. IFLC
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("2. Ионосферно-свободная линейная комбинация (IFLC)"))
    s.append(HR())

    s.append(H2("2.1 Зависимость ионосферной задержки от частоты"))
    s.append(P("Первичная ионосферная задержка обратно пропорциональна квадрату частоты:"))
    s.append(formula_block(r"I_j^s = \frac{f_1^2}{f_j^2}\,I_1^s"))
    s.append(P("Это позволяет устранить ионосферу линейной комбинацией. Коэффициенты:"))
    s.append(formula_block(
        r"\alpha = \frac{f_1^2}{f_1^2 - f_2^2},\qquad \beta = -\frac{f_2^2}{f_1^2 - f_2^2}"))

    s.append(H2("2.2 Ионосферно-свободные наблюдения"))
    s.append(formula_block(
        r"P_{IF} = \alpha P_1 + \beta P_2"
        r" = \rho^s + c\,\delta t_r - c\,\delta t^s + T^s + \varepsilon_{IF}"))
    s.append(formula_block(
        r"L_{IF} = \alpha L_1 + \beta L_2"
        r" = \rho^s + c\,\delta t_r - c\,\delta t^s + T^s + N_{IF}^s + \varepsilon_{IF}"))
    s.append(P("Нецелочисленная неоднозначность IF-комбинации:"))
    s.append(formula_block(
        r"N_{IF}^s = \alpha\,\lambda_1 N_1^s + \beta\,\lambda_2 N_2^s"))

    s.append(H2("2.3 Числовые значения для GPS"))
    s.append(var_table([
        ['f₁ = 1575.42 МГц', 'Несущая L1'],
        ['f₂ = 1227.60 МГц', 'Несущая L2'],
        ['α ≈ +2.546',       'Коэффициент для L1 в IFLC'],
        ['β ≈ −1.546',       'Коэффициент для L2 в IFLC'],
    ]))
    s += [SP(6)]
    s.append(P(
        "⚠ IFLC увеличивает дисперсию шума примерно в √(α²+β²) ≈ 3 раза. "
        "Поэтому в коде (<i>varerr()</i>) дисперсия умножается на <b>3.0</b> "
        "при режиме IONOOPT_IFLC."))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 3. ТРОПОСФЕРА
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("3. Тропосферная задержка"))
    s.append(HR())

    s.append(H2("3.1 Общая модель задержки"))
    s.append(formula_block(
        r"T^s = m_h(El)\cdot ZHD + m_w(El)\cdot ZWD"))
    s.append(var_table([
        ['ZHD', 'Зенитная гидростатическая задержка (вычисляется из давления)'],
        ['ZWD', 'Зенитная влажная задержка (оценивается фильтром Калмана)'],
        ['m_h(El)', 'Гидростатическая отображающая функция (mapping function)'],
        ['m_w(El)', 'Влажная отображающая функция'],
    ]))

    s.append(H2("3.2 Формула Саастамойнена для ZHD"))
    s.append(formula_block(
        r"ZHD = \frac{0.0022768\cdot P}"
        r"{1 - 0.00266\cos(2\varphi) - 2.8\times10^{-7}\,h}"))
    s.append(P("P — давление (мбар), φ — геодезическая широта, h — высота (м)."))

    s.append(H2("3.3 Оценка ZTD фильтром (TROPOPT_EST)"))
    s.append(P("Состояние x_trp — отклонение ZTD от модели:"))
    s.append(formula_block(
        r"T^s = m_h\cdot ZHD + m_w\cdot (ZTD_0 + x_{trp})"))
    s.append(P("Частная производная для строки матрицы H:"))
    s.append(formula_block(
        r"\frac{\partial T^s}{\partial x_{trp}} = m_w(El)"))

    s.append(H2("3.4 Градиентная модель (TROPOPT_ESTG)"))
    s.append(P("Три состояния: x_trp, ∇_N, ∇_E:"))
    s.append(formula_block(
        r"T^s = m_h\cdot ZHD + m_w\cdot x_{trp}"
        r" + m_w \frac{\cos(El)}{\sin(El)}"
        r"\left(\nabla_N\cos(Az) + \nabla_E\sin(Az)\right)"))
    s.append(P("Частные производные:"))
    s.append(formula_block(
        r"\frac{\partial T^s}{\partial \nabla_N}"
        r" = m_w\frac{\cos(El)}{\sin(El)}\cos(Az)"))
    s.append(formula_block(
        r"\frac{\partial T^s}{\partial \nabla_E}"
        r" = m_w\frac{\cos(El)}{\sin(El)}\sin(Az)"))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 4. ИОНОСФЕРА IONOOPT_EST
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("4. Ионосферная задержка в режиме IONOOPT_EST"))
    s.append(HR())

    s.append(P("Для каждого спутника вводится состояние I^s (вертикальная задержка). "
               "Наклонная задержка:"))
    s.append(formula_block(
        r"I_{slant}^s = \mathrm{IF}(El)\cdot I^s"))
    s.append(P("Отображающая функция ионосферного слоя (одно-слойная модель):"))
    s.append(formula_block(
        r"\mathrm{IF}(El) = \frac{1}{\cos\!\left("
        r"\arcsin\!\left(\frac{R_E}{R_E + H_{ion}}\cos(El)\right)"
        r"\right)}"))
    s.append(var_table([
        ['R_E = 6 378 137 м', 'Радиус Земли (WGS-84)'],
        ['H_ion = 350 000 м', 'Высота ионосферного слоя'],
    ]))
    s += [SP(4)]
    s.append(P("Коэффициент ионосферного слагаемого для частоты j:"))
    s.append(formula_block(
        r"C_j^P = +\frac{f_1^2}{f_j^2}\cdot\mathrm{IF}(El)"
        r"\quad (\mathrm{pseudorange})"))
    s.append(formula_block(
        r"C_j^L = -\frac{f_1^2}{f_j^2}\cdot\mathrm{IF}(El)"
        r"\quad (\mathrm{phase})"))
    s.append(P("Частная производная для матрицы H: ∂/∂I^s = C_j."))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 5. НЕВЯЗКИ И МАТРИЦА H
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("5. Вектор невязок и матрица частных производных H"))
    s.append(HR())

    s.append(H2("5.1 Скорректированные наблюдения"))
    s.append(P("Перед фильтрацией применяются поправки DCB/OSB из продуктов IGS/MGEX:"))
    s.append(P("Для псевдодальности:"))
    s.append(formula_block(
        r"y_j^s = P_j^s - b_{r,j}^{DCB} - b_j^{s,\mathrm{OSB}}"))
    s.append(P("Для фазы (в метрах):"))
    s.append(formula_block(
        r"y_j^s = L_j^s\cdot\lambda_j - b_{r,j}^{\phi}"))

    s.append(H2("5.2 Вектор невязок"))
    s.append(formula_block(
        r"v_j^s = y_j^s - \hat{\rho}^s - c\,\hat{\delta t}_r"
        r" - \widehat{T}^s - C_j\,\hat{I}^s - \hat{N}_j^s + c\,\delta t^s"))

    s.append(H2("5.3 Единичный вектор направления на спутник"))
    s.append(formula_block(
        r"\mathbf{e}^s = \frac{\mathbf{r}^s - \mathbf{r}_r}{\|\mathbf{r}^s - \mathbf{r}_r\|}"))

    s.append(H2("5.4 Структура строки матрицы H"))
    s += [SP(4)]
    avail = W - 2*MARGIN
    c1, c2 = 4.0*cm, 2.8*cm
    s.append(wide_table([
        ['Положение (x,y,z)',      'IC = 0..2',  '−eˢ'],
        ['Скорость',               '3..5',       '0'],
        ['Ускорение',              '6..8',       '0'],
        ['Часы приёмника (сист. k)', 'IC(k)',     '+1'],
        ['ZTD (тропосфера)',       'IT',          'm_w(El)'],
        ['Градиент ∇_N',          'IT+1',        'm_w · cos(El)/sin(El) · cos(Az)'],
        ['Градиент ∇_E',          'IT+2',        'm_w · cos(El)/sin(El) · sin(Az)'],
        ['Ионосфера спутника s',   'II(s)',       'C_j  (знак зависит от типа набл.)'],
        ['DCB L5',                 'ID',          '±1 (для L5)'],
        ['Неоднозначность N_j^s',  'IB(s,j)',    '+1'],
    ], headers=['Блок состояний', 'Индекс', 'H[i]'],
       col_widths=[c1, c2, avail-c1-c2]))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 6. МАТРИЦА R
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("6. Матрица шума измерений R"))
    s.append(HR())

    s.append(P("Матрица R диагональна. Дисперсия i-го измерения "
               "(<i>varerr()</i> в ppp.c):"))
    s.append(formula_block(
        r"\sigma_i^2 = \sigma_0^2 + \sigma_{trp}^2 + C_j^2\,\sigma_{ion}^2 + \sigma_{sat}^2"))

    s.append(H2("6.1 Базовая дисперсия"))
    s.append(formula_block(
        r"\sigma_0^2 = f_{sys}^2 \cdot f_{IF} \cdot"
        r"\left(a_\sigma^2 + \frac{b_\sigma^2}{\sin^2(El)}\right) \cdot \eta_j^2"))
    s.append(var_table([
        ['a_σ = err[1]',     'Постоянная составляющая погрешности, м'],
        ['b_σ = err[2]',     'Угломестнозависимая составляющая, м'],
        ['f_sys',            'GPS=1.0, ГЛОНАСС=1.5, Galileo=1.0, BeiDou=1.0'],
        ['f_IF = 3.0',       'Для IFLC; иначе 1.0 (учитывает усиление шума)'],
        ['η_j = eratio[j]',  '~100 для псевдодальности, ~1 для фазы'],
    ]))

    s.append(H2("6.2 Дополнительные составляющие дисперсии"))
    s.append(formula_block(
        r"\sigma_{trp}^2 = \mathrm{var\_trop}(El),\quad"
        r"\sigma_{ion}^2 = \mathrm{var\_iono}(El),\quad"
        r"\sigma_{sat}^2 = \mathrm{var\_rs}^2"))
    s.append(P("Для ГЛОНАСС (FDMA) дополнительная межчастотная погрешность:"))
    s.append(formula_block(
        r"\sigma_{GLO,IFB}^2 = (0.6)^2 = 0.36\;\mathrm{m}^2"
        r"\quad (\mathrm{GLONASS\ pseudorange\ only})"))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 7. ФИЛЬТР КАЛМАНА
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("7. Фильтр Калмана (общий вид)"))
    s.append(HR())

    s.append(H2("7.1 Вектор состояния"))
    s.append(P("Полный вектор состояния x объединяет все неизвестные. "
               "Суммарная размерность NX = NR + NB, где:"))
    s.append(var_table([
        ['NP = 3 (или 9)',    'Положение (+ скорость + ускорение при dynamics=1)'],
        ['NC = N_sys',        'Часы приёмника — по одному на каждую НС'],
        ['NT = 1 (или 3)',    'ZTD (+ горизонтальные градиенты при ESTG)'],
        ['NI = N_sat',        'Ионосфера — по одному на каждый спутник (при IONOOPT_EST)'],
        ['ND = 1',            'DCB L5 (при nf ≥ 3)'],
        ['NB = N_f × N_sat',  'Фазовые неоднозначности — на каждую частоту и спутник'],
    ]))

    s += [SP(6)]
    s.append(P("Структура вектора (слева направо):"))
    s.append(formula_block(
        r"\mathbf{x} = [\,\mathbf{r}_{xyz},\;"
        r"\mathbf{v}_{xyz},\;\mathbf{a}_{xyz},\;"
        r"\delta t_{r,k},\;ZTD,\nabla_N,\nabla_E,\;"
        r"I^s,\;DCB_{L5},\;N_j^s\,]^T", fontsize=12))

    s.append(H2("7.2 Уравнения прогноза (predict)"))
    s.append(formula_block(
        r"\mathbf{x}_{k|k-1} = \mathbf{F}_k\,\mathbf{x}_{k-1}"))
    s.append(formula_block(
        r"\mathbf{P}_{k|k-1} = \mathbf{F}_k\,\mathbf{P}_{k-1}\,\mathbf{F}_k^T + \mathbf{Q}_k"))

    s.append(H2("7.3 Уравнения обновления (update)"))
    s.append(formula_block(
        r"\mathbf{K}_k = \mathbf{P}_{k|k-1}\,\mathbf{H}_k^T"
        r"\left(\mathbf{H}_k\,\mathbf{P}_{k|k-1}\,\mathbf{H}_k^T + \mathbf{R}_k\right)^{-1}"))
    s.append(formula_block(
        r"\mathbf{x}_k = \mathbf{x}_{k|k-1} + \mathbf{K}_k\,\mathbf{v}_k"))
    s.append(formula_block(
        r"\mathbf{P}_k = \left(\mathbf{I} - \mathbf{K}_k\,\mathbf{H}_k\right)"
        r"\mathbf{P}_{k|k-1}"))
    s.append(var_table([
        ['F_k',                          'Матрица перехода состояния'],
        ['Q_k',                          'Матрица шума процесса (блочно-диагональная)'],
        ['H_k',                          'Матрица частных производных (Якобиан)'],
        ['R_k',                          'Диагональная матрица шума измерений'],
        ['v_k = y_k − H_k x_{k|k-1}',   'Вектор невязок (innovation)'],
        ['K_k',                          'Матрица усиления Калмана'],
    ]))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 8. ДИНАМИЧЕСКИЕ МОДЕЛИ
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("8. Динамические модели состояний (матрицы F и Q)"))
    s.append(HR())

    s.append(H2("8.1 Кинематическая модель с динамикой (opt.dynamics = 1)"))
    s.append(P("Блок положение–скорость–ускорение, матрица перехода <b>F</b> (9×9), "
               "Δt — интервал между эпохами. Строки: x,y,z, vx,vy,vz, ax,ay,az:"))
    # Matrix as table (pmatrix not supported in mpl mathtext)
    dt = 'Δt';  dt2 = 'Δt²/2'
    row0 = ['1','0','0', dt, '0','0', dt2,'0','0']
    row1 = ['0','1','0', '0',dt, '0', '0',dt2,'0']
    row2 = ['0','0','1', '0','0', dt, '0','0',dt2]
    row3 = ['0','0','0', '1','0','0', dt, '0','0']
    row4 = ['0','0','0', '0','1','0', '0',dt, '0']
    row5 = ['0','0','0', '0','0','1', '0','0',dt ]
    row6 = ['0','0','0', '0','0','0', '1','0','0']
    row7 = ['0','0','0', '0','0','0', '0','1','0']
    row8 = ['0','0','0', '0','0','0', '0','0','1']
    mat_rows = [row0,row1,row2,row3,row4,row5,row6,row7,row8]
    avail_m = W - 2*MARGIN
    cw = avail_m / 9
    mat_tbl = Table(mat_rows, colWidths=[cw]*9)
    mat_tbl.setStyle(TableStyle([
        ('FONTSIZE',      (0,0),(-1,-1), 8),
        ('ALIGN',         (0,0),(-1,-1), 'CENTER'),
        ('TOPPADDING',    (0,0),(-1,-1), 2),
        ('BOTTOMPADDING', (0,0),(-1,-1), 2),
        ('GRID',          (0,0),(-1,-1), 0.5, colors.HexColor('#aaaaaa')),
        ('BACKGROUND',    (0,0),(-1,-1), colors.HexColor('#f0f4ff')),
        ('FONTNAME',      (0,0),(-1,-1), 'Helvetica'),
        # Highlight diagonal blocks
        ('BACKGROUND',    (0,0),(2,2), colors.HexColor('#d0e8ff')),  # pos
        ('BACKGROUND',    (3,3),(5,5), colors.HexColor('#d0e8ff')),  # vel
        ('BACKGROUND',    (6,6),(8,8), colors.HexColor('#d0e8ff')),  # acc
        ('FONTNAME',      (3,0),(5,2), 'Helvetica-Bold'),
        ('FONTNAME',      (6,0),(8,2), 'Helvetica-Bold'),
    ]))
    s.append(SP(4))
    s.append(mat_tbl)
    s.append(SP(4))
    s.append(P("Голубые диагональные блоки — единичные (I). "
               "Надиагональные: Δt связывает pos↔vel, Δt²/2 связывает pos↔acc, "
               "Δt связывает vel↔acc."))

    s.append(P("Шум процесса добавляется в блок ускорений (3×3):"))
    s.append(formula_block(
        r"\mathbf{Q}_{acc} = \mathrm{diag}\!\left("
        r"q_h^2|\Delta t|,\;q_h^2|\Delta t|,\;q_v^2|\Delta t|\right)"))
    s.append(P("q_h = prn[3] — горизонтальный шум, q_v = prn[4] — вертикальный."))

    s.append(H2("8.2 Часы приёмника — белый шум"))
    s.append(P("Сбрасываются каждую эпоху (нет корреляции между эпохами):"))
    s.append(formula_block(
        r"P_{clk} \leftarrow \sigma_{clk}^2 = 60^2 = 3600\;\mathrm{m}^2"))

    s.append(H2("8.3 Тропосфера — случайное блуждание"))
    s.append(formula_block(
        r"P_{ZTD} \leftarrow P_{ZTD} + q_{ZTD}^2\cdot|\Delta t|"
        r",\quad q_{ZTD} = \mathrm{prn}[2]"))
    s.append(formula_block(
        r"P_{\nabla} \leftarrow P_{\nabla} + (0.1\,q_{ZTD})^2\cdot|\Delta t|"
        r"\quad (\mathrm{gradients})"))

    s.append(H2("8.4 Ионосфера — случайное блуждание с угломестной зависимостью"))
    s.append(formula_block(
        r"P_{I^s} \leftarrow P_{I^s}"
        r"+ \left(\frac{q_{ion}}{\sin(\max(El,\,5^{\circ}))}\right)^{2}\cdot|\Delta t|"
        r",\quad q_{ion} = \mathrm{prn}[1]"))

    s.append(H2("8.5 Фазовые неоднозначности"))
    s.append(P("Без срыва цикла — случайное блуждание:"))
    s.append(formula_block(
        r"P_{N_j^s} \leftarrow P_{N_j^s} + q_{bias}^2\cdot|\Delta t|"
        r",\quad q_{bias} = \mathrm{prn}[0]"))
    s.append(P("При срыве цикла — полный сброс:"))
    s.append(formula_block(
        r"P_{N_j^s} \leftarrow \sigma_{bias}^2 = 60^2 = 3600\;\mathrm{m}^2"))

    s.append(H2("8.6 Начальные дисперсии состояний"))
    s.append(wide_table([
        ['Положение',         '60 м',   '3 600 м²'],
        ['Часы',              '60 м',   '3 600 м²'],
        ['ZTD',               '0.6 м',  '0.36 м²'],
        ['Градиенты ∇',       '0.01 м', '10⁻⁴ м²'],
        ['Ионосфера I^s',     '60 м',   '3 600 м²'],
        ['Неоднозначность N', '60 м',   '3 600 м²'],
        ['DCB',               '60 м',   '3 600 м²'],
    ], headers=['Состояние', 'σ₀', 'σ₀²'],
       col_widths=[5*cm, 2.5*cm, avail-7.5*cm]))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 9. СРЫВЫ ЦИКЛА
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("9. Обнаружение срывов цикла"))
    s.append(HR())

    s.append(H2("9.1 Флаг LLI"))
    s.append(P("Приёмник выставляет бит LLI (Loss of Lock Indicator). "
               "При <b>LLI &amp; 1 = 1</b> — срыв зафиксирован аппаратурно."))

    s.append(H2("9.2 Геометрически-свободная комбинация (GF)"))
    s.append(P("Устраняет геометрию, тропосферу, часы — остаётся ионосфера + неоднозначности:"))
    s.append(formula_block(
        r"GF = L_1 - L_2 = 2I_1^s + (\lambda_1 N_1 - \lambda_2 N_2)"))
    s.append(P("Срыв, если:"))
    s.append(formula_block(
        r"|GF_k - GF_{k-1}| > \varepsilon_{GF}(El)"))

    s.append(H2("9.3 Комбинация Мельбурна–Вюббены (MW)"))
    s.append(P("Не зависит от ионосферы, тропосферы, геометрии:"))
    s.append(formula_block(
        r"MW = \frac{f_1 L_1 - f_2 L_2}{f_1 - f_2}"
        r" - \frac{f_1 P_1 + f_2 P_2}{f_1 + f_2}"))
    s.append(P("В стационарном состоянии MW ≈ λ_WL · N_WL. "
               "Широкополосная длина волны (GPS):"))
    s.append(formula_block(
        r"\lambda_{WL} = \frac{c}{f_1 - f_2} \approx 0.862\;\mathrm{m}"))
    s.append(P("Порог обнаружения срыва: |ΔMW / λ_WL| > 10 циклов."))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 10. PPP-AR
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("10. Разрешение фазовых неоднозначностей (PPP-AR)"))
    s.append(HR())

    s.append(H2("10.1 Wide-Lane неоднозначность"))
    s.append(P("После применения OSB-поправок из продуктов IGS:"))
    s.append(formula_block(
        r"\hat{N}_{WL}^s = \frac{MW^s}{\lambda_{WL}} = N_{WL}^s + \varepsilon_{WL}"))
    s.append(formula_block(
        r"\check{N}_{WL}^s = \mathrm{round}\!\left(\hat{N}_{WL}^s\right)"))

    s.append(H2("10.2 Narrow-Lane неоднозначность"))
    s.append(P("После фиксации WL, из IF-неоднозначности:"))
    s.append(formula_block(
        r"\hat{N}_{NL}^s = \frac{\hat{N}_{IF}^s - \beta\,\lambda_2\,\check{N}_{WL}^s}{\alpha\,\lambda_1}"))
    s.append(P("Длина волны NL для GPS:"))
    s.append(formula_block(
        r"\lambda_{NL} = \frac{c}{f_1 + f_2} \approx 0.107\;\mathrm{m}"))

    s.append(H2("10.3 Межспутниковые одиночные разности (SD)"))
    s.append(P("Для исключения приёмных биасов — разности между опорным r и тестируемым i:"))
    s.append(formula_block(
        r"\nabla\hat{N}_{NL}^{ri} = \hat{N}_{NL}^r - \hat{N}_{NL}^i"))

    s.append(H2("10.4 Метод LAMBDA (Integer Least Squares)"))
    s.append(formula_block(
        r"\check{\mathbf{N}} = \arg\min_{\mathbf{N}\in\mathbf{Z}^n}"
        r"\;(\hat{\mathbf{N}} - \mathbf{N})^T \mathbf{Q}_N^{-1} (\hat{\mathbf{N}} - \mathbf{N})"))
    s.append(P("Решается через Z-преобразование для декорреляции матрицы Q_N."))

    s.append(H2("10.5 Ratio Test"))
    s.append(formula_block(
        r"R = \frac{Q(\check{\mathbf{N}}_2)}{Q(\check{\mathbf{N}}_1)} > \mathrm{thresh}_{AR}"))
    s.append(P("Ň₁ — лучшее, Ň₂ — второе лучшее целочисленное решение. "
               "Типичный порог: thresh_AR = 3.0."))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 11. ГЛОНАСС IFB
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("11. ГЛОНАСС: межчастотные поправки IFB"))
    s.append(HR())

    s.append(P("ГЛОНАСС использует FDMA: каждый спутник передаёт на своей частоте:"))
    s.append(formula_block(
        r"f_k = f_0 + k\cdot\Delta f,\quad"
        r"f_0 = 1602\;\mathrm{MHz},\;\Delta f = 562.5\;\mathrm{kHz}"))
    s.append(P("Приборная задержка приёмника зависит от частотного номера k:"))
    s.append(formula_block(
        r"\delta_{WL,r}^{GLO}(k) = \delta_{WL,r} + k\cdot\Delta\delta_{WL}"))
    s.append(P("При однократном дифференцировании между спутниками i и j:"))
    s.append(formula_block(
        r"\nabla N_{WL}^{ij} = N_{WL}^{ij} + (k_i - k_j)\cdot\Delta\delta_{WL}"))
    s.append(P("Оценка IFB методом наименьших квадратов по всем парам:"))
    s.append(formula_block(
        r"\Delta\hat{\delta}_{WL}"
        r"= \frac{\sum_i \Delta k_i\cdot\mathrm{frac}(N_{WL}^{ri})}{\sum_i \Delta k_i^2}"))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 12. CPP
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("12. CPP — Collaborative Precise Positioning"))
    s.append(HR())

    s.append(H2("12.1 Концепция"))
    s.append(P("CPP ускоряет сходимость PPP путём инъекции RTK-решения от партнёрского приёмника. "
               "RTK сходится за секунды, но требует базовой станции вблизи. "
               "PPP не требует базы, но сходится 20–40 минут. "
               "CPP объединяет оба подхода."))

    s.append(H2("12.2 Критерий качества сходимости PPP"))
    s.append(P("Геометрическое среднее дисперсий положения:"))
    s.append(formula_block(
        r"\sigma_{PPP}^{geo} = \sqrt[3]{P_{xx}\cdot P_{yy}\cdot P_{zz}}"))

    s.append(H2("12.3 Условия инъекции RTK-решения"))
    s.append(P("<b>Условие 1.</b> Дисперсия PPP хуже RTK:"))
    s.append(formula_block(
        r"\sigma_{PPP}^{geo} > \sigma_{RTK}^{geo}"))
    s.append(P("<b>Условие 2.</b> Расхождение в положении превышает 2σ:"))
    s.append(formula_block(
        r"|\hat{x}_{PPP,i} - \hat{x}_{RTK,i}| > 2\sqrt{P_{PPP,ii}}"
        r"\quad\forall\,i\in\{x,y,z\}"))

    s.append(H2("12.4 Процедура инъекции"))
    s.append(P("RTK-решение инжектируется как псевдоизмерение — "
               "эквивалент шага обновления Калмана с H = [I₃ₓ₃ 0] и R = P_RTK:"))
    s.append(formula_block(
        r"\mathbf{x}_r \leftarrow \mathbf{x}_{RTK},\qquad"
        r"\mathbf{P}_{r,pos} \leftarrow \mathbf{P}_{RTK} + \mathbf{Q}_{inject}"))

    s.append(H2("12.5 Логика ожидания стабилизации"))
    s.append(wide_table([
        ['t < 30 эпох с момента потери сходимости',  'Ждём (warmup)'],
        ['Улучшение σ > 2% за эпоху',                'Продолжаем ждать'],
        ['3 эпохи подряд без улучшения > 2%',        'Стабильно → проверяем инъекцию'],
    ], headers=['Условие', 'Действие'],
       col_widths=[(W-2*MARGIN)*0.62, (W-2*MARGIN)*0.38]))
    s += [SP(6)]
    s.append(var_table([
        ['CPP_WARMUP_EPOCHS = 30',   'Эпохи прогрева после потери сходимости'],
        ['CPP_STABLE_EPOCHS = 3',    'Число стабильных эпох для подтверждения'],
        ['CPP_IMPROVE_THRESH = 0.02','Порог относительного улучшения (2%)'],
    ]))
    s.append(PageBreak())

    # ══════════════════════════════════════════════════════════════════════════
    # 13. СПРАВОЧНЫЕ ТАБЛИЦЫ
    # ══════════════════════════════════════════════════════════════════════════
    s.append(H1("13. Справочные таблицы"))
    s.append(HR())

    s.append(H2("Физические константы"))
    s.append(wide_table([
        ['c = 299 792 458 м/с',         'Скорость света'],
        ['f₁ = 1 575.42 МГц',           'GPS L1'],
        ['f₂ = 1 227.60 МГц',           'GPS L2'],
        ['f₅ = 1 176.45 МГц',           'GPS L5'],
        ['Ω_E = 7.292115×10⁻⁵ рад/с',   'Угловая скорость Земли (WGS-84)'],
        ['R_E = 6 378 137 м',            'Радиус Земли (WGS-84)'],
        ['H_ion = 350 000 м',            'Высота ионосферного слоя'],
        ['μ = 3.986×10¹⁴ м³/с²',        'Гравитационный параметр Земли'],
    ], headers=['Константа', 'Описание'],
       col_widths=[5*cm, avail-5*cm]))

    s += [SP(8)]
    s.append(H2("Шумы процесса (prn[])"))
    s.append(wide_table([
        ['prn[0]', '~10⁻⁴',   'Шум фазовой неоднозначности, м/√с'],
        ['prn[1]', '~10⁻⁴',   'Шум ионосферы, м/√с'],
        ['prn[2]', '~10⁻⁵',   'Шум ZTD, м/√с'],
        ['prn[3]', '~10⁻⁴',   'Шум горизонтального ускорения, м/с^{3/2}'],
        ['prn[4]', '~10⁻⁴',   'Шум вертикального ускорения, м/с^{3/2}'],
        ['prn[5]', '~10⁻¹²',  'Шум часов приёмника, м/√с'],
    ], headers=['Параметр', 'Значение', 'Описание'],
       col_widths=[2.5*cm, 2.5*cm, avail-5*cm]))

    s += [SP(8)]
    s.append(H2("Ключевые файлы кода RTKLIB"))
    s.append(wide_table([
        ['src/ppp.c',    'Основной алгоритм PPP: динамика, невязки, фильтрация, CPP'],
        ['src/ppp_ar.c', 'Разрешение неоднозначностей: WL, NL, LAMBDA, GLONASS IFB'],
        ['src/rtkcmn.c', 'Ядро фильтра Калмана, тропосфера, ионосфера'],
        ['src/rtklib.h', 'Константы и макросы индексов: IC(), IT(), II(), IB(), ID()'],
    ], headers=['Файл', 'Содержимое'],
       col_widths=[3.5*cm, avail-3.5*cm]))

    return s


# ─── Build ────────────────────────────────────────────────────────────────────

OUTPUT = "/Users/sergeidolin/RTKLIB/PPP_CPP_Lecture.pdf"

doc = SimpleDocTemplate(
    OUTPUT, pagesize=A4,
    leftMargin=MARGIN, rightMargin=MARGIN,
    topMargin=MARGIN, bottomMargin=MARGIN,
    title="PPP/CPP Лекция — RTKLIB",
    author="RTKLIB",
)

print("Building story…")
story = build_story()
print(f"  {len(story)} elements")
print("Rendering PDF…")
doc.build(story)
print(f"Done → {OUTPUT}")
