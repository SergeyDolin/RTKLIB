# RTKLIB 2.4.3 d4

## Новые возможности / New Features

### PPP-AR (Precise Point Positioning with Ambiguity Resolution)

В версии d4 добавлена поддержка **PPP-AR** — Precise Point Positioning с разрешением целочисленных неоднозначностей (Ambiguity Resolution). Это позволяет достигать сантиметровой точности позиционирования в режиме real-time без использования опорных станций.

**Version d4 introduces PPP-AR (Precise Point Positioning with Ambiguity Resolution)** — a new mode enabling centimeter-level positioning accuracy in real-time without base station networks by resolving integer ambiguities.
 
## Возможности / Features
 
### Multi-System GNSS Processing in PPP
- Обработка измерений нескольких систем ГНСС в режиме PPP с учётом дифференциальных задержек кода / Support for multi-system GNSS measurements in PPP mode with differential code delay correction
- L1/L5 комбинация поддерживается / L1/L5 combination supported
- Поддержка GPS, ГЛОНАСС, Galileo, BDS / Support for GPS, GLONASS, Galileo, BDS
 
### Differential Code Bias (DCB)
- **Post-processing:** используйте BSX файлы от CAS (https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/gnss_differential_code_bias_product.html)
- **Real-time:** любой SSR-поток с DCB или SSRA00CNE1 от IGS
 
### Advanced Kalman Filter
- Реализован алгоритм робастного адаптивного фильтра Калмана на основе вариационного Байеса (VBBRA)
- Variational Bayesian-based robust adaptive Kalman filter algorithm
- В конфигурационном файле можно выбрать базовый алгоритм или VBBRA (pos1-kalman)
 
## Улучшения в d4 / Improvements in d4
- ✅ Исправлены ошибки в решателе PPP-AR
- ✅ Улучшена обработка BDS сигналов
- ✅ Исправлены проблемы с парсингом DCB файлов
- ✅ Оптимизирована производительность вычислений