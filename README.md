# Optimización de la señal lipídica en resonancia magnética orientada a la sarcopenia

Este repositorio contiene los códigos desarrollados para el Trabajo de Fin de Grado titulado **"Optimización de la señal lipídica en resonancia magnética orientada a la sarcopenia"**.

El objetivo del trabajo es estudiar distintas estrategias de supresión de grasa en resonancia magnética mediante simulaciones computacionales. En concreto, se analizan técnicas clásicas como STIR y CHESS, junto con pulsos espectral-espaciales (SPSP), evaluando su comportamiento en condiciones controladas y ante inhomogeneidades simuladas del campo magnético principal B0.

## Estructura del repositorio

- `codigo/CHESS/`: simulaciones relacionadas con pulsos CHESS y su respuesta espectral.
- `codigo/STIR/`: simulaciones de recuperación longitudinal y supresión de grasa mediante STIR.
- `codigo/SPSP/`: diseño y evaluación de pulsos espectral-espaciales mediante la aproximación STA.
- `codigo/phantom_B0/`: simulaciones sobre phantom digital e incorporación de mapas sintéticos de B0.
- `documentacion/descripcion_codigos.md`: explicación detallada de cada script.

## Herramientas utilizadas

Los códigos se han desarrollado principalmente en **MATLAB R2024b**.

También se han empleado funciones externas relacionadas con simulaciones de Bloch y herramientas de reconstrucción/diseño de pulsos, como MIRT, cuando ha sido necesario para construir matrices del sistema o resolver problemas de optimización.

## Técnicas estudiadas

- **STIR**: supresión de grasa basada en diferencias de tiempos de relajación T1.
- **CHESS**: saturación selectiva en frecuencia de la señal grasa.
- **SPSP**: pulsos espectral-espaciales que combinan selectividad en frecuencia y en posición espacial.

## Nota

Los resultados incluidos proceden de simulaciones computacionales. No representan una validación experimental en un sistema de resonancia magnética real.
