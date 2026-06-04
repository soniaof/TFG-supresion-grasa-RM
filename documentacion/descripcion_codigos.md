# Descripción de los códigos

Este documento describe la función principal de los scripts incluidos en el repositorio.

La organización se ha dividido en tres carpetas principales dentro de `Código TFG`: `Código CHESS`, `Código STIR` y `Código SPSP`. Cada una contiene los códigos empleados para simular o analizar una técnica de supresión de grasa en resonancia magnética.

## Código CHESS

Los scripts de la carpeta `Código TFG/Código CHESS/` permiten simular la técnica CHESS, basada en la saturación selectiva en frecuencia de la señal grasa.

En general, estos códigos:
- definen los parámetros de un sistema de resonancia magnética de 1,5 T,
- generan un pulso de radiofrecuencia centrado en la frecuencia de resonancia de la grasa,
- simulan la respuesta de la magnetización mediante las ecuaciones de Bloch,
- representan la respuesta espectral del pulso,
- permiten analizar la evolución de la magnetización longitudinal y transversal.

Los resultados principales son gráficas de la respuesta espectral del pulso y de la evolución de la magnetización durante la aplicación del pulso RF.

## Código STIR

Los scripts de la carpeta `Código TFG/Código STIR/` permiten simular la técnica STIR, basada en recuperación con inversión.

En general, estos códigos:
- definen valores de T1 para distintos tejidos, como grasa, músculo y agua,
- simulan la recuperación longitudinal de la magnetización tras un pulso de inversión,
- calculan el tiempo de inversión necesario para anular la señal de la grasa,
- aplican la supresión sobre un phantom digital,
- comparan la imagen original, la imagen con supresión y la señal eliminada.

Los resultados principales son curvas de recuperación longitudinal y figuras que muestran la supresión de la región grasa en el phantom.

## Código SPSP

Los scripts de la carpeta `Código TFG/Código SPSP/` permiten diseñar y evaluar pulsos espectral-espaciales.

En general, estos códigos:
- definen un target de magnetización en frecuencia y posición espacial,
- construyen una región de interés para agua y grasa,
- generan trayectorias en el espacio-k,
- calculan el pulso RF mediante la aproximación STA,
- comparan distintas trayectorias, como VDS, espiral de Arquímedes y espiral de densidad variable,
- analizan el efecto de parámetros como beta y FA,
- evalúan el ajuste mediante el NRMSE.

Los resultados principales son figuras de comparación entre el target y la magnetización obtenida, tablas de NRMSE y gráficas del pulso RF diseñado.

## Requisitos generales

Los códigos se han desarrollado en MATLAB R2024b.

Algunos scripts pueden requerir funciones externas relacionadas con simulación de Bloch o con el toolbox MIRT. En caso de utilizarse, estas funciones deben añadirse al path de MATLAB antes de ejecutar los scripts correspondientes.

## Nota

Los códigos se han desarrollado con fines académicos y los resultados proceden de simulaciones computacionales. No se han validado experimentalmente en un sistema real de resonancia magnética.
