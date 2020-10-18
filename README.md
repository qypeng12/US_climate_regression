# Regression of climate anomalies in continental U.S.
**Inspiration**
Climate change makes extreme rainfall more likely. An understanding of the variation of temperature influence upon precipitation is crucial for long term forecasting and climate modeling, and crucial for prediction and preparedness of droughts or floods. The purpose of this project is to determine the regression of rainfall on El Niño, as defined by temperature anomaly using Niño 3.4 index, and whether that varies across continental US as well as across seasons.

**What it does**
I explored the teleconnection of Niño 3.4 (temperature anomaly proxy) with seasonal rainfall totals with linear regression and correlation analysis. An ensemble of eight simulations for the 66-year period (1948–2014) is completed to investigate the role of temperature in modulating rainfall over the US.

**How I built it**
The dataset was downloaded from http://www.esrl.noaa.gov/psd/data/gridded/data.noaa.ersst.html Monthly precipitation data are aggregated into four seasons: Spring (Mar-Apr-May), Summer (Jun-Jul-Aug), Autumn (Sep-Oct-Nov) and Winter (Dec-Jan-Feb). For each season, spatial linear regression and correlation is performed and made into maps of regression/correlation coefficients for convenience of comparison. Finally, Student’s t-test is applied to determine if the correlations are statistically significant in order to attain a robust result.

**Challenges I ran into**
There is no spatially coherent ENSO-precipitation correlations identiﬁed except for summertime. The relationships are more complicated and vague between precipitation and Niño 3.4 index in all other seasons with opposite correlation signs in different locations. The lack of robust relationships may be because regional precipitation itself has complex inherent variability, and there are various climate factors that have varying inﬂuences across the region, making it difﬁcult to recognize one single inﬂuence over the whole domain.

**Accomplishments that I'm proud of**
Developed MATLAB code for analyzing geospatial dataset, and made cool maps that reveal the regression/correlation pattern of US precipitation on temperature indexes.
