# cal_brd (archived coursework, 2022–2023)

> ⚠️ **This script contains known errors and should not be used as is.**

Known problems:

- GM is written as 3.986005; the correct value is 3.986005e14 m³/s².
- The longitude of the ascending node uses the argument of perigee (omega) where the Earth rotation rate (7.2921151467e-5 rad/s) is required.
- Harmonic corrections use 2*omega+v instead of 2*(omega+v); the true anomaly needs atan2 and the -e belongs in the denominator.

A corrected and unit-tested version (validated against pyproj and real IGS products)
is part of **gnsskit**, a consolidated GNSS/geodesy toolkit by the same author
(repository `xyz2plh`, currently private). The original code here is kept unchanged
for history.
