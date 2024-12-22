# XYZgeomag

[![Build Status](https://github.com/nhz2/XYZgeomag/workflows/test/badge.svg)](https://github.com/nhz2/XYZgeomag/actions)

Lightweight C++ header-only library for calculating the magnetic field on earth given geocentric cartesian coordinates using the [World Magnetic Model(WMM)](https://www.ncei.noaa.gov/products/world-magnetic-model). Compatible with Arduino.

The main function `geomag::GeoMag` calculates the magnetic field around earth in the International Terrestrial Reference System(ITRS) and uses units of decimal year, meter, and tesla.

Unlike most WMM software, which uses latitude, longitude, and altitude inputs to calculate the North East Down components of the magnetic field, `geomag::GeoMag` uses geocentric cartesian coordinates as input, and outputs the magnetic field in the same geocentric cartesian coordinate system as the inputs.

If you want to provide geodetic latitude, longitude, and height, and receive the local North East Down components of the magnetic field and the magnetic declination: 
see the `geomag::geodetic2ecef` and `geomag::magField2Elements` example below.
Note that latitude and longitude are in units of degrees, and the seven magnetic elements are in units of nanotesla and degrees.

## Error

XYZgeomag is within 0.5 nT of the official WMM software.

For more information on the limitations of the WMM model, see:
<https://www.ncei.noaa.gov/products/world-magnetic-model/accuracy-limitations-error-model>

## Performance

XYZgeomag uses single precision floating points. It's designed to minimize ram usage for embedded systems.

| Device      | Speed    |
|-------------|----------|
| Arduino Uno R3 | 52 ms    |
| Raspberry Pi Pico RP2040 | 4.5 ms |
| Teensy 3.6  |  83 µs |
| Teensy 4.0  |  21 µs |

## Using XYZgeomag

Just download [XYZgeomag.hpp](https://github.com/nhz2/XYZgeomag/releases/download/v2.0.0/XYZgeomag.hpp) and include it.
Here is an example Arduino sketch:

~~~cpp
#include "XYZgeomag.hpp"
void setup() {
  // put your setup code here, to run once:
  pinMode(1,INPUT);
  Serial.begin(9600);
}

void loop() {
  // put your main code here, to run repeatedly:
  int val= digitalRead(1);
  geomag::Vector in;
  in.x=val+1128529.6885767058f;
  in.y=val+0.0;
  in.z=val+6358023.736329913f;
  geomag::Vector out;
  int starttime=micros();
  int starttimemil=millis();
  out=geomag::GeoMag(2022.5,in,geomag::WMM2020);
  int endtime=micros();
  int endtimemil=millis();
  Serial.print(out.x*1E9);
  Serial.println(" nT x");
  Serial.print(out.y*1E9);
  Serial.println(" nT y");
  Serial.print(out.z*1E9);
  Serial.println(" nT z");
  Serial.print("time in micro seconds: ");
  Serial.println(endtime-starttime);
  Serial.print("time in milli seconds: ");
  Serial.println(endtimemil-starttimemil);
  delay(2000);
}
~~~

If you have a position in latitude, longitude, and height, 
you can convert it to geocentric cartesian coordinates 
with `geodetic2ecef`. Note that `geodetic2ecef` uses 
single precision floats, so it will only be accurate to about 1 meter.
You can also convert the magnetic field to 
the [seven magnetic elements](https://www.ncei.noaa.gov/products/world-magnetic-model) 
in units of nanotesla and degrees.
~~~cpp
#include "XYZgeomag.hpp"
void setup() {
  // put your setup code here, to run once:
  pinMode(1,INPUT);
  Serial.begin(9600);
}

void loop() {
  // put your main code here, to run repeatedly:
  int val= digitalRead(1);
  float lat = val + 43.0f; // latitude in degrees
  float lon = val + 75.0f; // longitude in degrees
  float height = val + 305; // height above WGS84 ellipsoid in meters
  geomag::Vector position = geomag::geodetic2ecef(lat,lon,height);
  geomag::Vector mag_field = geomag::GeoMag(2022.5,position,geomag::WMM2020);
  geomag::Elements out = geomag::magField2Elements(mag_field, lat, lon);
  Serial.print(out.north);
  Serial.println(" nT north");
  Serial.print(out.east);
  Serial.println(" nT east");
  Serial.print(out.down);
  Serial.println(" nT down");
  Serial.print(out.horizontal);
  Serial.println(" nT horizontal");
  Serial.print(out.total);
  Serial.println(" nT total");
  Serial.print(out.inclination);
  Serial.println(" deg inclination");
  Serial.print(out.declination);
  Serial.println(" deg declination");
  Serial.println();
  delay(2000);
}
~~~



## Adding New Coefficents

To add new coefficents, download the new `.COF` file from <https://www.ncei.noaa.gov/products/world-magnetic-model>

Add the .COF file to the `extras` directory.

Then run for example
`python3 wmmcodeupdate.py -f WMM2015.COF -f WMM2015v2.COF -f WMM2020.COF -f WMM2025.COF -o ../src/XYZgeomag.hpp -n 12` from the `extras` directory.

In this example, `WMM2015.COF` ,  `WMM2015v2.COF`, and  `WMM2020.COF` are the `.COF` files to use in `src/XYZgeomag.hpp`.

## Run Tests

In the `extras` directory.

Compile `geomag_test.cpp` for example with the command `g++ geomag_test.cpp -std=c++14`

Run the tests for example with the command `./a.out`

To add new models to the test update `wmmtestgen.py` and run it.

## References

Using spherical harmonics algorithm, described in sections 3.2.4 and 3.2.5:

  Satellite Orbits Models, Methods and Applications,
    by Oliver Montenbruck and Eberhard Gill 2000
    
Using geodetic2ecef algorithm from https://geographiclib.sourceforge.io/

Using coefficients and test points from:

NOAA NCEI Geomagnetic Modeling Team; British Geological Survey. 2024: World Magnetic Model 2025. NOAA National Centers for Environmental Information. https://doi.org/10.25921/aqfd-sd83. Accessed [22 DEC 2024].

NCEI Geomagnetic Modeling Team and British Geological Survey. 2019. World Magnetic Model 2020. NOAA National Centers for Environmental Information. doi: 10.25921/11v3-da71, 2020, [10 DEC 2019].

Chulliat, A., W. Brown, P. Alken, S. Macmillan, M. Nair, C. Beggan, A. Woods, B. Hamilton, B. Meyer and R. Redmon, 2019, Out-of-Cycle Update of the US/UK World Magnetic Model for 2015-2020: Technical Note, National Centers for Environmental Information, NOAA. doi: 10.25921/xhr3-0t19.

Chulliat, A., S. Macmillan, P. Alken, C. Beggan, M. Nair, B. Hamilton, A. Woods, V. Ridley, S. Maus and A. Thomson, 2015. The US/UK World Magnetic Model for 2015-2020: Technical Report, NOAA National Geophysical Data Center, Boulder, CO, doi: 10.7289/V5TB14V7.
