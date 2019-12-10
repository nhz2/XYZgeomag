# XYZgeomag
Lightweight C++ header-only library for calculating the magnetic field on earth given geocentric cartesian coordinates instead of latitude, longitude, and altitude. Compatible with CUDA and Arduino.

XYZgeomag calculates the magnetic field around earth in the International Terrestrial Reference System(ITRS).

Unlike most World Magnetic Model(WMM) software, which uses latitude, longitude, and altitude inputs to calculate the North East Down(NED) components of the magnetic field, XYZgeomag uses geocentric cartesian coordinates as input, and outputs the magnetic field in the same geocentric cartesian coordinate system as the inputs.

## Error

XYZgeomag is within 10nT of the official WMM software.

## Performance

XYZgeomag uses single precision floating points. It's designed to minimize ram usage for embedded systems.

| Device      | Speed    |
|-------------|----------|
| Arduino Uno | 50 ms    |
| Teensy 3.6  |  80 µs |
| GeForce GT 750M  |  262 µs   |

## Using XYZgeomag

Just download `geomag.hpp` and include it.
Here is an example Arduino sketch:

~~~cpp
#include "geomag.hpp"
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
  out=geomag::GeoMag(2017.5,in,geomag::WMM2015);
  int endtime=micros();
  int endtimemil=millis();
  Serial.println(out.x*1E9);
  Serial.println(out.y*1E9);
  Serial.println(out.z*1E9);
  Serial.print("time in micro seconds: ");
  Serial.println(endtime-starttime);
  Serial.print("time in milli seconds: ");
  Serial.println(endtimemil-starttimemil);
  delay(2000);
}
~~~


Here is an example CUDA program:

~~~cpp
#include <iostream>
#include "geomag.hpp"

__constant__ const geomag::ConstModel WMM = geomag::WMM2015;

__global__
void mag(float3* result)
{
  int index = threadIdx.x+blockDim.x*blockIdx.x;
  geomag::Vector in;
  in.x=1128529.6885767058f;
  in.y=0.0f + index*100.0f;
  in.z=6358023.736329913f;
  geomag::Vector out;
  out=geomag::GeoMag(2017.5,in,WMM);
  result[index].x= out.x*1E9f;
  result[index].y= out.y*1E9f;
  result[index].z= out.z*1E9f;
}

int main(void)
{
  float3 *magv;
  int blocks= 32;
  int threadsperblock= 64;
  // Allocate Unified Memory – accessible from CPU or GPU
  cudaMallocManaged(&magv, sizeof(float3)*threadsperblock*blocks);
  mag<<<blocks, threadsperblock>>>(magv);
  // Wait for GPU to finish before accessing on host
  cudaDeviceSynchronize();
  // Print magnetic field
  for (int i=0; i<threadsperblock*blocks; i++){
    std::cout << "Mag field: " << magv[i].x << magv[i].y << magv[i].z << std::endl;
  }
  // Free memory
  cudaFree(magv);
  return 0;
}
~~~

## Adding New Coefficents

To add new coefficents, download the new `.COF` file from [https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml](https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml)

Add the .COF file to the `test_codegen` directory.

Then run for example
`python wmmcodeupdate.py -f WMM2015.COF -f WMM2015v2.COF -f WMM2020.COF -o ../geomag.hpp -n 12` from the `test_codegen` directory.

In this example, `WMM2015.COF` ,  `WMM2015v2.COF`, and  `WMM2020.COF` are the `.COF` files to use in `geomag.hpp`.

## Run Tests

In the `test_codegen` directory.

Compile `geomag_test.cpp` for example with the command `g++ geomag_test.cpp -std=c++14`

Run the tests for example with the command `./a.out`

To add new models to the test update `wmmtestgen.py` and run it.

## References

Using spherical harmonics algorithm, described in sections 3.2.4 and 3.2.5:

  Satellite Orbits Models, Methods and Applications,
    by Oliver Montenbruck and Eberhard Gill 2000

Using coefficients and test points from:

NCEI Geomagnetic Modeling Team and British Geological Survey. 2019. World Magnetic Model 2020. NOAA National Centers for Environmental Information. doi: 10.25921/11v3-da71, 2020, [10 DEC 2019].

Chulliat, A., W. Brown, P. Alken, S. Macmillan, M. Nair, C. Beggan, A. Woods, B. Hamilton, B. Meyer and R. Redmon, 2019, Out-of-Cycle Update of the US/UK World Magnetic Model for 2015-2020: Technical Note, National Centers for Environmental Information, NOAA. doi: 10.25921/xhr3-0t19.

Chulliat, A., S. Macmillan, P. Alken, C. Beggan, M. Nair, B. Hamilton, A. Woods, V. Ridley, S. Maus and A. Thomson, 2015. The US/UK World Magnetic Model for 2015-2020: Technical Report, NOAA National Geophysical Data Center, Boulder, CO, doi: 10.7289/V5TB14V7.
