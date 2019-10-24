# geomag
Lightweight C++ header-only library for calculating the magnetic field on earth. Compatible with CUDA and Arduino.

Calculates the magnetic field around earth in the International Terrestrial Reference System(ITRS).

Unlike most World Magnetic Model(WMM) software, which uses latitude, longitude, and altitude inputs to calculate the North East Down(NED) components of the magnetic field, this library uses geocentric cartesian coordinates as input, and outputs the magnetic field in the same geocentric cartesian coordinate system as the inputs.

## Error

geomag is within 10nT of the official WMM software.

## Performance

geomag uses single precision floating points. It's designed to minimize ram usage for embedded systems.

| Device      | Speed    |
|-------------|----------|
| Arduino Uno | 50 ms    |
| Teensy 3.6  |  80 µs |
|             |          |

## Using geomag

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
  int threadsperblock= 128;
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
`python wmmcodeupdate.py -f WMM2015.COF -f WMM2015v2.COF -o ../geomag.hpp -n 12` from the `test_codegen` directory.

In this example, `WMM2015.COF` and `WMM2015v2.COF` are the `.COF` files to use in `geomag.hpp`.
