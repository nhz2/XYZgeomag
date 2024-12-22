#include <XYZgeomag.hpp>

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
  geomag::Vector mag_field = geomag::GeoMag(2027.5,position,geomag::WMM2025);
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