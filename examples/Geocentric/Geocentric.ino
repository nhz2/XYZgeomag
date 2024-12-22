#include <XYZgeomag.hpp>

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
  out=geomag::GeoMag(2027.5,in,geomag::WMM2025);
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