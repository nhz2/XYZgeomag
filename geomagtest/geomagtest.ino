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
  out=geomag::GeoMag(2015.0,in,geomag::WMM2015);
  int endtime=micros();
  int endtimemil=millis();
  Serial.println(out.x);
  Serial.println(out.y);
  Serial.println(out.z);
  Serial.print("time in micro seconds: ");
  Serial.println(endtime-starttime);
  Serial.print("time in milli seconds: ");
  Serial.println(endtimemil-starttimemil);
  delay(2000);
}
