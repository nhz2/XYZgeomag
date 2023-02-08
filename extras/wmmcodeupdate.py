# MIT License
#
# Copyright (c) 2019 Nathan Zimmerberg
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#23 OCT 2019
#wmmcodeupdate.py
#This is a simple script to parse WMM.COF into a c++ header file
#To run, use command "python wmmcodeupdate.py", add -h flag for help
#it un Schmidt semi-normalizes the coefficents and stores them by:
#    index=((2*maxdegree-m+1)*m)/2+n
import math

def main(infilenames, headerfilename, maxdegree):
    """parse infilenames into headerfilename c++ header file

    The coefficents are un Schmidt semi-normalized.
    indexing is by ((2*maxdegree-m+1)*m)/2+n
    Args:
        infilenames(list of filenames): the .COF files that contains the
            WMM coefficents, download this from https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml
        headerfilename(string ending in .hpp): the c++ header file.
            Stored in a ConstModel, index=((2*maxdegree-m+1)*m)/2+n
                typedef struct {
                        const float epoch;
                        const float Main_Field_Coeff_C[NUMCOF];
                        const float Main_Field_Coeff_S[NUMCOF];
                        const float Secular_Var_Coeff_C[NUMCOF];
                        const float Secular_Var_Coeff_S[NUMCOF];
                } ConstModel;
        maxdegree(positive integer): maximum degree"""
    outstr= header_file_header(headerfilename,maxdegree)
    for infilename in infilenames:
        data= parseescof(infilename,maxdegree)
        dyear= data[0]
        cofs= data[1:]
        c_cofs=[0]*(((maxdegree+1)*(maxdegree+2))//2)
        s_cofs=[0]*(((maxdegree+1)*(maxdegree+2))//2)
        c_secvars=[0]*(((maxdegree+1)*(maxdegree+2))//2)
        s_secvars=[0]*(((maxdegree+1)*(maxdegree+2))//2)
        for cof in cofs:
            n= cof[0]
            m= cof[1]
            g= cof[2]
            h= cof[3]
            gsec= cof[4]
            hsec= cof[5]
            if (m==0):
                unnorm= 1.0
            else:
                unnorm= math.sqrt(2.0*float(math.factorial(n-m))/float(math.factorial(n+m)))
            c_cofs[((2*maxdegree-m+1)*m)//2+n]= g*unnorm
            s_cofs[((2*maxdegree-m+1)*m)//2+n]= h*unnorm
            c_secvars[((2*maxdegree-m+1)*m)//2+n]= gsec*unnorm
            s_secvars[((2*maxdegree-m+1)*m)//2+n]= hsec*unnorm
        outstr = outstr + header_file_model_code(infilename[:-4],dyear,c_cofs,s_cofs,c_secvars,s_secvars)
    outstr = outstr + "}\n#endif /* GEOMAG_HPP */"
    with open(headerfilename,'w') as f:
        f.write(outstr)


def header_file_header(headerfilename,maxdegree):
    """returns the header of the header file

    Args:
        headerfilename(string ending in .hpp): the c++ header file with the coefficents
            Stored in a ConstModel, index=((2*maxdegree-m+1)*m)/2+n
                typedef struct {
                        const float epoch;
                        const float Main_Field_Coeff_C[NUMCOF];
                        const float Main_Field_Coeff_S[NUMCOF];
                        const float Secular_Var_Coeff_C[NUMCOF];
                        const float Secular_Var_Coeff_S[NUMCOF];
                } ConstModel;
        maxdegree(positive integer): maximum degree"""
    head="""/*
MIT License

Copyright (c) 2019 Nathan Zimmerberg

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


// %s Generated by python script wmmcodeupdate.py
/** \\file
 * \\author Nathan Zimmerberg (nhz2@cornell.edu)
 * \\date 07 AUG 2021
 * \\brief Header-only library to calculate the magnetic field in the International Terrestrial Reference System(ITRS).
 * \\details Designed to minimize ram usage for embedded systems.

The data from WMM models is not subject to copyright protection.
Modifications are:
  using ITRS coordinates,
  conversion from nT to T,
  Using unnormalized coefficents genrated by the python script wmmcodeupdate.py
  using a different spherical harmonics calculation, described in sections 3.2.4 and 3.2.5:
    Satellite Orbits Models, Methods and Applications,
      by Oliver Montenbruck and Eberhard Gill 2000
*/
#ifndef GEOMAG_HPP
#define GEOMAG_HPP

#include <math.h>

namespace geomag
{
constexpr int NMAX= %d;//order of the Model
constexpr int NUMCOF= (NMAX+1)*(NMAX+2)/2;//number of coefficents
struct ConstModel{
    float epoch;//decimal year
    float Main_Field_Coeff_C[NUMCOF];
    float Main_Field_Coeff_S[NUMCOF];
    float Secular_Var_Coeff_C[NUMCOF];
    float Secular_Var_Coeff_S[NUMCOF];
    /** Function for indexing the C spherical component n,m at dyear time.*/
    inline float C(int n, int m, float dyear) const{
      int index= (m*(2*NMAX-m+1))/2+n;
      #ifdef PROGMEM
        return pgm_read_float_near(Main_Field_Coeff_C+index)+(dyear-epoch)*pgm_read_float_near(Secular_Var_Coeff_C+index);
      #endif /* PROGMEM */
      return Main_Field_Coeff_C[index]+(dyear-epoch)*Secular_Var_Coeff_C[index];
    }
    /** Function for indexing the S spherical component n,m at dyear time.*/
    inline float S(int n, int m, float dyear) const{
      int index= (m*(2*NMAX-m+1))/2+n;
      #ifdef PROGMEM
        return pgm_read_float_near(Main_Field_Coeff_S+index)+(dyear-epoch)*pgm_read_float_near(Secular_Var_Coeff_S+index);
      #endif /* PROGMEM */
      return Main_Field_Coeff_S[index]+(dyear-epoch)*Secular_Var_Coeff_S[index];
    }
};
//mean radius of  ellipsoid in meters from section 1.2 of the WMM2015 Technical report
constexpr float EARTH_R= 6371200.0f;

typedef struct {
    float x;
    float y;
    float z;
} Vector;


typedef struct {
    float north;// local north magnetic field (nT)
    float east;// local east magnetic field (nT)
    float down;// local down magnetic field (nT)
    float horizontal;// local horizontal magnetic field intensity (nT)
    float total;// local total magnetic field intensity (nT)
    float inclination;// also called the dip angle, 
    // the angle measured from the horizontal plane to the 
    // magnetic field vector; a downward field is positive (deg)
    float declination;// also called the magnetic variation,
    // the angle between true north and the horizontal component 
    // of the field, a eastward magnetic field of true North is positive (deg)
} Elements;

/** Return a struct containing the 7 magnetic elements.
See https://www.geomag.nrcan.gc.ca/mag_fld/comp-en.php and
https://www.ngdc.noaa.gov/geomag/icons/faqelems.gif for more info.
 INPUT:
    mag_field_itrs: local magnetic field in the itrs coordinate system (T)
    lat: latitude in degrees, -90 at the south pole, 90 at the north pole.
    lon: longitude in degrees.
**/
inline Elements magField2Elements(Vector mag_field_itrs, float lat, float lon){
    float x = mag_field_itrs.x*1E9f;
    float y = mag_field_itrs.y*1E9f;
    float z = mag_field_itrs.z*1E9f;
    float phi = lat*((float)(M_PI/180.0));
    float lam = lon*((float)(M_PI/180.0));
    float sphi = sinf(phi);
    float cphi = cosf(phi);
    float slam = sinf(lam);
    float clam = cosf(lam);
    float x1 = clam*x + slam*y;
    float north = -sphi*x1 + cphi*z;
    float east = -slam*x + clam*y;
    float down = -cphi*x1 + -sphi*z;
    float horizontal = sqrtf(north*north + east*east);
    float total = sqrtf(horizontal*horizontal + down*down);
    float inclination = atan2f(down, horizontal)*((float)(180.0/M_PI));
    float declination = atan2f(east, north)*((float)(180.0/M_PI));
    return {north, east, down, horizontal, total, inclination, declination};
}


/** Return the position in International Terrestrial Reference System coordinates, units meters.
Using the WGS 84 ellipsoid and the algorithm from https://geographiclib.sourceforge.io/
 INPUT:
    lat: Geodetic latitude in degrees, -90 at the south pole, 90 at the north pole.
    lon: Geodetic longitude in degrees.
    h: Height above the WGS 84 ellipsoid in meters.
**/
inline Vector geodetic2ecef(float lat, float lon, float h){
    // Convert to radians
    float phi = lat*((float)(M_PI/180.0));
    float lam = lon*((float)(M_PI/180.0));
    // WGS 84 constants
    const float a = 6378137;
    // const float f = 1.0/298.257223563;
    const float e2 = 0.0066943799901413165;//f*(2-f);
    const float e2m = 0.9933056200098587;//(1-f)*(1-f);
    float sphi = sinf(phi);
    float cphi = cosf(phi);
    float slam = sinf(lam);
    float clam = cosf(lam);
    float n = a/sqrtf(1.0f - e2*(sphi*sphi));
    float z = (e2m*n + h) * sphi;
    float r = (n + h) * cphi;
    return {r*clam, r*slam, z};
}


/** Return the magnetic field in International Terrestrial Reference System coordinates, units Tesla.
 INPUT:
    position_itrs(Above the surface of earth): The location where the field is predicted, units m.
    dyear(should be around the epoch of the model): The decimal year, for example 2015.0
    WMM(): Magnetic field model to use.
 */
inline Vector GeoMag(float dyear,Vector position_itrs, const ConstModel& WMM){
    float x= position_itrs.x;
    float y= position_itrs.y;
    float z= position_itrs.z;
    float px= 0;
    float py= 0;
    float pz= 0;
    float rsqrd= x*x+y*y+z*z;
    float temp= EARTH_R/rsqrd;
    float a= x*temp;
    float b= y*temp;
    float f= z*temp;
    float g= EARTH_R*temp;

    int n,m;
    //first m==0 row, just solve for the Vs
    float Vtop= EARTH_R/sqrtf(rsqrd);//V0,0
    float Wtop= 0;//W0,0
    float Vprev= 0;
    float Wprev= 0;
    float Vnm= Vtop;
    float Wnm= Wtop;

    //iterate through all ms
    for ( m = 0; m <= NMAX+1; m++)
    {
        // iterate through all ns
        for (n = m; n <= NMAX+1; n++)
        {
            if (n==m){
                if(m!=0){
                    temp= Vtop;
                    Vtop= (2*m-1)*(a*Vtop-b*Wtop);
                    Wtop= (2*m-1)*(a*Wtop+b*temp);
                    Vprev= 0;
                    Wprev= 0;
                    Vnm= Vtop;
                    Wnm= Wtop;
                }
            }
            else{
                temp= Vnm;
                float invs_temp=1.0f/((float)(n-m));
                Vnm= ((2*n-1)*f*Vnm - (n+m-1)*g*Vprev)*invs_temp;
                Vprev= temp;
                temp= Wnm;
                Wnm= ((2*n-1)*f*Wnm - (n+m-1)*g*Wprev)*invs_temp;
                Wprev= temp;
            }
            if (m<NMAX && n>=m+2){
                px+= 0.5f*(n-m)*(n-m-1)*(WMM.C(n-1,m+1,dyear)*Vnm+WMM.S(n-1,m+1,dyear)*Wnm);
                py+= 0.5f*(n-m)*(n-m-1)*(-WMM.C(n-1,m+1,dyear)*Wnm+WMM.S(n-1,m+1,dyear)*Vnm);
            }
            if (n>=2 && m>=2){
                px+= 0.5f*(-WMM.C(n-1,m-1,dyear)*Vnm-WMM.S(n-1,m-1,dyear)*Wnm);
                py+= 0.5f*(-WMM.C(n-1,m-1,dyear)*Wnm+WMM.S(n-1,m-1,dyear)*Vnm);
            }
            if (m==1 && n>=2){
                px+= -WMM.C(n-1,0,dyear)*Vnm;
                py+= -WMM.C(n-1,0,dyear)*Wnm;
            }
            if (n>=2 && n>m){
                pz+= (n-m)*(-WMM.C(n-1,m,dyear)*Vnm-WMM.S(n-1,m,dyear)*Wnm);
            }
        }
    }
    return {-px*1.0E-9f,-py*1.0E-9f,-pz*1E-9f};
}
// Model parameters\n"""%(headerfilename,maxdegree)
    return head

def header_file_model_code(modelname, dyear,c_cofs,s_cofs,c_secvars,s_secvars):
    """return the code defining the WMM coefficents and model
                Stored in a ConstModel, index=((2*maxdegree-m+1)*m)/2+n
                    typedef struct {
                            const float epoch;
                            const float Main_Field_Coeff_C[NUMCOF];
                            const float Main_Field_Coeff_S[NUMCOF];
                            const float Secular_Var_Coeff_C[NUMCOF];
                            const float Secular_Var_Coeff_S[NUMCOF];
                    } ConstModel;

    Args:
        modelname(str, a valid C++ name): name of the model
        dyear(positive float): the year of the magnetic model, ex 2015.0
        c_cofs,s_cofs,c_secvars,s_secvars(list of floats): coefficents of the model"""

    head="""constexpr
#ifdef PROGMEM
    PROGMEM
#endif /* PROGMEM */
ConstModel %s = {%f,\n"""%(modelname,dyear)
    modeltail= '};\n\n'
    cs='{'
    ss='{'
    csec='{'
    ssec='{'
    for i in range(len(c_cofs)):
        cs+= repr(c_cofs[i])+','
        ss+= repr(s_cofs[i])+','
        csec+= repr(c_secvars[i])+','
        ssec+= repr(s_secvars[i])+','
    cs= cs[:-1]+'}'
    ss= ss[:-1]+'}'
    csec= csec[:-1]+'}'
    ssec= ssec[:-1]+'}'
    return head+cs+',\n'+ss+',\n'+csec+',\n'+ssec+modeltail


def parseescof(infilename, maxdegree):
    """return a list of lists from the infilename cof data file
    dyear,
    n,m,G,H,SecG,SecH,
    ...

    Args:
        infilename(string ending in .COF): the .COF file that contains the
            WMM coefficents, download this from https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml
        maxdegree(positive integer): maximum degree"""
    with open(infilename,'r') as f:
        s= f.read()
        a= s.split('\n')
        dyear= float(a[0].split()[0])
        l=[dyear]
        l.append([0,0,0.0,0.0,0.0,0.0])
        for i in range(1,((maxdegree+1)*(maxdegree+2))//2):
            rowa= a[i].split()
            l.append([int(rowa[0]),int(rowa[1]),float(rowa[2]),float(rowa[3]),float(rowa[4]),float(rowa[5])])
        return l



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="""
    #Nathan Zimmerberg
    #07 AUG 2021
    #wmmcodeupdate.py
    #This is a simple script to parse WMM.COF into a c++ header file
    #it un Schmidt semi-normalizes the coefficents and stores them by:
    #    index=((2*maxdegree-m+1)*m)/2+n
    """)
    parser.add_argument('-f',action='append',help="""the .COF files that contains the
        WMM coefficents, download from https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml.""")
    parser.add_argument('-o',type=str,default='../src/XYZgeomag.hpp',help='the c++ header filename to write the coefficents and model')
    parser.add_argument('-n',type=int,default=12,help='maximum number of degrees to use')
    arg=parser.parse_args()

    main(arg.f,arg.o,arg.n)
