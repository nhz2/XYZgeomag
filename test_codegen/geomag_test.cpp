
// geomag_test.cpp Generated by python script wmmcodeupdate.py
/** \file
 * \author Nathan Zimmerberg (nhz2@cornell.edu)
 * \date 24 OCT 2019
 * \brief c++ catch2 tests for magnetic field header only library.
 * \details Compile with g++ geomag_test.cpp -std=c++1z
 */
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../geomag.hpp"


TEST_CASE( "geomag test 0 of WMM2015 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= 1111164.8708100126;
    in.y= 0.0;
    in.z= 6259542.961028692;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -1.5978489161206863e-05;
    truth.y= -4.459e-07;
    truth.z= -5.24545672160699e-05;
    out= geomag::GeoMag(2015.0,in,geomag::WMM2015);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 1 of WMM2015 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -3189068.4999999986;
    in.y= 5523628.670817469;
    in.z= 0.0;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -5.966461381146904e-06;
    truth.y= 9.54841425354402e-06;
    truth.z= 3.95182e-05;
    out= geomag::GeoMag(2015.0,in,geomag::WMM2015);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 2 of WMM2015 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -555582.4354050067;
    in.y= -962297.0059143244;
    in.z= -6259542.961028692;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= 6.2002473589310396e-06;
    truth.y= -2.078305655483668e-05;
    truth.z= -5.110844938203984e-05;
    out= geomag::GeoMag(2015.0,in,geomag::WMM2015);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 3 of WMM2015 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= 1128529.6885767058;
    in.y= 0.0;
    in.z= 6358023.736329913;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -1.5294927111859905e-05;
    truth.y= -4.7160000000000007e-07;
    truth.z= -5.0379237600155215e-05;
    out= geomag::GeoMag(2015.0,in,geomag::WMM2015);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 4 of WMM2015 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -3239068.4999999986;
    in.y= 5610231.211195913;
    in.z= 0.0;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -5.702279657139047e-06;
    truth.y= 9.147838085131275e-06;
    truth.z= 3.75356e-05;
    out= geomag::GeoMag(2015.0,in,geomag::WMM2015);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 5 of WMM2015 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -564264.8442883533;
    in.y= -977335.3792323681;
    in.z= -6358023.736329913;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= 5.671826519155492e-06;
    truth.y= -1.9759108297106144e-05;
    truth.z= -4.8638531279838574e-05;
    out= geomag::GeoMag(2015.0,in,geomag::WMM2015);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 6 of WMM2015 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= 1111164.8708100126;
    in.y= 0.0;
    in.z= 6259542.961028692;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -1.5955881122427664e-05;
    truth.y= -3.1710000000000005e-07;
    truth.z= -5.2485868599147296e-05;
    out= geomag::GeoMag(2017.5,in,geomag::WMM2015);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 7 of WMM2015 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -3189068.4999999986;
    in.y= 5523628.670817469;
    in.z= 0.0;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -5.707740652342035e-06;
    truth.y= 9.44109680628274e-06;
    truth.z= 3.9571400000000005e-05;
    out= geomag::GeoMag(2017.5,in,geomag::WMM2015);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 8 of WMM2015 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -555582.4354050067;
    in.y= -962297.0059143244;
    in.z= -6259542.961028692;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= 6.200232507413443e-06;
    truth.y= -2.082368227841973e-05;
    truth.z= -5.086747774395191e-05;
    out= geomag::GeoMag(2017.5,in,geomag::WMM2015);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 9 of WMM2015 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= 1128529.6885767058;
    in.y= 0.0;
    in.z= 6358023.736329913;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -1.5275465230606787e-05;
    truth.y= -3.4850000000000003e-07;
    truth.z= -5.040592252432766e-05;
    out= geomag::GeoMag(2017.5,in,geomag::WMM2015);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 10 of WMM2015 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -3239068.4999999986;
    in.y= 5610231.211195913;
    in.z= 0.0;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -5.463532322092838e-06;
    truth.y= 9.044115570659571e-06;
    truth.z= 3.75855e-05;
    out= geomag::GeoMag(2017.5,in,geomag::WMM2015);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 11 of WMM2015 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -564264.8442883533;
    in.y= -977335.3792323681;
    in.z= -6358023.736329913;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= 5.670862799287426e-06;
    truth.y= -1.979537750888189e-05;
    truth.z= -4.841398189658139e-05;
    out= geomag::GeoMag(2017.5,in,geomag::WMM2015);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 12 of WMM2015v2 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= 1111164.8708100126;
    in.y= 0.0;
    in.z= 6259542.961028692;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -1.5983781467503073e-05;
    truth.y= -4.519e-07;
    truth.z= -5.242987305696158e-05;
    out= geomag::GeoMag(2015.0,in,geomag::WMM2015v2);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 13 of WMM2015v2 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -3189068.4999999986;
    in.y= 5523628.670817469;
    in.z= 0.0;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -5.94149779500938e-06;
    truth.y= 9.535576054014707e-06;
    truth.z= 3.9521100000000004e-05;
    out= geomag::GeoMag(2015.0,in,geomag::WMM2015v2);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 14 of WMM2015v2 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -555582.4354050067;
    in.y= -962297.0059143244;
    in.z= -6259542.961028692;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= 6.19831311928931e-06;
    truth.y= -2.0782406756170174e-05;
    truth.z= -5.111650149224161e-05;
    out= geomag::GeoMag(2015.0,in,geomag::WMM2015v2);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 15 of WMM2015v2 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= 1128529.6885767058;
    in.y= 0.0;
    in.z= 6358023.736329913;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -1.530029434513461e-05;
    truth.y= -4.776000000000001e-07;
    truth.z= -5.035727188125109e-05;
    out= geomag::GeoMag(2015.0,in,geomag::WMM2015v2);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 16 of WMM2015v2 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -3239068.4999999986;
    in.y= 5610231.211195913;
    in.z= 0.0;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -5.679611519268714e-06;
    truth.y= 9.135175718626881e-06;
    truth.z= 3.75381e-05;
    out= geomag::GeoMag(2015.0,in,geomag::WMM2015v2);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 17 of WMM2015v2 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -564264.8442883533;
    in.y= -977335.3792323681;
    in.z= -6358023.736329913;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= 5.66973929331642e-06;
    truth.y= -1.9758323478306285e-05;
    truth.z= -4.864577817902017e-05;
    out= geomag::GeoMag(2015.0,in,geomag::WMM2015v2);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 18 of WMM2015v2 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= 1111164.8708100126;
    in.y= 0.0;
    in.z= 6259542.961028692;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -1.5969771836563245e-05;
    truth.y= -2.987e-07;
    truth.z= -5.253124588488371e-05;
    out= geomag::GeoMag(2017.5,in,geomag::WMM2015v2);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 19 of WMM2015v2 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -3189068.4999999986;
    in.y= 5523628.670817469;
    in.z= 0.0;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -5.752448209374811e-06;
    truth.y= 9.458932566545789e-06;
    truth.z= 3.95694e-05;
    out= geomag::GeoMag(2017.5,in,geomag::WMM2015v2);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 20 of WMM2015v2 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -555582.4354050067;
    in.y= -962297.0059143244;
    in.z= -6259542.961028692;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= 6.188200185175059e-06;
    truth.y= -2.0809922871869646e-05;
    truth.z= -5.0886998808291256e-05;
    out= geomag::GeoMag(2017.5,in,geomag::WMM2015v2);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 21 of WMM2015v2 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= 1128529.6885767058;
    in.y= 0.0;
    in.z= 6358023.736329913;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -1.5287039032881014e-05;
    truth.y= -3.3110000000000005e-07;
    truth.z= -5.044967749091339e-05;
    out= geomag::GeoMag(2017.5,in,geomag::WMM2015v2);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 22 of WMM2015v2 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -3239068.4999999986;
    in.y= 5610231.211195913;
    in.z= 0.0;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -5.50437218767199e-06;
    truth.y= 9.062452292816946e-06;
    truth.z= 3.7584400000000006e-05;
    out= geomag::GeoMag(2017.5,in,geomag::WMM2015v2);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 23 of WMM2015v2 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -564264.8442883533;
    in.y= -977335.3792323681;
    in.z= -6358023.736329913;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= 5.660068276320213e-06;
    truth.y= -1.978267417110458e-05;
    truth.z= -4.8431724598834037e-05;
    out= geomag::GeoMag(2017.5,in,geomag::WMM2015v2);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 24 of WMM2020 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= 1111164.8708100126;
    in.y= 0.0;
    in.z= 6259542.961028692;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -1.5952813250071817e-05;
    truth.y= -1.4630000000000003e-07;
    truth.z= -5.263547417444184e-05;
    out= geomag::GeoMag(2020.0,in,geomag::WMM2020);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 25 of WMM2020 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -3189068.4999999986;
    in.y= 5523628.670817469;
    in.z= 0.0;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -5.561426191875908e-06;
    truth.y= 9.412872726873378e-06;
    truth.z= 3.9624300000000006e-05;
    out= geomag::GeoMag(2020.0,in,geomag::WMM2020);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 26 of WMM2020 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -555582.4354050067;
    in.y= -962297.0059143244;
    in.z= -6259542.961028692;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= 6.177267161005052e-06;
    truth.y= -2.084485942521248e-05;
    truth.z= -5.0651924360034926e-05;
    out= geomag::GeoMag(2020.0,in,geomag::WMM2020);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 27 of WMM2020 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= 1128529.6885767058;
    in.y= 0.0;
    in.z= 6358023.736329913;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -1.5270886859529107e-05;
    truth.y= -1.8550000000000001e-07;
    truth.z= -5.054523400453757e-05;
    out= geomag::GeoMag(2020.0,in,geomag::WMM2020);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 28 of WMM2020 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -3239068.4999999986;
    in.y= 5610231.211195913;
    in.z= 0.0;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -5.328246064856985e-06;
    truth.y= 9.01899289956124e-06;
    truth.z= 3.76367e-05;
    out= geomag::GeoMag(2020.0,in,geomag::WMM2020);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 29 of WMM2020 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -564264.8442883533;
    in.y= -977335.3792323681;
    in.z= -6358023.736329913;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= 5.6493843086129155e-06;
    truth.y= -1.9813979346000043e-05;
    truth.z= -4.821266111748948e-05;
    out= geomag::GeoMag(2020.0,in,geomag::WMM2020);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 30 of WMM2020 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= 1111164.8708100126;
    in.y= 0.0;
    in.z= 6259542.961028692;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -1.5931578350356247e-05;
    truth.y= 1.1000000000000001e-09;
    truth.z= -5.274827527831086e-05;
    out= geomag::GeoMag(2022.5,in,geomag::WMM2020);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 31 of WMM2020 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -3189068.4999999986;
    in.y= 5523628.670817469;
    in.z= 0.0;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -5.368203727960294e-06;
    truth.y= 9.382401602207892e-06;
    truth.z= 3.9684699999999996e-05;
    out= geomag::GeoMag(2022.5,in,geomag::WMM2020);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 32 of WMM2020 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -555582.4354050067;
    in.y= -962297.0059143244;
    in.z= -6259542.961028692;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= 6.163777504796277e-06;
    truth.y= -2.0877424195142713e-05;
    truth.z= -5.04130265263596e-05;
    out= geomag::GeoMag(2022.5,in,geomag::WMM2020);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 33 of WMM2020 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= 1128529.6885767058;
    in.y= 0.0;
    in.z= 6358023.736329913;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -1.5250661283058837e-05;
    truth.y= -4.45e-08;
    truth.z= -5.0648210584673275e-05;
    out= geomag::GeoMag(2022.5,in,geomag::WMM2020);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 34 of WMM2020 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -3239068.4999999986;
    in.y= 5610231.211195913;
    in.z= 0.0;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= -5.150429303246407e-06;
    truth.y= 8.991405234014355e-06;
    truth.z= 3.7694e-05;
    out= geomag::GeoMag(2022.5,in,geomag::WMM2020);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        
TEST_CASE( "geomag test 35 of WMM2020 model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= -564264.8442883533;
    in.y= -977335.3792323681;
    in.z= -6358023.736329913;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= 5.636486923202327e-06;
    truth.y= -1.984331827281598e-05;
    truth.z= -4.798964104031512e-05;
    out= geomag::GeoMag(2022.5,in,geomag::WMM2020);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(5) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(5) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(5) );
}


        