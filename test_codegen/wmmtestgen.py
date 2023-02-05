#Nathan Zimmerberg
#24 OCT 2019
#wmmtestgen.py
#This is a simple script to generate c++ catch2 tests.
import pymap3d as pm
import numpy as np


def add_geodetic_2_ecef_tests(outstr, num_tests, seed):
    """Return outstr append with geodetic 2 ecef conversion tests
    
    Args:
        outstr(str): the c++ tests to add to.
        num_tests(int): number of test points to use.
        seed(int): Seed to use to generate test points."""
    rng = np.random.default_rng(seed)
    for i in range(num_tests):
        lat = rng.uniform(-90,90) # test lat (deg)
        lon = rng.uniform(0,360) # test lon (deg)
        height = rng.uniform(-1E6,1E6) # test height (m) with respect to wgs 84 ellipsoid
        x, y, z = pm.geodetic2ecef(lat,lon,height)
        margin = 1.0 # error margin in meters
        testcase="""
TEST_CASE( "geodetic 2 ecef test %s model", "[Geodetic]" ) {
    geomag::Vector out = geomag::geodetic2ecef(%s, %s, %s);
    geomag::Vector truth;
    truth.x= %s;
    truth.y= %s;
    truth.z= %s;
    CHECK( out.x == Approx(truth.x).margin(%s) );
    CHECK( out.y == Approx(truth.y).margin(%s) );
    CHECK( out.z == Approx(truth.z).margin(%s) );
}


        """%(i,repr(lat),repr(lon),repr(height),repr(x),repr(y),repr(z),margin,margin,margin)
        outstr=outstr+testcase
    return outstr


def add_full_tests(
        outstr, #(str): the c++ tests to add to.
        modelnames, #(list of str): The list of model names for each test
        dates, #(numpy array): decimal test years
        heights, #(numpy array): test heights (km) with respect to wgs 84 ellipsoid
        lats, #(numpy array): test lat (deg)
        lons, #(numpy array): test lon (deg)
        norths, #(numpy array): test X(local north) magnetic field (nT)
        easts, #(numpy array): test Y(local east) magnetic field (nT)
        downs, #(numpy array): test Z(local down) magnetic field (nT)
        horizontal_intensities, #(numpy array): test H (nT)
        total_intensities, #(numpy array): test F(nT)
        inclinations, #(numpy array): test I, the angle measured from the horizontal plane to the magnetic field vector; downward  field is positive (deg)
        declinations, #(numpy array): test D, the angle between true north and the horizontal component of the field, positive eastward of true North (deg)
        margin_nT, #(float or int): Acceptable error in each field component (nT)
        margin_deg, #(float or int): Acceptable error in angular components (deg)
    ):
    """Return outstr append full tests."""
    for i in range(len(modelnames)):
        testcase=f"""
TEST_CASE( "full geomag test {i} of {modelnames[i]} model", "[Full GeoMag]" ) {{
    geomag::Elements out = geomag::magField2Elements(
        geomag::GeoMag(
            {dates[i]!r},
            geomag::geodetic2ecef({lats[i]!r}, {lons[i]!r}, {heights[i]*1000.0!r}),
            geomag::{modelnames[i]}
        ),
        {lats[i]!r},
        {lons[i]!r}
    );
    CHECK( out.north                == Approx({norths[i]!r}).margin({margin_nT!r}) );
    CHECK( out.east                 == Approx({easts[i]!r}).margin({margin_nT!r}) );
    CHECK( out.down                 == Approx({downs[i]!r}).margin({margin_nT!r}) );
    CHECK( out.horizontal_intensity == Approx({horizontal_intensities[i]!r}).margin({margin_nT!r}) );
    CHECK( out.total_intensity      == Approx({total_intensities[i]!r}).margin({margin_nT!r}) );
    CHECK( out.inclination          == Approx({inclinations[i]!r}).margin({margin_deg!r}) );
    CHECK( out.declination          == Approx({declinations[i]!r}).margin({margin_deg!r}) );
}}


        """
        outstr=outstr+testcase
    return outstr

def main(header_to_test,modelnames,dates,heights,lats,lons,bns,bes,bds,bhs,bfs,bincs,bdecs,margin_nT,margin_deg):
    """create and write catch2 tests as a c++ source file to test header_to_test
        defining the WMM coefficents and model

     The file created has the name header_to_test[:-4]+'_test.cpp'

    Args:
        header_to_test(str ending in .hpp): the c++ header file with the model and coefficents
        modelnames(list of str): The list of model names for each test
        dates(numpy array): decimal test years
        heights(numpy array): test heights (km) with respect to wgs 84 ellipsoid
        lats(numpy array): test lat (deg)
        lons(numpy array): test lon (deg)
        bns(numpy array): test X(local north) magnetic field (nT)
        bes(numpy array): test Y(local east) magnetic field (nT)
        bds(numpy array): test Z(local down) magnetic field (nT)
        bhs(numpy array): test H(local horizontal) magnetic field intensity (nT)
        bfs(numpy array): test F magnetic field intensity (nT)
        bincs(numpy array): test magnetic field inclination (deg)
        bdecs(numpy array): test magnetic field declination (deg)
        margin_nT(float or int): Acceptable error in each component (nT)
        margin_deg(float or int): Acceptable error in bincs and bdecs (deg)"""
    xs,ys,zs=pm.geodetic2ecef(lats,lons,heights*1000)
    bxs,bys,bzs=pm.enu2uvw(bes*1E-9,bns*1E-9,-bds*1E-9,lats,lons)
    file_name=header_to_test[:-4]+'_test.cpp'
    with open(file_name,'w') as f:
        outstr="""
// %s Generated by python script wmmtestgen.py
/** \\file
 * \\author Nathan Zimmerberg (nhz2@cornell.edu)
 * \\date 24 OCT 2019
 * \\brief c++ catch2 tests for magnetic field header only library.
 * \\details Compile with g++ %s -std=c++1z
 */
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../%s"

"""%(file_name,file_name,header_to_test)

        outstr = add_geodetic_2_ecef_tests(outstr, 100, 1234)

        for i in range(len(xs)):
            testcase="""
TEST_CASE( "geomag test %s of %s model", "[GeoMag]" ) {
    geomag::Vector in;
    in.x= %s;
    in.y= %s;
    in.z= %s;
    geomag::Vector out;
    geomag::Vector truth;
    truth.x= %s;
    truth.y= %s;
    truth.z= %s;
    out= geomag::GeoMag(%s,in,geomag::%s);
    CHECK( out.x*1E9 == Approx(truth.x*1E9).margin(%s) );
    CHECK( out.y*1E9 == Approx(truth.y*1E9).margin(%s) );
    CHECK( out.z*1E9 == Approx(truth.z*1E9).margin(%s) );
}


        """%(i,modelnames[i],repr(xs[i]),repr(ys[i]),repr(zs[i]),repr(bxs[i]),repr(bys[i]),repr(bzs[i]),dates[i],modelnames[i],margin_nT,margin_nT,margin_nT)
            outstr=outstr+testcase

        outstr = add_full_tests(outstr, modelnames, dates, heights, lats, lons, bns, bes, bds, bhs, bfs, bincs, bdecs, margin_nT, margin_deg)
        f.write(outstr)

if __name__ == '__main__':
    header_to_test='geomag.hpp'
    dates=([2015]*6+[2017.5]*6)*2 + [2020]*6 + [2022.5]*6
    modelnames= ['WMM2015']*12 + ['WMM2015v2']*12 + ['WMM2020']*12
    heights=([0]*3+[100]*3+[0]*3+[100]*3)*3
    lats=[80,0,-80]*4*3
    lons=[0,120,240]*4*3

    bns=[6627.1,39518.2,5797.3,6314.3,37535.6,5613.1,6599.4,39571.4,5873.8,6290.5,37585.5,5683.5]
    bes=[-445.9,392.9,15761.1,-471.6,364.4,14791.5,-317.1,222.5,15781.4,-348.5,209.5,14808.8]
    bds=[54432.3,-11252.4,-52919.1,52269.8,-10773.4,-50378.6,54459.2,-11030.1,-52687.9,52292.7,-10564.2,-50163.0]
    bhs=[6642.1,39520.2,16793.5,6331.9,37537.3,15820.7,6607.0,39572.0,16839.1,6300.1,37586.1,15862.0]
    bfs=[54836.0,41090.9,55519.8,52652.0,39052.7,52804.4,54858.5,41080.5,55313.4,52670.9,39042.5,52611.1]
    bincs=[83.04,-15.89,-72.39,83.09,-16.01,-72.57,83.08,-15.57,-72.28,83.13,-15.70,-72.45]
    bdecs=[-3.85,0.57,69.81,-4.27,0.56,69.22,-2.75,0.32,69.58,-3.17,0.32,69.00]

    bns+=[6636.6,39521.1,5796.3,6323.4,37538.1,5612.2,6605.2,39569.4,5864.6,6294.3,37584.4,5674.9]
    bes+=[-451.9,377.7,15759.1,-477.6,351.1,14789.3,-298.7,252.3,15764.1,-331.1,235.7,14793.1]
    bds+=[54408.9,-11228.8,-52927.1,52249.1,-10751.1,-50385.8,54506.3,-11067.9,-52706.1,52337.8,-10600.5,-50179.5]
    bhs+=[6651.9,39522.9,16791.2,6341.4,37539.7,15818.3,6612.0,39570.2,16819.7,6303.0,37585.1,15844.2]
    bfs+=[54814.0,41087.1,55526.8,52632.5,39048.9,52810.5,54905.9,41088.9,55324.8,52716.0,39051.4,52621.5]
    bincs+=[83.03,-15.86,-72.40,83.08,-15.98,-72.57,83.08,-15.63,-72.30,83.13,-15.75,-72.48]
    bdecs+=[-3.90,0.55,69.81,-4.32,0.54,69.22,-2.59,0.37,69.59,-3.01,0.36,69.01]
    

    bns+=[6570.4,39624.3,5940.6,6261.8,37636.7,5744.9,6529.9,39684.7,6016.5,6224.0,37694.0,5815.0]
    bes+=[-146.3,109.9,15772.1,-185.5,104.9,14799.5,1.1,-42.2,15776.7,-44.5,-35.3,14803.0]
    bds+=[54606.0,-10932.5,-52480.8,52429.1,-10474.8,-49969.4,54713.4,-10809.5,-52251.6,52527.0,-10362.0,-49755.3]
    bhs+=[6572.0,39624.4,16853.8,6264.5,37636.9,15875.4,6529.9,39684.7,16885.0,6224.2,37694.1,15904.1]
    bfs+=[55000.1,41104.9,55120.6,52802.0,39067.3,52430.6,55101.7,41130.5,54912.1,52894.5,39092.4,52235.4]
    bincs+=[83.14,-15.42,-72.20,83.19,-15.55,-72.37,83.19,-15.24,-72.09,83.24,-15.37,-72.27]
    bdecs+=[-1.28,0.16,69.36,-1.70,0.16,68.78,0.01,-0.06,69.13,-0.41,-0.05,68.55]

    margin_nT= 0.1#nT error acceptable
    margin_deg = 1E-2#degrees acceptable
    main(
        header_to_test,
        modelnames,
        np.array(dates),
        np.array(heights),
        np.array(lats),
        np.array(lons),
        np.array(bns),
        np.array(bes),
        np.array(bds),
        np.array(bhs),
        np.array(bfs),
        np.array(bincs),
        np.array(bdecs),
        margin_nT,
        margin_deg,
    )
