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
        horizontal, #(numpy array): test H (nT)
        total, #(numpy array): test F(nT)
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
    CHECK( out.north       == Approx({norths[i]!r}).margin({margin_nT!r}) );
    CHECK( out.east        == Approx({easts[i]!r}).margin({margin_nT!r}) );
    CHECK( out.down        == Approx({downs[i]!r}).margin({margin_nT!r}) );
    CHECK( out.horizontal  == Approx({horizontal[i]!r}).margin({margin_nT!r}) );
    CHECK( out.total       == Approx({total[i]!r}).margin({margin_nT!r}) );
    CHECK( out.inclination == Approx({inclinations[i]!r}).margin({margin_deg!r}) );
    CHECK( out.declination == Approx({declinations[i]!r}).margin({margin_deg!r}) );
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
    file_name='geomag_test.cpp'
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
    header_to_test='src/XYZgeomag.hpp'
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

    dates+= [2025]*10 + [2025.5]*10 + [2026]*10 + [2026.5]*10 + [2027]*10 + [2027.5]*10 + [2028]*10 + [2028.5]*10 + [2029]*10 + [2029.5]*10
    modelnames+=['WMM2025']*100
    heights+=[28, 48, 54, 65, 51, 39, 3, 94, 66, 18, 6, 63, 69, 50, 8, 8, 22, 40, 44, 50, 74, 46, 69, 33, 47, 62, 83, 82, 34, 56, 14, 12, 44, 43, 64, 12, 12, 14, 19, 86, 37, 67, 44, 54, 57, 44, 12, 38, 61, 67, 8, 77, 98, 34, 60, 73, 96, 0, 16, 72, 49, 71, 95, 86, 30, 75, 21, 1, 45, 11, 28, 68, 39, 27, 11, 72, 55, 59, 65, 95, 95, 95, 50, 58, 57, 38, 49, 90, 41, 19, 31, 93, 51, 64, 26, 66, 18, 63, 33, 77]
    lats+=[89, 80, 82, 43, -33, -59, -50, -29, 14, 0, -36, 26, 38, -70, -52, -66, -37, -12, 33, -81, -57, -24, 23, -3, -72, -14, 86, -64, -19, -81, 0, -82, -46, 17, 10, 33, -79, -33, 29, -11, -66, 72, 22, 54, -43, -43, -63, 27, 59, -47, 62, -68, -5, -29, 27, -72, -46, -13, 66, -87, 20, 5, 14, -85, -36, 79, 6, -76, -46, -22, 54, -58, -65, -23, 34, -62, 86, 32, 48, 30, -60, -70, 87, 32, 34, -76, -50, -55, 42, 46, 13, -2, -76, 22, -65, -21, 9, 88, 17, -18]
    lons+=[-121, -96, 87, 93, 109, -8, -103, -110, 143, 21, -137, 81, -144, -133, -75, 17, 140, -129, -118, -67, 3, -122, 63, -147, -22, 99, -46, 87, 43, 40, 80, -68, -42, 52, 78, -145, 115, -114, 66, 167, -5, -115, 174, 178, 50, -111, 178, -169, -77, -32, 53, -7, 159, -107, 65, 95, -85, -59, -178, 38, 167, -13, 65, -79, -64, 125, -32, -75, -41, -21, -120, 156, -88, 81, 0, 65, 70, 163, 148, 28, -59, 42, -154, 19, -13, 49, -179, -171, -19, -22, -132, 158, 40, -132, 55, 32, -172, 26, 5, 138]
    bns+=[-255.388723, 1875.982280, 1324.336929, 24299.852822, 21737.778822, 14358.095523, 19526.532799, 23275.471080, 35003.441364, 29274.811882, 23781.930678, 34802.613433, 22510.804524, 9021.847823, 19331.857134, 15201.438652, 21404.761773, 28594.451238, 23235.835850, 16132.682326, 13268.119649, 25846.118652, 34558.605297, 30514.086409, 18280.800055, 33437.828630, 2582.099403, 2007.074870, 19525.694050, 9000.536024, 39431.415992, 16052.819558, 13863.176808, 35994.907968, 39427.724480, 24136.839639, -9613.650418, 23393.133584, 32725.000691, 32578.377044, 16390.771268, 4883.287833, 28683.972510, 20615.871559, 11203.748204, 20471.522067, 6233.774935, 25821.061601, 10437.775689, 12450.832543, 12260.425578, 16578.099495, 33544.960996, 23540.396213, 33062.159421, -2912.418045, 18947.658009, 21365.849306, 13821.326215, 6926.051596, 30131.436134, 28142.647915, 36932.466767, 12713.115116, 17341.333537, 2447.595900, 27567.094260, 16993.742946, 13609.973828, 12430.606447, 14735.145798, 6944.959164, 17946.905687, 24940.672257, 29078.234597, 6567.891205, 901.741754, 28217.104175, 23365.167010, 29692.133298, 17893.018948, 10426.249753, 257.877420, 29359.395781, 28241.863242, 7991.121323, 15315.042731, 12821.507915, 24439.977735, 22522.142232, 27795.620707, 33806.520595, 10263.650094, 25601.266218, 8347.198524, 15578.335534, 30520.999126, 2041.140972, 34021.976056, 31751.497580]
    bes+=[-1482.460628, -1079.269389, 1883.428620, 210.517066, -2090.274098, -4049.107540, 10362.990980, 6559.046509, -116.624818, 659.800118, 8786.698927, 308.431881, 5167.636206, 14001.865232, 5147.774727, -9925.501123, 3498.897430, 5432.201489, 4558.379362, 8623.494405, -5498.179626, 6448.863284, 708.078839, 5221.840663, -2024.105320, -835.919473, -1527.839829, -13829.362571, -5224.017527, -15446.900890, -2133.456168, 9190.416754, -2785.834358, 747.876551, -1049.971869, 5112.225422, -8785.714482, 7653.110952, 1278.343380, 5882.794932, -5079.626554, 1192.857435, 3249.857362, 225.041793, -12563.437494, 9245.505620, 9925.455386, 3851.701313, -3087.391846, -2994.148743, 4315.497527, -4813.676173, 4589.735713, 6592.363667, 1068.555172, -12984.118939, 6129.590452, -6732.560620, 89.244328, -15153.941754, 2689.876669, -3192.970202, -329.910612, 11085.480728, -1411.102609, -823.235047, -7045.034525, 9761.147808, -2812.623735, -5337.987778, 4066.959949, 6159.591995, 10131.662696, -5883.907569, 798.434048, -16150.440331, 2192.139033, 75.200106, -3931.047029, 2367.965026, 2700.105868, -14920.796381, -869.484879, 2108.188210, -933.225495, -16588.167977, 9610.272521, 10252.398258, -1762.893879, -2227.393290, 4421.061344, 4207.156292, -15412.758102, 4619.955329, -16728.868464, -4066.430831, 4966.941695, 1511.567287, 531.446418, 2472.257962]
    bds+=[56194.288771, 55623.044051, 56740.772059, 50037.923998, -52710.003920, -24389.086374, -31437.562789, -19063.605287, 7966.315182, -14316.722540, -32577.518648, 30332.056989, 35525.990264, -51084.838301, -23532.664273, -30881.123197, -55397.587760, -8052.525092, 37727.715752, -44412.283868, -23576.062921, -18080.399376, 25043.405885, -1146.978999, -33397.486160, -33100.922253, 54279.308437, -53663.514288, -26182.862330, -45267.890393, -12188.838466, -45940.146795, -19744.304022, 15990.895993, 5214.811025, 32162.934508, -58104.306533, -23854.293436, 33958.454002, -20368.224196, -28608.243575, 55689.615237, 17961.198575, 45149.631876, -33221.366617, -29347.556906, -61075.730989, 24058.365347, 54397.713552, -20475.298079, 54849.810218, -29680.383821, -14525.780052, -18723.097084, 30667.425573, -55399.217725, -21631.434316, -6112.313520, 54092.646409, -48295.635435, 15295.611788, -9017.928693, 11601.180703, -46988.258092, -14639.967305, 57308.751756, -4352.046687, -42479.847250, -19816.196211, -21373.837958, 52393.678137, -62264.415929, -35982.872979, -41948.733657, 30945.577813, -44267.327811, 55926.154052, 26405.862611, 44177.553125, 29039.808439, -26011.845842, -38237.047006, 55992.308183, 30508.465542, 28997.458189, -44151.303941, -53522.651338, -52995.143972, 36929.897501, 40651.688191, 17187.983590, -10964.647547, -42018.541264, 24914.378368, -41431.528035, -24494.877642, 8779.472279, 55286.620082, 8341.137769, -34817.395113]
    bhs+=[1504.298146, 2164.285547, 2302.427342, 24300.764692, 21838.046477, 14918.115796, 22106.041373, 24181.990098, 35003.635649, 29282.246275, 25353.230658, 34803.980117, 23096.337032, 16656.709403, 20005.506364, 18154.870136, 21688.847590, 29105.866326, 23678.743422, 18292.842720, 14362.206593, 26638.500090, 34565.858527, 30957.666082, 18392.516222, 33448.275663, 3000.255301, 13974.248411, 20212.448819, 17877.818542, 39489.089663, 18497.480257, 14140.316272, 36002.676553, 39441.702532, 24672.289649, 13023.480845, 24613.183584, 32749.959268, 33105.255278, 17159.836500, 5026.868699, 28867.487799, 20617.099795, 16833.417225, 22462.470699, 11720.691727, 26106.758229, 10884.812802, 12805.786103, 12997.751893, 17262.817301, 33857.496691, 24446.053109, 33079.422542, 13306.747292, 19914.457641, 22401.493010, 13821.614337, 16661.696834, 30251.262453, 28323.200567, 36933.940251, 16867.459172, 17398.651081, 2582.332595, 28453.070088, 19597.635210, 13897.562372, 13528.270036, 15286.094495, 9282.943032, 20609.270068, 25625.329284, 29089.194286, 17434.847799, 2370.361097, 28217.204381, 23693.546804, 29786.406936, 18095.598879, 18202.660480, 906.920459, 29434.989012, 28257.277809, 18412.640681, 18080.593790, 16416.538468, 24503.475397, 22632.016516, 28145.022897, 34067.301020, 18517.441118, 26014.780783, 18695.741849, 16100.322907, 30922.514410, 2539.900024, 34026.126581, 31847.600506]
    bfs+=[56214.419888, 55665.134163, 56787.466800, 55626.621348, 57054.752538, 28589.818346, 38431.724126, 30792.688931, 35898.700342, 32594.761714, 41280.516301, 46166.554053, 42373.774537, 53731.803175, 30886.996822, 35822.382382, 59492.006517, 30199.248583, 44542.826874, 48032.062762, 27606.226129, 32194.883578, 42684.549360, 30978.906535, 38127.112857, 47058.030120, 54362.163830, 55453.154865, 33076.961273, 48670.301997, 41327.424134, 49524.275496, 24285.511845, 39394.180708, 39784.948821, 40536.110231, 59545.961164, 34275.882505, 47177.711160, 38869.300019, 33360.029814, 55916.032175, 33999.066253, 49634.202547, 37242.759503, 36957.295441, 62190.188377, 35501.658671, 55476.034370, 24150.072239, 56368.814385, 34335.550744, 36841.937629, 30792.269761, 45108.083389, 56974.931751, 29402.458634, 23220.406233, 55830.559897, 51088.947370, 33898.298187, 29724.177504, 38713.089985, 49924.018042, 22738.551012, 57366.902212, 28783.980055, 46782.525885, 24203.798713, 25295.356080, 54578.037650, 62952.605366, 41466.964689, 49156.421313, 42471.284539, 47576.992646, 55976.363930, 38645.571587, 50130.233993, 41599.765773, 31687.013474, 42348.655378, 55999.652503, 42393.219362, 40488.595068, 47836.837024, 56494.088877, 55479.618058, 44319.720622, 46527.066578, 32978.312476, 35788.329028, 45917.898858, 36020.758857, 45454.397791, 29312.444941, 32144.689001, 55344.931586, 35033.582023, 47186.021876]
    bincs+=[88.47, 87.77, 87.68, 64.10, -67.50, -58.55, -54.89, -38.25, 12.82, -26.06, -52.11, 41.07, 56.97, -71.94, -49.63, -59.55, -68.62, -15.46, 57.89, -67.61, -58.65, -34.17, 35.92, -2.12, -61.16, -44.70, 86.84, -75.40, -52.33, -68.45, -17.15, -68.07, -54.39, 23.95, 7.53, 52.51, -77.37, -44.10, 46.04, -31.60, -59.04, 84.84, 31.89, 65.46, -63.13, -52.57, -79.14, 42.66, 78.68, -57.98, 76.67, -59.82, -23.22, -37.45, 42.83, -76.49, -47.37, -15.26, 75.67, -70.97, 26.82, -17.66, 17.44, -70.25, -40.08, 87.42, -8.70, -65.23, -54.96, -57.67, 73.74, -81.52, -60.20, -58.58, 46.77, -68.50, 87.57, 43.10, 61.79, 44.27, -55.17, -64.54, 89.07, 46.03, 45.74, -67.36, -71.33, -72.79, 56.44, 60.89, 31.41, -17.84, -66.22, 43.76, -65.71, -56.68, 15.85, 87.37, 13.77, -47.55]
    bdecs+=[-99.77, -29.91, 54.89, 0.50, -5.49, -15.75, 27.96, 15.74, -0.19, 1.29, 20.28, 0.51, 12.93, 57.21, 14.91, -33.14, 9.28, 10.76, 11.10, 28.13, -22.51, 14.01, 1.17, 9.71, -6.32, -1.43, -30.61, -81.74, -14.98, -59.77, -3.10, 29.79, -11.36, 1.19, -1.53, 11.96, -137.58, 18.12, 2.24, 10.24, -17.22, 13.73, 6.46, 0.63, -48.27, 24.31, 57.87, 8.48, -16.48, -13.52, 19.39, -16.19, 7.79, 15.64, 1.85, -102.64, 17.93, -17.49, 0.37, -65.44, 5.10, -6.47, -0.51, 41.09, -4.65, -18.59, -14.34, 29.87, -11.68, -23.24, 15.43, 41.57, 29.45, -13.27, 1.57, -67.87, 67.64, 0.15, -9.55, 4.56, 8.58, -55.06, -73.48, 4.11, -1.89, -64.28, 32.11, 38.65, -4.13, -5.65, 9.04, 7.09, -56.34, 10.23, -63.48, -14.63, 9.24, 36.52, 0.89, 4.45]

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
