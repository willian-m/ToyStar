#include "vec3.h"
#define BOOST_TEST_MODULE Vec3Test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructor_vec3_test ){
    Vec3<double> test_vec3(1.,2.,3.);

    BOOST_CHECK( test_vec3.x == 1. );
    BOOST_CHECK( test_vec3.y == 2. );
    BOOST_CHECK( test_vec3.z == 3. );
}

BOOST_AUTO_TEST_CASE( update_vec3_test ){
    Vec3<double> test_vec3(1.,2.,3.);
    
    test_vec3.x = 10.;
    test_vec3.y = 20.;
    test_vec3.z = 30.;

    BOOST_CHECK( test_vec3.x == 10. );
    BOOST_CHECK( test_vec3.y == 20. );
    BOOST_CHECK( test_vec3.z == 30. );
}

BOOST_AUTO_TEST_CASE( multiply_vec3_by_scalar ){
    Vec3<double> test_vec3(1.,2.,3.);
    
    Vec3<double> res = test_vec3*10.;

    BOOST_CHECK( res.x == 10. );
    BOOST_CHECK( res.y == 20. );
    BOOST_CHECK( res.z == 30. );
}

BOOST_AUTO_TEST_CASE( sum_vec3 ){
    Vec3<double> test_vec3_A(1.,2.,3.);
    Vec3<double> test_vec3_B(3.,4.,5.);
    
    Vec3<double> res = test_vec3_A+test_vec3_B;

    BOOST_CHECK( res.x == 4. );
    BOOST_CHECK( res.y == 6. );
    BOOST_CHECK( res.z == 8. );
}