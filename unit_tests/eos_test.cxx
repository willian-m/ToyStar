#include "eos_polytropic.h"
#define BOOST_TEST_MODULE Vec3Test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( get_pressure_test ){

    EOSPolytropic eos(3.,2.);

    BOOST_CHECK( abs(eos.get_pressure(4.) - 3.*pow(4,1.5) ) < 1E-15 );

}