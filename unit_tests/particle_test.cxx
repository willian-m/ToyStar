#include "particle.h"
#define BOOST_TEST_MODULE ParticleTest
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructor_particle_test ){
    
    Vec3<double> r(1., 2., 3.);
    Vec3<double> v(.1, .2, .3);
    Vec3<double> a(10., 12., 13.);
    Particle particle1_test(1., 2., 3., .1, .2, .3, 10., 12., 13., 3.14);
    Particle particle2_test(r, v, a, 6.28);

    BOOST_CHECK( particle1_test.get_x() == 1. );
    BOOST_CHECK( particle1_test.get_y() == 2. );
    BOOST_CHECK( particle1_test.get_z() == 3. );

    BOOST_CHECK( particle1_test.get_vx() == .1 );
    BOOST_CHECK( particle1_test.get_vy() == .2 );
    BOOST_CHECK( particle1_test.get_vz() == .3 );

    BOOST_CHECK( particle1_test.get_ax() == 10. );
    BOOST_CHECK( particle1_test.get_ay() == 12. );
    BOOST_CHECK( particle1_test.get_az() == 13. );

    BOOST_CHECK( particle2_test.get_position().x == 1. );
    BOOST_CHECK( particle2_test.get_position().y == 2. );
    BOOST_CHECK( particle2_test.get_position().z == 3. );

    BOOST_CHECK( particle2_test.get_velocity().x == .1 );
    BOOST_CHECK( particle2_test.get_velocity().y == .2 );
    BOOST_CHECK( particle2_test.get_velocity().z == .3 );

    BOOST_CHECK( particle2_test.get_acceleration().x == 10. );
    BOOST_CHECK( particle2_test.get_acceleration().y == 12. );
    BOOST_CHECK( particle2_test.get_acceleration().z == 13. );

    BOOST_CHECK( particle2_test.get_mass() == 6.28 );
   
}

BOOST_AUTO_TEST_CASE( update_particle_test ){
    
    Vec3<double> r(1., 2., 3.);
    Vec3<double> v(.1, .2, .3);
    Vec3<double> a(10., 12., 13.);
    Particle particle1_test(1., 2., 3., .1, .2, .3, 10., 12., 13., 3.14);
    Particle particle2_test(r, v, a, 6.28);
    Particle particle3_test(1., 2., 3., .1, .2, .3, 10., 12., 13., 3.14);
   
    particle1_test.set_x(1.1);
    particle1_test.set_y(2.1);
    particle1_test.set_z(3.1);

    particle1_test.set_vx(11.1);
    particle1_test.set_vy(12.1);
    particle1_test.set_vz(13.1);

    particle1_test.set_ax(21.1);
    particle1_test.set_ay(22.1);
    particle1_test.set_az(23.1);

    BOOST_CHECK( particle1_test.get_x() == 1.1 );
    BOOST_CHECK( particle1_test.get_y() == 2.1 );
    BOOST_CHECK( particle1_test.get_z() == 3.1 );
    
    BOOST_CHECK( particle1_test.get_vx() == 11.1 );
    BOOST_CHECK( particle1_test.get_vy() == 12.1 );
    BOOST_CHECK( particle1_test.get_vz() == 13.1 );

    BOOST_CHECK( particle1_test.get_ax() == 21.1 );
    BOOST_CHECK( particle1_test.get_ay() == 22.1 );
    BOOST_CHECK( particle1_test.get_az() == 23.1 );

    r.x = 1.1;
    r.y = 2.1;
    r.z = 3.1;

    v.x = 11.1;
    v.y = 12.1;
    v.z = 13.1;

    a.x = 21.1;
    a.y = 22.1;
    a.z = 23.1;

    particle2_test.set_position(r);
    particle2_test.set_velocity(v);
    particle2_test.set_acceleration(a);

    BOOST_CHECK( particle2_test.get_x() == 1.1 );
    BOOST_CHECK( particle2_test.get_y() == 2.1 );
    BOOST_CHECK( particle2_test.get_z() == 3.1 );

    BOOST_CHECK( particle2_test.get_vx() == 11.1 );
    BOOST_CHECK( particle2_test.get_vy() == 12.1 );
    BOOST_CHECK( particle2_test.get_vz() == 13.1 );

    BOOST_CHECK( particle2_test.get_ax() == 21.1 );
    BOOST_CHECK( particle2_test.get_ay() == 22.1 );
    BOOST_CHECK( particle2_test.get_az() == 23.1 );

    particle3_test.set_position(1.1, 2.1, 3.1);
    particle3_test.set_velocity(11.1, 12.1, 13.1);
    particle3_test.set_acceleration(21.1, 22.1, 23.1);

    BOOST_CHECK( particle3_test.get_x() == 1.1 );
    BOOST_CHECK( particle3_test.get_y() == 2.1 );
    BOOST_CHECK( particle3_test.get_z() == 3.1 );

    BOOST_CHECK( particle3_test.get_vx() == 11.1 );
    BOOST_CHECK( particle3_test.get_vy() == 12.1 );
    BOOST_CHECK( particle3_test.get_vz() == 13.1 );

    BOOST_CHECK( particle3_test.get_ax() == 21.1 );
    BOOST_CHECK( particle3_test.get_ay() == 22.1 );
    BOOST_CHECK( particle3_test.get_az() == 23.1 );

    
}