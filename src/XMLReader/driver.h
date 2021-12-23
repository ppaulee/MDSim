//
// Created by ethan on 12/21/2021.
//

#ifndef TESTX_DRIVER_H
#define TESTX_DRIVER_H


#include "library.h"
#include "MolSimImpl.h"
#include <iostream>

library::molsim generateFromXML(std::string xmlInput){
    library::molsim m;
    try {
        // Instantiate individual parsers.
        //
        ::molsim_pimpl molsim_p;
        ::xml_schema::string_pimpl string_p;
        ::time_pimpl time_p;
        ::xml_schema::int_pimpl int_p;
        ::epsilon_pimpl epsilon_p;
        ::sigma_pimpl sigma_p;
        ::xml_schema::double_pimpl double_p;
        ::algorithm_pimpl algorithm_p;
        ::simulationContainer_pimpl simulationContainer_p;
        ::boundaryConditions_pimpl boundaryConditions_p;
        ::dimension_pimpl dimension_p;
        ::containerAlgorithm_pimpl containerAlgorithm_p;
        ::particles_pimpl particles_p;
        ::Cube_pimpl Cube_p;
        ::point_pimpl point_p;
        ::velocity_pimpl velocity_p;
        ::Sphere_pimpl Sphere_p;
        ::thermostats_pimpl thermostats_p;
        ::benchmark_pimpl benchmark_p;

        // Connect the parsers together.
        //
        molsim_p.parsers(string_p,
                         time_p,
                         time_p,
                         int_p,
                         epsilon_p,
                         sigma_p,
                         double_p,
                         double_p,
                         algorithm_p,
                         simulationContainer_p,
                         particles_p,
                         thermostats_p,
                         benchmark_p);

        simulationContainer_p.parsers(boundaryConditions_p,
                                      dimension_p,
                                      double_p,
                                      double_p,
                                      containerAlgorithm_p);

        dimension_p.parsers(int_p,
                            int_p,
                            int_p);

        particles_p.parsers(Cube_p,
                            Sphere_p);

        Cube_p.parsers(dimension_p,
                       point_p,
                       double_p,
                       double_p,
                       velocity_p,
                       epsilon_p,
                       sigma_p);

        point_p.parsers(double_p,
                        double_p,
                        double_p);

        velocity_p.parsers(double_p,
                           double_p,
                           double_p);

        Sphere_p.parsers(point_p,
                         double_p,
                         double_p,
                         double_p,
                         velocity_p,
                         epsilon_p,
                         sigma_p);

        thermostats_p.parsers(double_p,
                              double_p,
                              double_p,
                              int_p);

        // Parse the XML document.
        //
        ::xml_schema::document doc_p(molsim_p, "molsim");

        molsim_p.pre();
        doc_p.parse(xmlInput);
        m = molsim_p.post_molsim();
    }
    catch (const ::xml_schema::exception_ &e) {
        std::cerr << e << std::endl;
    }
    catch (const std::ios_base::failure &) {
        std::cerr << xmlInput << ": error: io failure" << std::endl;
    }
    return m;
}


#endif //TESTX_DRIVER_H
