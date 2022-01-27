//
// Created by ethan on 12/15/2021.
//

#ifndef PSEMOLDYN_GROUPB_LIBRARY_H
#define PSEMOLDYN_GROUPB_LIBRARY_H

#include <string>
#include <vector>

namespace library
{
    //
    //
    typedef double time;

    //
    //
    typedef double epsilon;

    //
    //
    typedef double sigma;

    //
    //
    typedef std::string algorithm;

    //
    //
    typedef std::string containerAlgorithm;

    //
    //
    typedef int boundary;

    //
    //
    typedef std::string benchmark;

    //
    //
    struct dimension
    {
    public:
        const int&
        x () const
        {
            return x_;
        }
        void x (const int& x)
        {
            x_ = x;
        }

        const int&
        y () const
        {
            return y_;
        }
        void y (const int& y)
        {
            y_ = y;
        }

        const int&
        z () const
        {
            return z_;
        }
        void z (const int& z)
        {
            z_ = z;
        }
    private:
        int x_;
        int y_;
        int z_;
    };

    //
    //
    struct point
    {
    public:
        const double&
        x () const
        {
            return x_;
        }
        void x (const double& x)
        {
            x_ = x;
        }

        const double&
        y () const
        {
            return y_;
        }
        void y (const double& y)
        {
            y_ = y;
        }

        const double&
        z () const
        {
            return z_;
        }
        void z (const double& z)
        {
            z_ = z;
        }
    private:
        double x_;
        double y_;
        double z_;
    };

    //
    //
    struct velocity
    {
    public:
        const double&
        x () const
        {
            return x_;
        }
        void x (const double& x)
        {
            x_ = x;
        }

        const double&
        y () const
        {
            return y_;
        }
        void y (const double& y)
        {
            y_ = y;
        }

        const double&
        z () const
        {
            return z_;
        }
        void z (const double& z)
        {
            z_ = z;
        }
    private:
        double x_;
        double y_;
        double z_;
    };

    //
    //
    struct boundaryConditions
    {
    public:
        const library::boundary
        x () const
        {
            return xBoundary_;
        }
        void x (const library::boundary& x)
        {
            xBoundary_ = x;
        }

        const library::boundary
        y () const
        {
            return yBoundary_;
        }
        void y (const library::boundary& y)
        {
            yBoundary_ = y;
        }

        const library::boundary
        z () const
        {
            return zBoundary_;
        }
        void z (const library::boundary& z)
        {
            zBoundary_ = z;
        }
    private:
        library::boundary xBoundary_;
        library::boundary yBoundary_;
        library::boundary zBoundary_;
    };

    //
    //
    struct simulationContainer
    {
    public:
        library::boundaryConditions
        boundaryConditions () const
        {
            return boundaryConditions_;
        }
        void
        boundaryConditions(const library::boundaryConditions& boundaryConditions)
        {
            boundaryConditions_=boundaryConditions;
        }

        library::dimension
        dimension () const
        {
            return dimension_;
        }
        void
        dimension(const library::dimension& dimension)
        {
            dimension_=dimension;
        }

        library::containerAlgorithm
        containerAlgorithm()const
        {
            return cAlgorithm_;
        }
        void
        containerAlgorithm(const library::containerAlgorithm& c)
        {
            cAlgorithm_=c;
        }


        const double&
        mesh () const
        {
            return mesh_;
        }
        void
        mesh(const double& mesh)
        {
            mesh_ = mesh;
        }

        const double&
        cutOff () const
        {
            return cutOff_;
        }
        void
        cutOff(const double& cutOff)
        {
            cutOff_=cutOff;
        }
    private:
        library::boundaryConditions boundaryConditions_;
        library::dimension dimension_;
        library::containerAlgorithm cAlgorithm_;
        double mesh_;
        double cutOff_;
    };

    //
    //
    struct Cube
    {
    public:
        library::dimension
        dimension () const
        {
            return dimension_;
        }
        void
        dimension(const library::dimension& dimension)
        {
            dimension_ = dimension;
        }

        library::point
        startPoint () const
        {
            return startPoint_;
        }
        void
        startPoint (const library::point& startPoint)
        {
            startPoint_=startPoint;
        }

        const double&
        h () const
        {
            return h_;
        }
        void
        h(const double& h)
        {
            h_=h;
        }

        const double&
        m () const
        {
            return m_;
        }
        void
        m(const double& m)
        {
            m_=m;
        }

        library::velocity
        velocity () const
        {
            return velocity_;
        }
        void
        velocity(const library::velocity& velocity)
        {
            velocity_=velocity;
        }

        library::epsilon
        epsilon () const
        {
            return epsilon_;
        }
        void
        epsilon(const library::epsilon& epsilon)
        {
            epsilon_=epsilon;
        }

        library::sigma
        sigma () const
        {
            return sigma_;
        }
        void
        sigma(const library::sigma& sigma)
        {
            sigma_=sigma;
        }
    private:
        library::dimension dimension_;
        library::point startPoint_;
        double h_;
        double m_;
        library::velocity velocity_;
        library::epsilon epsilon_;
        library::sigma sigma_;
    };

    //
    //
    struct Sphere
    {
    public:
        const double&
        radius () const
        {
            return radius_;
        }
        void
        radius(const double& radius)
        {
            radius_=radius;
        }

        library::point
        center () const
        {
            return center_;
        }
        void
        center (const library::point& center)
        {
            center_=center;
        }

        const double&
        h () const
        {
            return h_;
        }
        void
        h(const double& h)
        {
            h_=h;
        }

        const double&
        m () const
        {
            return m_;
        }
        void
        m(const double& m)
        {
            m_=m;
        }

        library::velocity
        velocity () const
        {
            return velocity_;
        }
        void
        velocity(const library::velocity& velocity)
        {
            velocity_=velocity;
        }

        library::epsilon
        epsilon () const
        {
            return epsilon_;
        }
        void
        epsilon(const library::epsilon& epsilon)
        {
            epsilon_=epsilon;
        }

        library::sigma
        sigma () const
        {
            return sigma_;
        }
        void
        sigma(const library::sigma& sigma)
        {
            sigma_=sigma;
        }
    private:
        double radius_;
        library::point center_;
        double h_;
        double m_;
        library::velocity velocity_;
        library::epsilon epsilon_;
        library::sigma sigma_;
    };

    //
    //
    struct particles
    {
    public:
        //
        //
        typedef std::vector<library::Cube> cubes;

        //
        //
        typedef std::vector<library::Sphere> spheres;

        const cubes&
        cube () const {
            return cubes_;
        }
        cubes&
        cube ()
        {
            return cubes_;
        }

        const spheres&
        sphere () const
        {
            return spheres_;
        }

        spheres&
        sphere ()
        {
            return spheres_;
        }

    private:
        cubes cubes_;
        spheres spheres_;
    };

    //
    //
    struct thermostats
    {
    public:
        const double&
        initialTemperature () const
        {
            return initialTempterature_;
        }
        void
        initialTemperature(const double& initialTemperature)
        {
            initialTempterature_=initialTemperature;
        }

        const double&
        targetTemperature () const
        {
            return targetTemperature_;
        }
        void targetTemperature(const double& targetTemperature)
        {
            targetTemperature_=targetTemperature;
        }

        const double&
        maxDelta()const
        {
            return maxDelta_;
        }
        void
        maxDelta(const double& maxDelta)
        {
            maxDelta_=maxDelta;
        }

        const int&
        stepSize () const
        {
            return stepSize_;
        }
        void
        stepSize(const int& stepSize)
        {
            stepSize_=stepSize;
        }

    private:
        double initialTempterature_;
        double targetTemperature_;
        double maxDelta_;
        int stepSize_;
    };

    //
    //
    typedef int parallelizationStrategy;



    struct molsim
    {
    public:
        const std::string&
        input()const{
            return input_;
        }
        void
        input(const std::string& input){
            input_=input;
        }

        library::time
        delta_t()const
        {
            return deltaT;
        }
        void
        delta_t(library::time t)
        {
            deltaT=t;
        }

        library::time
        endtime () const
        {
            return endTime;
        }
        void
        endtime(library::time t)
        {
            endTime=t;
        }

        const int&
        outputStep () const
        {
            return outputStep_;
        }
        void
        outputStep(const int& step)
        {
            outputStep_=step;
        }

        library::epsilon
        epsilon () const
        {
            return epsilon_;
        }
        void
        epsilon(const library::epsilon& epsilon)
        {
            epsilon_=epsilon;
        }

        library::sigma
        sigma () const
        {
            return sigma_;
        }
        void
        sigma(const library::sigma& sigma)
        {
            sigma_=sigma;
        }

        const double&
        gravity()const
        {
            return gravity_;
        }
        void
        gravity(const double& g)
        {
            gravity_=g;
        }

        const double&
        averageV()const
        {
            return averageV_;
        }
        void
        averageV(const double& v)
        {
            averageV_ = v;
        }

        std::string
        algorithm () const
        {
            return algorithm_;
        }
        void
        algorithm(const std::string& algo)
        {
            algorithm_=algo;
        }

        library::simulationContainer
        simu()const{
            return simulationContainer_;
        }
        void
        simu(const library::simulationContainer& c)
        {
            simulationContainer_=c;
        }

        library::particles
        particles()const{
            return particles_;
        }
        void
        particles(const library::particles& p)
        {
            particles_=p;
        }

        library::thermostats
        thermostats()const
        {
            return thermostats_;
        }
        void
        thermostats(const library::thermostats& t){
            thermostats_=t;
        }

        library::parallelizationStrategy
        parStrat()const
        {
            return parallelizationStrategy_;
        }
        void
        parStrat(const library::parallelizationStrategy& p){
            parallelizationStrategy_=p;
        }

        library::benchmark
        benchmark()const
        {
            return benchmark_;
        }
        void
        benchmark(const library::benchmark& b){
            benchmark_=b;
        }
    private:
        std::string input_;
        library::time deltaT;
        library::time endTime;
        int outputStep_;
        library::epsilon epsilon_;
        library::sigma sigma_;
        double gravity_;
        double averageV_;
        std::string algorithm_;
        library::simulationContainer simulationContainer_;
        library::particles particles_;
        library::thermostats thermostats_;
        library::parallelizationStrategy parallelizationStrategy_;
        library::benchmark benchmark_;
    };

}
#endif //PSEMOLDYN_GROUPB_LIBRARY_H
