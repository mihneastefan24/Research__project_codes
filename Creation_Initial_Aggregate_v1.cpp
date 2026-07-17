// Main code for creating spherical particle aggregate

// =======================================




// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Alessandro Tasora
// =============================================================================
//
// Demo code about
// - using the ChParticleEmitter to create a cluster of random shapes
// - applying custom force field to particles
// - using Irrlicht to display objects
//
// =============================================================================

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/particlefactory/ChParticleEmitter.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/physics/ChContactContainer.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono/utils/ChUtilsCreators.h"
#include <fstream> // Include for file output
#include <random>
#include <string>


using namespace chrono;
using namespace chrono::particlefactory;
using namespace chrono::irrlicht;
using namespace chrono::utils;

//     A callback executed at each particle creation can be attached to the emitter.
//     For example, we need that new particles will be bound to Irrlicht visualization:
class MyCreatorForAll : public ChRandomShapeCreator::AddBodyCallback {
public:
    virtual void OnAddBody(std::shared_ptr<ChBody> mbody,
        ChCoordsys<> mcoords,
        ChRandomShapeCreator& mcreator) override {
        // Bind visual model to the visual system
        mbody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
        vis->BindItem(mbody);

        // Bind the collision model to the collision system
        if (mbody->GetCollisionModel())
            coll->Add(mbody->GetCollisionModel());

        // Disable gyroscopic forces for increased integrator stability
        mbody->SetUseGyroTorque(false);
    }
    ChVisualSystem* vis;
    ChCollisionSystem* coll;
};
class ContactBodyReporter : public chrono::ChContactContainer::ReportContactCallback {
public:
    virtual bool OnReportContact(
        const ChVector3<>& pA,            // contact point on object A
        const ChVector3<>& pB,            // contact point on object B
        const ChMatrix33<>& plane_coord, // contact plane coords (normal, U, V)
        const double& distance,          // penetration distance
        const double& eff_radius,        // effective radius of curvature
        const ChVector3<>& react_forces,  // forces in contact plane
        const ChVector3<>& react_torques, // torques in contact plane
        chrono::ChContactable* objA,     // contactable object A
        chrono::ChContactable* objB      // contactable object B
    ) override {
        // Attempt to cast contactable objects to ChBody
        auto bodyA = dynamic_cast<chrono::ChBody*>(objA);
        auto bodyB = dynamic_cast<chrono::ChBody*>(objB);

        if (bodyA && bodyB) {
            std::cout << "Contact between Body A and Body B:" << std::endl;
            std::cout << " - Body A ID: " << bodyA->GetIdentifier() << std::endl;
            std::cout << " - Body B ID: " << bodyB->GetIdentifier() << std::endl;
        }
        else {
            std::cout << "Contact involves non-body objects." << std::endl;
        }
        return true; // Continue reporting contacts
    }
};

int main(int argc, char* argv[]) {
    std::cout << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << std::endl;

    // Create a Chrono physical system
    ChSystemNSC sys;
    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    // Simulation parameters
    int velcase = 2;
    double initialSpin = 0.0000;
    double maxT = 10.0;
    double time_step = 0.5;
    double density = 2800;
    double restitution = 0.1;
    double adhesion = 1;
    double adhesionMult = 0.0;
    double Kn = 200000;
    double Kt = 200000;
    double Gn = 40;
    double Gt = 20;

    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Particle emitter spherical fixed particles");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddTypicalLights();
    vis->AddCamera(ChVector3d(0, 50, -35));


    auto particle_mat = chrono_types::make_shared<ChContactMaterialNSC>();
    particle_mat->SetFriction(0.6f);
    particle_mat->SetRestitution(restitution);
    // This part is if you change in MSC
//    particle_mat->SetAdhesion(adhesion);
//    particle_mat->SetGn(Gn);
//    particle_mat->SetGt(Gt);
//    particle_mat->SetKn(Kn);
//    particle_mat->SetKt(Kt);

//  This particle is left, if a core of much larger size is desired R_c >> R_particle
//  In this case it is used as any other particle
    auto sphere = chrono_types::make_shared<ChBodyEasySphere>(2,
        density,
        true,
        true,
        particle_mat);
    sphere->SetPos(ChVector3d(0, 0, 0));
    sphere->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
    sys.Add(sphere);
    // Create an emitter:
    ChParticleEmitter emitter;
    
    emitter.ParticlesPerSecond() = 20000;

    emitter.SetUseParticleReservoir(true);
    emitter.ParticleReservoirAmount() = 9999;


    // Our ChParticleEmitter object, among the main settings, it requires
    // that you give him four 'randomizer' objects: one is in charge of
    // generating random shapes, one is in charge of generating
    // random positions, one for random alignements, and one for random velocities.
    // In the following we need to instance such objects. (There are many ready-to-use
    // randomizer objects already available in chrono, but note that you could also
    // inherit your own class from these randomizers if the choice is not enough).

    // ---Initialize the randomizer for POSITIONS: random points in a large cube
    auto emitter_positions = chrono_types::make_shared<ChRandomParticlePositionOnGeometry>();
    //emitter_positions->SetGeometry(chrono_types::make_shared<ChSphere>(500000), ChFrame<>());
    emitter_positions->SetGeometry(chrono_types::make_shared<ChBox>(250, 250, 250), ChFrame<>());
    emitter.SetParticlePositioner(emitter_positions);

    // ---Initialize the randomizer for ALIGNMENTS
    auto emitter_rotations = chrono_types::make_shared<ChRandomParticleAlignmentUniform>();
    emitter.SetParticleAligner(emitter_rotations);

    // ---Initialize the randomizer for VELOCITIES, with statistical distribution
    auto mvelo = chrono_types::make_shared<ChRandomParticleVelocityAnyDirection>();
    mvelo->SetModulusDistribution(chrono_types::make_shared<ChUniformDistribution>(0.0, 0.0000));
    emitter.SetParticleVelocity(mvelo);

    // ---Initialize the randomizer for ANGULAR VELOCITIES, with statistical distribution
    auto mangvelo = chrono_types::make_shared<ChRandomParticleVelocityAnyDirection>();
    mangvelo->SetModulusDistribution(chrono_types::make_shared<ChUniformDistribution>(initialSpin, initialSpin));
    emitter.SetParticleAngularVelocity(mangvelo);

    //Idea: To speed up the process, create first a bulk of particles that reach equilibrium, over this add many other particles that are not in 
    //equilibrium and make the simulation like this. As the core is already in equilibrium , it will act as one and attract the others around it,
    //until all are in equilibrium. Check if it does make sens to have all the particles in contact in the initial part, or it should be some in contact
    //most of them not, as it is g = 0
        // ---Initialize the randomizer for CREATED SHAPES, with statistical distribution


        // Create a ChRandomShapeCreator object (ex. here for sphere particles)
    auto mcreator_spheres = chrono_types::make_shared<ChRandomShapeCreatorSpheres>();
    mcreator_spheres->SetDiameterDistribution(chrono_types::make_shared<ChConstantDistribution>(2));
    mcreator_spheres->SetDensityDistribution(chrono_types::make_shared<ChConstantDistribution>(density));
    mcreator_spheres->SetAddCollisionShape(true);
    emitter.SetParticleCreator(mcreator_spheres);

    // --- Optional: what to do by default on ALL newly created particles?
    //     A callback executed at each particle creation can be attached to the emitter.
    //     For example, we need that new particles will be bound to Irrlicht visualization:

    // a- define a class that implement your custom OnAddBody method (see top of source file)
    // b- create the callback object...
    auto mcreation_callback = chrono_types::make_shared<MyCreatorForAll>();
    // c- set callback own data that he might need...
    mcreation_callback->vis = vis.get();
    mcreation_callback->coll = sys.GetCollisionSystem().get();
    // d- attach the callback to the emitter!
    emitter.RegisterAddBodyCallback(mcreation_callback);
    
    // Bind all existing visual shapes to the visualization system
    vis->AttachSystem(&sys);

    // Modify some setting of the physical system for the simulation, if you want
    sys.SetSolverType(ChSolver::Type::PSOR);
    sys.GetSolver()->AsIterative()->SetMaxIterations(50);

    // Turn off default -9.8 downward gravity
    sys.SetGravitationalAcceleration(ChVector3d(0, 0, 0));

    // Write the documents
    std::ofstream pos("Position_Fixed.txt");
    std::ofstream vel("Velocity_Fixed.txt");
    std::ofstream mass("Mass_Fixed.txt");

    // Simulation loop
    double timestep = 0.5; // Set a timestep between 0.05 and 0.5, depending on the initial sparsity of the map
    int out_time = 0;
    double out_step = 1.0 / 20;
    while (vis->Run()) {

        // Create particle flow
        emitter.EmitParticles(sys, timestep);

        // Apply custom forcefield (brute force approach..)
        // A) reset 'user forces accumulators':
        unsigned int iter = 0;
        for (auto body : sys.GetBodies()) {
            body->EmptyAccumulators();
        }

        // B) store user computed force:
        // double G_constant = 6.674e-11; // gravitational constant
        double G_constant = 6.674e-9;  // gravitational constant - HACK to speed up simulation
        for (unsigned int i = 0; i < sys.GetBodies().size(); i++) {
            auto abodyA = sys.GetBodies()[i];
            for (unsigned int j = i + 1; j < sys.GetBodies().size(); j++) {
                auto abodyB = sys.GetBodies()[j];
                ChVector3d D_attract = abodyB->GetPos() - abodyA->GetPos();
                double r_attract = D_attract.Length();
                double f_attract = G_constant * (abodyA->GetMass() * abodyB->GetMass()) / (std::pow(r_attract, 2));
                ChVector3d F_attract = (D_attract / r_attract) * f_attract;

                abodyA->AccumulateForce(F_attract, abodyA->GetPos(), false);
                abodyB->AccumulateForce(-F_attract, abodyB->GetPos(), false);
            }
        }

        // Perform the integration timestep
        sys.DoStepDynamics(timestep);
        std::cout << sys.GetChTime() << std::endl;

        double time = sys.GetChTime();
        if (time >= out_time) {
            vis->BeginScene();
            vis->Render();
            vis->EndScene();
            out_time += out_step;
        }


    }
    // Take the values of the particles at the end when the body is aggregated
    unsigned int iterafter = 0;
    for (auto body : sys.GetBodies()) {
        body->EmptyAccumulators();
        pos << iterafter << "\t" << body->GetPos() << std::endl;
        vel << iterafter << "\t" << body->GetPosDt() << std::endl;
        mass << iterafter << "\t" << body->GetMass() << std::endl;
        iterafter++;
    }
    mass.close();
    pos.close();
    vel.close();
    return 0;
}