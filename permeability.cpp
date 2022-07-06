#include "palabos3D.h"
#include "palabos3D.hh"
#include "ctime"

using namespace plb;

#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

typedef double T;

// Smagorinsky constant for LES model.
// const T cSmago = 0.14;

const T rhoEmpty = T(1);
const T rho = T(1325);

plint writeImagesIter = 500;
plint getStatisticsIter = 20;

plint maxIter;
plint N;
plint nx, ny, nz;
T delta_t, delta_x;
Array<T, 3> externalForce;
T nuPhys, nuLB, tau, omega, surfaceTensionPhys, surfaceTensionLB, contactAngle;

std::string fNameIn, fNameOut;

void readGeometry(std::string fNameIn, std::string fNameOut, MultiScalarField3D<int>& geometry)
{
    const plint nx = geometry.getNx();
    const plint ny = geometry.getNy();
    const plint nz = geometry.getNz();

    Box3D sliceBox(0, 0, 0, ny - 1, 0, nz - 1);
    std::unique_ptr<MultiScalarField3D<int> > slice = generateMultiScalarField<int>(geometry, sliceBox);
    plb_ifstream geometryFile(fNameIn.c_str());
    for (plint iX = 0; iX < nx; ++iX) {
        if (!geometryFile.is_open()) {
            pcout << "Error: could not open geometry file " << fNameIn << std::endl;
            exit(EXIT_FAILURE);
        }
        geometryFile >> *slice;
        copy(*slice, slice->getBoundingBox(), geometry, Box3D(iX, iX, 0, ny - 1, 0, nz - 1));
    }
}

void setupParameters() {
    delta_x = 0.001 / N;
    // Gravity in lattice units.
    T gLB = 9.8 * delta_t * delta_t / delta_x;
    externalForce = Array<T, 3>(0., 0., -gLB); // Gravity in negative penetration-direction
    tau = (nuPhys * DESCRIPTOR<T>::invCs2 * delta_t) / (delta_x * delta_x) + 0.5;
    omega = 1. / tau;
    nuLB = (tau - 0.5) * DESCRIPTOR<T>::cs2; // Viscosity in lattice units.
    surfaceTensionLB = (rhoEmpty / rho) * (delta_t * delta_t) / (delta_x * delta_x * delta_x) * surfaceTensionPhys;
}
// inlet definition
bool insideFluid(T x, T y, T z)
{
    Array<T, 3> pos(x, y, z);
    Array<T, 3> center(nx / 2.0, ny / 2.0, 0); 
    T r = norm(pos - center);
    if (r <= 1) {  // r <= nz / 25
        return true;
    }
    return false;
}
//
int initialFluidFlags(plint iX, plint iY, plint iZ) {
    if (insideFluid(iX, iY, iZ)) {
        return freeSurfaceFlag::fluid;
    }
    else {
        return freeSurfaceFlag::empty;
    }
}
void porousMediaSetup(FreeSurfaceFields3D<T, DESCRIPTOR>& fields, FreeSurfaceFields3D<T, DESCRIPTOR>& fields2,
    MultiScalarField3D<int>& geometry)
{
    setToConstant(fields.flag, fields.flag.getBoundingBox(), (int)freeSurfaceFlag::wall);
    setToFunction(fields.flag, fields.flag.getBoundingBox().enlarge(-1), initialFluidFlags); // getBoundingBox().enlarge(-1)
    setToConstant(fields.flag, geometry, 1,
        fields.flag.getBoundingBox().enlarge(-1), (int)freeSurfaceFlag::wall);
    setToConstant(fields.flag, Box3D(0, nx - 1, 0, ny - 1, 0, 0), (int)freeSurfaceFlag::fluid);
    // No periodicity (top != bottom ; left != right)
    //pcout << "5" << std::endl;
    fields.periodicityToggle(2, false);
    fields.periodicityToggle(1, false);
    fields.periodicityToggle(0, false);
    fields.defaultInitialize();
}

void writeResults(MultiBlockLattice3D<T, DESCRIPTOR>& lattice, MultiScalarField3D<T>& volumeFraction, plint iT)
{
    static const plint nx = lattice.getNx();
    static const plint ny = lattice.getNy();
    static const plint nz = lattice.getNz();

    const plint imSize = 600;
    Box3D slice(0, nx - 1, ny / 2, ny / 2, 0, nz - 1);
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iT, 6),
        *computeVelocityNorm(lattice, slice),
        imSize, imSize);
    imageWriter.writeScaledGif(createFileName("rho", iT, 6),
        *computeDensity(lattice, slice),
        imSize, imSize);
    //VtkImageOutput3D<T> vtkOut(createFileName("rho", iT, 6), 1.);
    //vtkOut.writeData<float>(*computeDensity(lattice, slice), "density", 1.);
    

    imageWriter.writeScaledGif(createFileName("volumeFraction", iT, 6), *extractSubDomain(volumeFraction, slice),
        imSize, imSize);
    // Set vtk-image outputs to compare with experiment datas
    //VtkImageOutput3D<T> vtkOut(createFileName("u", iT, 6), 1.);
    //vtkOut.writeData<float>(*computeVelocityNorm(lattice, slice), "velocityNorm", 1.);
    //vtkOut.writeData<3, float>(*computeVelocity(lattice, slice), "velocity", 1.);

    VtkImageOutput3D<T> vtkOut(createFileName("volumeFraction", iT, 6), 1.);
    vtkOut.writeData<float>(volumeFraction, "vf", 1.);
}

void writeStatistics(FreeSurfaceFields3D<T, DESCRIPTOR>& fields) {
    pcout << " -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- " << std::endl;
    T averageMass = freeSurfaceAverageMass<T, DESCRIPTOR>(fields.freeSurfaceArgs, fields.lattice.getBoundingBox());
    pcout << "Average Mass: " << averageMass << std::endl;

    T averageDensity = freeSurfaceAverageDensity<T, DESCRIPTOR>(fields.freeSurfaceArgs, fields.lattice.getBoundingBox());
    pcout << "Average Density: " << std::setprecision(12) << averageDensity << std::endl;

    T averageVolumeFraction = freeSurfaceAverageVolumeFraction<T, DESCRIPTOR>(fields.freeSurfaceArgs, fields.lattice.getBoundingBox());
    pcout << "Average Volume-Fraction: " << std::setprecision(12) << averageVolumeFraction << std::endl;
    pcout << " -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- " << std::endl;
}

int main(int argc, char** argv)
{

// record the time of the program
    clock_t start , end;
    start =  clock();

    plbInit(&argc, &argv);

    if (global::argc() < 15) {
        pcout << "Error missing some input parameter\n";
    }

    try {
        global::argv(1).read(fNameOut);
        global::directories().setOutputDir(fNameOut + "/");

        global::argv(2).read(nuPhys);
        global::argv(3).read(surfaceTensionPhys);
        global::argv(4).read(contactAngle);
        global::argv(5).read(N);
        global::argv(6).read(delta_t);
        global::argv(7).read(maxIter);
        global::argv(8).read(nx);
        global::argv(9).read(ny);
        global::argv(10).read(nz);
        global::argv(11).read(fNameIn);
    }
    catch (PlbIOException& except) {
        pcout << except.what() << std::endl;
        pcout << "The parameters for this program are:\n";
        pcout << "1. Output directory name.\n";
        pcout << "2. Kinematic viscosity in physical Units (m^2/s) of layer 1.\n";
        pcout << "3. Surface tension in physical units (N/m) of layer 1.\n";
        pcout << "4. Contact angle (in degrees) of layer 1.\n";
        pcout << "8. Number of lattice nodes for lz (suggestion: 1000).\n";
        pcout << "9. delta_t (suggestion: 1.e-8).\n";
        pcout << "10. Maximum Iteration cycle.\n";
        pcout << "11. nx nz ny.\n";
        pcout << "12. Name of geometry input data.\n";
        pcout << "Reasonable parameters on a parallel machine are: " << (std::string)global::argv(0) << " tmp 1.01e-6 0.0728 150.0 1000 1.e-8 80000 501(x) 501(z) 501(y) Raw Data.dat\n";
        exit(EXIT_FAILURE);
    }

    setupParameters();

    pcout << "Reading the geometry file." << std::endl;
    MultiScalarField3D<int> geometry(nx, ny, nz);
    readGeometry(fNameIn, fNameOut, geometry);
    MultiScalarField3D<int> inputgeometry(nx, ny, nz / 2 + 1);
    //MultiScalarField3D<int> inputgeometry2(nx, ny, nz / 2 + 1);
    //copy(geometry, Box3D(0, nx / 2, 0, ny - 1, 0, nz / 2 + 1), inputgeometry, inputgeometry.getBoundingBox());
    copy(geometry, Box3D(0, nx - 1, 0, ny - 1, 0, nz / 2 + 1), inputgeometry, inputgeometry.getBoundingBox());
    //copy(geometry, Box3D(0, nx - 1, 0, ny / 2, 0, nz / 2 + 1), inputgeometry2, inputgeometry2.getBoundingBox());
    
    pcout << "delta_t = " << delta_t << std::endl;
    pcout << "delta_x = " << delta_x << std::endl;
    pcout << "delta_t*delta_t/delta_x = " << delta_t * delta_t / delta_x << std::endl;
    pcout << "externalForce = " << externalForce[2] << std::endl;
    pcout << "surfaceTensionLB L1 = " << surfaceTensionLB << std::endl;
    pcout << "relaxation time L1 = " << tau << std::endl;
    pcout << "omega L1 = " << omega << std::endl;
    pcout << "kinematic viscosity physical units L1 = " << nuPhys << std::endl;
    pcout << "kinematic viscosity lattice units L1 = " << nuLB << std::endl;

    global::timer("initialization").start();
    SparseBlockStructure3D blockStructure2(createRegularDistribution3D(2*nx, ny, nz / 2 + 1));
    SparseBlockStructure3D blockStructure(createRegularDistribution3D(nx, ny, nz / 2 + 1));

    /*Dynamics<T,DESCRIPTOR>* dynamics
        = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega, cSmago);*/
    Dynamics<T, DESCRIPTOR>* dynamics = new IncBGKdynamics<T, DESCRIPTOR>(omega);

    FreeSurfaceFields3D<T, DESCRIPTOR> fields(blockStructure, dynamics->clone(), rhoEmpty,
        surfaceTensionLB, contactAngle, externalForce);
    
    std::vector<MultiBlock3D*> freeSurfaceArgs = aggregateFreeSurfaceParams(fields.lattice, fields.rhoBar, fields.j, fields.mass, fields.volumeFraction,
        fields.flag, fields.normal, fields.helperLists, fields.curvature, fields.outsideDensity);

    Array<T, 3> injectionVelocity((T)0., (T)0., (T)0.1);

    pcout << "Time spent for setting up lattices: "
        << global::timer("initialization").stop() << std::endl;
    T lastIterationTime = T();

    T tmp_mass = T(); //def a temporary number
    tmp_mass = 0;

    for (plint iT = 0; iT <= maxIter; ++iT) {
        global::timer("iteration").restart();

        T sum_of_mass_matrix = T();
        T lost_mass = T();
        if (iT % getStatisticsIter == 0) {
            pcout << std::endl;
            pcout << "ITERATION = " << iT << std::endl;
            pcout << "Time of last iteration is " << lastIterationTime << " seconds" << std::endl;
            writeStatistics(fields);
            sum_of_mass_matrix = fields.lattice.getInternalStatistics().getSum(0);
            pcout << "Sum of mass matrix: " << sum_of_mass_matrix << std::endl;
            lost_mass = fields.lattice.getInternalStatistics().getSum(1);
            pcout << "Lost mass: " << lost_mass << std::endl;
            pcout << "Total mass: " << sum_of_mass_matrix + lost_mass << std::endl;
            pcout << "Interface cells: " << fields.lattice.getInternalStatistics().getIntSum(0) << std::endl;
        }



        if (iT % writeImagesIter == 0) {
            global::timer("images").start();
            writeResults(fields.lattice, fields.volumeFraction, iT);
            pcout << "Total time spent for writing images: "
                << global::timer("images").stop() << std::endl;
        }
        
        // This includes the collision-streaming cycle, plus all free-surface operations.
        fields.lattice.executeInternalProcessors();
        /*applyProcessingFunctional (
            new PouringLiquid3D<T,DESCRIPTOR>(dynamics->clone(), injectionVelocity),
                                   Box3D(4*nx/9, 5*nx/9, 4*ny/9, 5*ny/9, 0, 0), freeSurfaceArgs );*/
        applyProcessingFunctional(
            new PouringLiquid3D<T, DESCRIPTOR>(dynamics->clone(), injectionVelocity),
            //Box3D(9 * nx / 19, 10 * nx / 19, 9 * ny / 19, 10 * ny / 19, 0, 0), freeSurfaceArgs);
            Box3D(0, nx - 1, 0, ny - 1, 0, 0), freeSurfaceArgs);
        fields.lattice.evaluateStatistics();
        fields.lattice.incrementTime();

        lastIterationTime = global::timer("iteration").stop();

        if(tmp_mass == sum_of_mass_matrix + lost_mass){
            maxIter = iT;
            break;
        }else{
            tmp_mass = sum_of_mass_matrix + lost_mass;
        }


    }

    end = clock();
    std::cout << (double)( end - start )/CLOCKS_PER_SEC;
}
