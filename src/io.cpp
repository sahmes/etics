#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <fstream>
#include <sstream>
// #include <sstream>
#include "common.hpp"

#ifndef NOBOOST
    #include <boost/property_tree/ptree.hpp>
    #include <boost/property_tree/ini_parser.hpp>
#endif

using namespace std;

void ReadICs(string Filename, int N, int Skip, Particle **P_h) {

    *P_h = new Particle[N];

    string Line;
    ifstream InputFile(Filename.c_str());
    int LineNum;
    for (LineNum = 0; LineNum < Skip; LineNum++) getline(InputFile, Line);
    int i = 0;
    if (InputFile.is_open()) {
        while (getline(InputFile, Line))
        {
            Real m, x, y, z, vx, vy, vz;
            Particle p;
// #ifdef ETICS_DOUBLE_PRECISION
//             int Count = sscanf(Line.c_str(), "%*d %*f %lf %lf %lf %lf %lf %lf", &x, &y, &z, &vx, &vy, &vz);
// #else
//             int Count = sscanf(Line.c_str(), "%*d %*f %f %f %f %f %f %f", &x, &y, &z, &vx, &vy, &vz);
// #endif
//             if (Count != 6) {
#ifdef ETICS_DOUBLE_PRECISION
            int Count = sscanf(Line.c_str(), "%*d %lf %lf %lf", &x, &y, &z);
#else
            int Count = sscanf(Line.c_str(), "%*d %f %f %f", &x, &y, &z);
#endif
            if (Count != 3) {
                cerr << "Problem in input file, line " << LineNum + 1 << "." << endl;
                exit(1);
            }
//             p.m = m;
            p.m = 1.0/N;
            p.pos = vec3(x, y, z);
            p.vel = vec3(vx, vy, vz);
            p.ID = i;
            p.Status = 0;
            p.CalculateR2();
            (*P_h)[i] = p;
            i++;
            LineNum++;
            if (i == N) break;
        }
        InputFile.close();
    } else {
        cerr << "Can't open file" << Filename << "." << endl;
        exit(1);
    }
    if (i < N) {
        cerr << "Was only able to read " << i << " particles while " << N << " were requested." << endl;
        exit(1);
    }
    // Transfer everything to the device; the host memory is implicitly freed
    // when we leave the current scope (probably).
//     CmCorrection(&P_h[0]);
}

void Dooooo(string Prefix, int SnapNumber, Particle *P_h, int N, Real T);
void WriteSnapshot(string Prefix, int SnapNumber, Particle *P_h, int N, Real T) {
    Dooooo(Prefix, SnapNumber, P_h, N, T);
    return;
    char S[512];
    sprintf(S, "%s%04d", Prefix.c_str(), SnapNumber);
    ofstream SnapshotFile;
    SnapshotFile.open(S);
    sprintf(S, "%06d\n", SnapNumber); SnapshotFile << S;
    sprintf(S, "%06d\n", N); SnapshotFile << S;
    sprintf(S, "%.16E\n", T); SnapshotFile << S;
    for (int i = 0; i < N; i++) {
        sprintf(S, "%06d%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e\n", P_h[i].ID, P_h[i].m, P_h[i].pos.x, P_h[i].pos.y, P_h[i].pos.z, P_h[i].vel.x, P_h[i].vel.y, P_h[i].vel.z); SnapshotFile << S;
    }
    SnapshotFile.close();
}

#define THROW_EXCEPTION(str, stat) {cerr << str << endl; exit((stat));}
void ParseInput(int argc, char *argv[], ParametersStruct *Params) {
    ParametersStruct P;
#ifndef NOBOOST
    if (argc < 2) THROW_EXCEPTION("Usage: etics [FILE]...", 1)
    ifstream TestFileObject(argv[1]);
    if (!TestFileObject) THROW_EXCEPTION("Problem reading file...", 1)
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(argv[1], pt); // what if file doesn't exist?

    P.N = pt.get<int>("N", 0);
    if (P.N <= 0) THROW_EXCEPTION("Could not read number of particles (N) from ini file.", 1)

    P.Tcrit = pt.get<Real>("Tcrit", -1);
    if (P.Tcrit < 0) THROW_EXCEPTION("Could not read finish time (Tcrit) from ini file.", 1)

    P.dT1 = pt.get<Real>("dT1", 0);
    if (P.dT1 <= 0) THROW_EXCEPTION("Could output time interval (dT1) from ini file.", 1)

    P.dT2 = pt.get<Real>("dT2", 0);
    if (P.dT2 <= 0) THROW_EXCEPTION("Could snapshot time interval (dT2) from ini file.", 1)

    P.ConstantStep = pt.get<Real>("StepSize", -1);
    if (P.ConstantStep < 0) THROW_EXCEPTION("Could not read step size (StepSize) from ini file.", 1)

    P.Filename = pt.get<string>("Filename", "\n");
    if (P.Filename == "\n") THROW_EXCEPTION("Could not read initial condition file name (Filename) from ini file.", 1)

    P.Skip = pt.get<int>("Skip", 0);
    P.Prefix = pt.get<string>("Prefix", "");
    P.DeviceID = pt.get<int>("device", -1);
    if (argc >= 3) { // If there is a argument after the file, it must be either --device or -d
        int Version1 = strncmp(argv[2], "--device", 8);
        int Version2 = strncmp(argv[2], "-d", 2);
        if ((Version1!=0) && (Version2!=0)) THROW_EXCEPTION("Commandline argument after input file must be device number.", 1)
        char *DeviceIDStr = strchr(argv[2], '=');
        if (DeviceIDStr!=NULL) DeviceIDStr = DeviceIDStr + 1;
        else {
            if (argc <= 3) THROW_EXCEPTION("Device number must follow.", 1)
            DeviceIDStr = argv[3];
        }
        P.DeviceID = atoi(DeviceIDStr);
        if ((P.DeviceID==0) && (strcmp(DeviceIDStr, "0")!=0)) THROW_EXCEPTION("Error understanding device number.", 1)
    }
#else
    // If no BOOST, we include a file with the parameters and compile it.
    int N, Skip = 0, DeviceID = -1;
    Real dT1, dT2, Tcrit, StepSize;
    std::string Filename, Prefix = "";
    #include "noboost.inc"
    P.N = N;
    P.Filename = Filename;
    P.Skip = Skip;
    P.dT1 = dT1;
    P.dT2 = dT2;
    P.Tcrit = Tcrit;
    P.ConstantStep = StepSize;
    P.DeviceID = DeviceID;
#endif
    P.Seed = INT_MIN;
    if ((P.Filename[0]=='_') && (P.Filename.rfind("_") > 0)) {
        int Break = P.Filename.rfind("_") + 1;
        string SeedStr = P.Filename.substr(Break, P.Filename.length() - Break);
        int SeedInt;
        istringstream iss(SeedStr);
        iss >> ws >> P.Seed >> ws;
        if(!iss.eof()) THROW_EXCEPTION("Could not understand random seed (in Filename).", 1)
        P.Filename = P.Filename.substr(0, Break);
    }
    if (P.Seed == INT_MIN) P.Seed = (int)time(NULL);
    *Params = P;
}

#define ETICS_HDF5
// make safety: if single precision, fail to compile
#ifdef ETICS_HDF5
#include "H5Cpp.h"
void Dooooo(string Prefix, int SnapNumber, Particle *P_h, int N, Real T) {
    
    char Filename[512];
    sprintf(Filename, "%s.h5part", Prefix.c_str());

    
//     WriteSnapshot(string Prefix, int SnapNumber, Particle *P_h, int N, Real T)
    double *X, *Y, *Z, *VX, *VY, *VZ, *AX, *AY, *AZ;
    X  = new double[N];
    Y  = new double[N];
    Z  = new double[N];
    VX = new double[N];
    VY = new double[N];
    VZ = new double[N];
    AX = new double[N];
    AY = new double[N];
    AZ = new double[N];
    int *ID;
    ID = new int[N];

    for (int i = 0; i < N; i++) {
        Particle p = P_h[i];
        ID[i] = p.ID;
        X[i]  = p.pos.x;
        Y[i]  = p.pos.y;
        Z[i]  = p.pos.z;
        VX[i] = p.vel.x;
        VY[i] = p.vel.y;
        VZ[i] = p.vel.z;
        AX[i] = p.acc.x;
        AY[i] = p.acc.y;
        AZ[i] = p.acc.z;
    }

    H5::H5File file;
    std::ifstream infile(Filename);
    if (!infile.good()) file = H5::H5File(Filename, H5F_ACC_TRUNC);
    else file = H5::H5File(Filename, H5F_ACC_RDWR);
#warning be more sophosticated when checking if file exists. Open it and if the snapshot we are trying to write already exist, then fail.
    
#warning also make it so when we read the ini file, if the ic file ends with .h5part, then we treat it as continue.
    
#warning you forgot to include the mass in the parameters

    char GroupName[64];
    sprintf(GroupName, "/Step#%d", SnapNumber);

    H5::Group group(file.createGroup(GroupName));

    H5::DataSpace attr_dataspace = H5::DataSpace(H5S_SCALAR);

    H5::Attribute attribute = group.createAttribute("Time", H5::PredType::NATIVE_DOUBLE, attr_dataspace);
    attribute.write( H5::PredType::NATIVE_DOUBLE, &T);

    attribute = group.createAttribute("TotalN", H5::PredType::NATIVE_UINT32, attr_dataspace);
    attribute.write( H5::PredType::NATIVE_UINT32, &N);

    hsize_t dimsxxx = N;
    H5::DataSpace dataspacexxx(1, &dimsxxx);
    H5::DataSet dataset3;
    dataset3 = H5::DataSet(group.createDataSet("ID", H5::PredType::NATIVE_UINT32, dataspacexxx));
    dataset3.write(ID, H5::PredType::NATIVE_UINT32);
    dataset3 = H5::DataSet(group.createDataSet("X",  H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(X, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("Y",  H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(Y, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("Z",  H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(Z, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("VX", H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(VX, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("VY", H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(VY, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("VZ", H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(VZ, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("AX", H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(AX, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("AY", H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(AY, H5::PredType::NATIVE_DOUBLE);
    dataset3 = H5::DataSet(group.createDataSet("AZ", H5::PredType::NATIVE_DOUBLE, dataspacexxx));
    dataset3.write(AZ, H5::PredType::NATIVE_DOUBLE);

    
    dataset3.close();
    dataspacexxx.close();

    attr_dataspace.close();
    attribute.close();
    group.close();
    file.close();

}

#endif
