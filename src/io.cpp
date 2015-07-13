#include <iostream>
#include <cstring>
#include <cstdlib>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#define __host__
#define __device__
#include "common.hpp"

using namespace std;

void ReadICs(string FileName, int N, int Skip, Particle **P_h) {

    *P_h = new Particle[N];

    string Line;
    ifstream InputFile(FileName.c_str());
    int LineNum;
    for (LineNum = 0; LineNum < Skip; LineNum++) getline(InputFile, Line);
    int i = 0;
    if (InputFile.is_open()) {
        while (getline(InputFile, Line))
        {
            Real m, x, y, z, vx, vy, vz;
            Particle p;
// #ifdef DOUBLE_PRECISION
//             int Count = sscanf(Line.c_str(), "%*d %*f %lf %lf %lf %lf %lf %lf", &x, &y, &z, &vx, &vy, &vz);
// #else
//             int Count = sscanf(Line.c_str(), "%*d %*f %f %f %f %f %f %f", &x, &y, &z, &vx, &vy, &vz);
// #endif
//             if (Count != 6) {
#ifdef DOUBLE_PRECISION
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
        cerr << "Can't open file" << FileName << "." << endl;
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

void WriteSnapshot(string Prefix, int SnapNumber, Particle *P_h, int N, Real T) {
    char S[512];
    sprintf(S, "%s%04d", Prefix.c_str(), SnapNumber);
    ofstream SnapshotFile;
    SnapshotFile.open(S);
    SnapshotFile << S << endl;
    sprintf(S, "%d\n", N); SnapshotFile << S;
    sprintf(S, "%.16E\n", T); SnapshotFile << S;
    for (int i = 0; i < N; i++) {
#warning Snapshot has fake mass
        sprintf(S, "%d 1E-6 %E %E %E %E %E %E\n", P_h[i].ID, P_h[i].pos.x, P_h[i].pos.y, P_h[i].pos.z, P_h[i].vel.x, P_h[i].vel.y, P_h[i].vel.z); SnapshotFile << S;
    }
    SnapshotFile.close();
}

#define THROW_EXCEPTION(str, stat) {cerr << str << endl; exit((stat));}
void ParseInput(int argc, char *argv[], ParametersStruct *Params) {
    ParametersStruct P;
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

    P.FileName = pt.get<string>("Filename", "\n");
    if (P.FileName == "\n") THROW_EXCEPTION("Could not read initial condition file name (Filename) from ini file.", 1)
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
    *Params = P;
}
