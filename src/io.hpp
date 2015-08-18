#pragma once
void ReadICs(std::string Filename, int N, Particle **P_h, int *FileSnapshotNum, Real *FileTime);
void WriteSnapshot(std::string Prefix, int SnapNumber, Particle *hostP, int N, Real T);
void ParseInput(int argc, char *argv[], ParametersStruct *Params);
void GenerateRandomData(int N, Particle *hostP);