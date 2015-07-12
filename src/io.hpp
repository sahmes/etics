#pragma once
void ReadICs(std::string FileName, int N, int Skip, Particle **hostP);
void WriteSnapshot(std::string Prefix, int SnapNumber, Particle *hostP, int N, Real T);
void ParseInput(int argc, char *argv[], ParametersStruct *Params);
void GenerateRandomData(int N, Particle *hostP);