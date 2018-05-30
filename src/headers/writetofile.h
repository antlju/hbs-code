#pragma once
std::string set_fname(std::string base, std::string end, Int step);
Int openStats(const std::string fname);
Int writeToFile(std::string fname, const Mesh &u, const Int vi, const Pencil &x, const Pencil &y, const Pencil &z);

Int writeToFile_1DArr(const std::string fname, const Mesh &u, const Int vi, const Grid &grid);
Int writeStatsToFile(const std::string fname, const MeshContainer &meshCntr, const Stats &stats, const Grid &grid,const Int t);
