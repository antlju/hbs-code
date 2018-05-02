#pragma once
std::string set_fname(std::string base, std::string end, Int step);
Int writeToFile(std::string fname, const Mesh &u, const Int vi, const Pencil &x, const Pencil &y, const Pencil &z);
